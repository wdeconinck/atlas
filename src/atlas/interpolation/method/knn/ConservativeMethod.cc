/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iomanip>
#include <vector>

#include "eckit/log/ProgressTimer.h"

#include "atlas/grid.h"
#include "atlas/interpolation/method/MethodFactory.h"
#include "atlas/interpolation/method/knn/ConservativeMethod.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildNode2CellConnectivity.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Topology.h"

#define USE_EDGE_CONNECTIVITY 1

#if USE_EDGE_CONNECTIVITY
#include "atlas/mesh/actions/BuildEdges.h"
#endif

namespace atlas {
namespace interpolation {
namespace method {

using CSPolygon      = util::ConvexSphericalPolygon;
using CSPolygonArray = ConservativeMethod::CSPolygonArray;

namespace {
MethodBuilder<ConservativeMethod> __builder("conservative");
}

ConservativeMethod::ConservativeMethod(const Config& config): Method(config) {
    config.get("order", order_ = 1);
    config.get("normalise_intersections", normalise_intersections_ = 1);
    config.get("field_value_type", fvtype_ = 0);
    config.get("matrix_free", matrix_free_ = true);
    config.get("src_cell_data", src_cell_data_ = true);
    config.get("tgt_cell_data", tgt_cell_data_ = true);
}

// get cyclically sorted edges from a cell
std::vector<idx_t> ConservativeMethod::sort_cell_edges(Mesh& mesh, idx_t cell_id) const {
    const auto& cell2edge = mesh.cells().edge_connectivity();
    const auto& cell2node = mesh.cells().node_connectivity();
    const auto& edge2node = mesh.edges().node_connectivity();
    const int nnodes      = cell2node.cols(cell_id);
    const int nedges      = cell2edge.cols(cell_id);
    std::vector<idx_t> edges;
    edges.resize(nedges);
    idx_t ii = 0;
    for (int inode = 0; inode < nnodes; ++inode) {
        idx_t node  = cell2node(cell_id, inode);
        idx_t nnode = cell2node(cell_id, (inode != nnodes - 1 ? inode + 1 : 0));
        for (int iedge = 0; iedge < nedges; ++iedge) {
            const idx_t edge  = cell2edge(cell_id, iedge);
            const idx_t node0 = edge2node(edge, 0);
            const idx_t node1 = edge2node(edge, 1);
            if ((node0 == node && node1 == nnode) or (node0 == nnode && node1 == node)) {
                edges[ii++] = edge;
                break;
            }
        }
    }
    return edges;
}

#if USE_EDGE_CONNECTIVITY

// get cyclically sorted edges from a node
std::vector<idx_t> ConservativeMethod::sort_node_edges(Mesh& mesh, idx_t node_id) const {
    const auto& node2edge = mesh.nodes().edge_connectivity();
    const auto& edge2cell = mesh.edges().cell_connectivity();
    const int nedges      = node2edge.cols(node_id);
    std::vector<idx_t> edges;
    edges.resize(nedges);
    if (nedges <= 2) {
        for (int i = 0; i < nedges; i++) {
            edges[i] = node2edge.row(node_id)(i);
        }
        return edges;
    }
    std::vector<std::array<idx_t, 3>> ecc;
    ecc.resize(nedges);
    idx_t count = 0;
    for (int iedge = 0; iedge < nedges; ++iedge) {
        edges[iedge] = -1;
        idx_t edge   = node2edge(node_id, iedge);
        ecc[count++] = std::array<idx_t, 3>({edge2cell(edge, 0), edge2cell(edge, 1), edge});
    }
    count          = 0;
    edges[count++] = ecc[0][2];
    idx_t prev     = ecc[0][0];  // previous cell
    idx_t find     = ecc[0][1];  // search cell
    idx_t ctrl     = 0;
    for (; ctrl < nedges; ++ctrl) {
        for (int i = 1; i < nedges; ++i) {
            if (ecc[i][0] == find and ecc[i][1] != prev) {
                edges[count++] = ecc[i][2];
                prev           = find;
                find           = ecc[i][1];
                continue;
            }
            if (ecc[i][1] == find and ecc[i][0] != prev) {
                edges[count++] = ecc[i][2];
                prev           = find;
                find           = ecc[i][0];
                continue;
            }
        }
    }
    return edges;
}

#endif

// get cyclically sorted neighbours of a cell
std::vector<idx_t> ConservativeMethod::get_cell_neighbours(Mesh& mesh, idx_t cell_id) const {
    const auto& cell2edge  = mesh.cells().edge_connectivity();
    const auto& edge2cell  = mesh.edges().cell_connectivity();
    auto c2e_missval       = cell2edge.missing_value();
    const auto& edges_sort = sort_cell_edges(mesh, cell_id);
    const idx_t nedges     = cell2edge.cols(cell_id);
    std::vector<idx_t> nbr_cells;
    nbr_cells.reserve(nedges);

    for (idx_t iedge = 0; iedge < nedges; ++iedge) {
        const idx_t edge  = edges_sort[iedge];
        const idx_t c1_id = edge2cell(edge, 0);
        if (c1_id != c2e_missval && c1_id != cell_id) {
            nbr_cells.emplace_back(c1_id);
            continue;
        }
        const idx_t c2_id = edge2cell(edge, 1);
        if (c2_id != c2e_missval && c2_id != cell_id) {
            nbr_cells.emplace_back(c2_id);
            continue;
        }
    }
    return nbr_cells;
}

#if USE_EDGE_CONNECTIVITY

// get cyclically sorted node neighbours using edge connectivity
std::vector<idx_t> ConservativeMethod::get_node_neighbours(Mesh& mesh, idx_t node_id) const {
    const auto& node2edge  = mesh.nodes().edge_connectivity();
    const auto& edge2node  = mesh.edges().node_connectivity();
    auto n2e_missval       = node2edge.missing_value();
    const auto& edges_sort = sort_node_edges(mesh, node_id);
    const int nedges       = node2edge.cols(node_id);
    std::vector<idx_t> nbr_nodes;
    nbr_nodes.reserve(nedges);
    for (idx_t iedge = 0; iedge < nedges; ++iedge) {
        const idx_t edge  = edges_sort[iedge];
        const idx_t n1_id = edge2node(edge, 0);
        if (n1_id != n2e_missval && n1_id != node_id) {
            nbr_nodes.emplace_back(n1_id);
            continue;
        }
        const idx_t n2_id = edge2node(edge, 1);
        if (n2_id != n2e_missval && n2_id != node_id) {
            nbr_nodes.emplace_back(n2_id);
            continue;
        }
    }
    ATLAS_ASSERT(nedges == nbr_nodes.size());
    return nbr_nodes;
}

#else 

// get cyclically sorted node neighbours without using edge connectivity
std::vector<idx_t> ConservativeMethod::get_node_neighbours(Mesh& mesh, idx_t node_id) const {
    const auto& node2cell  = mesh.nodes().cell_connectivity();
    const auto& cell2node  = mesh.cells().node_connectivity();
    std::vector<idx_t> nbr_nodes;
    std::vector<idx_t> nbr_nodes_od;
    const int ncells       = node2cell.cols( node_id );
	ATLAS_ASSERT( ncells > 0 );
    idx_t cnodes[ncells][2];
    nbr_nodes.reserve( ncells + 1 );
    nbr_nodes_od.reserve( ncells + 1 );
    for (idx_t icell = 0; icell < ncells; ++icell) {
        const idx_t cell = node2cell( node_id, icell );
		const int nnodes = cell2node.cols( cell );
		idx_t cnode = 0;
    	for (; cnode < nnodes; ++cnode) {
			if ( node_id == cell2node( cell, cnode ) ) {
				break;
			}
		}
		cnodes[icell][0] = cell2node( cell, (cnode!=0 ? cnode-1 : nnodes-1) );
		cnodes[icell][1] = cell2node( cell, (cnode!=nnodes-1 ? cnode+1 : 0) );
    }
	if ( ncells == 1 ) {
		nbr_nodes.emplace_back( cnodes[0][0] );
		nbr_nodes.emplace_back( cnodes[0][1] );
		return nbr_nodes;
	}
	// cycle one direction
	idx_t find = cnodes[0][1];
	idx_t prev = cnodes[0][0];
	nbr_nodes.emplace_back( prev );
	nbr_nodes.emplace_back( find );
    for (idx_t icycle = 0; nbr_nodes[0] != find; ) {
		idx_t jcell = 0;
    	for ( ; jcell < ncells; ++jcell) {
			idx_t ocell = (icycle + jcell + 1)%ncells;
			idx_t cand0 = cnodes[ocell][0];
			idx_t cand1 = cnodes[ocell][1];
			if ( find == cand0 && prev != cand1 ) {
				if ( cand1 == nbr_nodes[0] ) {
					return nbr_nodes;
				}
				nbr_nodes.emplace_back( cand1 );
				prev = find;
				find = cand1;
				break;
			}
			if ( find == cand1 && prev != cand0 ) {
				if ( cand0 == nbr_nodes[0] ) {
					return nbr_nodes;
				}
				nbr_nodes.emplace_back( cand0 );
				prev = find;
				find = cand0;
				break;
			}
		}
		if ( jcell == ncells ) { // not found
			if ( nbr_nodes[0] != find && find != nbr_nodes[ nbr_nodes.size()-1 ] ) {
				nbr_nodes.emplace_back( find );
			}
			break;
		}
		else {
			icycle++;
		}
	}
	if ( nbr_nodes[0] == find ) {
		return nbr_nodes;
	}
	// cycle the oposite direction
	find = cnodes[0][0];
	prev = cnodes[0][1];
	nbr_nodes_od.emplace_back( prev );
	nbr_nodes_od.emplace_back( find );
    for (idx_t icycle = 0; nbr_nodes_od[0] != find; ) {
		idx_t jcell = 0;
    	for ( ; jcell < ncells; ++jcell) {
			idx_t ocell = (icycle + jcell + 1)%ncells;
			if ( find == cnodes[ocell][0] && prev != cnodes[ocell][1] ) {
				nbr_nodes_od.emplace_back( cnodes[ocell][1] );
				prev = find;
				find = cnodes[ocell][1];
				break;
			}
			if ( find == cnodes[ocell][1] && prev != cnodes[ocell][0] ) {
				nbr_nodes_od.emplace_back( cnodes[ocell][0] );
				prev = find;
				find = cnodes[ocell][0];
				break;
			}
		}
		if ( jcell == ncells ) {
			if ( find != nbr_nodes_od[ nbr_nodes_od.size()-1 ] ) {
				nbr_nodes_od.emplace_back( find );
			}
			break;
		}
		icycle++;
	}
	// put together
	int ow_size = nbr_nodes_od.size();
	for( int i = 0; i < ow_size-2; i++ ) {
		nbr_nodes.emplace_back( nbr_nodes_od[ ow_size - 1 - i] );
	}
    return nbr_nodes;
}

#endif

// Create polygons for cell-centred data. Here, the polygons are mesh cells
CSPolygonArray ConservativeMethod::get_polygons_celldata(Mesh& mesh) const {
    CSPolygonArray cspolygons;
    const idx_t n_cells = mesh.cells().size();
    cspolygons.resize(n_cells);
    const auto& cell2node  = mesh.cells().node_connectivity();
    const auto lonlat      = array::make_view<double, 2>(mesh.nodes().lonlat());
    const auto cell_halo   = array::make_view<int, 1>(mesh.cells().halo());
    const auto& cell_flags = array::make_view<int, 1>(mesh.cells().flags());
    const auto& cell_part  = array::make_view<int, 1>(mesh.cells().partition());
    std::vector<PointLonLat> pts_ll;
    for (idx_t icell = 0; icell < n_cells; ++icell) {
        const idx_t n_nodes = cell2node.cols(icell);
        pts_ll.clear();
        pts_ll.resize(n_nodes);
        for (idx_t jnode = 0; jnode < n_nodes; ++jnode) {
            idx_t inode   = cell2node(icell, jnode);
            pts_ll[jnode] = PointLonLat{lonlat(inode, 0), lonlat(inode, 1)};
        }
        std::get<0>(cspolygons[icell]) = CSPolygon(pts_ll);
        int halo_type                  = cell_halo(icell);
        if (util::Bitflags::view(cell_flags(icell)).check(util::Topology::PERIODIC)) {
            halo_type = -1;
        }
        std::get<1>(cspolygons[icell]) = halo_type;
    }
    return cspolygons;
}

// Create polygons for cell-vertex data. Here, the polygons are subcells of mesh cells created as
// 	 (cell_centre, edge_centre, cell_vertex, edge_centre)
// additionally, subcell-to-node and node-to-subcells mapping are computed
CSPolygonArray ConservativeMethod::get_polygons_nodedata(Mesh& mesh, std::vector<idx_t>& csp2node,
                                                         std::vector<std::vector<idx_t>>& node2csp) const {
    CSPolygonArray cspolygons;
    csp2node.clear();
    node2csp.clear();
    node2csp.resize(mesh.nodes().size());
    const auto xy         = array::make_view<double, 2>(mesh.nodes().xy());
    const auto nodes_ll   = array::make_view<double, 2>(mesh.nodes().lonlat());
    auto edge_flags       = array::make_view<int, 1>(mesh.edges().flags());
    const auto& cell2edge = mesh.cells().edge_connectivity();
    const auto& cell2node = mesh.cells().node_connectivity();
    const auto& edge2node = mesh.edges().node_connectivity();
    const auto cell_halo  = array::make_view<int, 1>(mesh.cells().halo());
    const auto cell_flags = array::make_view<int, 1>(mesh.cells().flags());
    const auto cell_part  = array::make_view<int, 1>(mesh.cells().partition());
    const auto cell_gidx  = array::make_view<gidx_t, 1>(mesh.cells().global_index());
    const auto nodes_gidx = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
    auto cell_patch       = [&cell_flags](idx_t e) {
        using Topology = atlas::mesh::Nodes::Topology;
        return Topology::check(cell_flags(e), Topology::PATCH);
    };
    auto xyz2ll = [](atlas::PointXYZ& p_xyz) {
        PointLonLat p_ll;
        eckit::geometry::Sphere::convertCartesianToSpherical(1., p_xyz, p_ll);
        return p_ll;
    };
    idx_t cspol_id = 0; // subpolygon enumeration
    for (idx_t cell = 0; cell < mesh.cells().size(); ++cell) {
        ATLAS_ASSERT(cell < cell2edge.rows());
        ATLAS_ASSERT(cell < cell2node.rows());
        const idx_t n_nodes = cell2node.cols(cell);
        PointXYZ cell_mid(0., 0., 0.);	// cell centre
        for (idx_t inode = 0; inode < n_nodes; ++inode) {
            idx_t node0             = cell2node(cell, inode);
            idx_t node1             = cell2node(cell, inode!=n_nodes-1 ? inode+1 : 0);
            const PointLonLat p0_ll = PointLonLat{nodes_ll(node0, 0), nodes_ll(node0, 1)};
            const PointLonLat p1_ll = PointLonLat{nodes_ll(node1, 0), nodes_ll(node1, 1)};
            PointXYZ p0, p1;
            eckit::geometry::Sphere::convertSphericalToCartesian(1., p0_ll, p0);
            eckit::geometry::Sphere::convertSphericalToCartesian(1., p1_ll, p1);
            if (PointXYZ::norm(p0 - p1) < 1e-14) {
                continue;  // skip this edge, it is a pole point
            }
            cell_mid = cell_mid + p0;
            cell_mid = cell_mid + p1;
        }
        cell_mid = PointXYZ::div(cell_mid, PointXYZ::norm(cell_mid));
        PointLonLat cell_ll;
        eckit::geometry::Sphere::convertCartesianToSpherical(1., cell_mid, cell_ll);
        // get CSPolygon for each valid edge
        for (idx_t inode = 0; inode < n_nodes; ++inode) {
            idx_t node0              = cell2node(cell, inode);
            idx_t node1              = cell2node(cell, inode!=n_nodes-1 ? inode+1 : 0);
            const PointLonLat pi0_ll = PointLonLat{nodes_ll(node0, 0), nodes_ll(node0, 1)};
            const PointLonLat pi1_ll = PointLonLat{nodes_ll(node1, 0), nodes_ll(node1, 1)};
            PointXYZ pi0, pi1;
            eckit::geometry::Sphere::convertSphericalToCartesian(1., pi0_ll, pi0);
            eckit::geometry::Sphere::convertSphericalToCartesian(1., pi1_ll, pi1);
            if (PointXYZ::norm(pi0 - pi1) < 1e-14) {
                continue;  // skip this edge, it is a pole point
            }
            PointXYZ iedge_mid = pi0 + pi1;
            iedge_mid          = PointXYZ::div(iedge_mid, PointXYZ::norm(iedge_mid));
			csp2node.emplace_back(node1);
			node2csp[node0].emplace_back(cspol_id);
			idx_t node2 = cell2node(cell, inode<n_nodes-2 ? inode+2 : inode+2-n_nodes);	// the end point of the other real edge touching pi1
			auto pi2_ll = PointLonLat{nodes_ll(node2, 0), nodes_ll(node2, 1)};
			PointXYZ pi2;
			eckit::geometry::Sphere::convertSphericalToCartesian(1., pi2_ll, pi2);
			if ( PointXYZ::norm( pi1 - pi2 ) < 1e-14 ) { // we need real edge [pi1,pi2]
				node2 = cell2node(cell, inode<n_nodes-3 ? inode+3 : inode+3-n_nodes); 
				pi2_ll = PointLonLat{nodes_ll(node2, 0), nodes_ll(node2, 1)};
				eckit::geometry::Sphere::convertSphericalToCartesian(1., pi2_ll, pi2);
			}
			if ( PointXYZ::norm( pi1 - pi2 ) < 1e-14 ) {
				ATLAS_THROW_EXCEPTION("Three cell vertices on a same great arc!");
			}
            PointXYZ jedge_mid;
			jedge_mid = pi1 + pi2;
            jedge_mid = PointXYZ::div(jedge_mid, PointXYZ::norm(jedge_mid));	
            std::vector<PointLonLat> pts_ll(4);
            pts_ll[0]     = cell_ll;
            pts_ll[1]     = xyz2ll(iedge_mid);
            pts_ll[2]     = pi1_ll;
            pts_ll[3]     = xyz2ll(jedge_mid);
            int halo_type = cell_halo(cell);
            if (util::Bitflags::view(cell_flags(cell)).check(util::Topology::PERIODIC)) {
                halo_type = -1;
            }
            cspolygons.emplace_back(CSPolygon(pts_ll), halo_type);
            cspol_id++;
        }
    }
    Log::info() << "ConservativeMethod::get_polygons_nodedata : Created " << cspolygons.size() << " CSPolygons from "
                << mesh.cells().size() << " mesh cells\n";
    return cspolygons;
}

void ConservativeMethod::do_setup(const Grid& src_grid, const Grid& tgt_grid) {
    ATLAS_TRACE("ConservativeMethod::do_setup( Grid, Grid )");
    ATLAS_ASSERT(src_grid);
    ATLAS_ASSERT(tgt_grid);
    auto src_mesh_config = src_grid.meshgenerator();
    auto tgt_mesh_config = tgt_grid.meshgenerator();
    tgt_mesh_            = MeshGenerator(tgt_mesh_config).generate(tgt_grid);
    functionspace::NodeColumns tmp_tgt_fs(tgt_mesh_, option::halo(0));
    if (mpi::size() > 1) {
        src_mesh_ = MeshGenerator(src_mesh_config).generate(src_grid, grid::MatchingPartitioner(tgt_mesh_));
    }
    else {
        src_mesh_ = MeshGenerator(src_mesh_config).generate(src_grid);
    }
    functionspace::NodeColumns tmp_src_fs(src_mesh_, option::halo(2));
    mesh::actions::build_edges(src_mesh_, util::Config("pole_edges", false));
    if (not src_cell_data_) {
        mesh::actions::build_node_to_edge_connectivity(src_mesh_);
    }
    mesh::actions::build_edges(tgt_mesh_, util::Config("pole_edges", false));
    if (src_cell_data_) {
        functionspace::CellColumns src_fs(src_mesh_, option::halo(2));
        src_fs_ = src_fs;
    }
    else {
        functionspace::NodeColumns src_fs(src_mesh_, option::halo(2));
        src_fs_ = src_fs;
    }
    if (tgt_cell_data_) {
        functionspace::CellColumns tgt_fs(tgt_mesh_, option::halo(0));
        tgt_fs_ = tgt_fs;
    }
    else {
        functionspace::NodeColumns tgt_fs(tgt_mesh_, option::halo(0));
        tgt_fs_ = tgt_fs;
    }
    do_setup(src_fs_, tgt_fs_);
}

void ConservativeMethod::do_setup(const FunctionSpace& src_fs, const FunctionSpace& tgt_fs) {
    ATLAS_TRACE("ConservativeMethod::do_setup( FunctionSpace, FunctionSpace )");
    ATLAS_ASSERT(src_fs);
    ATLAS_ASSERT(tgt_fs);
    if (functionspace::CellColumns(src_fs)) {
        src_cell_data_ = true;
        src_fs_        = functionspace::CellColumns(src_fs);
        src_mesh_      = functionspace::CellColumns(src_fs).mesh();
    }
    else if (functionspace::NodeColumns(src_fs)) {
        src_cell_data_ = false;
        src_fs_        = functionspace::NodeColumns(src_fs);
        src_mesh_      = functionspace::NodeColumns(src_fs).mesh();
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
    if (functionspace::CellColumns(tgt_fs)) {
        tgt_cell_data_ = true;
        tgt_fs_        = functionspace::CellColumns(tgt_fs);
        tgt_mesh_      = functionspace::CellColumns(tgt_fs).mesh();
    }
    else if (functionspace::NodeColumns(tgt_fs)) {
        tgt_cell_data_ = false;
        tgt_fs_        = functionspace::NodeColumns(tgt_fs);
        tgt_mesh_      = functionspace::NodeColumns(tgt_fs).mesh();
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }

    {
        // TODO: Check if this is still required, and if so, work needs to be done to make it disappear
        functionspace::NodeColumns tmp_tgt_fs(tgt_mesh_, option::halo(0));
    }

    auto src_grid        = src_mesh_.grid();
    auto src_mesh_config = src_grid.meshgenerator();
    if (mpi::size() > 1) {
        src_mesh_ = MeshGenerator(src_mesh_config).generate(src_grid, grid::MatchingPartitioner(tgt_mesh_));
    }
    else {
        src_mesh_ = MeshGenerator(src_mesh_config).generate(src_grid);
    }
    functionspace::NodeColumns tmp_src_fs(src_mesh_, option::halo(2));

    {
        // TODO: make everything in this scope unnecessary, relying only on cells and nodes.
        mesh::actions::build_edges(src_mesh_, util::Config("pole_edges", false));
        if (not src_cell_data_) {
            mesh::actions::build_node_to_edge_connectivity(src_mesh_);
        }
        mesh::actions::build_edges(tgt_mesh_, util::Config("pole_edges", false));
    }

	{
		// todo: maybe do not need to build this connectivity in all cases
		if ( not src_cell_data_ ) {
    		mesh::actions::build_node_to_cell_connectivity(src_mesh_);
		}
		if ( not tgt_cell_data_ ) {
    		mesh::actions::build_node_to_cell_connectivity(tgt_mesh_);
		}
	}

    CSPolygonArray src_csp;
    CSPolygonArray tgt_csp;
    {
        ATLAS_TRACE("Get source polygons");
        if (src_cell_data_) {
            functionspace::CellColumns src_fs(src_mesh_, option::halo(2));
            src_fs_ = src_fs;
            src_csp = get_polygons_celldata(src_mesh_);
        }
        else {
            functionspace::NodeColumns src_fs(src_mesh_, option::halo(2));
            src_fs_ = src_fs;
            src_csp = get_polygons_nodedata(src_mesh_, src_csp2node_, src_node2csp_);
        }
    }
    {
        ATLAS_TRACE("Get target polygons");
        if (tgt_cell_data_) {
            functionspace::CellColumns tgt_fs(tgt_mesh_, option::halo(0));
            tgt_fs_ = tgt_fs;
            tgt_csp = get_polygons_celldata(tgt_mesh_);
        }
        else {
            functionspace::NodeColumns tgt_fs(tgt_mesh_, option::halo(0));
            tgt_fs_ = tgt_fs;
            tgt_csp = get_polygons_nodedata(tgt_mesh_, tgt_csp2node_, tgt_node2csp_);
        }
    }
    intersect_polygons(src_csp, tgt_csp);

    n_spoints_ = (src_cell_data_ ? src_mesh_.cells().size() : src_mesh_.nodes().size());
    n_tpoints_ = (tgt_cell_data_ ? tgt_mesh_.cells().size() : tgt_mesh_.nodes().size());
    src_points_.resize(n_spoints_);
    tgt_points_.resize(n_tpoints_);
    src_areas_       = src_fs_.createField<double>();
    auto src_areas_v = array::make_view<double, 1>(src_areas_);
    if (src_cell_data_) {
        for (idx_t spt = 0; spt < n_spoints_; ++spt) {
            const auto& s_csp = std::get<0>(src_csp[spt]);
            src_points_[spt]  = s_csp.centroid();
            src_areas_v(spt)  = s_csp.area();
        }
    }
    else {
        const auto lonlat = array::make_view<double, 2>(src_mesh_.nodes().lonlat());
        for (idx_t spt = 0; spt < n_spoints_; ++spt) {
            auto p = PointLonLat{lonlat(spt, 0), lonlat(spt, 1)};
            eckit::geometry::Sphere::convertSphericalToCartesian(1., p, src_points_[spt]);
            src_points_[spt] = PointXYZ{0., 0., 0.};
            src_areas_v(spt) = 0.;
            for (idx_t isubcell = 0; isubcell < src_node2csp_[spt].size(); ++isubcell) {
                idx_t subcell     = src_node2csp_[spt][isubcell];
                const auto& s_csp = std::get<0>(src_csp[subcell]);
                src_areas_v(spt) += s_csp.area();
                src_points_[spt] = src_points_[spt] + PointXYZ::mul(s_csp.centroid(), s_csp.area());
            }
            double src_point_norm = PointXYZ::norm(src_points_[spt]);
            src_points_[spt]      = PointXYZ::div(src_points_[spt], (src_point_norm > 1e-16 ? src_point_norm : 1.));
        }
    }
    tgt_areas_       = tgt_fs_.createField<double>();
    auto tgt_areas_v = array::make_view<double, 1>(tgt_areas_);
    if (tgt_cell_data_) {
        for (idx_t tpt = 0; tpt < n_tpoints_; ++tpt) {
            const auto& t_csp = std::get<0>(tgt_csp[tpt]);
            tgt_points_[tpt]  = t_csp.centroid();
            tgt_areas_v(tpt)  = t_csp.area();
        }
    }
    else {
        const auto lonlat = array::make_view<double, 2>(tgt_mesh_.nodes().lonlat());
        for (idx_t tpt = 0; tpt < n_tpoints_; ++tpt) {
            auto p = PointLonLat{lonlat(tpt, 0), lonlat(tpt, 1)};
            eckit::geometry::Sphere::convertSphericalToCartesian(1., p, tgt_points_[tpt]);
            tgt_points_[tpt] = PointXYZ{0., 0., 0.};
            tgt_areas_v(tpt) = 0.;
            for (idx_t isubcell = 0; isubcell < tgt_node2csp_[tpt].size(); ++isubcell) {
                idx_t subcell = tgt_node2csp_[tpt][isubcell];
                if (std::get<1>(tgt_csp[subcell]) == 0) {
                    const auto& t_csp = std::get<0>(tgt_csp[subcell]);
                    tgt_areas_v(tpt) += t_csp.area();
                    tgt_points_[tpt] = tgt_points_[tpt] + PointXYZ::mul(t_csp.centroid(), t_csp.area());
                }
            }
            double tgt_point_norm = PointXYZ::norm(tgt_points_[tpt]);
            tgt_points_[tpt]      = PointXYZ::div(tgt_points_[tpt], (tgt_point_norm > 1e-16 ? tgt_point_norm : 1.));
        }
    }
    //src_areas_.set_dirty( true );
    //src_areas_.haloExchange();
    //tgt_areas_.set_dirty( true );
    //tgt_areas_.haloExchange();
    setup_1st_order_matrix();
    setup_2nd_order_matrix();
}


void ConservativeMethod::intersect_polygons(const CSPolygonArray& src_csp, const CSPolygonArray& tgt_csp) {
    ATLAS_TRACE();
    util::KDTree<idx_t> kdt_search;
    kdt_search.reserve(tgt_csp.size());

    double max_tgtcell_rad = 0.;
    for (idx_t jcell = 0; jcell < tgt_csp.size(); ++jcell) {
        if (std::get<1>(tgt_csp[jcell]) == 0) {
            const auto& t_csp = std::get<0>(tgt_csp[jcell]);
            kdt_search.insert(t_csp.centroid(), jcell);
            max_tgtcell_rad = std::max(max_tgtcell_rad, t_csp.cell_radius());
        }
    }
    kdt_search.build();

    size_t nonintersect        = 0;
    size_t n_icsp              = 0;
    double src_area_notcovered = 0.;
    iparam_.resize(src_csp.size());
    eckit::Channel blackhole;
    eckit::ProgressTimer progress("Intersecting polygons ", src_csp.size(), " cell", double(10),
                                  src_csp.size() > 50 ? Log::info() : blackhole);
    for (idx_t scell = 0; scell < src_csp.size(); ++scell, ++progress) {
        if (std::get<1>(src_csp[scell]) == -1) {
            // skip periodic cells
            continue;
        }
        const auto& s_csp   = std::get<0>(src_csp[scell]);
        double covered_area = 0.;
        auto tgt_cells = kdt_search.closestPointsWithinRadius(s_csp.centroid(), s_csp.cell_radius() + max_tgtcell_rad);
        for (idx_t ttcell = 0; ttcell < tgt_cells.size(); ++ttcell) {
            auto tcell        = tgt_cells[ttcell].payload();
            const auto& t_csp = std::get<0>(tgt_csp[tcell]);
            CSPolygon csp_i   = s_csp.intersect(t_csp);
            if (csp_i.area() > 0.) {
                iparam_[scell].tcell_id.emplace_back(tcell);
                iparam_[scell].weights.emplace_back(csp_i.area());
                iparam_[scell].sweights.emplace_back(csp_i.area() / t_csp.area());
                iparam_[scell].centroids.emplace_back(csp_i.centroid());
                covered_area += csp_i.area();
            }
        }
        const double loc_csp_error = std::abs(s_csp.area() - covered_area) / s_csp.area();
        src_area_notcovered += loc_csp_error;
        if (iparam_[scell].tcell_id.size() == 0.) {
            ++nonintersect;
        }
        if (normalise_intersections_ && loc_csp_error < 1e-4) {
            double wfactor = s_csp.area() / (covered_area > 1e-10 ? covered_area : 1.);
            for (idx_t i = 0; i < iparam_[scell].weights.size(); i++) {
                iparam_[scell].weights[i] *= wfactor;
                iparam_[scell].sweights[i] *= wfactor;
            }
        }
        if (false && loc_csp_error > 1e-5) {
            Log::info() << "* src cell area NOT covered: " << loc_csp_error << "\n";
            //dump_intersection( s_csp, tgt_csp, tgt_cells );
            //ATLAS_ASSERT( false );
        }
        n_icsp += iparam_[scell].weights.size();
    }
    Log::info() << "ConservativeMethod::intersect_polygons : size of src_grid, tgt_grid, supergrid: " << src_csp.size()
                << " " << tgt_csp.size() << " " << n_icsp << "\n";
    Log::info() << "ConservativeMethod::intersect_polygons : " << nonintersect
                << " source mesh polygons do NOT intersect any other polygon.\n";
    Log::info() << "ConservativeMethod::intersect_polygons : " << src_area_notcovered
                << " area of source mesh NOT covered by target mesh.\n";
    geo_err_intsc_l1_   = 0.;
    geo_err_intsc_linf_ = 0.;
    size_t no_iplg      = 0;
    for (idx_t scell = 0; scell < src_csp.size(); ++scell) {
        const int cell_flag = std::get<1>(src_csp[scell]);
        if (cell_flag == -1 or cell_flag > 0) {
            // skip periodic cells
            continue;
        }
        double diff_cell = std::get<0>(src_csp[scell]).area();
        for (idx_t icell = 0; icell < iparam_[scell].weights.size(); ++icell) {
            diff_cell -= iparam_[scell].weights[icell];
        }
        no_iplg += iparam_[scell].weights.size();
        geo_err_intsc_l1_ += std::abs(diff_cell);
        geo_err_intsc_linf_ = std::max(geo_err_intsc_linf_, std::abs(diff_cell));
    }
    geo_err_intsc_l1_ *= 0.25 * M_1_PI;
    Log::info() << "ConservativeMethod::intersect_polygons : cons err in polygon intersect  : (L1) "
                << geo_err_intsc_l1_ << " (Lmax) " << geo_err_intsc_linf_ << "\n";
}

void ConservativeMethod::setup_1st_order_matrix() {
    if (order_ != 1 or matrix_free_) {
        return;
    }
    ATLAS_TRACE("ConservativeMethod::setup: build cons-1 interpolant matrix");
    Triplets triplets;
    size_t triplets_size = 0;
    // determine the size of array of triplets used to define the sparse matrix
    if (src_cell_data_) {
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            triplets_size += iparam_[scell].centroids.size();
        }
    }
    else {
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                idx_t subcell = src_node2csp_[snode][isubcell];
                triplets_size += iparam_[subcell].sweights.size();
            }
        }
    }
    triplets.reserve(triplets_size);
    // assemble triplets to define the sparse matrix
    const auto src_areas_v = array::make_view<double, 1>(src_areas_);
    const auto tgt_areas_v = array::make_view<double, 1>(tgt_areas_);
    if (src_cell_data_ && tgt_cell_data_) {
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            const auto& iparam = iparam_[scell];
            for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                idx_t tcell = iparam.tcell_id[icell];
                triplets.emplace_back(tcell, scell, iparam.sweights[icell]);
            }
        }
    }
    else if (not src_cell_data_ && tgt_cell_data_) {
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                const idx_t subcell = src_node2csp_[snode][isubcell];
                const auto& iparam  = iparam_[subcell];
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    idx_t tcell = iparam.tcell_id[icell];
                    triplets.emplace_back(tcell, snode, iparam.sweights[icell]);
                }
            }
        }
    }
    else if (src_cell_data_ && not tgt_cell_data_) {
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            const auto& iparam = iparam_[scell];
            for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                idx_t tcell = iparam.tcell_id[icell];
                idx_t tnode = tgt_csp2node_[tcell];
                triplets.emplace_back(tnode, scell, iparam.weights[icell] / tgt_areas_v(tnode));
            }
        }
    }
    else if (not src_cell_data_ && not tgt_cell_data_) {
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                const idx_t subcell = src_node2csp_[snode][isubcell];
                const auto& iparam  = iparam_[subcell];
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    idx_t tcell = iparam.tcell_id[icell];
                    idx_t tnode = tgt_csp2node_[tcell];
                    triplets.emplace_back(tnode, snode, iparam.weights[icell] / tgt_areas_v(tnode));
                }
            }
        }
    }
    std::sort(std::begin(triplets), std::end(triplets), [](const Triplet& t1, const Triplet& t2) {
        return (t1.row() < t2.row() or (t1.row() == t2.row() and (t1.col() < t2.col())));
    });
    Matrix A(n_tpoints_, n_spoints_, triplets);
    matrix_shared_->swap(A);
}

void ConservativeMethod::setup_2nd_order_matrix() {
    if (order_ != 2 or matrix_free_) {
        return;
    }
    ATLAS_TRACE("ConservativeMethod::setup: build cons-2 interpolant matrix");
    Triplets triplets;
    size_t triplets_size   = 0;
    const auto tgt_areas_v = array::make_view<double, 1>(tgt_areas_);
    if (src_cell_data_) {
        const auto halo = array::make_view<int, 1>(src_mesh_.cells().halo());
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            if (halo(scell)) {
                continue;
            }
            const auto nb_cells = get_cell_neighbours(src_mesh_, scell);
            triplets_size += (2 * nb_cells.size() + 1) * iparam_[scell].centroids.size();
        }
        triplets.reserve(triplets_size);
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            const auto nb_cells = get_cell_neighbours(src_mesh_, scell);
            const auto& iparam  = iparam_[scell];
            if (iparam.centroids.size() == 0 && not halo(scell)) {
                Log::info() << " WARNING source cell " << scell << " not covered"
                            << "\n";
                continue;
            }
            /* // better conservation after Kritsikis et al. (2017)
            PointXYZ Cs = {0., 0., 0.};
            for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
                Cs = Cs + PointXYZ::mul( iparam.centroids[icell], iparam.weights[icell] );
            }
            const double Cs_norm = PointXYZ::norm( Cs );
            Cs = PointXYZ::div( Cs, Cs_norm );
            ATLAS_ASSERT( Cs_norm > 0. );
			*/
            const PointXYZ& Cs = src_points_[scell];
            // compute gradient from cells
            double dual_area_inv = 0.;
            std::vector<PointXYZ> Rsj;
            Rsj.resize(nb_cells.size());
            for (idx_t j = 0; j < nb_cells.size(); ++j) {
                idx_t nj         = (j != nb_cells.size() - 1) ? j + 1 : 0;
                idx_t sj         = nb_cells[j];
                idx_t nsj        = nb_cells[nj];
                const auto& Csj  = src_points_[sj];
                const auto& Cnsj = src_points_[nsj];
                if (CSPolygon::leftOf(Cnsj, Cs, Csj, 1e-16, 0)) {
                    Rsj[j] = PointXYZ::cross(Cnsj, Csj);
                    dual_area_inv += CSPolygon({Cs, Csj, Cnsj}).area();
                }
                else {
                    Rsj[j] = PointXYZ::cross(Csj, Cnsj);
                    dual_area_inv += CSPolygon({Cs, Cnsj, Csj}).area();
                }
            }
            dual_area_inv = (dual_area_inv > 0.) ? 1. / dual_area_inv : 1.;
            PointXYZ Rs   = {0., 0., 0.};
            for (idx_t j = 0; j < nb_cells.size(); ++j) {
                Rs = Rs + Rsj[j];
            }
            // assemble the matrix
            std::vector<PointXYZ> Aik;
            Aik.resize(iparam.centroids.size());
            for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                const PointXYZ& Csk   = iparam.centroids[icell];
                const PointXYZ Csk_Cs = Csk - Cs;
                Aik[icell]            = Csk_Cs - PointXYZ::mul(Cs, PointXYZ::dot(Cs, Csk_Cs));
                Aik[icell]            = PointXYZ::mul(Aik[icell], iparam.sweights[icell] * dual_area_inv);
            }
            if (tgt_cell_data_) {
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    const idx_t tcell = iparam.tcell_id[icell];
                    for (idx_t j = 0; j < nb_cells.size(); ++j) {
                        idx_t nj  = (j != nb_cells.size() - 1) ? j + 1 : 0;
                        idx_t sj  = nb_cells[j];
                        idx_t nsj = nb_cells[nj];
                        triplets.emplace_back(tcell, sj, 0.5 * PointXYZ::dot(Rsj[j], Aik[icell]));
                        triplets.emplace_back(tcell, nsj, 0.5 * PointXYZ::dot(Rsj[j], Aik[icell]));
                    }
                    triplets.emplace_back(tcell, scell, iparam.sweights[icell] - PointXYZ::dot(Rs, Aik[icell]));
                }
            }
            else {
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    idx_t tcell                = iparam.tcell_id[icell];
                    idx_t tnode                = tgt_csp2node_[tcell];
                    const double csp2node_coef = iparam.weights[icell] / iparam.sweights[icell] / tgt_areas_v(tnode);
                    for (idx_t j = 0; j < nb_cells.size(); ++j) {
                        idx_t nj  = (j != nb_cells.size() - 1) ? j + 1 : 0;
                        idx_t sj  = nb_cells[j];
                        idx_t nsj = nb_cells[nj];
                        triplets.emplace_back(tnode, sj, (0.5 * PointXYZ::dot(Rsj[j], Aik[icell])) * csp2node_coef);
                        triplets.emplace_back(tnode, nsj, (0.5 * PointXYZ::dot(Rsj[j], Aik[icell])) * csp2node_coef);
                    }
                    triplets.emplace_back(tnode, scell,
                                          (iparam.sweights[icell] - PointXYZ::dot(Rs, Aik[icell])) * csp2node_coef);
                }
            }
        }
    }
    else {  // if ( not src_cell_data_ )
        const auto src_halo = array::make_view<int, 1>(src_mesh_.nodes().halo());
        //       const auto glidx = array::make_view<gidx_t, 1>( src_mesh_.nodes().global_index() );
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            const auto nb_nodes = get_node_neighbours(src_mesh_, snode);
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                idx_t subcell = src_node2csp_[snode][isubcell];
                triplets_size += (2 * nb_nodes.size() + 1) * iparam_[subcell].centroids.size();
            }
        }
        triplets.reserve(triplets_size);
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            const auto nb_nodes = get_node_neighbours(src_mesh_, snode);
            // get the barycentre of the dual cell
            /* // better conservation
            PointXYZ Cs = {0., 0., 0.};
            for ( idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell ) {
                idx_t subcell      = src_node2csp_[snode][isubcell];
                const auto& iparam = iparam_[subcell];
                for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
                    Cs = Cs + PointXYZ::mul( iparam.centroids[icell], iparam.weights[icell] );
                }
            }
            const double Cs_norm = PointXYZ::norm( Cs );
            ATLAS_ASSERT( Cs_norm > 0. );
            Cs = PointXYZ::div( Cs, Cs_norm );
*/
            const PointXYZ& Cs = src_points_[snode];
            // compute gradient from nodes
            double dual_area_inv = 0.;
            std::vector<PointXYZ> Rsj;
            Rsj.resize(nb_nodes.size());
            const auto& Ns = src_points_[snode];
            for (idx_t j = 0; j < nb_nodes.size(); ++j) {
                idx_t nj         = (j != nb_nodes.size() - 1) ? j + 1 : 0;
                idx_t sj         = nb_nodes[j];
                idx_t snj        = nb_nodes[nj];
                const auto& Nsj  = src_points_[sj];
                const auto& Nsnj = src_points_[snj];
                if (CSPolygon::leftOf(Nsnj, Ns, Nsj, 1e-16, 0)) {
                    Rsj[j] = PointXYZ::cross(Nsnj, Nsj);
                    dual_area_inv += CSPolygon({Ns, Nsj, Nsnj}).area();
                }
                else {
                    Rsj[j] = PointXYZ::cross(Nsj, Nsnj);
                    dual_area_inv += CSPolygon({Ns, Nsnj, Nsj}).area();
                }
            }
            dual_area_inv = (dual_area_inv > 0.) ? 1. / dual_area_inv : 1.;
            PointXYZ Rs   = {0., 0., 0.};
            for (idx_t j = 0; j < nb_nodes.size(); ++j) {
                Rs = Rs + Rsj[j];
            }
            // assemble the matrix
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                idx_t subcell      = src_node2csp_[snode][isubcell];
                const auto& iparam = iparam_[subcell];
                if (iparam.centroids.size() == 0 and not src_halo(snode)) {
                    Log::info() << " WARNING source subcell around node " << snode << " not covered "
                                << "\n";
                    continue;
                }
                std::vector<PointXYZ> Aik;
                Aik.resize(iparam.centroids.size());
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    const PointXYZ& Csk   = iparam.centroids[icell];
                    const PointXYZ Csk_Cs = Csk - Cs;
                    Aik[icell]            = Csk_Cs - PointXYZ::mul(Cs, PointXYZ::dot(Cs, Csk_Cs));
                    Aik[icell]            = PointXYZ::mul(Aik[icell], iparam.sweights[icell] * dual_area_inv);
                }
                if (tgt_cell_data_) {
                    for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                        const idx_t tcell = iparam.tcell_id[icell];
                        for (idx_t j = 0; j < nb_nodes.size(); ++j) {
                            idx_t nj  = (j != nb_nodes.size() - 1) ? j + 1 : 0;
                            idx_t sj  = nb_nodes[j];
                            idx_t snj = nb_nodes[nj];
                            triplets.emplace_back(tcell, sj, 0.5 * PointXYZ::dot(Rsj[j], Aik[icell]));
                            triplets.emplace_back(tcell, snj, 0.5 * PointXYZ::dot(Rsj[j], Aik[icell]));
                        }
                        triplets.emplace_back(tcell, snode, iparam.sweights[icell] - PointXYZ::dot(Rs, Aik[icell]));
                    }
                }
                else {
                    for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                        idx_t tcell = iparam.tcell_id[icell];
                        idx_t tnode = tgt_csp2node_[tcell];
                        const double csp2node_coef =
                            iparam.weights[icell] / iparam.sweights[icell] / tgt_areas_v(tnode);
                        for (idx_t j = 0; j < nb_nodes.size(); ++j) {
                            idx_t nj  = (j != nb_nodes.size() - 1) ? j + 1 : 0;
                            idx_t sj  = nb_nodes[j];
                            idx_t snj = nb_nodes[nj];
                            triplets.emplace_back(tnode, sj, (0.5 * PointXYZ::dot(Rsj[j], Aik[icell])) * csp2node_coef);
                            triplets.emplace_back(tnode, snj,
                                                  (0.5 * PointXYZ::dot(Rsj[j], Aik[icell])) * csp2node_coef);
                        }
                        triplets.emplace_back(tnode, snode,
                                              (iparam.sweights[icell] - PointXYZ::dot(Rs, Aik[icell])) * csp2node_coef);
                    }
                }
            }
        }
    }
    std::sort(std::begin(triplets), std::end(triplets), [](const Triplet& t1, const Triplet& t2) {
        return (t1.row() < t2.row() or (t1.row() == t2.row() and t1.col() < t2.col()));
    });
    Matrix A(n_tpoints_, n_spoints_, triplets);
    matrix_shared_->swap(A);
}

void ConservativeMethod::do_execute(const Field& src_field, Field& tgt_field) {
    ATLAS_TRACE("ConservativeMethod::do_execute()");
    {
        ATLAS_TRACE("halo exchange source");
        src_field.set_dirty(true);
        src_field.haloExchange();
    }
    const auto tgt_areas_v = array::make_view<double, 1>(tgt_areas_);

    if (order_ == 1) {
        ATLAS_TRACE("order 1");
        if (matrix_free_) {
            if (not src_cell_data_ or not tgt_cell_data_) {
                ATLAS_NOTIMPLEMENTED;
            }

            const auto src_vals = array::make_view<double, 1>(src_field);
            auto tgt_vals       = array::make_view<double, 1>(tgt_field);

            for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                tgt_vals(tcell) = 0.;
            }
            for (idx_t scell = 0; scell < src_vals.size(); ++scell) {
                const auto& iparam = iparam_[scell];
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    tgt_vals(iparam.tcell_id[icell]) += iparam.weights[icell] * src_vals(scell);
                }
            }
            for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                tgt_vals(tcell) /= tgt_areas_v(tcell);
            }
        }
        else {
            Method::do_execute(src_field, tgt_field);
        }
    }
    else if (order_ == 2) {
        ATLAS_TRACE("order 2");
        if (matrix_free_) {
            if (not src_cell_data_ or not tgt_cell_data_) {
                ATLAS_NOTIMPLEMENTED;
            }
            const auto src_vals       = array::make_view<double, 1>(src_field);
            auto tgt_vals             = array::make_view<double, 1>(tgt_field);
            const auto& src_cell2edge = src_mesh_.cells().edge_connectivity();
            const auto& src_edge2node = src_mesh_.edges().node_connectivity();
            const auto halo           = array::make_view<int, 1>(src_mesh_.cells().halo());
            for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                tgt_vals(tcell) = 0.;
            }
            for (idx_t scell = 0; scell < src_vals.size(); ++scell) {
                if (halo(scell)) {
                    continue;
                }
                const auto& iparam       = iparam_[scell];
                const PointXYZ& P        = src_points_[scell];
                PointXYZ grad            = {0., 0., 0.};
                PointXYZ src_barycenter  = {0., 0., 0.};
                auto src_neighbour_cells = get_cell_neighbours(src_mesh_, scell);
                double dual_area         = 0.;
                for (idx_t nb_id = 0; nb_id < src_neighbour_cells.size(); ++nb_id) {
                    idx_t nnb_id    = (nb_id != src_neighbour_cells.size() - 1) ? nb_id + 1 : 0;
                    idx_t ncell     = src_neighbour_cells[nb_id];
                    idx_t nncell    = src_neighbour_cells[nnb_id];
                    const auto& Pn  = src_points_[ncell];
                    const auto& Pnn = src_points_[nncell];
                    if (ncell != scell && nncell != scell) {
                        double val = 0.5 * (src_vals(ncell) + src_vals(nncell)) - src_vals(scell);
                        auto csp   = CSPolygon({Pn, Pnn, P});
                        if (csp.area() < std::numeric_limits<double>::epsilon()) {
                            csp = CSPolygon({Pn, P, Pnn});
                        }
                        val *= (csp.leftOf(Pnn, P, Pn, 1e-16, 0) ? -1 : 1);
                        dual_area += std::abs(csp.area());
                        grad = grad + PointXYZ::mul(PointXYZ::cross(Pn, Pnn), val);
                    }
                    else if (ncell != scell) {
                        ATLAS_NOTIMPLEMENTED;
                        //double val = 0.5 * ( src_vals( ncell ) - src_vals( scell ) );
                        //grad = grad + PointXYZ::mul( PointXYZ::cross( Pn, P ), val );
                    }
                    else if (nncell != scell) {
                        ATLAS_NOTIMPLEMENTED;
                        //double val = 0.5 * ( src_vals( nncell ) - src_vals( scell ) );
                        //grad = grad + PointXYZ::mul( PointXYZ::cross( P, Pnn ), val );
                    }
                }
                if (dual_area > std::numeric_limits<double>::epsilon()) {
                    grad = PointXYZ::div(grad, dual_area);
                }
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    src_barycenter = src_barycenter + PointXYZ::mul(iparam.centroids[icell], iparam.weights[icell]);
                }
                src_barycenter = PointXYZ::div(src_barycenter, PointXYZ::norm(src_barycenter));
                grad           = grad - PointXYZ::mul(src_barycenter, PointXYZ::dot(grad, src_barycenter));
                ATLAS_ASSERT(std::abs(PointXYZ::dot(grad, src_barycenter)) < 1e-14);
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    tgt_vals(iparam.tcell_id[icell]) +=
                        iparam.weights[icell] *
                        (src_vals(scell) + PointXYZ::dot(grad, iparam.centroids[icell] - src_barycenter));
                }
            }
            for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                tgt_vals(tcell) /= tgt_areas_v(tcell);
            }
        }
        else {
            Method::do_execute(src_field, tgt_field);
        }
    }

    {
        ATLAS_TRACE("halo exchange target");
        tgt_field.set_dirty(true);
        tgt_field.haloExchange();
    }
}

void ConservativeMethod::setup_stat(double& geo_create_err) const {
    const auto& src_cell_halo  = array::make_view<int, 1>(src_mesh_.cells().halo());
    const auto& src_node_ghost = array::make_view<int, 1>(src_mesh_.nodes().ghost());
    const auto src_areas_v     = array::make_view<double, 1>(src_areas_);
    const auto tgt_areas_v     = array::make_view<double, 1>(tgt_areas_);
    double src_sum             = 0.;
    if (src_cell_data_) {
        for (idx_t spt = 0; spt < src_areas_v.size(); ++spt) {
            if (not src_cell_halo(spt)) {
                src_sum += src_areas_v(spt);
            }
        }
    }
    else {
        for (idx_t src = 0; src < src_areas_v.size(); ++src) {
            if (not src_node_ghost(src)) {
                src_sum += src_areas_v(src);
            }
        }
    }
    const auto& tgt_cell_halo  = array::make_view<int, 1>(tgt_mesh_.cells().halo());
    const auto& tgt_node_ghost = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
    double tgt_sum             = 0.;
    if (tgt_cell_data_) {
        for (idx_t tpt = 0; tpt < tgt_areas_v.size(); ++tpt) {
            if (not tgt_cell_halo(tpt)) {
                tgt_sum += tgt_areas_v(tpt);
            }
        }
    }
    else {
        for (idx_t tpt = 0; tpt < tgt_areas_v.size(); ++tpt) {
            if (not tgt_node_ghost(tpt)) {
                tgt_sum += tgt_areas_v(tpt);
            }
        }
    }
    geo_create_err = std::abs(src_sum - tgt_sum) * 0.25 * M_1_PI;
    Log::info() << " ConservativeMethod::stat : global error in polygon create   : " << geo_create_err << "\n";
}

void ConservativeMethod::remap_stat(const FieldArray& src_vals, const FieldArray& tgt_vals, FieldArray& diff_vals,
                                    double& global_cons_err, double func(const PointLonLat&), double& remap_error_l2,
                                    double& remap_error_linf) const {
    const auto& src_cell_halo  = array::make_view<int, 1>(src_mesh_.cells().halo());
    const auto& src_node_ghost = array::make_view<int, 1>(src_mesh_.nodes().ghost());
    const auto& tgt_cell_halo  = array::make_view<int, 1>(tgt_mesh_.cells().halo());
    const auto& tgt_node_ghost = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
    const auto src_areas_v     = array::make_view<double, 1>(src_areas_);
    const auto tgt_areas_v     = array::make_view<double, 1>(tgt_areas_);
    global_cons_err            = 0;
    if (src_cell_data_) {
        for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
            if (src_cell_halo(spt)) {
                continue;
            }
            double diff = src_vals(spt) * src_areas_v(spt);
            global_cons_err += diff;
            const auto& iparam = iparam_[spt];
            if (tgt_cell_data_) {
                for (idx_t icell = 0; icell < iparam.weights.size(); ++icell) {
                    idx_t tcell = iparam.tcell_id[icell];
                    if (tgt_cell_halo(tcell)) {
                        diff -= tgt_vals(iparam.tcell_id[icell]) * iparam.weights[icell];
                    }
                }
            }
            else {
                for (idx_t icell = 0; icell < iparam.weights.size(); ++icell) {
                    idx_t tcell = iparam.tcell_id[icell];
                    idx_t tnode = tgt_csp2node_[tcell];
                    diff -= tgt_vals(tnode) * iparam.weights[icell];
                }
            }
            diff_vals(spt) = std::abs(diff) / src_areas_v(spt);
        }
    }
    else {
        for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
            if (src_node_ghost(spt) or src_areas_v(spt) < 1e-14) {
                continue;
            }
            double diff = src_vals(spt) * src_areas_v(spt);
            global_cons_err += diff;
            const auto& node2csp = src_node2csp_[spt];
            for (idx_t subcell = 0; subcell < node2csp.size(); ++subcell) {
                const auto& iparam = iparam_[node2csp[subcell]];
                if (tgt_cell_data_) {
                    for (idx_t icell = 0; icell < iparam.weights.size(); ++icell) {
                        diff -= tgt_vals(iparam.tcell_id[icell]) * iparam.weights[icell];
                    }
                }
                else {
                    for (idx_t icell = 0; icell < iparam.weights.size(); ++icell) {
                        idx_t tcell = iparam.tcell_id[icell];
                        idx_t tnode = tgt_csp2node_[tcell];
                        diff -= tgt_vals(tnode) * iparam.weights[icell];
                    }
                }
            }
            diff_vals(spt) = std::abs(diff) / src_areas_v(spt);
        }
    }
    remap_error_l2   = 0.;
    remap_error_linf = 0.;
    if (tgt_cell_data_) {
        for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
            if (tgt_cell_halo(tpt)) {
                continue;
            }
            global_cons_err -= tgt_vals(tpt) * tgt_areas_v(tpt);
            auto p = tgt_points_[tpt];
            PointLonLat pll;
            eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
            double err_l = std::abs(tgt_vals(tpt) - func(pll));
            remap_error_l2 += err_l * err_l * tgt_areas_v(tpt);
            remap_error_linf = std::max(remap_error_linf, err_l);
        }
    }
    else {
        for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
            if (tgt_node_ghost(tpt)) {
                continue;
            }
            global_cons_err -= tgt_vals(tpt) * tgt_areas_v(tpt);
            auto p = tgt_points_[tpt];
            PointLonLat pll;
            eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
            double err_l = std::abs(tgt_vals(tpt) - func(pll));
            remap_error_l2 += err_l * err_l * tgt_areas_v(tpt);
            remap_error_linf = std::max(remap_error_linf, err_l);
        }
    }
    remap_error_l2  = std::sqrt(remap_error_l2 * 0.25 * M_1_PI);
    global_cons_err = std::sqrt(std::abs(global_cons_err) * 0.25 * M_1_PI);
}

template <class TargetCellsIDs>
void ConservativeMethod::dump_intersection(const CSPolygon& s_csp, const CSPolygonArray& tgt_csp,
                                           const TargetCellsIDs& tgt_cells) const {
    Log::info().flush();
    Log::info() << "\n === DEBUG ===\n\n";
    Log::info() << "* src cell: " << std::setprecision(15) << s_csp << "\n";
    Log::info() << "* src area: " << s_csp.area() << "\n\n";
    double area_ncov = s_csp.area();
    for (int i = 0; i < tgt_cells.size(); ++i) {
        const auto tcell  = tgt_cells[i].payload();
        const auto& t_csp = std::get<0>(tgt_csp[tcell]);
        Log::info() << "* src cell: " << s_csp << "\n";
        Log::info() << "* tgt cell, tgt_cell_part, tgt_centroid : " << t_csp << " " << std::get<1>(tgt_csp[tcell])
                    << " " << t_csp.centroid() << "\n";
        auto iplg          = s_csp.intersect(t_csp);
        auto jplg          = t_csp.intersect(s_csp);
        const double darea = std::abs(iplg.area() - jplg.area());
        Log::info() << "* src ^ tgt      : " << iplg << "\n";
        Log::info() << "* src ^ tgt area : " << iplg.area() << "\n";
        if (darea > 5e-8) {
            s_csp.intersect(t_csp, 1);
            Log::info() << "* (!!) intersect comm area diff: " << darea << "\n";
            Log::info() << "* (!!) tgt ^ src      : " << jplg << "\n";
            Log::info() << "* (!!) tgt ^ src area : " << jplg.area() << "\n";
            t_csp.intersect(s_csp, 1);
            ATLAS_ASSERT(false);
        }
        Log::info() << "\n";
        area_ncov -= iplg.area();
    }
    Log::info() << "\n=== END DEBUG ===\n\n";
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
