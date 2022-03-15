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

#define PLG_DEBUG 0

namespace atlas {
namespace interpolation {
namespace method {

using CSPolygon      = util::ConvexSphericalPolygon;
using CSPolygonArray = ConservativeMethod::CSPolygonArray;
using RemapStat      = ConservativeMethod::RemapStat;

namespace {
MethodBuilder<ConservativeMethod> __builder("conservative");
}

int inside_vertices( const CSPolygon& plg1, const CSPolygon& plg2, int& pout ) {
    int points_in = 0;
    pout = 0;
    for (int j = 0; j < plg2.size(); j++ ) {
        int i = 0;
        for (; i < plg1.size(); i++ ) {
            int in = (i!=plg1.size()-1) ? i+1 : 0;
            auto gss = CSPolygon::GreatCircleSegment(plg1[i], plg1[in]);
            if (not gss.inLeftHemisphere(plg2[j], -std::numeric_limits<double>::epsilon())) {
                pout++;
                break;
            };
        }
        if (i==plg1.size()) {
            points_in++;
        }
    }
    ATLAS_ASSERT( points_in + pout == plg2.size() );
    return points_in;
}

ConservativeMethod::ConservativeMethod(const Config& config): Method(config) {
    config.get("order", order_ = 1);
    config.get("normalise_intersections", normalise_intersections_ = 0);
    config.get("matrix_free", matrix_free_ = false);
    config.get("src_cell_data", src_cell_data_ = true);
    config.get("tgt_cell_data", tgt_cell_data_ = true);
    remap_stat_.setup_computed = false;
    remap_stat_.remap_computed = false;
}

int ConservativeMethod::next_index(int current_index, int size, int offset) const {
    ATLAS_ASSERT(current_index >= 0 && current_index < size);
    ATLAS_ASSERT(offset >= 0 && offset <= size);
    return (current_index < size - offset) ? current_index + offset : current_index + offset - size;
}
int ConservativeMethod::prev_index(int current_index, int size, int offset) const {
    ATLAS_ASSERT(current_index >= 0 && current_index < size);
    ATLAS_ASSERT(offset >= 0 && offset <= size);
    return (current_index >= offset) ? current_index - offset : current_index - offset + size;
}

RemapStat& ConservativeMethod::remap_stat() const {
    if (not remap_stat_.setup_computed or not remap_stat_.remap_computed) {
        std::cerr << "WARNING RemapStat not computed before accessed.\n";
    }
    return remap_stat_;
}

void ConservativeMethod::set_order(int order) {
    if (matrix_free_ && (not src_cell_data_ or not tgt_cell_data_)) {
        ATLAS_NOTIMPLEMENTED;
    }
    if (order == 1) {
        if (not matrix_free_) {
            setup_1st_order_matrix();
        }
    }
    else if (order == 2) {
        if (not matrix_free_) {
            setup_2nd_order_matrix();
        }
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
    remap_stat_.remap_computed = false;
}

// get counter-clockwise sorted neighbours of a cell
std::vector<idx_t> ConservativeMethod::get_cell_neighbours(Mesh& mesh, idx_t cell) const {
    const auto& cell2node  = mesh.cells().node_connectivity();
    const auto& nodes_ll   = array::make_view<double, 2>(mesh.nodes().lonlat());
    const idx_t n_nodes    = cell2node.cols(cell);
    std::vector<idx_t> nbr_cells;
    nbr_cells.reserve(n_nodes);
    if (mesh.nodes().cell_connectivity().rows() == 0) {
        mesh::actions::build_node_to_cell_connectivity(mesh);
    }
    const auto& node2cell  = mesh.nodes().cell_connectivity();
    const auto n2c_missval = node2cell.missing_value();

    for (idx_t inode = 0; inode < n_nodes; ++inode) {
		idx_t node0             = cell2node(cell, inode);
        idx_t node1             = cell2node(cell, next_index(inode, n_nodes));
		const PointLonLat p0_ll = PointLonLat{nodes_ll(node0, 0), nodes_ll(node0, 1)};
		const PointLonLat p1_ll = PointLonLat{nodes_ll(node1, 0), nodes_ll(node1, 1)};
		PointXYZ p0, p1;
		eckit::geometry::Sphere::convertSphericalToCartesian(1., p0_ll, p0);
		eckit::geometry::Sphere::convertSphericalToCartesian(1., p1_ll, p1);
		if (PointXYZ::norm(p0 - p1) < 1e-14) {
			continue;  // edge = point
		}
        bool still_search = true; // still search the cell having vertices node0 & node1, not having index "cell"
		int n_cells0 = node2cell.cols( node0 );
		int n_cells1 = node2cell.cols( node1 );
		for( int icell0 = 0; still_search && icell0 < n_cells0; icell0++ ) {
			int cell0 = node2cell( node0, icell0 );
			if ( cell0 == cell ) {
				continue;
			}
			for( int icell1 = 0; still_search && icell1 < n_cells1; icell1++ ) {
				int cell1 = node2cell( node1, icell1 );
                if ( cell0 == cell1 && cell0 != n2c_missval ) {
					nbr_cells.emplace_back( cell0 );
					still_search = false;
				}
			}
		}
    }
    return nbr_cells;
}

// get cyclically sorted node neighbours without using edge connectivity
std::vector<idx_t> ConservativeMethod::get_node_neighbours(Mesh& mesh, idx_t node_id) const {
    const auto& cell2node  = mesh.cells().node_connectivity();
    if (mesh.nodes().cell_connectivity().rows() == 0) {
        mesh::actions::build_node_to_cell_connectivity(mesh);
    }
    const auto& node2cell  = mesh.nodes().cell_connectivity();
    std::vector<idx_t> nbr_nodes;
    std::vector<idx_t> nbr_nodes_od;
    const int ncells       = node2cell.cols( node_id );
    ATLAS_ASSERT( ncells > 0, "There is a node which does not connect to any cell" );
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
        cnodes[icell][0] = cell2node( cell, prev_index(cnode, nnodes) );
        cnodes[icell][1] = cell2node( cell, next_index(cnode, nnodes) );
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
    for (idx_t cell = 0; cell < n_cells; ++cell) {
        int halo_type         = cell_halo(cell);
        const idx_t n_nodes   = cell2node.cols(cell);
        pts_ll.clear();
        pts_ll.resize(n_nodes);
        for (idx_t jnode = 0; jnode < n_nodes; ++jnode) {
            idx_t inode   = cell2node(cell, jnode);
            pts_ll[jnode] = PointLonLat{lonlat(inode, 0), lonlat(inode, 1)};
        }
        const auto& bitflag = util::Bitflags::view(cell_flags(cell));
        if (bitflag.check(util::Topology::PERIODIC) and mpi::rank()==cell_part(cell)) {
              halo_type = -1;
        }
        std::get<0>(cspolygons[cell]) = CSPolygon(pts_ll);
        std::get<1>(cspolygons[cell]) = halo_type;
    }
    return cspolygons;
}

// Create polygons for cell-vertex data. Here, the polygons are subcells of mesh cells created as
// 	 (cell_centre, edge_centre, cell_vertex, edge_centre)
// additionally, subcell-to-node and node-to-subcells mapping are computed
CSPolygonArray ConservativeMethod::get_polygons_nodedata(Mesh& mesh, std::vector<idx_t>& csp2node,
                                                         std::vector<std::vector<idx_t>>& node2csp,
                                                         std::array<double,2>& errors) const {
    CSPolygonArray cspolygons;
    csp2node.clear();
    node2csp.clear();
    node2csp.resize(mesh.nodes().size());
    const auto nodes_ll   = array::make_view<double, 2>(mesh.nodes().lonlat());
    const auto& cell2node = mesh.cells().node_connectivity();
    const auto cell_halo  = array::make_view<int, 1>(mesh.cells().halo());
    const auto cell_flags = array::make_view<int, 1>(mesh.cells().flags());
    const auto cell_part  = array::make_view<int, 1>(mesh.cells().partition());
    auto xyz2ll = [](const atlas::PointXYZ& p_xyz) {
        PointLonLat p_ll;
        eckit::geometry::Sphere::convertCartesianToSpherical(1., p_xyz, p_ll);
        return p_ll;
    };
    auto ll2xyz = [](const atlas::PointLonLat& p_ll) {
        PointXYZ p_xyz;
        eckit::geometry::Sphere::convertSphericalToCartesian(1., p_ll, p_xyz);
        return p_xyz;
    };
    idx_t cspol_id = 0; // subpolygon enumeration
    errors = {0., 0.}; // over/undershoots in creation of subpolygons
    for (idx_t cell = 0; cell < mesh.cells().size(); ++cell) {
        ATLAS_ASSERT(cell < cell2node.rows());
        const idx_t n_nodes = cell2node.cols(cell);
        ATLAS_ASSERT(n_nodes > 2);
        PointXYZ cell_mid(0., 0., 0.);	// cell centre
        std::vector<PointXYZ> pts_xyz;
        std::vector<PointLonLat> pts_ll;
        std::vector<int> pts_idx;
        pts_xyz.reserve( n_nodes );
        pts_ll.reserve( n_nodes );
        pts_idx.reserve( n_nodes );
        for (idx_t inode = 0; inode < n_nodes; ++inode) {
            idx_t node0             = cell2node(cell, inode);
            idx_t node1             = cell2node(cell, next_index(inode, n_nodes));
            const PointLonLat p0_ll = PointLonLat{nodes_ll(node0, 0), nodes_ll(node0, 1)};
            const PointLonLat p1_ll = PointLonLat{nodes_ll(node1, 0), nodes_ll(node1, 1)};
            PointXYZ p0 = ll2xyz( p0_ll );
            PointXYZ p1 = ll2xyz( p1_ll );
            if (PointXYZ::norm(p0 - p1) < 1e-14) {
                continue;  // skip this edge = a pole point
            }
            pts_xyz.emplace_back( p0 );
            pts_ll.emplace_back( p0_ll );
            pts_idx.emplace_back( inode );
            cell_mid = cell_mid + p0;
            cell_mid = cell_mid + p1;
        }
        cell_mid = PointXYZ::div(cell_mid, PointXYZ::norm(cell_mid));
        PointLonLat cell_ll = xyz2ll( cell_mid );
        double loc_csp_area_shoot = CSPolygon( pts_ll ).area();
        // get CSPolygon for each valid edge
        for ( int inode = 0; inode < pts_idx.size(); inode++ ) {
            int inode_n              = next_index(inode, pts_idx.size());
            idx_t node               = cell2node(cell, inode);
            idx_t node_n             = cell2node(cell, inode_n);
            PointXYZ iedge_mid = pts_xyz[inode] + pts_xyz[inode_n];
            iedge_mid          = PointXYZ::div(iedge_mid, PointXYZ::norm(iedge_mid));
			csp2node.emplace_back( node_n );
			node2csp[node_n].emplace_back( cspol_id );
            int inode_nn             = next_index(inode_n, pts_idx.size());
			if ( PointXYZ::norm( pts_xyz[inode_nn] - pts_xyz[inode_n] ) < 1e-14 ) {
				ATLAS_THROW_EXCEPTION("Three cell vertices on a same great arc!");
			}
            PointXYZ jedge_mid;
			jedge_mid = pts_xyz[inode_nn] + pts_xyz[inode_n];
            jedge_mid = PointXYZ::div(jedge_mid, PointXYZ::norm(jedge_mid));	
            std::vector<PointLonLat> subpol_pts_ll(4);
            subpol_pts_ll[0]     = cell_ll;
            subpol_pts_ll[1]     = xyz2ll(iedge_mid);
            subpol_pts_ll[2]     = pts_ll[inode_n];
            subpol_pts_ll[3]     = xyz2ll(jedge_mid);
            int halo_type = cell_halo(cell);
            if (util::Bitflags::view(cell_flags(cell)).check(util::Topology::PERIODIC)
                and cell_part(cell)==mpi::rank()) {
                halo_type = -1;
            }
            CSPolygon cspi( subpol_pts_ll );
            loc_csp_area_shoot -= cspi.area();
            cspolygons.emplace_back(cspi, halo_type);
            cspol_id++;
        }
        errors[0] += std::abs(loc_csp_area_shoot);
        errors[1] = std::max( std::abs(loc_csp_area_shoot), errors[1] );
    }
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(&errors[0], 1, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&errors[1], 1, eckit::mpi::max());
    }
    return cspolygons;
}

void ConservativeMethod::do_setup(const Grid& src_grid, const Grid& tgt_grid) {
    ATLAS_TRACE("ConservativeMethod::do_setup( Grid, Grid )");
    ATLAS_ASSERT(src_grid);
    ATLAS_ASSERT(tgt_grid);
    auto src_mesh_config = src_grid.meshgenerator() | option::halo(2);
    auto tgt_mesh_config = tgt_grid.meshgenerator() | option::halo(0);
    tgt_mesh_            = MeshGenerator(tgt_mesh_config).generate(tgt_grid);
    if (mpi::size() > 1) {
        src_mesh_ = MeshGenerator(src_mesh_config).generate(src_grid, grid::MatchingPartitioner(tgt_mesh_));
    }
    else {
        src_mesh_ = MeshGenerator(src_mesh_config).generate(src_grid);
    }
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
        ATLAS_THROW_EXCEPTION("ConservativeMethod: source function space invalid");
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
        ATLAS_THROW_EXCEPTION("ConservativeMethod: target function space invalid");
    }
    {
        // we need src_halo_size >= 2, whereas tgt_halo_size >= 0 is enough
        int src_halo_size = 0;
        src_mesh_.metadata().get("halo", src_halo_size);
        ATLAS_ASSERT(src_halo_size > 1);
    }
    CSPolygonArray src_csp;
    CSPolygonArray tgt_csp;
    std::array<double,2> errors = {0., 0.};
    {
        ATLAS_TRACE("Get source polygons");
        if (src_cell_data_) {
            src_csp = get_polygons_celldata(src_mesh_);
        }
        else {
            src_csp = get_polygons_nodedata(src_mesh_, src_csp2node_, src_node2csp_, errors);
        }
    }
    remap_stat_.errors[RemapStat::Errors::SRC_PLG_L1] = errors[0];
    remap_stat_.errors[RemapStat::Errors::SRC_PLG_LINF] = errors[1];
    {
        ATLAS_TRACE("Get target polygons");
        if (tgt_cell_data_) {
            tgt_csp = get_polygons_celldata(tgt_mesh_);
        }
        else {
            tgt_csp = get_polygons_nodedata(tgt_mesh_, tgt_csp2node_, tgt_node2csp_, errors);
        }
    }
    remap_stat_.counts[RemapStat::Counts::SRC_PLG] = src_csp.size();
    remap_stat_.counts[RemapStat::Counts::TGT_PLG] = tgt_csp.size();
    remap_stat_.errors[RemapStat::Errors::TGT_PLG_L1] = errors[0];
    remap_stat_.errors[RemapStat::Errors::TGT_PLG_LINF] = errors[1];

    intersect_polygons(src_csp, tgt_csp);

    n_spoints_ = src_fs_.size();
    n_tpoints_ = tgt_fs_.size();
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
            if (src_node2csp_[spt].size() == 0) {
                // this is a node to which no subpolygon is associated
                // maximal twice per mesh we end here, and that is only when mesh has nodes on poles
                auto p = PointLonLat{lonlat(spt, 0), lonlat(spt, 1)};
                eckit::geometry::Sphere::convertSphericalToCartesian(1., p, src_points_[spt]);
            }
            else {
                // .. in the other case, start computing the barycentre
                src_points_[spt] = PointXYZ{0., 0., 0.};
            }
            src_areas_v(spt) = 0.;
            for (idx_t isubcell = 0; isubcell < src_node2csp_[spt].size(); ++isubcell) {
                idx_t subcell     = src_node2csp_[spt][isubcell];
                const auto& s_csp = std::get<0>(src_csp[subcell]);
                src_areas_v(spt) += s_csp.area();
                src_points_[spt] = src_points_[spt] + PointXYZ::mul(s_csp.centroid(), s_csp.area());
            }
            double src_point_norm = PointXYZ::norm(src_points_[spt]);
            ATLAS_ASSERT(src_point_norm > 0.);
            src_points_[spt]      = PointXYZ::div(src_points_[spt], src_point_norm);
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
            if (tgt_node2csp_[tpt].size() == 0) {
                // this is a node to which no subpolygon is associated
                // maximal twice per mesh we end here, and that is only when mesh has nodes on poles
                auto p = PointLonLat{lonlat(tpt, 0), lonlat(tpt, 1)};
                eckit::geometry::Sphere::convertSphericalToCartesian(1., p, tgt_points_[tpt]);
            }
            else {
                // .. in the other case, start computing the barycentre
                tgt_points_[tpt] = PointXYZ{0., 0., 0.};
            }
            tgt_areas_v(tpt) = 0.;
            for (idx_t isubcell = 0; isubcell < tgt_node2csp_[tpt].size(); ++isubcell) {
                idx_t subcell = tgt_node2csp_[tpt][isubcell];
                const auto& t_csp = std::get<0>(tgt_csp[subcell]);
                tgt_areas_v(tpt) += t_csp.area();
                tgt_points_[tpt] = tgt_points_[tpt] + PointXYZ::mul(t_csp.centroid(), t_csp.area());
            }
            double tgt_point_norm = PointXYZ::norm(tgt_points_[tpt]);
            ATLAS_ASSERT(tgt_point_norm > 0.);
            tgt_points_[tpt]      = PointXYZ::div(tgt_points_[tpt], tgt_point_norm);
        }
    }
    if (not matrix_free_) {
        setup_1st_order_matrix();
        setup_2nd_order_matrix();
    }
}

// needed for intersect_polygons only
struct ComparePointXYZ {
    bool operator()(const PointXYZ& f, const PointXYZ& s) const {
        // eps = ConvexSphericalPolygon::EPS which is the threshold when two points are "same"
        double eps = 1e4 * std::numeric_limits<double>::epsilon();
        if (f[0] < s[0] - eps) {
            return true;
        }
        else if ( std::abs(f[0]-s[0]) < eps ) {
            if ( f[1] < s[1] - eps ) {
                return true;
            }
            else if ( std::abs(f[1]-s[1]) < eps ) {
                if ( f[2] < s[2] - eps ) {
                    return true;
                }
            }
        }
        return false;
    }
};

void ConservativeMethod::intersect_polygons(const CSPolygonArray& src_csp, const CSPolygonArray& tgt_csp) {
    ATLAS_TRACE();
    util::KDTree<idx_t> kdt_search;
    kdt_search.reserve(tgt_csp.size());
    double max_tgtcell_rad = 0.;
    for (idx_t jcell = 0; jcell < tgt_csp.size(); ++jcell) {
        if (std::get<1>(tgt_csp[jcell]) == 0) {
            const auto& t_csp = std::get<0>(tgt_csp[jcell]);
            kdt_search.insert(t_csp.centroid(), jcell);
            max_tgtcell_rad = std::max(max_tgtcell_rad, t_csp.radius());
        }
    }
    kdt_search.build();

    std::set<PointXYZ, ComparePointXYZ> src_cent;
    auto polygon_point = [](const CSPolygon& pol) {
        PointXYZ p{0.,0.,0.};
        for ( int i =0; i < pol.size() ; i++ ) {
            p = p + pol[i];
        }
        p /= pol.size();
        return p;
    };
    auto src_already_in = [&](const PointXYZ& halo_cent) {
          if( src_cent.find(halo_cent) == src_cent.end() ) {
            src_cent.insert(halo_cent);
            return false;
          }
        return true;
    };

    enum MeshSizeId {SRC, TGT, SRC_TGT_INTERSECT, SRC_NONINTERSECT};
    std::array<size_t,4> num_pol{0,0,0,0};
    enum AreaCoverageId {TOTAL_SRC, MAX_SRC};
    std::array<double,2> area_coverage{0.,0.};
    src_iparam_.resize(src_csp.size());
#if PLG_DEBUG
    std::vector<InterpolationParameters> tgt_iparam;
    tgt_iparam.resize(tgt_csp.size());
#endif
    eckit::Channel blackhole;
    eckit::ProgressTimer progress("Intersecting polygons ", src_csp.size(), " cell", double(10),
                                  src_csp.size() > 50 ? Log::info() : blackhole);
    for (idx_t scell = 0; scell < src_csp.size(); ++scell, ++progress) {
        if (src_already_in( polygon_point(std::get<0>(src_csp[scell])) )) {
            continue;
        }
        const auto& s_csp   = std::get<0>(src_csp[scell]);
		const double s_csp_area = s_csp.area();
        double src_cover_area = 0.;
        auto tgt_cells = kdt_search.closestPointsWithinRadius(s_csp.centroid(), s_csp.radius() + max_tgtcell_rad);
        for (idx_t ttcell = 0; ttcell < tgt_cells.size(); ++ttcell) {
            auto tcell        = tgt_cells[ttcell].payload();
            const auto& t_csp = std::get<0>(tgt_csp[tcell]);
            CSPolygon csp_i   = s_csp.intersect(t_csp);
            double csp_i_area = csp_i.area();
#if PLG_DEBUG
            // check zero area intersections with inside_vertices
            int pout;
            if ( inside_vertices(s_csp, t_csp, pout) > 2 && csp_i.area() < 3e-16 ) {
                dump_intersection( s_csp, tgt_csp, tgt_cells );
            }
#endif
            if (csp_i_area > 0.) {
#if PLG_DEBUG
                tgt_iparam[tcell].cell_idx.emplace_back(scell);
                tgt_iparam[tcell].tgt_weights.emplace_back(csp_i_area);
#endif
                src_iparam_[scell].cell_idx.emplace_back(tcell);
                src_iparam_[scell].src_weights.emplace_back(csp_i_area);
                double target_weight = csp_i_area / t_csp.area();
                src_iparam_[scell].tgt_weights.emplace_back( target_weight );
                src_iparam_[scell].centroids.emplace_back(csp_i.centroid());
                src_cover_area += csp_i_area;
                ATLAS_ASSERT( target_weight < 1.1 );
                ATLAS_ASSERT( csp_i_area / s_csp_area < 1.1 );
            }
        }
        const double src_cover_err = std::abs(s_csp_area - src_cover_area);
        const double src_cover_err_percent = 100. * src_cover_err / s_csp_area;
        if (src_cover_err_percent > 10. and std::get<1>(src_csp[scell]) == 0) {
            // HACK: source cell at process boundary will not be covered by target cells, skip them
            // TODO: mark these source cells beforehand and compute error in them among the processes
#if PLG_DEBUG
            if ( mpi::size() == 1) {
                Log::info() << "WARNING src cell covering error : " << src_cover_err_percent << "%\n";
                dump_intersection( s_csp, tgt_csp, tgt_cells );
            }
#endif
            area_coverage[TOTAL_SRC] += src_cover_err;
            area_coverage[MAX_SRC]   = std::max( area_coverage[MAX_SRC], src_cover_err );
        }
        if (src_iparam_[scell].cell_idx.size() == 0) {
            num_pol[SRC_NONINTERSECT]++;
        }
        if (normalise_intersections_ && src_cover_err_percent < 1.) {
            double wfactor = s_csp.area() / (src_cover_area > 0. ? src_cover_area : 1.);
            for (idx_t i = 0; i < src_iparam_[scell].src_weights.size(); i++) {
                src_iparam_[scell].src_weights[i] *= wfactor;
                src_iparam_[scell].tgt_weights[i] *= wfactor;
            }
        }
        num_pol[SRC_TGT_INTERSECT] += src_iparam_[scell].src_weights.size();
    }
    num_pol[SRC] = src_csp.size();
    num_pol[TGT] = tgt_csp.size();
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(&num_pol[0], 4, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&area_coverage[0], 2, eckit::mpi::max());
    }
    remap_stat_.counts[RemapStat::Counts::INT_PLG] = num_pol[SRC_TGT_INTERSECT];
    remap_stat_.counts[RemapStat::Counts::UNCVR_SRC] = num_pol[SRC_NONINTERSECT];
    remap_stat_.errors[RemapStat::Errors::GEO_L1] = area_coverage[TOTAL_SRC];
    remap_stat_.errors[RemapStat::Errors::GEO_LINF] = area_coverage[MAX_SRC];

    double geo_err_l1    =  0.;
    double geo_err_linf  =  0.;
    for (idx_t scell = 0; scell < src_csp.size(); ++scell) {
        const int cell_flag = std::get<1>(src_csp[scell]);
        if (cell_flag == -1 or cell_flag > 0) {
            // skip periodic & halo cells
            continue;
        }
        double diff_cell = std::get<0>(src_csp[scell]).area();
        for (idx_t icell = 0; icell < src_iparam_[scell].src_weights.size(); ++icell) {
            diff_cell -= src_iparam_[scell].src_weights[icell];
        }
        geo_err_l1    += std::abs(diff_cell);
        geo_err_linf  = std::max(geo_err_linf, std::abs(diff_cell));
    }
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(&geo_err_l1, 1, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&geo_err_linf, 1, eckit::mpi::max());
    }
    remap_stat_.errors[RemapStat::Errors::GEO_L1] = 0.25 * M_1_PI * geo_err_l1;
    remap_stat_.errors[RemapStat::Errors::GEO_LINF] = geo_err_linf;
    remap_stat_.setup_computed = true;
#if PLG_DEBUG
    for (idx_t tcell = 0; tcell < tgt_csp.size(); ++tcell) {
        const auto& t_csp = std::get<0>( tgt_csp[tcell] );
        double tgt_cover_area = 0.;
        const auto& tiparam = tgt_iparam[tcell];
        for (idx_t icell = 0; icell < tiparam.cell_idx.size(); ++icell) {
            tgt_cover_area += tiparam.tgt_weights[icell];
        }
        const double tgt_cover_err_percent = 100. * std::abs(t_csp.area() - tgt_cover_area) / t_csp.area();
        if (tgt_cover_err_percent > 1e-4 and std::get<1>(tgt_csp[tcell]) == 0) {
            Log::info() << "WARNING tgt cell covering error : " << tgt_cover_err_percent << " %\n";
            dump_intersection( t_csp, src_csp, tiparam.cell_idx );
        }
    }
#endif
}

void ConservativeMethod::setup_1st_order_matrix() {
    ATLAS_TRACE("ConservativeMethod::setup: build cons-1 interpolant matrix");
	order_ = 1;
	ATLAS_ASSERT(not matrix_free_);
    Triplets triplets;
    size_t triplets_size = 0;
    // determine the size of array of triplets used to define the sparse matrix
    if (src_cell_data_) {
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            triplets_size += src_iparam_[scell].centroids.size();
        }
    }
    else {
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                idx_t subcell = src_node2csp_[snode][isubcell];
                triplets_size += src_iparam_[subcell].tgt_weights.size();
            }
        }
    }
    triplets.reserve(triplets_size);
    // assemble triplets to define the sparse matrix
    const auto src_areas_v = array::make_view<double, 1>(src_areas_);
    const auto tgt_areas_v = array::make_view<double, 1>(tgt_areas_);
    if (src_cell_data_ && tgt_cell_data_) {
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            const auto& iparam = src_iparam_[scell];
            for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                idx_t tcell = iparam.cell_idx[icell];
                triplets.emplace_back(tcell, scell, iparam.tgt_weights[icell]);
            }
        }
    }
    else if (not src_cell_data_ && tgt_cell_data_) {
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                const idx_t subcell = src_node2csp_[snode][isubcell];
                const auto& iparam  = src_iparam_[subcell];
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    idx_t tcell = iparam.cell_idx[icell];
                    triplets.emplace_back(tcell, snode, iparam.tgt_weights[icell]);
                }
            }
        }
    }
    else if (src_cell_data_ && not tgt_cell_data_) {
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            const auto& iparam = src_iparam_[scell];
            for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                idx_t tcell = iparam.cell_idx[icell];
                idx_t tnode = tgt_csp2node_[tcell];
                double inv_node_weight = (tgt_areas_v(tnode) > 0. ? 1./tgt_areas_v(tnode) : 0.);
                triplets.emplace_back(tnode, scell, iparam.src_weights[icell] * inv_node_weight);
            }
        }
    }
    else if (not src_cell_data_ && not tgt_cell_data_) {
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                const idx_t subcell = src_node2csp_[snode][isubcell];
                const auto& iparam  = src_iparam_[subcell];
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    idx_t tcell = iparam.cell_idx[icell];
                    idx_t tnode = tgt_csp2node_[tcell];
                    double inv_node_weight = (tgt_areas_v(tnode) > 0. ? 1./tgt_areas_v(tnode) : 0.);
                    triplets.emplace_back(tnode, snode, iparam.src_weights[icell] * inv_node_weight);
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
    ATLAS_TRACE("ConservativeMethod::setup: build cons-2 interpolant matrix");
	order_ = 2;
	ATLAS_ASSERT(not matrix_free_);
    Triplets triplets;
    size_t triplets_size   = 0;
    const auto tgt_areas_v = array::make_view<double, 1>(tgt_areas_);
    if (src_cell_data_) {
        const auto src_halo = array::make_view<int, 1>(src_mesh_.cells().halo());
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            const auto nb_cells = get_cell_neighbours(src_mesh_, scell);
            triplets_size += (2 * nb_cells.size() + 1) * src_iparam_[scell].centroids.size();
        }
        triplets.reserve(triplets_size);
        for (idx_t scell = 0; scell < n_spoints_; ++scell) {
            const auto nb_cells = get_cell_neighbours(src_mesh_, scell);
            const auto& iparam  = src_iparam_[scell];
            if (iparam.centroids.size() == 0 && not src_halo(scell)) {
                continue;
            }
            /* // better conservation after Kritsikis et al. (2017)
            PointXYZ Cs = {0., 0., 0.};
            for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
                Cs = Cs + PointXYZ::mul( iparam.centroids[icell], iparam.src_weights[icell] );
            }
            const double Cs_norm = PointXYZ::norm( Cs );
            ATLAS_ASSERT( Cs_norm > 0. );
            Cs = PointXYZ::div( Cs, Cs_norm );
			*/
            const PointXYZ& Cs = src_points_[scell];
            // compute gradient from cells
            double dual_area_inv = 0.;
            std::vector<PointXYZ> Rsj;
            Rsj.resize(nb_cells.size());
            for (idx_t j = 0; j < nb_cells.size(); ++j) {
                idx_t nj         = next_index(j, nb_cells.size());
                idx_t sj         = nb_cells[j];
                idx_t nsj        = nb_cells[nj];
                const auto& Csj  = src_points_[sj];
                const auto& Cnsj = src_points_[nsj];
                if (CSPolygon::GreatCircleSegment( Cs, Csj ).inLeftHemisphere( Cnsj,
-1e-16 )) {;
                    Rsj[j] = PointXYZ::cross(Cnsj, Csj);
                    dual_area_inv += CSPolygon({Cs, Csj, Cnsj}).area();
                }
                else {
                    Rsj[j] = PointXYZ::cross(Csj, Cnsj);
                    dual_area_inv += CSPolygon({Cs, Cnsj, Csj}).area();
                }
            }
            dual_area_inv = (dual_area_inv > 0.) ? 1. / dual_area_inv : 0.;
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
                Aik[icell]            = PointXYZ::mul(Aik[icell], iparam.tgt_weights[icell] * dual_area_inv);
            }
            if (tgt_cell_data_) {
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    const idx_t tcell = iparam.cell_idx[icell];
                    for (idx_t j = 0; j < nb_cells.size(); ++j) {
                        idx_t nj  = next_index(j, nb_cells.size());
                        idx_t sj  = nb_cells[j];
                        idx_t nsj = nb_cells[nj];
                        triplets.emplace_back(tcell, sj, 0.5 * PointXYZ::dot(Rsj[j], Aik[icell]));
                        triplets.emplace_back(tcell, nsj, 0.5 * PointXYZ::dot(Rsj[j], Aik[icell]));
                    }
                    triplets.emplace_back(tcell, scell, iparam.tgt_weights[icell] - PointXYZ::dot(Rs, Aik[icell]));
                }
            }
            else {
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    idx_t tcell                = iparam.cell_idx[icell];
                    idx_t tnode                = tgt_csp2node_[tcell];
                    double inv_node_weight = (tgt_areas_v(tnode) > 0.) ? 1./tgt_areas_v(tnode) : 0.;
                    double csp2node_coef = iparam.src_weights[icell] / iparam.tgt_weights[icell] * inv_node_weight;
                    for (idx_t j = 0; j < nb_cells.size(); ++j) {
                        idx_t nj  = next_index(j, nb_cells.size());
                        idx_t sj  = nb_cells[j];
                        idx_t nsj = nb_cells[nj];
                        triplets.emplace_back(tnode, sj, (0.5 * PointXYZ::dot(Rsj[j], Aik[icell])) * csp2node_coef);
                        triplets.emplace_back(tnode, nsj, (0.5 * PointXYZ::dot(Rsj[j], Aik[icell])) * csp2node_coef);
                    }
                    triplets.emplace_back(tnode, scell,
                                          (iparam.tgt_weights[icell] - PointXYZ::dot(Rs, Aik[icell])) * csp2node_coef);
                }
            }
        }
    }
    else {  // if ( not src_cell_data_ )
        const auto src_halo = array::make_view<int, 1>(src_mesh_.nodes().halo());
        for (idx_t snode = 0; snode < n_spoints_; ++snode) {
            const auto nb_nodes = get_node_neighbours(src_mesh_, snode);
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                idx_t subcell = src_node2csp_[snode][isubcell];
                triplets_size += (2 * nb_nodes.size() + 1) * src_iparam_[subcell].centroids.size();
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
                const auto& iparam = src_iparam_[subcell];
                for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
                    Cs = Cs + PointXYZ::mul( iparam.centroids[icell], iparam.src_weights[icell] );
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
                idx_t nj         = next_index(j, nb_nodes.size());
                idx_t sj         = nb_nodes[j];
                idx_t snj        = nb_nodes[nj];
                const auto& Nsj  = src_points_[sj];
                const auto& Nsnj = src_points_[snj];
                if (CSPolygon::GreatCircleSegment( Ns, Nsj ).inLeftHemisphere(
Nsnj,-1e-16)) {
                    Rsj[j] = PointXYZ::cross(Nsnj, Nsj);
                    dual_area_inv += CSPolygon({Ns, Nsj, Nsnj}).area();
                }
                else {
                    Rsj[j] = PointXYZ::cross(Nsj, Nsnj);
                    dual_area_inv += CSPolygon({Ns, Nsnj, Nsj}).area();
                }
            }
            dual_area_inv = (dual_area_inv > 0.) ? 1. / dual_area_inv : 0.;
            PointXYZ Rs   = {0., 0., 0.};
            for (idx_t j = 0; j < nb_nodes.size(); ++j) {
                Rs = Rs + Rsj[j];
            }
            // assemble the matrix
            for (idx_t isubcell = 0; isubcell < src_node2csp_[snode].size(); ++isubcell) {
                idx_t subcell      = src_node2csp_[snode][isubcell];
                const auto& iparam = src_iparam_[subcell];
                if (iparam.centroids.size() == 0) {
                    continue;
                }
                std::vector<PointXYZ> Aik;
                Aik.resize(iparam.centroids.size());
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    const PointXYZ& Csk   = iparam.centroids[icell];
                    const PointXYZ Csk_Cs = Csk - Cs;
                    Aik[icell]            = Csk_Cs - PointXYZ::mul(Cs, PointXYZ::dot(Cs, Csk_Cs));
                    Aik[icell]            = PointXYZ::mul(Aik[icell], iparam.tgt_weights[icell] * dual_area_inv);
                }
                if (tgt_cell_data_) {
                    for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                        const idx_t tcell = iparam.cell_idx[icell];
                        for (idx_t j = 0; j < nb_nodes.size(); ++j) {
                            idx_t nj  = next_index(j, nb_nodes.size());
                            idx_t sj  = nb_nodes[j];
                            idx_t snj = nb_nodes[nj];
                            triplets.emplace_back(tcell, sj, 0.5 * PointXYZ::dot(Rsj[j], Aik[icell]));
                            triplets.emplace_back(tcell, snj, 0.5 * PointXYZ::dot(Rsj[j], Aik[icell]));
                        }
                        triplets.emplace_back(tcell, snode, iparam.tgt_weights[icell] - PointXYZ::dot(Rs, Aik[icell]));
                    }
                }
                else {
                    for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                        idx_t tcell = iparam.cell_idx[icell];
                        idx_t tnode = tgt_csp2node_[tcell];
                        double inv_node_weight = (tgt_areas_v(tnode) > 1e-15) ? 1./tgt_areas_v(tnode) : 0.;
                        double csp2node_coef = iparam.src_weights[icell] / iparam.tgt_weights[icell] * inv_node_weight;
                        for (idx_t j = 0; j < nb_nodes.size(); ++j) {
                            idx_t nj  = next_index(j, nb_nodes.size());
                            idx_t sj  = nb_nodes[j];
                            idx_t snj = nb_nodes[nj];
                            triplets.emplace_back(tnode, sj, (0.5 * PointXYZ::dot(Rsj[j], Aik[icell])) * csp2node_coef);
                            triplets.emplace_back(tnode, snj,
                                                  (0.5 * PointXYZ::dot(Rsj[j], Aik[icell])) * csp2node_coef);
                        }
                        triplets.emplace_back(tnode, snode,
                                              (iparam.tgt_weights[icell] - PointXYZ::dot(Rs, Aik[icell])) * csp2node_coef);
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
        if (matrix_free_) {
            ATLAS_TRACE("matrix_free_order_1");
            if (not src_cell_data_ or not tgt_cell_data_) {
                ATLAS_NOTIMPLEMENTED;
            }
            const auto src_vals = array::make_view<double, 1>(src_field);
            auto tgt_vals       = array::make_view<double, 1>(tgt_field);
            for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                tgt_vals(tcell) = 0.;
            }
            for (idx_t scell = 0; scell < src_vals.size(); ++scell) {
                const auto& iparam = src_iparam_[scell];
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    tgt_vals(iparam.cell_idx[icell]) += iparam.src_weights[icell] * src_vals(scell);
                }
            }
            for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                tgt_vals(tcell) /= tgt_areas_v(tcell);
            }
        }
        else {
            ATLAS_TRACE("matrix_order_1");
            Method::do_execute(src_field, tgt_field);
        }
    }
    else if (order_ == 2) {
        if (matrix_free_) {
            ATLAS_TRACE("matrix_free_order_2");
            if (not src_cell_data_ or not tgt_cell_data_) {
                ATLAS_NOTIMPLEMENTED;
            }
            const auto src_vals       = array::make_view<double, 1>(src_field);
            auto tgt_vals             = array::make_view<double, 1>(tgt_field);
            const auto halo           = array::make_view<int, 1>(src_mesh_.cells().halo());
            for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                tgt_vals(tcell) = 0.;
            }
            for (idx_t scell = 0; scell < src_vals.size(); ++scell) {
                if (halo(scell)) {
                    continue;
                }
                const auto& iparam       = src_iparam_[scell];
                const PointXYZ& P        = src_points_[scell];
                PointXYZ grad            = {0., 0., 0.};
                PointXYZ src_barycenter  = {0., 0., 0.};
                auto src_neighbour_cells = get_cell_neighbours(src_mesh_, scell);
                double dual_area         = 0.;
                for (idx_t nb_id = 0; nb_id < src_neighbour_cells.size(); ++nb_id) {
                    idx_t nnb_id    = next_index(nb_id, src_neighbour_cells.size());
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
                        auto NsNsj = CSPolygon::GreatCircleSegment( P, Pn );
                        val *= (NsNsj.inLeftHemisphere(Pnn, -1e-16) ? -1 : 1);
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
                    src_barycenter = src_barycenter + PointXYZ::mul(iparam.centroids[icell], iparam.src_weights[icell]);
                }
                src_barycenter = PointXYZ::div(src_barycenter, PointXYZ::norm(src_barycenter));
                grad           = grad - PointXYZ::mul(src_barycenter, PointXYZ::dot(grad, src_barycenter));
                ATLAS_ASSERT(std::abs(PointXYZ::dot(grad, src_barycenter)) < 1e-14);
                for (idx_t icell = 0; icell < iparam.centroids.size(); ++icell) {
                    tgt_vals(iparam.cell_idx[icell]) +=
                        iparam.src_weights[icell] *
                        (src_vals(scell) + PointXYZ::dot(grad, iparam.centroids[icell] - src_barycenter));
                }
            }
            for (idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell) {
                tgt_vals(tcell) /= tgt_areas_v(tcell);
            }
        }
        else {
            ATLAS_TRACE("matrix_order_2");
            Method::do_execute(src_field, tgt_field);
        }
    }
    {
        ATLAS_TRACE("halo exchange target");
        tgt_field.set_dirty(true);
        tgt_field.haloExchange();
    }
}

void ConservativeMethod::setup_stat() const {
    const auto& src_cell_halo  = array::make_view<int, 1>(src_mesh_.cells().halo());
    const auto& src_node_ghost = array::make_view<int, 1>(src_mesh_.nodes().ghost());
    const auto src_areas_v     = array::make_view<double, 1>(src_areas_);
    const auto tgt_areas_v     = array::make_view<double, 1>(tgt_areas_);
    double geo_create_err = 0.;
    double src_tgt_sums[2]          = {0., 0.};
    if (src_cell_data_) {
        for (idx_t spt = 0; spt < src_areas_v.size(); ++spt) {
            if (not src_cell_halo(spt)) {
                src_tgt_sums[0] += src_areas_v(spt);
            }
        }
    }
    else {
        for (idx_t src = 0; src < src_areas_v.size(); ++src) {
            if (not src_node_ghost(src)) {
                src_tgt_sums[0] += src_areas_v(src);
            }
        }
    }
    const auto& tgt_cell_halo  = array::make_view<int, 1>(tgt_mesh_.cells().halo());
    const auto& tgt_node_ghost = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
    if (tgt_cell_data_) {
        for (idx_t tpt = 0; tpt < tgt_areas_v.size(); ++tpt) {
            if (not tgt_cell_halo(tpt)) {
                src_tgt_sums[1] += tgt_areas_v(tpt);
            }
        }
    }
    else {
        for (idx_t tpt = 0; tpt < tgt_areas_v.size(); ++tpt) {
            if (not tgt_node_ghost(tpt)) {
                src_tgt_sums[1] += tgt_areas_v(tpt);
            }
        }
    }
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(src_tgt_sums, 2, eckit::mpi::sum());
    }
    geo_create_err = std::abs(src_tgt_sums[0] - src_tgt_sums[1]) * 0.25 * M_1_PI;
    remap_stat_.errors[RemapStat::Errors::GEO_DIFF] = geo_create_err;
}

void ConservativeMethod::remap_stat(const FieldArray& src_vals, const FieldArray& tgt_vals,
                                    FieldArray* diff_vals, double func(const PointLonLat&)) const {
    const auto& src_cell_halo  = array::make_view<int, 1>(src_mesh_.cells().halo());
    const auto& src_node_ghost = array::make_view<int, 1>(src_mesh_.nodes().ghost());
    const auto& src_node_halo  = array::make_view<int, 1>(src_mesh_.nodes().halo());
    const auto& tgt_cell_halo  = array::make_view<int, 1>(tgt_mesh_.cells().halo());
    const auto& tgt_node_ghost = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
    const auto& tgt_node_halo  = array::make_view<int, 1>(tgt_mesh_.nodes().halo());
    const auto src_areas_v     = array::make_view<double, 1>(src_areas_);
    const auto tgt_areas_v     = array::make_view<double, 1>(tgt_areas_);
    double err_remap_cons  = 0.;
    double err_remap_l2    = 0.;
    double err_remap_linf  = 0.;
    if (src_cell_data_) {
        for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
            if (src_cell_halo(spt)) {
                continue;
            }
            double diff = src_vals(spt) * src_areas_v(spt);
            err_remap_cons += diff;
            const auto& iparam = src_iparam_[spt];
            if (tgt_cell_data_) {
                for (idx_t icell = 0; icell < iparam.src_weights.size(); ++icell) {
                    idx_t tcell = iparam.cell_idx[icell];
                    if (tgt_cell_halo(tcell) < 1) {
                        diff -= tgt_vals(iparam.cell_idx[icell]) * iparam.src_weights[icell];
                    }
                }
            }
            else {
                for (idx_t icell = 0; icell < iparam.src_weights.size(); ++icell) {
                    idx_t tcell = iparam.cell_idx[icell];
                    idx_t tnode = tgt_csp2node_[tcell];
                    if (tgt_node_halo(tnode) < 1) {
                        diff -= tgt_vals(tnode) * iparam.src_weights[icell];
                    }
                }
            }
            if (diff_vals) {
                (*diff_vals)(spt) = std::abs(diff) / src_areas_v(spt);
            }
        }
    }
    else {
        for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
            if (src_node_ghost(spt) or src_areas_v(spt) < 1e-14) {
                if (diff_vals) {
                    (*diff_vals)(spt) = 0.;
                }
                continue;
            }
            double diff = src_vals(spt) * src_areas_v(spt);
            err_remap_cons += diff;
            const auto& node2csp = src_node2csp_[spt];
            for (idx_t subcell = 0; subcell < node2csp.size(); ++subcell) {
                const auto& iparam = src_iparam_[node2csp[subcell]];
                if (tgt_cell_data_) {
                    for (idx_t icell = 0; icell < iparam.src_weights.size(); ++icell) {
                        diff -= tgt_vals(iparam.cell_idx[icell]) * iparam.src_weights[icell];
                    }
                }
                else {
                    for (idx_t icell = 0; icell < iparam.src_weights.size(); ++icell) {
                        idx_t tcell = iparam.cell_idx[icell];
                        idx_t tnode = tgt_csp2node_[tcell];
                        diff -= tgt_vals(tnode) * iparam.src_weights[icell];
                    }
                }
            }
            if (diff_vals) {
                (*diff_vals)(spt) = std::abs(diff) / src_areas_v(spt);
            }
        }
    }
    if (tgt_cell_data_) {
        for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
            if (tgt_cell_halo(tpt)) {
                continue;
            }
            err_remap_cons -= tgt_vals(tpt) * tgt_areas_v(tpt);
            auto p = tgt_points_[tpt];
            PointLonLat pll;
            eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
            double err_l = std::abs(tgt_vals(tpt) - func(pll));
            err_remap_l2 += err_l * err_l * tgt_areas_v(tpt);
            err_remap_linf = std::max(err_remap_linf, err_l);
        }
    }
    else {
        for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
            if (tgt_node_ghost(tpt)) {
                continue;
            }
            err_remap_cons -= tgt_vals(tpt) * tgt_areas_v(tpt);
            auto p = tgt_points_[tpt];
            PointLonLat pll;
            eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
            double err_l = std::abs(tgt_vals(tpt) - func(pll));
            err_remap_l2 += err_l * err_l * tgt_areas_v(tpt);
            err_remap_linf = std::max(err_remap_linf, err_l);
        }
    }
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(&err_remap_cons, 1, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&err_remap_l2, 1, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&err_remap_linf, 1, eckit::mpi::max());
    }
    remap_stat_.errors[RemapStat::Errors::REMAP_L2]    = std::sqrt(err_remap_l2 * 0.25 * M_1_PI);
    remap_stat_.errors[RemapStat::Errors::REMAP_LINF]  = err_remap_linf;
    remap_stat_.errors[RemapStat::Errors::REMAP_CONS]  = std::sqrt(std::abs(err_remap_cons) * 0.25 * M_1_PI);
    remap_stat_.remap_computed = true;
}

auto debug_polygons = [](const CSPolygon& plg_1, const CSPolygon& plg_2) {
    auto iplg          = plg_1.intersect( plg_2 );
    auto jplg          = plg_2.intersect( plg_1 );
    const double intersection_comm_err = std::abs(iplg.area() - jplg.area()) / (plg_1.area() > 0 ? plg_1.area() : 1.);
    Log::info().indent();
    if (intersection_comm_err > 1e-6) {
        Log::info() << "PLG_1       : " << std::setprecision(10) << plg_1 << "\n";
        Log::info() << "area(PLG_1) : " << plg_1.area() << "\n";
        Log::info() << "PLG_2      :" << plg_2 << "\n";
        Log::info() << "PLG_12 = PLG_1 ^ PLG_2  : " << iplg << "\n";
        Log::info() << "PLG_21 = PLG_2 ^ PLG_1  : " << jplg << "\n";
        Log::info() << "area(PLG_12 - PLG_21)   : " << intersection_comm_err << "\n";
        Log::info() << "area(PLG_21)            : " << jplg.area() << "\n";
        //ATLAS_ASSERT( false, "SRC.intersect.TGT =/= TGT.intersect.SRC.");
    }
    int pout;
    int pin = inside_vertices(plg_1, plg_2, pout);
    if ( pin > 2 && iplg.area() < 3e-16 ) {
        Log::info() << " pin : " << pin << ", pout :" << pout 
                    << ", total vertices : " << plg_2.size()<< "\n";
        Log::info() << "PLG_2                   :" << plg_2 << "\n";
        Log::info() << "PLG_12 = PLG_1 ^ PLG_2  : " << iplg << "\n";
        Log::info() << "area(PLG_12)            : " << iplg.area() << "\n\n";
        //ATLAS_ASSERT( false, "SRC must intersect TGT." );
    }
    Log::info().unindent();
};

void ConservativeMethod::dump_intersection(const CSPolygon& plg_1, const CSPolygonArray& plg_2_array,
                                           const std::vector<idx_t>& plg_2_idx_array) const {
    for (int i = 0; i < plg_2_idx_array.size(); ++i) {
        const auto plg_2_idx  = plg_2_idx_array[i];
        const auto& plg_2     = std::get<0>(plg_2_array[plg_2_idx]);
        debug_polygons(plg_1, plg_2);
    }
}

template <class TargetCellsIDs>
void ConservativeMethod::dump_intersection(const CSPolygon& plg_1, const CSPolygonArray& plg_2_array,
                                           const TargetCellsIDs& plg_2_idx_array) const {
    for (int i = 0; i < plg_2_idx_array.size(); ++i) {
        const auto plg_2_idx  = plg_2_idx_array[i].payload();
        const auto& plg_2     = std::get<0>(plg_2_array[plg_2_idx]);
        debug_polygons(plg_1, plg_2);
    }
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
