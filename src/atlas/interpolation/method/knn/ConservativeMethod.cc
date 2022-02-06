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

// get counter-clockwise sorted neighbours of a cell
std::vector<idx_t> ConservativeMethod::get_cell_neighbours(Mesh& mesh, idx_t cell) const {
    const auto& cell2node  = mesh.cells().node_connectivity();
    const auto& node2cell  = mesh.nodes().cell_connectivity();
    const auto& nodes_ll   = array::make_view<double, 2>(mesh.nodes().lonlat());
    const auto n2c_missval = node2cell.missing_value();
    const idx_t n_nodes    = cell2node.cols(cell);
    std::vector<idx_t> nbr_cells;
    nbr_cells.reserve(n_nodes);

    for (idx_t inode = 0; inode < n_nodes; ++inode) {
		idx_t node0             = cell2node(cell, inode);
		idx_t node1             = cell2node(cell, inode!=n_nodes-1 ? inode+1 : 0);
		const PointLonLat p0_ll = PointLonLat{nodes_ll(node0, 0), nodes_ll(node0, 1)};
		const PointLonLat p1_ll = PointLonLat{nodes_ll(node1, 0), nodes_ll(node1, 1)};
		PointXYZ p0, p1;
		eckit::geometry::Sphere::convertSphericalToCartesian(1., p0_ll, p0);
		eckit::geometry::Sphere::convertSphericalToCartesian(1., p1_ll, p1);
		if (PointXYZ::norm(p0 - p1) < 1e-14) {
			continue;  // edge = point
		}
		bool still_search = true; // still search the cell having vertices node0 & node1, not havin index "cell"
		int n_cells0 = node2cell.cols( node0 );
		int n_cells1 = node2cell.cols( node1 );
		for( int icell0 = 0; still_search && icell0 < n_cells0; icell0++ ) {
			int cell0 = node2cell( node0, icell0 );
			if ( cell0 == cell ) {
				continue;
			}
			for( int icell1 = 0; still_search && icell1 < n_cells1; icell1++ ) {
				int cell1 = node2cell( node1, icell1 );
				if ( cell0 == cell1 && cell0 != n2c_missval && cell0 != cell ) {
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
    const auto nodes_ll   = array::make_view<double, 2>(mesh.nodes().lonlat());
    const auto& cell2node = mesh.cells().node_connectivity();
    const auto cell_halo  = array::make_view<int, 1>(mesh.cells().halo());
    const auto cell_flags = array::make_view<int, 1>(mesh.cells().flags());
    const auto cell_part  = array::make_view<int, 1>(mesh.cells().partition());
    auto xyz2ll = [](atlas::PointXYZ& p_xyz) {
        PointLonLat p_ll;
        eckit::geometry::Sphere::convertCartesianToSpherical(1., p_xyz, p_ll);
        return p_ll;
    };
    idx_t cspol_id = 0; // subpolygon enumeration
    double total_csp_area_shoots = 0; // total over/undershoots in creation of subpolygons
    double max_csp_area_shoots = 0; // max over/undershoots in creation of subpolygons
    for (idx_t cell = 0; cell < mesh.cells().size(); ++cell) {
        ATLAS_ASSERT(cell < cell2node.rows());
        const idx_t n_nodes = cell2node.cols(cell);
        PointXYZ cell_mid(0., 0., 0.);	// cell centre
        std::vector<PointLonLat> pts_ll;
        pts_ll.reserve( n_nodes );
        for (idx_t inode = 0; inode < n_nodes; ++inode) {
            idx_t node0             = cell2node(cell, inode);
            idx_t node1             = cell2node(cell, inode!=n_nodes-1 ? inode+1 : 0);
            const PointLonLat p0_ll = PointLonLat{nodes_ll(node0, 0), nodes_ll(node0, 1)};
            const PointLonLat p1_ll = PointLonLat{nodes_ll(node1, 0), nodes_ll(node1, 1)};
            pts_ll.emplace_back( p0_ll );
            PointXYZ p0, p1;
            eckit::geometry::Sphere::convertSphericalToCartesian(1., p0_ll, p0);
            eckit::geometry::Sphere::convertSphericalToCartesian(1., p1_ll, p1);
            if (PointXYZ::norm(p0 - p1) < 1e-14) {
                continue;  // skip this edge = a pole point
            }
            cell_mid = cell_mid + p0;
            cell_mid = cell_mid + p1;
        }
        CSPolygon csp( pts_ll );
        double loc_csp_area_shoot = csp.area();
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
                continue;  // skip this edge = a pole point
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
            CSPolygon cspi(pts_ll);
            loc_csp_area_shoot -= cspi.area();
            cspolygons.emplace_back(cspi, halo_type);
            cspol_id++;
        }
        total_csp_area_shoots = std::abs(loc_csp_area_shoot) + total_csp_area_shoots;
        max_csp_area_shoots = std::max( std::abs(loc_csp_area_shoot), max_csp_area_shoots );
    }
    Log::info() << "Created " << cspolygons.size() << " polygons from "
                << mesh.cells().size() << " mesh cells.\n";
    Log::info() << "Total sum of over/undershoots " << total_csp_area_shoots << ", max over/undershoots per cell"
                << max_csp_area_shoots << "\n";
    return cspolygons;
}

void ConservativeMethod::do_setup(const Grid& src_grid, const Grid& tgt_grid) {
    ATLAS_TRACE("ConservativeMethod::do_setup( Grid, Grid )");
    ATLAS_ASSERT(src_grid);
    ATLAS_ASSERT(tgt_grid);
    auto src_mesh_config = src_grid.meshgenerator();
    auto tgt_mesh_config = tgt_grid.meshgenerator();
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
    mesh::actions::build_node_to_cell_connectivity(src_mesh_);
    mesh::actions::build_node_to_cell_connectivity(tgt_mesh_);

    CSPolygonArray src_csp;
    CSPolygonArray tgt_csp;
    {
        ATLAS_TRACE("Get source polygons");
        if (src_cell_data_) {
            src_csp = get_polygons_celldata(src_mesh_);
        }
        else {
            src_csp = get_polygons_nodedata(src_mesh_, src_csp2node_, src_node2csp_);
        }
    }
    {
        ATLAS_TRACE("Get target polygons");
        if (tgt_cell_data_) {
            tgt_csp = get_polygons_celldata(tgt_mesh_);
        }
        else {
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
            max_tgtcell_rad = std::max(max_tgtcell_rad, t_csp.radius());
        }
    }
    kdt_search.build();

    enum MeshSizeId {SRC, TGT, SRC_TGT_INTERSECT, SRC_NONINTERSECT};
    std::array<size_t,4> num_pol{0,0,0,0};
    enum AreaCoverageId {TOTAL_SRC, MAX_SRC};
    std::array<double,2> area_coverage{0.,0.};
    iparam_.resize(src_csp.size());
    eckit::Channel blackhole;
    eckit::ProgressTimer progress("Intersecting polygons ", src_csp.size(), " cell", double(10),
                                  src_csp.size() > 50 ? Log::info() : blackhole);
    for (idx_t scell = 0; scell < src_csp.size(); ++scell, ++progress) {
        const int cell_flag = std::get<1>(src_csp[scell]);
        if (cell_flag == -1) { // skip periodic cells
            continue;
        }
        const auto& s_csp   = std::get<0>(src_csp[scell]);
		const double s_csp_area = s_csp.area();
        double covered_area = 0.;
        auto tgt_cells = kdt_search.closestPointsWithinRadius(s_csp.centroid(), s_csp.radius() + max_tgtcell_rad);
        for (idx_t ttcell = 0; ttcell < tgt_cells.size(); ++ttcell) {
            auto tcell        = tgt_cells[ttcell].payload();
            const auto& t_csp = std::get<0>(tgt_csp[tcell]);
            CSPolygon csp_i   = s_csp.intersect(t_csp);
			double csp_i_area = csp_i.area();
            if (csp_i_area > 0.) {
                iparam_[scell].tcell_id.emplace_back(tcell);
                iparam_[scell].weights.emplace_back(csp_i_area);
                iparam_[scell].sweights.emplace_back(csp_i_area / t_csp.area());
                iparam_[scell].centroids.emplace_back(csp_i.centroid());
                covered_area += csp_i_area;
            }
        }
        const double loc_csp_error = s_csp_area - covered_area;
        if ( loc_csp_error > 1e-5 && cell_flag == 0) {
            Log::info() << "WARNING src cell area NOT covered: " << loc_csp_error << "\n";
            //dump_intersection( s_csp, tgt_csp, tgt_cells );
            //ATLAS_ASSERT( false );
        }
        if (cell_flag == 0) {
            area_coverage[TOTAL_SRC] += loc_csp_error;
            area_coverage[MAX_SRC]   = std::max( area_coverage[MAX_SRC], loc_csp_error );
        }
        if (iparam_[scell].tcell_id.size() == 0.) {
            num_pol[SRC_NONINTERSECT]++;
        }
        if (normalise_intersections_ && loc_csp_error < 1e-5) {
            double wfactor = s_csp.area() / (covered_area > 1e-10 ? covered_area : 1.);
            for (idx_t i = 0; i < iparam_[scell].weights.size(); i++) {
                iparam_[scell].weights[i] *= wfactor;
                iparam_[scell].sweights[i] *= wfactor;
            }
        }
        num_pol[SRC_TGT_INTERSECT] += iparam_[scell].weights.size();
    }
    num_pol[SRC] = src_csp.size();
    num_pol[TGT] = tgt_csp.size();
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(&num_pol[0], 4, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&area_coverage[0], 2, eckit::mpi::max());
    }
    Log::info() << "ConservativeMethod:: num_src_polygons, num_tgt_polygons, num_intersect_polygons: "
                << num_pol[0] << " " << num_pol[1] << " " << num_pol[2] << "\n";
    Log::info() << "ConservativeMethod::intersect_polygons : " << num_pol[SRC_NONINTERSECT]
                << " source mesh polygons do NOT intersect any other polygon.\n";
    Log::info() << "ConservativeMethod::intersect_polygons : " << area_coverage[TOTAL_SRC]
                << " total area of source polygons over/undercovered by target polygons.\n";
    Log::info() << "ConservativeMethod::intersect_polygons : " << area_coverage[MAX_SRC]
                << " maximal area of a source polygons NOT covered by target polygons.\n";
    geo_err_intsc_l1_   =  0.;
    geo_err_intsc_linf_ =  0.;
    for (idx_t scell = 0; scell < src_csp.size(); ++scell) {
        const int cell_flag = std::get<1>(src_csp[scell]);
        if (cell_flag == -1 or cell_flag > 0) {
            // skip periodic & halo cells
            continue;
        }
        double diff_cell = std::get<0>(src_csp[scell]).area();
        for (idx_t icell = 0; icell < iparam_[scell].weights.size(); ++icell) {
            diff_cell -= iparam_[scell].weights[icell];
        }
        geo_err_intsc_l1_   += std::abs(diff_cell);
        geo_err_intsc_linf_ = std::max(geo_err_intsc_linf_, std::abs(diff_cell));
    }
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(&geo_err_intsc_l1_, 1, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&geo_err_intsc_linf_, 1, eckit::mpi::sum());
    }
    geo_err_intsc_l1_ *= 0.25 * M_1_PI;
    Log::info() << "ConservativeMethod::intersect_polygons : cons err in polygon intersect  : (L1) "
                << geo_err_intsc_l1_ << " (Lmax) " << geo_err_intsc_linf_ << "\n";
}

void ConservativeMethod::setup_1st_order_matrix() {
    ATLAS_TRACE("ConservativeMethod::setup: build cons-1 interpolant matrix");
	order_ = 1;
	matrix_free_ = false;
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
					double weight = tgt_areas_v(tnode);
					weight = iparam.weights[icell] / (weight>0.? weight : 1.);
                    triplets.emplace_back(tnode, snode, weight);
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
	matrix_free_ = false;
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
            ATLAS_ASSERT( Cs_norm > 0. );
            Cs = PointXYZ::div( Cs, Cs_norm );
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
                if (CSPolygon::leftOf(Cnsj, Cs, Csj, 1e-16)) {
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
                if (CSPolygon::leftOf(Nsnj, Ns, Nsj, 1e-16)) {
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
                        double csp2node_coef = tgt_areas_v(tnode);
						csp2node_coef = ( csp2node_coef>0. ? csp2node_coef : 1. );
                        csp2node_coef = iparam.weights[icell] / iparam.sweights[icell] / csp2node_coef;
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
                        val *= (csp.leftOf(Pnn, P, Pn, 1e-16) ? -1 : 1);
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
    Log::info() << " ConservativeMethod::stat : global error in polygon create   : " << geo_create_err << "\n";
}

void ConservativeMethod::remap_stat(const FieldArray& src_vals, const FieldArray& tgt_vals, FieldArray& diff_vals,
                                    double func(const PointLonLat&), std::array<double,3>& errors ) const {
    const auto& src_cell_halo  = array::make_view<int, 1>(src_mesh_.cells().halo());
    const auto& src_node_ghost = array::make_view<int, 1>(src_mesh_.nodes().ghost());
    const auto& tgt_cell_halo  = array::make_view<int, 1>(tgt_mesh_.cells().halo());
    const auto& tgt_node_ghost = array::make_view<int, 1>(tgt_mesh_.nodes().ghost());
    const auto src_areas_v     = array::make_view<double, 1>(src_areas_);
    const auto tgt_areas_v     = array::make_view<double, 1>(tgt_areas_);
    errors = {{0., 0., 0.}};
    if (src_cell_data_) {
        for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
            if (src_cell_halo(spt)) {
                continue;
            }
            double diff = src_vals(spt) * src_areas_v(spt);
            errors[GLOBAL] += diff;
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
            errors[GLOBAL] += diff;
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
    if (tgt_cell_data_) {
        for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
            if (tgt_cell_halo(tpt)) {
                continue;
            }
            errors[GLOBAL] -= tgt_vals(tpt) * tgt_areas_v(tpt);
            auto p = tgt_points_[tpt];
            PointLonLat pll;
            eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
            double err_l = std::abs(tgt_vals(tpt) - func(pll));
            errors[L2] += err_l * err_l * tgt_areas_v(tpt);
            errors[LINF] = std::max(errors[LINF], err_l);
        }
    }
    else {
        for (idx_t tpt = 0; tpt < tgt_vals.size(); ++tpt) {
            if (tgt_node_ghost(tpt)) {
                continue;
            }
            errors[GLOBAL] -= tgt_vals(tpt) * tgt_areas_v(tpt);
            auto p = tgt_points_[tpt];
            PointLonLat pll;
            eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
            double err_l = std::abs(tgt_vals(tpt) - func(pll));
            errors[L2] += err_l * err_l * tgt_areas_v(tpt);
            errors[LINF] = std::max(errors[LINF], err_l);
        }
    }
    ATLAS_TRACE_MPI(ALLREDUCE) {
        mpi::comm().allReduceInPlace(&errors[0], 2, eckit::mpi::sum());
        mpi::comm().allReduceInPlace(&errors[2], 1, eckit::mpi::max());
    }
    errors[L2]     = std::sqrt(errors[L2] * 0.25 * M_1_PI);
    errors[GLOBAL] = std::sqrt(std::abs(errors[GLOBAL]) * 0.25 * M_1_PI);
}

template <class TargetCellsIDs>
void ConservativeMethod::dump_intersection(const CSPolygon& s_csp, const CSPolygonArray& tgt_csp,
                                           const TargetCellsIDs& tgt_cells) const {
    Log::info().flush();
    Log::info() << "\n === DEBUG ===\n\n";
    Log::info() << "* src cell: " << std::setprecision(10) << s_csp << "\n";
    Log::info() << "* src area: " << s_csp.area() << "\n\n";
    double area_ncov = s_csp.area();
    for (int i = 0; i < tgt_cells.size(); ++i) {
        const auto tcell  = tgt_cells[i].payload();
        const auto& t_csp = std::get<0>(tgt_csp[tcell]);
        Log::info() << "* src cell: " << s_csp << "\n";
        Log::info() << "* tgt cell: " << t_csp << "\n";
        auto iplg          = s_csp.intersect(t_csp);
        auto jplg          = t_csp.intersect(s_csp);
        const double darea = std::abs(iplg.area() - jplg.area());
        Log::info() << "* src ^ tgt      : " << iplg << "\n";
        Log::info() << "* src ^ tgt area : " << iplg.area() << "\n";
        if (darea > 5e-15) {
            s_csp.intersect(t_csp);
            Log::info() << "* (!!) intersect comm area diff: " << darea << "\n";
            Log::info() << "* (!!) tgt ^ src      : " << jplg << "\n";
            Log::info() << "* (!!) tgt ^ src area : " << jplg.area() << "\n";
            t_csp.intersect(s_csp);
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
