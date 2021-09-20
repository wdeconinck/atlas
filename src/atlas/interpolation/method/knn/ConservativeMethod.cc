/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#include <iomanip>
#include <vector>

#include "eckit/log/ProgressTimer.h"

#include "atlas/grid.h"
#include "atlas/interpolation/method/knn/ConservativeMethod.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/KDTree.h"


namespace atlas {
namespace interpolation {
namespace method {

using CSPolygon = util::ConvexSphericalPolygon;

ConservativeMethod::ConservativeMethod( const util::Config& config ) : Method( config ) {
    config.get( "order", order_ = 1 );
    config.get( "normalise_intersections", normalise_intersections_ = 1 );
    config.get( "field_value_type", fvtype_ = 0 );
    config.get( "matrix_free", matrix_free_ = true );
    config.get( "src_cell_data", src_cell_data_ = true );
    config.get( "tgt_cell_data", tgt_cell_data_ = true );
}

// get counter-clockwise sorted edges of a cell
std::vector<idx_t> ConservativeMethod::sort_cell_edges( Mesh& mesh, idx_t cell_id ) const {
    const auto& cell2edge = mesh.cells().edge_connectivity();
    const auto& cell2node = mesh.cells().node_connectivity();
    const auto& edge2node = mesh.edges().node_connectivity();
    const int nnodes      = cell2node.cols( cell_id );
    const int nedges      = cell2edge.cols( cell_id );
    std::vector<idx_t> edges;
    edges.resize( nedges );
    idx_t ii = 0;
    for ( int inode = 0; inode < nnodes; ++inode ) {
        idx_t node  = cell2node( cell_id, inode );
        idx_t nnode = cell2node( cell_id, ( inode != nnodes - 1 ? inode + 1 : 0 ) );
        for ( int iedge = 0; iedge < nedges; ++iedge ) {
            const idx_t edge  = cell2edge( cell_id, iedge );
            const idx_t node0 = edge2node( edge, 0 );
            const idx_t node1 = edge2node( edge, 1 );
            if ( ( node0 == node && node1 == nnode ) or ( node0 == nnode && node1 == node ) ) {
                edges[ii++] = edge;
                break;
            }
        }
    }
    return edges;
}

// get counter-clockwise sorted neighbours of a cell
std::vector<idx_t> ConservativeMethod::get_cell_neighbours( Mesh& mesh, idx_t cell_id ) const {
    const auto& cell2edge  = mesh.cells().edge_connectivity();
    const auto& edge2cell  = mesh.edges().cell_connectivity();
    const auto& edge2node  = mesh.edges().node_connectivity();
    auto c2e_missval       = cell2edge.missing_value();
    const auto& edges_sort = sort_cell_edges( mesh, cell_id );
    const idx_t nedges     = cell2edge.cols( cell_id );
    std::vector<idx_t> nbr_cells;
    nbr_cells.reserve( nedges );

    for ( idx_t iedge = 0, inbr = 0; iedge < nedges; ++iedge ) {
        const idx_t edge  = edges_sort[iedge];
        const idx_t c1_id = edge2cell( edge, 0 );
        if ( c1_id != c2e_missval && c1_id != cell_id ) {
            nbr_cells.emplace_back( c1_id );
            continue;
        }
        const idx_t c2_id = edge2cell( edge, 1 );
        if ( c2_id != c2e_missval && c2_id != cell_id ) {
            nbr_cells.emplace_back( c2_id );
            continue;
        }
    }
    return nbr_cells;
}


std::vector<CSPolygon> ConservativeMethod::get_polygons_celldata( Mesh& mesh ) const {
    std::vector<CSPolygon> cspolygons;
    // cell centred data -> polygons are mesh elements
    const idx_t n_cells = mesh.cells().size();
    cspolygons.resize( n_cells );
    const auto& cell2node = mesh.cells().node_connectivity();
    const auto lonlat     = array::make_view<double, 2>( mesh.nodes().lonlat() );
    std::vector<PointLonLat> pts_ll;
    for ( idx_t icell = 0; icell < n_cells; ++icell ) {
        const idx_t n_nodes = cell2node.cols( icell );
        pts_ll.clear();
        pts_ll.resize( n_nodes );
        for ( idx_t jnode = 0; jnode < n_nodes; ++jnode ) {
            idx_t inode   = cell2node( icell, jnode );
            pts_ll[jnode] = PointLonLat{lonlat( inode, 0 ), lonlat( inode, 1 )};
        }
        cspolygons[icell] = CSPolygon( pts_ll );
    }
    return cspolygons;
}

std::vector<CSPolygon> ConservativeMethod::get_polygons_nodedata( Mesh& mesh ) const {
    std::vector<CSPolygon> cspolygons;
    // cell vertex data -> polygons are (cell_centre, edge_centre, cell_vertex, edge_centre)
    mesh::actions::build_edges( mesh, util::Config( "pole_edges", false ) );
    const auto xy           = array::make_view<double, 2>( mesh.nodes().xy() );
    const auto nodes_ll     = array::make_view<double, 2>( mesh.nodes().lonlat() );
    auto edge_flags         = array::make_view<int, 1>( mesh.edges().flags() );
    const auto& cell2edge   = mesh.cells().edge_connectivity();
    const auto& edge2node   = mesh.edges().node_connectivity();
    const auto& field_flags = array::make_view<int, 1>( mesh.cells().flags() );

    // compute gradient estimate now
    //Field node_grad_field   = Field( "node_grad", array::make_datatype<double>(), array::make_shape( 3 ) );
    //		auto node_grad      = array::make_view<double, 3>( node_grad_field );

    /*      TODO would be nice but does not work with healpix!!
	auto patch = [&field_flags]( idx_t e ) {
		using Topology = atlas::mesh::Nodes::Topology;
		return Topology::check( field_flags( e ), Topology::PATCH );
	};
	auto is_pole_edge = [&edge_flags]( idx_t e ) { // does not work with H-grids
		using Topology = atlas::mesh::Nodes::Topology;
		return Topology::check( edge_flags( e ), Topology::POLE );
	};
*/
    auto xyz2ll = []( atlas::PointXYZ& p_xyz ) {
        PointLonLat p_ll;
        eckit::geometry::Sphere::convertCartesianToSpherical( 1., p_xyz, p_ll );
        return p_ll;
    };

    for ( idx_t icell = 0; icell < mesh.cells().size(); ++icell ) {
        // get cell centre
        PointXYZ cell_mid   = PointXYZ{0., 0., 0.};
        const idx_t n_edges = cell2edge.cols( icell );
        for ( idx_t iedge = 0; iedge < n_edges; ++iedge ) {
            idx_t edge              = cell2edge( icell, iedge );
            idx_t node0             = edge2node( edge, 0 );
            idx_t node1             = edge2node( edge, 1 );
            const PointLonLat p0_ll = PointLonLat{nodes_ll( node0, 0 ), nodes_ll( node0, 1 )};
            const PointLonLat p1_ll = PointLonLat{nodes_ll( node1, 0 ), nodes_ll( node1, 1 )};
            PointXYZ p0, p1;
            eckit::geometry::Sphere::convertSphericalToCartesian( 1., p0_ll, p0 );
            eckit::geometry::Sphere::convertSphericalToCartesian( 1., p1_ll, p1 );
            if ( PointXYZ::norm( p0 - p1 ) < 1e-14 ) {
                continue;  // pole edge
            }
            cell_mid = cell_mid + p0;
            cell_mid = cell_mid + p1;
        }
        cell_mid = PointXYZ::div( cell_mid, PointXYZ::norm( cell_mid ) );
        PointLonLat cell_ll;
        eckit::geometry::Sphere::convertCartesianToSpherical( 1., cell_mid, cell_ll );
        // get CSPolygon for each valid edge
        for ( idx_t iedge = 0; iedge < n_edges; ++iedge ) {
            idx_t edge               = cell2edge( icell, iedge );
            idx_t node0              = edge2node( edge, 0 );
            idx_t node1              = edge2node( edge, 1 );
            const PointLonLat pi0_ll = PointLonLat{nodes_ll( node0, 0 ), nodes_ll( node0, 1 )};
            const PointLonLat pi1_ll = PointLonLat{nodes_ll( node1, 0 ), nodes_ll( node1, 1 )};
            PointXYZ pi0, pi1;
            eckit::geometry::Sphere::convertSphericalToCartesian( 1., pi0_ll, pi0 );
            eckit::geometry::Sphere::convertSphericalToCartesian( 1., pi1_ll, pi1 );
            PointXYZ iedge_mid = pi0 + pi1;
            iedge_mid          = PointXYZ::div( iedge_mid, PointXYZ::norm( iedge_mid ) );
            // get edge centre
            if ( PointXYZ::norm( iedge_mid - pi0 ) < 1e-14 ) {
                continue;  // skip this edge, it is a pole point
            }
            PointLonLat iedge_mid_ll = xyz2ll( iedge_mid );
            PointXYZ third_point;
            PointLonLat third_point_ll;
            if ( util::ConvexSphericalPolygon::leftOf( pi0, cell_mid, iedge_mid, 1e-14 ) ) {
                third_point    = pi0;
                third_point_ll = pi0_ll;
            }
            else {
                third_point    = pi1;
                third_point_ll = pi1_ll;
            }
            // find the other valid edge touching pi0
            PointXYZ jedge_mid;
            for ( idx_t jedge = 0; jedge < n_edges; ++jedge ) {
                idx_t edge               = cell2edge( icell, jedge );
                idx_t j0                 = edge2node( edge, 0 );
                idx_t j1                 = edge2node( edge, 1 );
                const PointLonLat pj0_ll = PointLonLat{nodes_ll( j0, 0 ), nodes_ll( j0, 1 )};
                const PointLonLat pj1_ll = PointLonLat{nodes_ll( j1, 0 ), nodes_ll( j1, 1 )};
                PointXYZ pj0, pj1;
                eckit::geometry::Sphere::convertSphericalToCartesian( 1., pj0_ll, pj0 );
                eckit::geometry::Sphere::convertSphericalToCartesian( 1., pj1_ll, pj1 );
                if ( PointXYZ::norm( pj0 - pj1 ) < 1e-14 ) {
                    continue;  // pole edge
                }
                if ( jedge != iedge &&
                     ( PointXYZ::norm( third_point - pj0 ) < 1e-14 or PointXYZ::norm( third_point - pj1 ) < 1e-14 ) ) {
                    jedge_mid = pj0 + pj1;
                    jedge_mid = PointXYZ::div( jedge_mid, PointXYZ::norm( jedge_mid ) );
                    break;
                }
            }
            std::vector<PointLonLat> pts_ll( 4 );
            pts_ll[0] = cell_ll;
            pts_ll[1] = iedge_mid_ll;
            pts_ll[2] = third_point_ll;
            pts_ll[3] = xyz2ll( jedge_mid );
            cspolygons.emplace_back( CSPolygon( pts_ll ) );
        }
    }
    ( Log::info() << "Created " << cspolygons.size() << " CSPolygons from " << mesh.cells().size() << " mesh cells\n" )
        .flush();
    return cspolygons;
}

void ConservativeMethod::do_setup( const Grid& src_grid, const Grid& tgt_grid ) {
    ATLAS_TRACE( "ConservativeMethod::do_setup()" );
    ATLAS_ASSERT( src_grid );
    ATLAS_ASSERT( tgt_grid );
    if ( mpi::size() > 1 ) {
        ATLAS_NOTIMPLEMENTED;
    }
    auto src_mesh_config = src_grid.meshgenerator();
    src_mesh_config.set( "include_pole", true );
    src_mesh_ = MeshGenerator( src_mesh_config ).generate( src_grid );
    tgt_mesh_ = MeshGenerator( tgt_grid.meshgenerator() ).generate( tgt_grid );

    std::vector<CSPolygon> src_csp;
    std::vector<CSPolygon> tgt_csp;
    if ( src_cell_data_ ) {
        functionspace::CellColumns src_fs( src_mesh_ );
        src_fs_ = src_fs;
        src_csp = get_polygons_celldata( src_mesh_ );
    }
    else {
        functionspace::NodeColumns src_fs( src_mesh_ );
        src_fs_ = src_fs;
        src_csp = get_polygons_nodedata( src_mesh_ );
    }
    if ( tgt_cell_data_ ) {
        functionspace::CellColumns tgt_fs( tgt_mesh_ );
        tgt_fs_ = tgt_fs;
        tgt_csp = get_polygons_celldata( tgt_mesh_ );
    }
    else {
        functionspace::NodeColumns tgt_fs( tgt_mesh_ );
        tgt_fs_ = tgt_fs;
        tgt_csp = get_polygons_nodedata( tgt_mesh_ );
    }
    do_setup_with_polygons( src_csp, tgt_csp );
}


void ConservativeMethod::do_setup_with_polygons( const std::vector<CSPolygon>& src_csp,
                                                 const std::vector<CSPolygon>& tgt_csp ) {
    n_scells_ = src_csp.size();
    n_tcells_ = tgt_csp.size();

    src_centroids_.resize( n_scells_ );
    tgt_centroids_.resize( n_tcells_ );
    src_areas_.resize( n_scells_ );
    tgt_areas_.resize( n_tcells_ );
    util::KDTree<idx_t> kdt_search;
    kdt_search.reserve( n_tcells_ );

    double max_tgtcell_rad = 0.;
    for ( idx_t jcell = 0; jcell < n_tcells_; ++jcell ) {
        kdt_search.insert( tgt_csp[jcell].centroid(), jcell );
        max_tgtcell_rad = std::max( max_tgtcell_rad, tgt_csp[jcell].cell_radius() );
    }
    kdt_search.build();

    size_t nonintersect        = 0;
    double src_area_notcovered = 0.;
    iparam_.resize( n_scells_ );
    eckit::Channel blackhole;
    eckit::ProgressTimer progress( "Intersecting polygons ", n_scells_, " cell", double( 10 ),
                                   n_scells_ > 50 ? Log::info() : blackhole );
    for ( idx_t scell = 0; scell < n_scells_; ++scell, ++progress ) {
        const auto& s_csp     = src_csp[scell];
        src_centroids_[scell] = s_csp.centroid();
        src_areas_[scell]     = s_csp.area();
        double covered_area   = 0.;
        auto tgt_cells =
            kdt_search.closestPointsWithinRadius( s_csp.centroid(), s_csp.cell_radius() + max_tgtcell_rad );
        for ( idx_t ttcell = 0; ttcell < tgt_cells.size(); ++ttcell ) {
            auto tcell      = tgt_cells[ttcell].payload();
            CSPolygon csp_i = s_csp.intersect( tgt_csp[tcell] );
            if ( csp_i.area() > 0. ) {
                iparam_[scell].tcell_id.emplace_back( tcell );
                iparam_[scell].weights.emplace_back( csp_i.area() );
                iparam_[scell].sweights.emplace_back( csp_i.area() / tgt_csp[tcell].area() );
                iparam_[scell].centroids.emplace_back( csp_i.centroid() );
                covered_area += csp_i.area();
            }
        }
        const double loc_csp_error = std::abs( src_csp[scell].area() - covered_area );
        src_area_notcovered += loc_csp_error;
        if ( iparam_[scell].tcell_id.size() == 0. ) {
            ++nonintersect;
        }
        if ( normalise_intersections_ ) {
            double wfactor = src_csp[scell].area() / ( covered_area > 1e-10 ? covered_area : 1. );
            for ( idx_t i = 0; i < iparam_[scell].weights.size(); i++ ) {
                iparam_[scell].weights[i] *= wfactor;
                iparam_[scell].sweights[i] *= wfactor;
            }
        }
        if ( false && loc_csp_error > 1e-8 ) {
            Log::info() << "* src cell area NOT covered: " << loc_csp_error << "\n";
            dump_intersection( src_csp[scell], tgt_csp, tgt_cells );
            ATLAS_ASSERT( false );
        }
    }
    Log::info() << "WARNING " << nonintersect << " source mesh polygons do NOT intersect any other polygon.\n";
    Log::info() << "WARNING " << src_area_notcovered << " area of source mesh NOT covered by target mesh.\n";
    for ( idx_t tcell = 0; tcell < n_tcells_; ++tcell ) {
        tgt_centroids_[tcell] = tgt_csp[tcell].centroid();
        tgt_areas_[tcell]     = tgt_csp[tcell].area();
    }
    mesh::actions::build_halo( src_mesh_, 0 );
    setup_1st_order_matrix();
    setup_2nd_order_matrix();
}

void ConservativeMethod::setup_1st_order_matrix() {
    if ( order_ != 1 or matrix_free_ ) {
        return;
    }
    ATLAS_TRACE( "ConservativeMethod::setup: build cons-1 interpolant matrix" );
    Triplets triplets;
    size_t triplets_size = 0;
    for ( idx_t scell = 0; scell < n_scells_; ++scell ) {
        const auto& iparam = iparam_[scell];
        for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
            triplets_size++;
        }
    }
    triplets.reserve( triplets_size );
    for ( idx_t scell = 0; scell < n_scells_; ++scell ) {
        const auto& iparam = iparam_[scell];
        for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
            idx_t tcell = iparam.tcell_id[icell];
            triplets.emplace_back( tcell, scell, iparam.sweights[icell] );
        }
    }
    std::sort( std::begin( triplets ), std::end( triplets ), []( const Triplet& t1, const Triplet& t2 ) {
        return ( t1.row() < t2.row() or ( t1.row() == t2.row() and ( t1.col() < t2.col() ) ) );
    } );
    Matrix A( n_tcells_, n_scells_, triplets );
    matrix_shared_->swap( A );
}

void ConservativeMethod::setup_2nd_order_matrix() {
    if ( order_ != 2 or matrix_free_ ) {
        return;
    }
    ATLAS_TRACE( "ConservativeMethod::setup: build cons-2 interpolant matrix" );
    mesh::actions::build_edges( src_mesh_, util::Config( "pole_edges", false ) );
    //const auto halo     = array::make_view<int, 1>( src_mesh_.cells().halo() );
    Triplets triplets;
    size_t triplets_size = 0;
    for ( idx_t scell = 0; scell < n_scells_; ++scell ) {
        const auto nb_cells = get_cell_neighbours( src_mesh_, scell );
        triplets_size += ( 2 * nb_cells.size() + 1 ) * iparam_[scell].centroids.size();
    }
    triplets.reserve( triplets_size );
    for ( idx_t scell = 0; scell < n_scells_; ++scell ) {
        const auto& iparam = iparam_[scell];
        if ( iparam.centroids.size() == 0 ) {
            Log::info() << " WARNING source cell " << scell << " not covered "
                        << "\n";
            continue;
        }
        PointXYZ Ci = {0., 0., 0.};
        for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
            Ci = Ci + PointXYZ::mul( iparam.centroids[icell], iparam.weights[icell] );
        }
        const double Ci_norm = PointXYZ::norm( Ci );
        ATLAS_ASSERT( Ci_norm > 0. );
        Ci                   = PointXYZ::div( Ci, Ci_norm );
        double dual_area_inv = 0.;
        const auto nb_cells  = get_cell_neighbours( src_mesh_, scell );
        std::vector<PointXYZ> Rij;
        Rij.resize( nb_cells.size() );
        for ( idx_t nb_id = 0; nb_id < nb_cells.size(); ++nb_id ) {
            idx_t nnb_id = ( nb_id != nb_cells.size() - 1 ) ? nb_id + 1 : 0;
            idx_t ncell  = nb_cells[nb_id];
            idx_t nncell = nb_cells[nnb_id];
            if ( ncell != scell && nncell != scell ) {
                const auto& Cij1 = src_centroids_[ncell];
                const auto& Cij2 = src_centroids_[nncell];
                auto csp         = CSPolygon( {Cij1, Cij2, Ci} );
                if ( csp.area() < std::numeric_limits<double>::epsilon() ) {
                    csp = CSPolygon( {Cij1, Ci, Cij2} );
                }
                dual_area_inv += csp.area();
                if ( CSPolygon::leftOf( Cij2, Ci, Cij1, 1e-16, 0 ) ) {
                    Rij[nb_id] = PointXYZ::cross( Cij2, Cij1 );
                }
                else {
                    Rij[nb_id] = PointXYZ::cross( Cij1, Cij2 );
                }
            }
            else {
                Rij[nb_id] = {0., 0., 0.};
            }
        }
        dual_area_inv = ( dual_area_inv > 0. ) ? 1. / dual_area_inv : 1.;
        std::vector<PointXYZ> Aik;
        Aik.resize( iparam.centroids.size() );
        for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
            const PointXYZ& Cik   = iparam.centroids[icell];
            const PointXYZ Cik_Ci = Cik - Ci;
            Aik[icell]            = Cik_Ci - PointXYZ::mul( Ci, PointXYZ::dot( Ci, Cik_Ci ) );
            Aik[icell]            = PointXYZ::mul( Aik[icell], iparam.sweights[icell] * dual_area_inv );
        }
        PointXYZ Rij_sum = {0., 0., 0.};
        for ( idx_t nb_id = 0; nb_id < nb_cells.size(); ++nb_id ) {
            Rij_sum = Rij_sum + Rij[nb_id];
        }
        for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
            for ( idx_t nb_id = 0; nb_id < nb_cells.size(); ++nb_id ) {
                idx_t nnb_id = ( nb_id != nb_cells.size() - 1 ) ? nb_id + 1 : 0;
                idx_t ij1    = nb_cells[nb_id];
                idx_t ij2    = nb_cells[nnb_id];
                triplets.emplace_back( iparam.tcell_id[icell], ij1, 0.5 * PointXYZ::dot( Rij[nb_id], Aik[icell] ) );
                triplets.emplace_back( iparam.tcell_id[icell], ij2, 0.5 * PointXYZ::dot( Rij[nb_id], Aik[icell] ) );
            }
            triplets.emplace_back( iparam.tcell_id[icell], scell,
                                   iparam.sweights[icell] - PointXYZ::dot( Rij_sum, Aik[icell] ) );
        }
    }
    std::sort( std::begin( triplets ), std::end( triplets ), []( const Triplet& t1, const Triplet& t2 ) {
        return ( t1.row() < t2.row() or ( t1.row() == t2.row() and t1.col() < t2.col() ) );
    } );
    Matrix A( n_tcells_, n_scells_, triplets );
    matrix_shared_->swap( A );
}

void ConservativeMethod::do_execute( const Field& src_field, Field& tgt_field ) {
    ATLAS_TRACE( "ConservativeMethod::do_execute()" );

    const auto src_vals = array::make_view<double, 1>( src_field );
    auto tgt_vals       = array::make_view<double, 1>( tgt_field );
    //const auto halo     = array::make_view<int, 1>( src_mesh_.cells().halo() );

    const auto& src_cell2edge = src_mesh_.cells().edge_connectivity();
    const auto& src_edge2cell = src_mesh_.edges().cell_connectivity();
    const auto& src_edge2node = src_mesh_.edges().node_connectivity();

    if ( order_ == 1 ) {
        if ( matrix_free_ ) {
            for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
                tgt_vals( tcell ) = 0.;
            }
            for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
                const auto& iparam = iparam_[scell];
                for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
                    tgt_vals( iparam.tcell_id[icell] ) += iparam.weights[icell] * src_vals( scell );
                }
            }
            for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
                tgt_vals( tcell ) /= tgt_areas_[tcell];
            }
        }
        else {
            Method::do_execute( src_field, tgt_field );
        }
    }
    else if ( order_ == 2 ) {
        if ( matrix_free_ ) {
            for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
                tgt_vals( tcell ) = 0.;
            }
            for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
                //if ( halo( scell ) ) {
                //    continue;
                //}
                const auto& iparam       = iparam_[scell];
                const PointXYZ& P        = src_centroids_[scell];
                PointXYZ grad            = {0., 0., 0.};
                PointXYZ src_barycenter  = {0., 0., 0.};
                auto src_neighbour_cells = get_cell_neighbours( src_mesh_, scell );
                double dual_area         = 0.;
                for ( idx_t nb_id = 0; nb_id < src_neighbour_cells.size(); ++nb_id ) {
                    idx_t nnb_id    = ( nb_id != src_neighbour_cells.size() - 1 ) ? nb_id + 1 : 0;
                    idx_t ncell     = src_neighbour_cells[nb_id];
                    idx_t nncell    = src_neighbour_cells[nnb_id];
                    const auto& Pn  = src_centroids_[ncell];
                    const auto& Pnn = src_centroids_[nncell];
                    if ( ncell != scell && nncell != scell ) {
                        double val = 0.5 * ( src_vals( ncell ) + src_vals( nncell ) ) - src_vals( scell );
                        auto csp   = CSPolygon( {Pn, Pnn, P} );
                        if ( csp.area() < std::numeric_limits<double>::epsilon() ) {
                            csp = CSPolygon( {Pn, P, Pnn} );
                        }
                        //auto orientation = ( csp.leftOf( Pnn, P, Pn ) ? -1 : 1 );
                        //ATLAS_ASSERT( orientation == -1 ); // orientation changes !
                        val *= ( csp.leftOf( Pnn, P, Pn, 1e-16, 0 ) ? -1 : 1 );
                        dual_area += std::abs( csp.area() );
                        grad = grad + PointXYZ::mul( PointXYZ::cross( Pn, Pnn ), val );
                    }
                    else if ( ncell != scell ) {
                        ATLAS_NOTIMPLEMENTED;
                        //double val = 0.5 * ( src_vals( ncell ) - src_vals( scell ) );
                        //grad = grad + PointXYZ::mul( PointXYZ::cross( Pn, P ), val );
                    }
                    else if ( nncell != scell ) {
                        ATLAS_NOTIMPLEMENTED;
                        //double val = 0.5 * ( src_vals( nncell ) - src_vals( scell ) );
                        //grad = grad + PointXYZ::mul( PointXYZ::cross( P, Pnn ), val );
                    }
                }
                if ( dual_area > std::numeric_limits<double>::epsilon() ) {
                    grad = PointXYZ::div( grad, dual_area );
                }
                for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
                    src_barycenter = src_barycenter + PointXYZ::mul( iparam.centroids[icell], iparam.weights[icell] );
                }
                src_barycenter = PointXYZ::div( src_barycenter, PointXYZ::norm( src_barycenter ) );
                grad           = grad - PointXYZ::mul( src_barycenter, PointXYZ::dot( grad, src_barycenter ) );
                ATLAS_ASSERT( std::abs( PointXYZ::dot( grad, src_barycenter ) ) < 1e-14 );
                for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
                    tgt_vals( iparam.tcell_id[icell] ) +=
                        iparam.weights[icell] *
                        ( src_vals( scell ) + PointXYZ::dot( grad, iparam.centroids[icell] - src_barycenter ) );
                }
            }
            for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
                tgt_vals( tcell ) /= tgt_areas_[tcell];
            }
        }
        else {
            Method::do_execute( src_field, tgt_field );
        }
    }
}

template <class TargetCellsIDs>
void ConservativeMethod::dump_intersection( const CSPolygon& s_csp, const std::vector<CSPolygon>& tgt_csp,
                                            const TargetCellsIDs& tgt_cells ) const {
    Log::info().flush();
    Log::info() << "\n === DEBUG ===\n\n";
    Log::info() << "* src cell: " << std::setprecision( 15 ) << s_csp << "\n";
    Log::info() << "* src area: " << s_csp.area() << "\n\n";
    double area_ncov = s_csp.area();
    for ( int i = 0; i < tgt_cells.size(); ++i ) {
        const auto tcell  = tgt_cells[i].payload();
        const auto& t_csp = tgt_csp[tcell];
        Log::info() << "* src cell: " << s_csp << "\n";
        Log::info() << "* tgt cell: " << t_csp << "\n";
        auto iplg          = s_csp.intersect( t_csp );
        auto jplg          = t_csp.intersect( s_csp );
        const double darea = std::abs( iplg.area() - jplg.area() );
        Log::info() << "* src ^ tgt      : " << iplg << "\n";
        Log::info() << "* src ^ tgt area : " << iplg.area() << "\n";
        if ( darea > 5e-8 ) {
            s_csp.intersect( t_csp, 1 );
            Log::info() << "* (!!) intersect comm area diff: " << darea << "\n";
            Log::info() << "* (!!) tgt ^ src      : " << jplg << "\n";
            Log::info() << "* (!!) tgt ^ src area : " << jplg.area() << "\n";
            t_csp.intersect( s_csp, 1 );
            ATLAS_ASSERT( false );
        }
        Log::info() << "\n";
        area_ncov -= iplg.area();
    }
    Log::info() << "\n=== END DEBUG ===\n\n";
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
