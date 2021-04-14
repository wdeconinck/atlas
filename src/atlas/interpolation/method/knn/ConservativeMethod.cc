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


#define DEBUG_OUTPUT_DETAIL 1

namespace atlas {
namespace interpolation {
namespace method {

using CSPolygon = util::ConvexSphericalPolygon;

ConservativeMethod::ConservativeMethod( const util::Config& config ) : Method( config ) {
    config.get( "order", order_ = 2 );
    config.get( "normalise_intersections", normalise_intersections_ = 1 );
    config.get( "field_value_type", fvtype_ = 0 );
    config.get( "matrix_free", matrix_free_ = true );
}

// get cell neighbours
std::vector<idx_t> ConservativeMethod::get_neighbours( Mesh& mesh, idx_t jcell ) const {
    const auto& cell2edge = mesh.cells().edge_connectivity();
    const auto& edge2cell = mesh.edges().cell_connectivity();
    const auto& edge2node = mesh.edges().node_connectivity();
    auto c2e_missval      = cell2edge.missing_value();
    auto valid_nb_cell    = [jcell, c2e_missval]( idx_t cid1, idx_t cid2 ) {
        if ( cid1 != c2e_missval && cid1 != jcell ) {
            return cid1;
        }
        else if ( cid2 != c2e_missval && cid2 != jcell ) {
            return cid2;
        }
        return -1;
    };
    std::vector<idx_t> nb_cells;
    idx_t n_edges = cell2edge.cols( jcell );
    nb_cells.reserve( n_edges );
    std::vector<bool> edge_done;
    edge_done.resize( n_edges );
    std::vector<idx_t> loc_edge_id( n_edges );
    for ( int ledge = 0; ledge < n_edges; ++ledge ) {
        loc_edge_id[ledge] = cell2edge( jcell, ledge );
    }
    idx_t ledge = 0;
    idx_t iedge = cell2edge( jcell, ledge );
    idx_t nbid  = valid_nb_cell( edge2cell( iedge, 0 ), edge2cell( iedge, 1 ) );
    if ( nbid != -1 ) {
        nb_cells.emplace_back( nbid );
    }
    edge_done[ledge] = true;
    idx_t nedge_done = 1;
    auto last_node   = edge2node( iedge, 1 );  // take any end point

    for ( ledge = 0; nedge_done < n_edges; ++ledge ) {
        if ( edge_done[ledge] ) {
            ledge = ( ledge == n_edges - 1 ? -1 : ledge );
            continue;
        }
        idx_t node0 = edge2node( cell2edge( jcell, ledge ), 0 );
        idx_t node1 = edge2node( cell2edge( jcell, ledge ), 1 );
        if ( last_node == node0 or last_node == node1 ) {
            nbid =
                valid_nb_cell( edge2cell( cell2edge( jcell, ledge ), 0 ), edge2cell( cell2edge( jcell, ledge ), 1 ) );
            if ( nbid != -1 ) {
                nb_cells.emplace_back( nbid );
            }
            last_node        = ( last_node == node0 ? node1 : node0 );
            edge_done[ledge] = true;
            ++nedge_done;
        }
        ledge = ( ledge == n_edges - 1 ? -1 : ledge );
    }
    return nb_cells;
}

std::vector<CSPolygon> ConservativeMethod::get_polygons( Mesh& mesh ) const {
    std::vector<CSPolygon> src_csp;
    if ( fvtype_ == 0 ) {  // CellColumns
        const idx_t n_cells = mesh.cells().size();
        src_csp.resize( n_cells );
        const auto& node_connectivity = mesh.cells().node_connectivity();
        const auto lonlat             = array::make_view<double, 2>( mesh.nodes().lonlat() );
        std::vector<PointLonLat> pts_ll;
        for ( idx_t icell = 0; icell < n_cells; ++icell ) {
            const idx_t n_nodes = node_connectivity.cols( icell );
            pts_ll.clear();
            pts_ll.resize( n_nodes );
            for ( idx_t jnode = 0; jnode < n_nodes; ++jnode ) {
                idx_t inode   = node_connectivity( icell, jnode );
                pts_ll[jnode] = PointLonLat{lonlat( inode, 0 ), lonlat( inode, 1 )};
            }
            src_csp[icell] = CSPolygon( pts_ll );
        }
    }
    else {  // NodeColumns
        if ( !mesh.cells().has_field( "centroids_xy" ) ) {
            mesh.cells().add(
                Field( "centroids_xy", mesh::actions::build_centroids_xy( mesh.cells(), mesh.nodes().xy() ) ) );
        }
        if ( !mesh.edges().has_field( "centroids_xy" ) ) {
            mesh.edges().add(
                Field( "centroids_xy", mesh::actions::build_centroids_xy( mesh.edges(), mesh.nodes().xy() ) ) );
        }
        /*
		auto xy             = array::make_view<double, 2>( nodes.xy() );
		auto cell_centroids = array::make_view<double, 2>( cells.field( "centroids_xy" ) );
		auto edge_centroids = array::make_view<double, 2>( edges.field( "centroids_xy" ) );
		const mesh::HybridElements::Connectivity& cell_edge_connectivity = cells.edge_connectivity();
		const mesh::HybridElements::Connectivity& edge_node_connectivity = edges.node_connectivity();
		auto field_flags                                                 = array::make_view<int, 1>( cells.flags() );

		auto patch = [&field_flags]( idx_t e ) {
        	using Topology = atlas::mesh::Nodes::Topology;
        	return Topology::check( field_flags( e ), Topology::PATCH );
    	};

		// special ordering for bit-identical results
		idx_t nb_cells = cells.size();
		std::vector<Node> ordering( nb_cells );
		for ( idx_t jcell = 0; jcell < nb_cells; ++jcell ) {
			ordering[jcell] =
				Node( util::unique_lonlat( cell_centroids( jcell, XX ), cell_centroids( jcell, YY ) ), jcell );
		}
		std::sort( ordering.data(), ordering.data() + nb_cells );
        const idx_t n_nodes = mesh.nodes().size();
        src_csp.resize( n_nodes );
        const auto& node_edge_connectivity = mesh.nodes().edge_connectivity();
        const auto lonlat             = array::make_view<double, 2>( mesh.nodes().lonlat() );
        std::vector<PointLonLat> pts_ll;
        for ( idx_t inode = 0; inode < n_nodes; ++inode ) {
            const idx_t n_nodes = node_connectivity.cols( inode );
            pts_ll.clear();
            pts_ll.resize( n_nodes );
            for ( idx_t jnode = 0; jnode < n_nodes; ++jnode ) {
                idx_t inode   = node_connectivity( icell, jnode );
                pts_ll[jnode] = PointLonLat{lonlat( inode, 0 ), lonlat( inode, 1 )};
            }
            src_csp[icell] = CSPolygon( pts_ll );
        }
*/
        ATLAS_ASSERT( false );
    }
    return src_csp;
}

void ConservativeMethod::do_setup( const Grid& src_grid, const Grid& tgt_grid ) {
    ATLAS_TRACE( "ConservativeMethod::do_setup()" );
    ATLAS_ASSERT( src_grid );
    ATLAS_ASSERT( tgt_grid );
    if ( mpi::size() > 1 ) {
        ATLAS_NOTIMPLEMENTED;
    }
    order2_setup_ = false;
    src_mesh_     = MeshGenerator( src_grid.meshgenerator() ).generate( src_grid );
    tgt_mesh_     = MeshGenerator( tgt_grid.meshgenerator() ).generate( tgt_grid );
    functionspace::CellColumns src_fs( src_mesh_ );
    functionspace::CellColumns tgt_fs( tgt_mesh_ );
    source_ = src_fs;
    target_ = tgt_fs;

    src_centroids_.resize( src_mesh_.cells().size() );
    tgt_centroids_.resize( tgt_mesh_.cells().size() );
    src_areas_.resize( src_mesh_.cells().size() );
    tgt_areas_.resize( tgt_mesh_.cells().size() );
    util::KDTree<idx_t> kdt_search;
    kdt_search.reserve( tgt_mesh_.cells().size() );

    const idx_t src_nb_cells = src_mesh_.cells().size();
    const auto& src_csp      = get_polygons( src_mesh_ );
    const idx_t tgt_nb_cells = tgt_mesh_.cells().size();
    const auto& tgt_csp      = get_polygons( tgt_mesh_ );

    double max_tgtcell_rad = 0.;
    for ( idx_t jcell = 0; jcell < tgt_nb_cells; ++jcell ) {
        kdt_search.insert( tgt_csp[jcell].centroid(), jcell );
        max_tgtcell_rad = std::max( max_tgtcell_rad, tgt_csp[jcell].cell_radius() );
    }
    kdt_search.build();

    size_t nonintersect        = 0;
    double src_area_notcovered = 0.;
    iparam_.resize( src_nb_cells );
    eckit::Channel blackhole;
    eckit::ProgressTimer progress( "Intersecting polygons ", src_nb_cells, " cell", double( 10 ),
                                   src_nb_cells > 50 ? Log::info() : blackhole );
    n_weights_ = 0;
    for ( idx_t scell = 0; scell < src_nb_cells; ++scell, ++progress ) {
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
                iparam_[scell].cell_id.emplace_back( tcell );
                iparam_[scell].weights.emplace_back( csp_i.area() );
                iparam_[scell].centroids.emplace_back( csp_i.centroid() );
                covered_area += csp_i.area();
                n_weights_++;
            }
        }
        const double loc_csp_error = std::abs( src_csp[scell].area() - covered_area );
        src_area_notcovered += loc_csp_error;
        if ( iparam_[scell].cell_id.size() == 0. ) {
            ++nonintersect;
        }
        // normalise
        if ( normalise_intersections_ ) {
            double wfactor = src_csp[scell].area() / ( covered_area > 1e-10 ? covered_area : 1. );
            for ( idx_t i = 0; i < iparam_[scell].weights.size(); i++ ) {
                iparam_[scell].weights[i] *= wfactor;
            }
        }
        if ( false && loc_csp_error > 1e-8 ) {
            Log::info() << "* src cell area NOT covered: " << loc_csp_error << "\n";
            dump_intersection( src_csp[scell], tgt_csp, tgt_cells );
            ATLAS_ASSERT( false );
        }
    }
    if ( not matrix_free_ ) {
        weights_.resize( n_weights_ );
        scell_id_.resize( n_weights_ );
        tcell_id_.resize( n_weights_ );
        int cnt = 0;
        for ( idx_t scell = 0; scell < src_nb_cells; ++scell ) {
            const auto& iparam = iparam_[scell];
            for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
                scell_id_[cnt] = scell;
                tcell_id_[cnt] = iparam.cell_id[icell];
                weights_[cnt]  = iparam.weights[icell] / tgt_csp[tcell_id_[cnt]].area();
                cnt++;
            }
        }
        ATLAS_TRACE( "ConservativeMethod::setup: build interpolant matrix" );
        Triplets triplets;
        triplets.resize( n_weights_ );
        for ( idx_t i = 0; i < n_weights_; i++ ) {
            triplets[i] = Triplet( tcell_id_[i], scell_id_[i], weights_[i] );
        }
        std::sort( std::begin( triplets ), std::end( triplets ), []( const Triplet& t1, const Triplet& t2 ) {
            return ( t1.row() < t2.row() or ( t1.row() == t2.row() and ( t1.col() > t2.col() ) ) );
        } );
        Matrix A( tgt_nb_cells, src_nb_cells, triplets );
        matrix_shared_->swap( A );
    }
    Log::info() << "WARNING " << nonintersect << " source mesh polygons do NOT intersect any other polygon.\n";
    Log::info() << "WARNING " << src_area_notcovered << " area of source mesh NOT covered by target mesh.\n";
    for ( idx_t tcell = 0; tcell < tgt_nb_cells; ++tcell ) {
        tgt_centroids_[tcell] = tgt_csp[tcell].centroid();
        tgt_areas_[tcell]     = tgt_csp[tcell].area();
    }
    if ( order_ > 1 ) {
        mesh::actions::build_halo( src_mesh_, 0 );
        mesh::actions::build_edges( src_mesh_, util::Config( "pole_edges", false ) );
    }
    do_setup_2nd_order();
}

void ConservativeMethod::do_setup_2nd_order() {
    if ( order_ != 2 || matrix_free_ ) {
        return;
    }
    //const auto halo     = array::make_view<int, 1>( src_mesh_.cells().halo() );
    const auto& src_cell2edge = src_mesh_.cells().edge_connectivity();
    const auto& src_edge2cell = src_mesh_.edges().cell_connectivity();
    const auto& src_edge2node = src_mesh_.edges().node_connectivity();
    const idx_t src_nb_cells  = src_mesh_.cells().size();
    neighbours_.resize( src_nb_cells );
    std::vector<double> grad_area;
    grad_area.resize( src_nb_cells );

    int cnt = 0;
    std::vector<std::array<idx_t, 3>> nb_idx_h;
    std::vector<int> val_h;
    for ( idx_t scell = 0; scell < src_nb_cells; ++scell ) {
        //	if ( halo( scell ) ) {
        //		continue;
        //	}
        const auto nb_cells = get_neighbours( src_mesh_, scell );
        neighbours_[scell].resize( nb_cells.size() );
        for ( int i = 0; i < nb_cells.size(); i++ ) {
            neighbours_[scell][i] = nb_cells[i];
        }
        grad_area[scell] = 0.;
        const auto& P    = src_centroids_[scell];
        for ( idx_t nb_id = 0; nb_id < nb_cells.size(); ++nb_id ) {
            idx_t nnb_id    = ( nb_id != nb_cells.size() - 1 ) ? nb_id + 1 : 0;
            idx_t ncell     = nb_cells[nb_id];
            idx_t nncell    = nb_cells[nnb_id];
            const auto& Pn  = src_centroids_[ncell];
            const auto& Pnn = src_centroids_[nncell];
            auto csp        = CSPolygon( {Pn, Pnn, P} );
            if ( csp.area() < std::numeric_limits<double>::epsilon() ) {
                csp = CSPolygon( {Pn, P, Pnn} );
            }
            if ( ncell != scell && nncell != scell ) {
                grad_area[scell] += std::abs( csp.area() );
            }
            cnt++;
            nb_idx_h.emplace_back( std::array<idx_t, 3>( {scell, ncell, nncell} ) );
            val_h.emplace_back( ( CSPolygon::leftOf( Pnn, P, Pn, 1e-16, 0 ) ? -1 : 1 ) );
        }
    }
    nb_idx_.resize( cnt );
    grad_nb_prod_.resize( cnt );
    cnt = 0;
    for ( idx_t scell = 0; scell < src_nb_cells; ++scell ) {
        const double grad_area_inv = 1. / ( grad_area[scell] > 0. ? grad_area[scell] : 1. );
        for ( idx_t nb_id = 0; nb_id < neighbours_[scell].size(); ++nb_id ) {
            nb_idx_[cnt][0]    = nb_idx_h[cnt][0];
            nb_idx_[cnt][1]    = nb_idx_h[cnt][1];
            nb_idx_[cnt][2]    = nb_idx_h[cnt][2];
            double fct         = 0.5 * double( val_h[cnt] ) * grad_area_inv;
            grad_nb_prod_[cnt] = PointXYZ::mul(
                PointXYZ::cross( src_centroids_[nb_idx_[cnt][1]], src_centroids_[nb_idx_[cnt][2]] ), fct );
            cnt++;
        }
    }
    src_bary_.resize( src_nb_cells );
    src_dbary_.resize( n_weights_ );
    cnt = 0;
    for ( idx_t scell = 0; scell < src_nb_cells; ++scell ) {
        const auto& iparam      = iparam_[scell];
        PointXYZ src_barycenter = {0., 0., 0.};
        for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
            src_barycenter = src_barycenter + PointXYZ::mul( iparam.centroids[icell], iparam.weights[icell] );
        }
        double b_norm = PointXYZ::norm( src_barycenter );
        if ( b_norm > 0. ) {
            src_barycenter = PointXYZ::div( src_barycenter, b_norm );
        }
        src_bary_[scell] = src_barycenter;
        for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
            src_dbary_[cnt++] = iparam.centroids[icell] - src_barycenter;
        }
    }
    order2_setup_ = true;
}

void ConservativeMethod::do_execute( const Field& src_field, Field& tgt_field ) {
    ATLAS_TRACE( "ConservativeMethod::do_execute()" );

    const auto src_vals = array::make_view<double, 1>( src_field );
    auto tgt_vals       = array::make_view<double, 1>( tgt_field );
    //const auto halo     = array::make_view<int, 1>( src_mesh_.cells().halo() );

    const auto& src_cell2edge = src_mesh_.cells().edge_connectivity();
    const auto& src_edge2cell = src_mesh_.edges().cell_connectivity();
    const auto& src_edge2node = src_mesh_.edges().node_connectivity();

    for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
        tgt_vals( tcell ) = 0.;
    }
    if ( order_ == 1 ) {
        if ( matrix_free_ ) {
            for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
                const auto& iparam = iparam_[scell];
                for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
                    tgt_vals( iparam.cell_id[icell] ) += iparam.weights[icell] * src_vals( scell );
                }
            }
            for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
                tgt_vals( tcell ) /= tgt_areas_[tcell];
            }
        }
        else {
            //            for ( idx_t i = 0; i < n_weights_; ++i ) {
            //              tgt_vals( tcell_id_[i] ) += weights_[i] * src_vals( scell_id_[i] );
            //        }
            Method::do_execute( src_field, tgt_field );
        }
    }
    else if ( order_ == 2 ) {
        if ( matrix_free_ ) {
            for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
                //if ( halo( scell ) ) {
                //    continue;
                //}
                const auto& iparam      = iparam_[scell];
                const PointXYZ& P       = src_centroids_[scell];
                PointXYZ grad           = {0., 0., 0.};
                PointXYZ src_barycenter = {0., 0., 0.};
                if ( order_ == 2 ) {
                    auto c2e_missval   = src_cell2edge.missing_value();
                    auto valid_nb_cell = [scell, c2e_missval]( idx_t cid1, idx_t cid2 ) {
                        if ( cid1 != c2e_missval && cid1 != scell ) {
                            return cid1;
                        }
                        else if ( cid2 != c2e_missval && cid2 != scell ) {
                            return cid2;
                        }
                        return -1;
                    };
                    // get cell neighbours
                    auto src_neighbour_cells = get_neighbours( src_mesh_, scell );
                    // calculate gradient
                    double dual_area = 0.;
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
                        src_barycenter =
                            src_barycenter + PointXYZ::mul( iparam.centroids[icell], iparam.weights[icell] );
                    }
                    src_barycenter = PointXYZ::div( src_barycenter, PointXYZ::norm( src_barycenter ) );
                    grad           = grad - PointXYZ::mul( src_barycenter, PointXYZ::dot( grad, src_barycenter ) );
                    ATLAS_ASSERT( std::abs( PointXYZ::dot( grad, src_barycenter ) ) < 1e-14 );
                    for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
                        tgt_vals( iparam.cell_id[icell] ) +=
                            iparam.weights[icell] *
                            ( src_vals( scell ) + PointXYZ::dot( grad, iparam.centroids[icell] - src_barycenter ) );
                    }
                }
            }
            for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
                tgt_vals( tcell ) /= tgt_areas_[tcell];
            }
        }
        else {
            if ( not order2_setup_ ) {
                Log::info() << " 2nd order method not set up. order2_setup_ == false\n";
                ATLAS_ASSERT( false );
            }
            std::vector<PointXYZ> grad;
            grad.resize( src_vals.size() );
            for ( idx_t scell = 0; scell < grad.size(); ++scell ) {
                grad[scell] = {0., 0., 0.};
            }
            for ( idx_t cnt = 0; cnt < nb_idx_.size(); ++cnt ) {
                const auto& nbidx = nb_idx_[cnt];
                const auto& g     = grad_nb_prod_[cnt];
                auto& gr          = grad[nbidx[0]];
                double val        = src_vals( nbidx[1] ) + src_vals( nbidx[2] ) - 2 * src_vals( nbidx[0] );
                gr[0] += g[0] * val;
                gr[1] += g[1] * val;
                gr[2] += g[2] * val;
            }
            for ( idx_t scell = 0; scell < grad.size(); ++scell ) {
                const PointXYZ& b = src_bary_[scell];
                PointXYZ& g       = grad[scell];
                const double gb   = g[0] * b[0] + g[1] * b[1] + g[2] * b[2];
                g[0] -= gb * b[0];
                g[1] -= gb * b[1];
                g[2] -= gb * b[2];
            }
            for ( idx_t i = 0; i < n_weights_; ++i ) {
                const idx_t scell = scell_id_[i];
                const auto& g     = grad[scell];
                const auto& db    = src_dbary_[i];
                tgt_vals( tcell_id_[i] ) +=
                    weights_[i] * ( src_vals( scell ) + g[0] * db[0] + g[1] * db[1] + g[2] * db[2] );
            }
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
