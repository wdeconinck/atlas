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
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildHalo.h"
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

ConservativeMethod::ConservativeMethod( const util::Config& config ) {
    config.get( "order", order_ = 2 );
    config.get( "normalise_intersections", normalise_intersections_ = 1 );
}

void ConservativeMethod::do_setup( Mesh& src_mesh, const Mesh& tgt_mesh ) {
    ATLAS_TRACE( "ConservativeMethod::do_setup()" );

    if ( mpi::size() > 1 ) {
        ATLAS_NOTIMPLEMENTED;
    }

    src_centroids_.resize( src_mesh.cells().size() );
    tgt_centroids_.resize( tgt_mesh.cells().size() );
    src_areas_.resize( src_mesh.cells().size() );
    tgt_areas_.resize( tgt_mesh.cells().size() );
    util::KDTree<idx_t> kdt_search;
    kdt_search.reserve( tgt_mesh.cells().size() );
    double max_tgtcell_rad = 0.;

    ATLAS_ASSERT( src_mesh );
    ATLAS_ASSERT( tgt_mesh );

    std::vector<PointLonLat> pts_ll;
    const idx_t src_nb_cells          = src_mesh.cells().size();
    const auto& src_node_connectivity = src_mesh.cells().node_connectivity();
    std::vector<CSPolygon> src_csp;
    src_csp.resize( src_nb_cells );
    const auto src_lonlat = array::make_view<double, 2>( src_mesh.nodes().lonlat() );
    for ( idx_t jcell = 0; jcell < src_nb_cells; ++jcell ) {
        const idx_t nb_nodes = src_node_connectivity.cols( jcell );
        pts_ll.clear();
        pts_ll.resize( nb_nodes );
        for ( idx_t jnode = 0; jnode < nb_nodes; ++jnode ) {
            idx_t inode   = src_node_connectivity( jcell, jnode );
            pts_ll[jnode] = PointLonLat{src_lonlat( inode, 0 ), src_lonlat( inode, 1 )};
        }
        src_csp[jcell] = CSPolygon( pts_ll );
    }

    const idx_t tgt_nb_cells          = tgt_mesh.cells().size();
    const auto& tgt_node_connectivity = tgt_mesh.cells().node_connectivity();
    std::vector<CSPolygon> tgt_csp;
    tgt_csp.resize( tgt_nb_cells );
    const auto tgt_lonlat = array::make_view<double, 2>( tgt_mesh.nodes().lonlat() );
    for ( idx_t jcell = 0; jcell < tgt_nb_cells; ++jcell ) {
        const idx_t nb_nodes = tgt_node_connectivity.cols( jcell );
        pts_ll.clear();
        pts_ll.resize( nb_nodes );
        for ( idx_t jnode = 0; jnode < nb_nodes; ++jnode ) {
            idx_t inode   = tgt_node_connectivity( jcell, jnode );
            pts_ll[jnode] = PointLonLat{tgt_lonlat( inode, 0 ), tgt_lonlat( inode, 1 )};
        }
        tgt_csp[jcell] = CSPolygon( pts_ll );
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

        //#ifndef NDEBUG
        if ( false && loc_csp_error > 1e-5 ) {
            Log::info().flush();
            Log::info() << "\n === DEBUG ===\n\n";
            Log::info() << "* src cell area NOT covered: " << loc_csp_error << "\n";
            Log::info() << "* src cell: " << std::setprecision( 30 ) << s_csp << "\n";
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
                if ( darea > 1e-9 or iplg.area() > 0.1 ) {
                    //ATLAS_ASSERT( false );
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
            ATLAS_ASSERT( false );
        }
        //#endif
    }
    Log::info() << "WARNING " << nonintersect << " source mesh polygons do NOT intersect any other polygon.\n";
    Log::info() << "WARNING " << src_area_notcovered << " area of source mesh NOT covered by target mesh.\n";

    for ( idx_t tcell = 0; tcell < tgt_nb_cells; ++tcell ) {
        tgt_centroids_[tcell] = tgt_csp[tcell].centroid();
        tgt_areas_[tcell]     = tgt_csp[tcell].area();
    }

    if ( order_ > 1 ) {
        mesh::actions::build_halo( src_mesh, 0 );
        mesh::actions::build_edges( src_mesh, util::Config( "pole_edges", false ) );
    }
    src_mesh_ = src_mesh;
}

void ConservativeMethod::do_execute( const Field& src_field, Field& tgt_field ) const {
    ATLAS_TRACE( "ConservativeMethod::do_execute()" );

    auto src_vals = array::make_view<double, 1>( src_field );
    auto tgt_vals = array::make_view<double, 1>( tgt_field );
    auto halo     = array::make_view<int, 1>( src_mesh_.cells().halo() );

    const auto& src_cell2edge = src_mesh_.cells().edge_connectivity();
    const auto& src_edge2cell = src_mesh_.edges().cell_connectivity();
    const auto& src_edge2node = src_mesh_.edges().node_connectivity();

    for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
        tgt_vals( tcell ) = 0.;
    }
    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        //if ( halo( scell ) ) {
        //    continue;
        //}
        const auto& iparam      = iparam_[scell];
        const PointXYZ& P       = src_centroids_[scell];
        PointXYZ grad           = {0., 0., 0.};
        PointXYZ src_barycenter = {0., 0., 0.};
        if ( order_ > 1 ) {
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
            idx_t src_nb_edges = src_cell2edge.cols( scell );
            std::vector<idx_t> src_neighbour_cells;
            src_neighbour_cells.reserve( src_nb_edges );
            std::vector<bool> edge_done;
            edge_done.resize( src_nb_edges );
            std::vector<idx_t> loc_edge_id( src_nb_edges );
            for ( int ledge = 0; ledge < src_nb_edges; ++ledge ) {
                loc_edge_id[ledge] = src_cell2edge( scell, ledge );
            }
            idx_t ledge = 0;
            idx_t iedge = src_cell2edge( scell, ledge );
            idx_t nbid  = valid_nb_cell( src_edge2cell( iedge, 0 ), src_edge2cell( iedge, 1 ) );
            if ( nbid != -1 ) {
                src_neighbour_cells.emplace_back( nbid );
            }
            edge_done[ledge] = true;
            idx_t nedge_done = 1;
            auto last_node   = src_edge2node( iedge, 1 );  // take any end point

            for ( ledge = 0; nedge_done < src_nb_edges; ++ledge ) {
                if ( edge_done[ledge] ) {
                    ledge = ( ledge == src_nb_edges - 1 ? -1 : ledge );
                    continue;
                }
                idx_t node0 = src_edge2node( src_cell2edge( scell, ledge ), 0 );
                idx_t node1 = src_edge2node( src_cell2edge( scell, ledge ), 1 );
                if ( last_node == node0 or last_node == node1 ) {
                    nbid = valid_nb_cell( src_edge2cell( src_cell2edge( scell, ledge ), 0 ),
                                          src_edge2cell( src_cell2edge( scell, ledge ), 1 ) );
                    if ( nbid != -1 ) {
                        src_neighbour_cells.emplace_back( nbid );
                    }
                    last_node        = ( last_node == node0 ? node1 : node0 );
                    edge_done[ledge] = true;
                    ++nedge_done;
                }
                ledge = ( ledge == src_nb_edges - 1 ? -1 : ledge );
            }

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
                    //auto orientation = ( csp.leftOf( Pnn, P, Pn ) ? -1 : 1 );
                    //ATLAS_ASSERT( orientation == -1 ); // orientation changes !
                    val *= ( csp.leftOf( Pnn, P, Pn ) ? -1 : 1 );
                    dual_area += csp.area();
                    grad = grad + PointXYZ::mul( PointXYZ::cross( Pn, Pnn ), val );
                }
                else if ( ncell != scell ) {
                    double val = 0.5 * ( src_vals( ncell ) - src_vals( scell ) );
                    val *= -1;
                    //grad = grad + PointXYZ::mul( PointXYZ::cross( Pn, P ), val );
                }
                else if ( nncell != scell ) {
                    double val = 0.5 * ( src_vals( nncell ) - src_vals( scell ) );
                    val *= -1;
                    //grad = grad + PointXYZ::mul( PointXYZ::cross( P, Pnn ), val );
                }
            }
            grad = PointXYZ::div( grad, ( dual_area > 0. ? dual_area : 1. ) );
            for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
                src_barycenter = src_barycenter + PointXYZ::mul( iparam.centroids[icell], iparam.weights[icell] );
            }
            const double src_brc_norm = PointXYZ::norm( src_barycenter );
            //ATLAS_ASSERT( src_brc_norm > 1e-7 );
            if ( src_brc_norm < 1e-7 ) {
                src_barycenter = src_centroids_[scell];
            }
            else {
                src_barycenter = PointXYZ::div( src_barycenter, src_brc_norm );
            }
            grad = grad - PointXYZ::mul( src_barycenter, PointXYZ::dot( grad, src_barycenter ) );
        }
        //ATLAS_ASSERT( std::abs( PointXYZ::dot(grad, src_barycenter) ) < 1e-14 );
        for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
            tgt_vals( iparam.cell_id[icell] ) +=
                iparam.weights[icell] *
                ( src_vals( scell ) + PointXYZ::dot( grad, iparam.centroids[icell] - src_barycenter ) );
        }
    }
    for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
        tgt_vals( tcell ) /= tgt_areas_[tcell];
    }
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
