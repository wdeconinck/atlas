/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

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

#define DEBUG_OUTPUT_DETAIL 0

namespace atlas {
namespace interpolation {
namespace method {

using CSPolygon = util::ConvexSphericalPolygon;

ConservativeMethod::ConservativeMethod( const util::Config& config ) {
    config.get( "order", order_ = 2 );
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
    }

    // brute force (!) needs to be changed
    iparam_.resize( src_nb_cells );
    eckit::Channel blackhole;
    eckit::ProgressTimer progress( "Intersecting polygons ", src_nb_cells, " cell", double( 10 ),
                                   src_nb_cells > 50 ? Log::info() : blackhole );
    for ( idx_t scell = 0; scell < src_nb_cells; ++scell, ++progress ) {
        src_centroids_[scell] = src_csp[scell].centroid();
        src_areas_[scell]     = src_csp[scell].area();
        for ( idx_t tcell = 0; tcell < tgt_nb_cells; ++tcell ) {
            CSPolygon csp_i = src_csp[scell].intersect( tgt_csp[tcell] );
            if ( csp_i.area() > 0. ) {
                iparam_[scell].cell_id.emplace_back( tcell );
                iparam_[scell].weights.emplace_back( csp_i.area() );
                iparam_[scell].centroids.emplace_back( csp_i.centroid() );
            }
        }
    }

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

    for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
        tgt_vals( tcell ) = 0.;
    }
    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        if ( halo( scell ) ) {
            continue;
        }
        const auto& iparam      = iparam_[scell];
        const PointXYZ& P       = src_centroids_[scell];
        PointXYZ grad           = {0., 0., 0.};
        PointXYZ src_barycenter = {0., 0., 0.};
        if ( order_ > 1 ) {
            // get cell neighbours
            idx_t src_nb_edges = src_cell2edge.cols( scell );
            std::vector<idx_t> src_neighbour_cells;
            src_neighbour_cells.reserve( src_nb_edges );
            for ( idx_t sedge = 0; sedge < src_nb_edges; ++sedge ) {
                idx_t iedge = src_cell2edge( scell, sedge );
                ATLAS_ASSERT( iedge < src_mesh_.edges().size() );
                ATLAS_ASSERT( iedge < src_edge2cell.rows() );
                idx_t cell0 = src_edge2cell( iedge, 0 );
                idx_t cell1 = src_edge2cell( iedge, 1 );
#if DEBUG_OUTPUT_DETAIL
                Log::info() << " mv, cell0, cell1, scell: " << src_cell2edge.missing_value() << " " << cell0 << " "
                            << cell1 << " " << scell;
                Log::info().flush();
#endif
                if ( cell0 != src_cell2edge.missing_value() && cell0 != scell ) {
                    src_neighbour_cells.emplace_back( cell0 );
#if DEBUG_OUTPUT_DETAIL
                    Log::info() << ", add " << cell0 << "\n";
                    Log::info().flush();
#endif
                }
                else if ( cell1 != src_cell2edge.missing_value() && cell1 != scell ) {
                    src_neighbour_cells.emplace_back( cell1 );
#if DEBUG_OUTPUT_DETAIL
                    Log::info() << ", add " << cell1 << "\n";
                    Log::info().flush();
#endif
                }
                else {
                    // even for global meshes we come here, why?
                    src_neighbour_cells.emplace_back( scell );
#if DEBUG_OUTPUT_DETAIL
                    Log::info() << ", add " << scell << "\n";
                    Log::info().flush();
#endif
                }
            }
#if DEBUG_OUTPUT_DETAIL
            Log::info() << "\n";
            Log::info().flush();
#endif
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
            ATLAS_ASSERT( PointXYZ::norm( src_barycenter ) > 1e-7 );
            src_barycenter = PointXYZ::div( src_barycenter, PointXYZ::norm( src_barycenter ) );
            grad           = grad - PointXYZ::mul( src_barycenter, PointXYZ::dot( grad, src_barycenter ) );
        }
        //ATLAS_ASSERT( std::abs( PointXYZ::dot(grad, src_barycenter) ) < 1e-14 );
        for ( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
            tgt_vals( iparam.cell_id[icell] ) +=
                iparam.weights[icell] *
                ( src_vals( scell ) + PointXYZ::dot( grad, iparam.centroids[icell] - src_barycenter ) );
            //( src_vals( scell ) + PointXYZ::dot( grad, iparam.centroids[icell] - src_centroids_[scell] ) );
        }
    }
    for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
        tgt_vals( tcell ) /= tgt_areas_[tcell];
    }
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
