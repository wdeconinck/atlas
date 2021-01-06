/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */


#include "atlas/interpolation/method/knn/ConservativeMethod.h"

#include <algorithm>
#include <vector>

//#include "eckit/log/Plural.h"
//#include "eckit/log/ProgressTimer.h"
//#include "eckit/types/FloatCompare.h"

#include "atlas/grid.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

#include "eckit/log/ProgressTimer.h"

namespace atlas {
namespace interpolation {
namespace method {

using CSPolygon = util::ConvexSphericalPolygon;

ConservativeMethod::ConservativeMethod( const util::Config& config ) {
    config.get( "order", order_ = 2 );
}

//void ConservativeMethod::print( std::ostream& out ) const {
//    out << "ConservativeMethod[]";
//}


//void ConservativeMethod::do_setup( const FunctionSpace& /*source*/, const FunctionSpace& /*target*/ ) {
//    ATLAS_NOTIMPLEMENTED;
//}

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
        mesh::actions::build_edges( src_mesh );  // needed for gradient
    }
    src_mesh_ = src_mesh;
}

void ConservativeMethod::do_execute( const Field& src_field, Field& tgt_field ) const {
    ATLAS_TRACE( "ConservativeMethod::do_execute()" );

    auto src_vals = array::make_view<double, 1>( src_field );
    auto tgt_vals = array::make_view<double, 1>( tgt_field );

    const auto& src_cell2edge = src_mesh_.cells().edge_connectivity();
    const auto& src_edge2cell = src_mesh_.edges().cell_connectivity();

    for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
        tgt_vals( tcell ) = 0.;
    }
    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        PointXYZ grad = {0., 0., 0.};
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
                if ( cell0 != src_cell2edge.missing_value() && cell0 != scell ) {
                    src_neighbour_cells.emplace_back( cell0 );
                }
                else if ( cell1 != src_cell2edge.missing_value() ) {
                    src_neighbour_cells.emplace_back( cell1 );
                }
                else {
                    src_neighbour_cells.emplace_back( scell );
                }
            }
            // calculate gradient
            double dual_area = 0.;
            for ( idx_t nid = 0; nid < src_neighbour_cells.size(); ++nid ) {
                idx_t ncell  = src_neighbour_cells[nid];
                idx_t nncell = src_neighbour_cells[nid != src_neighbour_cells.size() - 1 ? nid + 1 : 0];
                if ( ncell != scell && nncell != scell ) {
                    double coeff = 0.5 * ( src_vals( ncell ) + src_vals( nncell ) ) - src_vals( scell );
                    dual_area +=
                        CSPolygon( {src_centroids_[ncell], src_centroids_[nncell], src_centroids_[scell]} ).area();
                    grad = PointXYZ::add(
                        grad,
                        PointXYZ::mul( PointXYZ::cross( src_centroids_[ncell], src_centroids_[nncell] ), coeff ) );
                }
            }
            grad = PointXYZ::mul( grad, ( dual_area > 0. ? 1. / dual_area : 1. ) );
        }
        for ( idx_t icell = 0; icell < iparam_[scell].centroids.size(); ++icell ) {
            // NOTE: check if this is correct barycenter of the source cell!!
            tgt_vals( iparam_[scell].cell_id[icell] ) +=
                iparam_[scell].weights[icell] *
                ( src_vals( scell ) + PointXYZ::dot( grad, iparam_[scell].centroids[icell] - src_centroids_[scell] ) );
        }
    }
    for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
        tgt_vals( tcell ) /= tgt_areas_[tcell];
    }
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
