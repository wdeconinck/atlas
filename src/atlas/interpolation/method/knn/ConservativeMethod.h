/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#pragma once

#include "atlas/interpolation/method/knn/KNearestNeighboursBase.h"

#include <forward_list>

#include "atlas/functionspace.h"
#include "atlas/util/Config.h"
#include "atlas/util/ConvexSphericalPolygon.h"


namespace atlas {
namespace interpolation {
namespace method {


class ConservativeMethod : public Method {
public:
    typedef util::ConvexSphericalPolygon CSPolygon;
    struct InterpolationParameters {
        std::vector<idx_t> tcell_id;
        std::vector<PointXYZ> centroids;
        std::vector<double> weights;
        std::vector<double> sweights;
    };

    ConservativeMethod( const util::Config& = util::NoConfig() );

    using Method::do_setup;
    void do_setup( const FunctionSpace& source, const FunctionSpace& target ) { ATLAS_NOTIMPLEMENTED; }
    void do_setup( const Grid& source, const Grid& target );
    void do_setup( const Grid& source, const Grid& target, const Cache& ) { ATLAS_NOTIMPLEMENTED; }

    void do_execute( const Field& src_field, Field& tgt_field );

    void print( std::ostream& out ) const { out << "ConservativeMethod[]"; }
    const FunctionSpace& source() const { return source_; }
    const FunctionSpace& target() const { return target_; }

    inline const std::vector<InterpolationParameters>& iparam() const { return iparam_; }
    inline const PointXYZ& src_centroid( size_t id ) const { return src_centroids_[id]; }
    inline const PointXYZ& tgt_centroid( size_t id ) const { return tgt_centroids_[id]; }
    inline const double& src_area( size_t id ) const { return src_areas_[id]; }
    inline const double& tgt_area( size_t id ) const { return tgt_areas_[id]; }
    void set_order( int order ) {
        order_ = order;
        setup_1st_order_matrix();
        setup_2nd_order_matrix();
    }
    int order() const { return order_; }
    Mesh src_mesh() const { return src_mesh_; }
    Mesh tgt_mesh() const { return tgt_mesh_; }
    void setup_1st_order_matrix();
    void setup_2nd_order_matrix();

private:
    template <class TargetCellsIDs>
    void dump_intersection( const util::ConvexSphericalPolygon& s_csp,
                            const std::vector<util::ConvexSphericalPolygon>& tgt_csp,
                            const TargetCellsIDs& tgt_cells ) const;

    std::vector<CSPolygon> get_polygons( Mesh& mesh ) const;
    std::vector<idx_t> get_cell_neighbours( Mesh& mesh, idx_t jcell ) const;

protected:
    FunctionSpace source_;
    FunctionSpace target_;
    Mesh src_mesh_;
    Mesh tgt_mesh_;
    int normalise_intersections_;
    int order_;
    int fvtype_;
    idx_t n_scells_;
    idx_t n_tcells_;
    bool matrix_free_;
    std::vector<PointXYZ> src_centroids_;
    std::vector<PointXYZ> tgt_centroids_;
    std::vector<double> src_areas_;
    std::vector<double> tgt_areas_;
    std::vector<InterpolationParameters> iparam_;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
