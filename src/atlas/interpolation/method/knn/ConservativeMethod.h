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


class ConservativeMethod {
public:
    static constexpr double tol = 1e-14;

    ConservativeMethod( const util::Config& = util::NoConfig() );

    struct InterpolationParameters {
        std::vector<int> cell_id;
        std::vector<PointXYZ> centroids;
        std::vector<double> weights;
        std::ostream& print( std::ostream& os ) {
            os << "centroids         : " << centroids << "\n"
               << "area              : " << weights << "\n";
			return os;
        }
    };

    /**
     * @brief setup ConvexSphericalPolygons for two grids
     * @param source functionspace containing source points
     * @param target functionspace containing target points
     */
    //virtual void do_setup( const FunctionSpace& source, const FunctionSpace& target ) override;
    void do_setup( Mesh& source, const Mesh& target );

    //void do_execute( const FieldSet& source, FieldSet& target ) const override;
    void do_execute( const Field& src_field, Field& tgt_field ) const;

    inline const std::vector<InterpolationParameters>& iparam() const { return iparam_; }
    inline const PointXYZ& src_centroid( size_t id ) const { return src_centroids_[id]; }
    inline const PointXYZ& tgt_centroid( size_t id ) const { return tgt_centroids_[id]; }
    inline const double& src_area( size_t id ) const { return src_areas_[id]; }
    inline const double& tgt_area( size_t id ) const { return tgt_areas_[id]; }

    void set_order( int order ) { order_ = order; }
    const int order() const { return order_; }

protected:
    Mesh src_mesh_;
    std::vector<PointXYZ> src_centroids_;
    std::vector<PointXYZ> tgt_centroids_;
    std::vector<double> src_areas_;
    std::vector<double> tgt_areas_;
    int order_;
    std::vector<InterpolationParameters> iparam_;
    int normalise_intersections_;
};


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
