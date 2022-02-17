/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include <cmath>

#include "eckit/geometry/Sphere.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/array.h"
#include "atlas/array/MakeView.h"
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/interpolation/method/knn/ConservativeMethod.h"
#include "atlas/mesh.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

using CSPolygon          = util::ConvexSphericalPolygon;
using ConservativeMethod = interpolation::method::ConservativeMethod;
using RemapStat          = ConservativeMethod::RemapStat;
using FieldArray         = array::ArrayView<double, 1>;

void compute_field_errors(const FieldArray& src_vals, const FieldArray& tgt_vals, 
                          ConservativeMethod& consMethod, double func(const PointLonLat&)) {
    consMethod.remap_stat(src_vals, tgt_vals, nullptr, func);
}

void do_remapping_test(Grid src_grid, Grid tgt_grid, double func(const PointLonLat&)) {
    // setup conservative remap: compute weights, polygon intersection, etc
    util::Config config;
    ConservativeMethod consMethod(config);
    consMethod.do_setup(src_grid, tgt_grid);

    // get errors in polygon intersections
    double geo_create_err;
    consMethod.setup_stat();

    // create source field from analytic function "func"
    const auto& src_fs   = consMethod.source();
    const auto& tgt_fs   = consMethod.target();
    auto src_field       = src_fs.createField<double>();
    auto tgt_field       = tgt_fs.createField<double>();
    auto src_vals        = array::make_view<double, 1>(src_field);
    auto tgt_vals        = array::make_view<double, 1>(tgt_field);
    ATLAS_ASSERT(src_vals.size() == src_fs.size());
    for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
        auto p = consMethod.src_points(spt);
        PointLonLat pll;
        eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
        src_vals(spt) = func(pll);
    }
  
    // project source field to target mesh in 1st order
    consMethod.set_order(1);
    consMethod.do_execute(src_field, tgt_field);
    compute_field_errors(src_vals, tgt_vals, consMethod, func);

    // project source field to target mesh in 2nd order
    consMethod.set_order(2);
    consMethod.do_execute(src_field, tgt_field);
    compute_field_errors(src_vals, tgt_vals, consMethod, func);
}

CASE("test_interpolation_conservative") {
    SECTION("analytic constfunc") {
        auto func = [](const PointLonLat& p) { return 1.; };
        do_remapping_test(Grid("H47"), Grid("H48"), func);
    }

    SECTION("analytic Y_2^2 as in Jones - scaling") {
        auto func = [](const PointLonLat& p) {
            double cos = std::cos(0.025 * p[0]);
            return 2. + cos * cos * std::cos(2 * 0.025 * p[1]);
        };
        do_remapping_test(Grid("H47"), Grid("H48"), func);
    }

    SECTION("analytic Hill as in Jones") {
        auto func = [](const PointLonLat& p) {
            PointXYZ c = {1., 0., 0.};
            PointXYZ p_sph;
            eckit::geometry::Sphere::convertSphericalToCartesian(1., p, p_sph);
            double r = PointXYZ::norm(p_sph - c);
            return 2. + std::cos(M_PI * r / 10.);
        };
        do_remapping_test(Grid("H47"), Grid("H48"), func);
    }
}

}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
