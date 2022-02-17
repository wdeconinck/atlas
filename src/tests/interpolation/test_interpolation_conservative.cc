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

void do_remapping_test(Grid src_grid, Grid tgt_grid, double func(const PointLonLat&),
                       RemapStat& remap_stat_1, RemapStat& remap_stat_2) {
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
    consMethod.remap_stat(src_vals, tgt_vals, nullptr, func);
    remap_stat_1 = consMethod.remap_stat();

    // project source field to target mesh in 2nd order
    consMethod.set_order(2);
    consMethod.do_execute(src_field, tgt_field);
    consMethod.remap_stat(src_vals, tgt_vals, nullptr, func);
    remap_stat_2 = consMethod.remap_stat();
}

void check(const RemapStat remap_stat_1, RemapStat remap_stat_2, std::array<double,6> tol) {
    auto improvement = [](double& e, double& r){ return (r-e)/r; };
    double err;
    // check polygon intersections
    err = remap_stat_1.errors[RemapStat::Errors::GEO_DIFF];
    Log::info() << "Polygon area computation improvement: "
                << improvement(err, tol[0]) << " %" << std::endl;
    EXPECT(err < tol[0]);
    err = remap_stat_1.errors[RemapStat::Errors::GEO_L1];
    Log::info() << "Polygon intersection improvement: "
                << improvement(err, tol[1]) << " %" << std::endl;
    EXPECT(err < tol[1]);

    // check remap accuracy
    err = remap_stat_1.errors[RemapStat::Errors::REMAP_L2];
    Log::info() << "1st order accuracy improvement: "
                << improvement(err, tol[2]) << " %" << std::endl;
    EXPECT(err < tol[2]);
    err = remap_stat_2.errors[RemapStat::Errors::REMAP_L2];
    Log::info() << "2nd order accuracy improvement: "
                << improvement(err, tol[3]) << " %" << std::endl;
    EXPECT(err < tol[3]);

    // check mass conservation
    err = remap_stat_1.errors[RemapStat::Errors::REMAP_CONS];
    Log::info() << "1st order conservation improvement: "
                << improvement(err, tol[4]) << " %" << std::endl;
    EXPECT(err < tol[4]);
    err = remap_stat_2.errors[RemapStat::Errors::REMAP_CONS];
    Log::info() << "2nd order conservation improvement: "
                << improvement(err, tol[5]) << " %" << std::endl;
Log::info() << err << " " << tol[5] <<std::endl;
    EXPECT(err < tol[5]);
}

CASE("test_interpolation_conservative") {
    SECTION("analytic constfunc") {
        auto func = [](const PointLonLat& p) { return 1.; };
        RemapStat remap_stat_1;
        RemapStat remap_stat_2;
        do_remapping_test(Grid("H47"), Grid("H48"), func, remap_stat_1, remap_stat_2);
        check(remap_stat_1, remap_stat_2, {1.e-13, 5.e-8, 2.5e-7, 2.5e-7, 1.5e-6, 1.5e-6});
    }

    SECTION("analytic Y_2^2 as in Jones - scaling") {
        auto func = [](const PointLonLat& p) {
            double cos = std::cos(0.025 * p[0]);
            return 2. + cos * cos * std::cos(2 * 0.025 * p[1]);
        };
        RemapStat remap_stat_1;
        RemapStat remap_stat_2;
        do_remapping_test(Grid("H47"), Grid("H48"), func, remap_stat_1, remap_stat_2);
        check(remap_stat_1, remap_stat_2, {1.e-13, 5.e-8, 4.8e-4, 1.1e-4, 5.8e-6, 5.3e-5});
    }
}

}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
