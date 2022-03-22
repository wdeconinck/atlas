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
#include "atlas/interpolation/Interpolation.h"
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

void do_remapping_test(Grid src_grid, Grid tgt_grid, double func(const PointLonLat&), RemapStat& remap_stat_1,
                       RemapStat& remap_stat_2) {
    Log::info().indent();
    // setup conservative remap: compute weights, polygon intersection, etc
    util::Config config("order", 1);
    ConservativeMethod consMethod(config);
    consMethod.setup(src_grid, tgt_grid);

    // create source field from analytic function "func"
    const auto& src_fs = consMethod.source();
    const auto& tgt_fs = consMethod.target();
    auto src_field     = src_fs.createField<double>();
    auto tgt_field     = tgt_fs.createField<double>();
    auto src_vals      = array::make_view<double, 1>(src_field);
    auto tgt_vals      = array::make_view<double, 1>(tgt_field);
    for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
        auto p = consMethod.src_points(spt);
        PointLonLat pll;
        eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
        src_vals(spt) = func(pll);
    }

    // project source field to target mesh in 1st order
    // consMethod.set_order(1);
    remap_stat_1 = RemapStat(consMethod.execute(src_field, tgt_field));
    remap_stat_1.compute(consMethod, src_vals, tgt_vals, nullptr, func);

    ATLAS_TRACE_SCOPE("test caching") {
        // We can create the interpolation without polygon intersections
        auto cache = consMethod.createCache();
        // cache = ConservativeMethod::Cache + MatrixCache (1st order)
        util::Config cfg(option::type("conservative"));
        {
            ATLAS_TRACE("cached -> 1st order using cached matrix");
            cfg.set("matrix_free", false);
            cfg.set("order", 1);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache);
            interpolation.execute(src_field, tgt_field);
            Log::info() << interpolation.source().type() << std::endl;
            Log::info() << interpolation.target().type() << std::endl;
        }
        {
            ATLAS_TRACE("cached -> 1st order constructing new matrix");
            cfg.set("matrix_free", false);
            cfg.set("order", 1);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, ConservativeMethod::Cache(cache));
            interpolation.execute(src_field, tgt_field);
        }
        {
            ATLAS_TRACE("cached -> 1st order matrix-free");
            cfg.set("matrix_free", true);
            cfg.set("order", 1);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache);
            interpolation.execute(src_field, tgt_field);
        }
        auto cache_2 = interpolation::Cache{};
        {
            ATLAS_TRACE("cached -> 2nd order constructing new matrix");
            cfg.set("matrix_free", false);
            cfg.set("order", 2);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, ConservativeMethod::Cache(cache));
            interpolation.execute(src_field, tgt_field);
            cache_2 = interpolation.createCache();
        }
        {
            ATLAS_TRACE("cached -> 2nd order matrix-free");
            cfg.set("matrix_free", true);
            cfg.set("order", 2);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache);
            interpolation.execute(src_field, tgt_field);
        }
        {
            ATLAS_TRACE("cached -> 2nd order using cached matrix");
            cfg.set("matrix_free", false);
            cfg.set("order", 2);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache_2);
            interpolation.execute(src_field, tgt_field);
        }
    }


    // TODO: We should not be allowed to change order like that. Rather a new interpolation object should
    //       be created of second order.
    {
        // project source field to target mesh in 2nd order
        consMethod.set_order(2);
        remap_stat_2 = RemapStat(consMethod.execute(src_field, tgt_field));
        remap_stat_2.compute(consMethod, src_vals, tgt_vals, nullptr, func);
    }
}

void check(const RemapStat remap_stat_1, RemapStat remap_stat_2, std::array<double, 6> tol) {
    auto improvement = [](double& e, double& r) { return 100. * (r - e) / r; };
    double err;
    // check polygon intersections
    err = remap_stat_1.errors[RemapStat::Errors::GEO_DIFF];
    Log::info() << "Polygon area computation improvement: " << improvement(err, tol[0]) << " %" << std::endl;
    EXPECT(err < tol[0]);
    err = remap_stat_1.errors[RemapStat::Errors::GEO_L1];
    Log::info() << "Polygon intersection improvement    : " << improvement(err, tol[1]) << " %" << std::endl;
    EXPECT(err < tol[1]);

    // check remap accuracy
    err = remap_stat_1.errors[RemapStat::Errors::REMAP_L2];
    Log::info() << "1st order accuracy improvement      : " << improvement(err, tol[2]) << " %" << std::endl;
    EXPECT(err < tol[2]);
    err = remap_stat_2.errors[RemapStat::Errors::REMAP_L2];
    Log::info() << "2nd order accuracy improvement      : " << improvement(err, tol[3]) << " %" << std::endl;
    EXPECT(err < tol[3]);

    // check mass conservation
    err = remap_stat_1.errors[RemapStat::Errors::REMAP_CONS];
    Log::info() << "1st order conservation improvement  : " << improvement(err, tol[4]) << " %" << std::endl;
    EXPECT(err < tol[4]);
    err = remap_stat_2.errors[RemapStat::Errors::REMAP_CONS];
    Log::info() << "2nd order conservation improvement  : " << improvement(err, tol[5]) << " %" << std::endl
                << std::endl;
    EXPECT(err < tol[5]);
    Log::info().unindent();
}

CASE("test_interpolation_conservative") {
#if 1
    SECTION("analytic constfunc") {
        auto func = [](const PointLonLat& p) { return 1.; };
        RemapStat remap_stat_1;
        RemapStat remap_stat_2;
        do_remapping_test(Grid("H47"), Grid("H48"), func, remap_stat_1, remap_stat_2);
        check(remap_stat_1, remap_stat_2, {1.e-13, 5.e-8, 2.9e-7, 2.9e-7, 5.5e-5, 5.5e-5});
    }

    SECTION("analytic Y_2^2 as in Jones(1998)") {
        auto func = [](const PointLonLat& p) {
            double cos = std::cos(0.025 * p[0]);
            return 2. + cos * cos * std::cos(2 * 0.025 * p[1]);
        };
        RemapStat remap_stat_1;
        RemapStat remap_stat_2;
        do_remapping_test(Grid("H47"), Grid("H48"), func, remap_stat_1, remap_stat_2);
        check(remap_stat_1, remap_stat_2, {1.e-13, 5.e-8, 4.8e-4, 1.1e-4, 8.9e-5, 1.1e-4});
    }
#endif
}

}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
