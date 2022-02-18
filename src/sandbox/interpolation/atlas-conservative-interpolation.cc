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
#include <fstream>

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

Grid localgrid(int nx, int ny) {
    util::Config gridspec;
    gridspec.set("type", "regional");
    gridspec.set("nx", nx);
    gridspec.set("ny", ny);
    gridspec.set("north", 80);
    gridspec.set("south", 0);
    gridspec.set("west", 0);
    gridspec.set("east", 90);
    return Grid{gridspec};
}

void print_remap_errors(ConservativeMethod& consMethod, std::ofstream& outfile) {
    const auto& remap_stat = consMethod.remap_stat();
    Log::info() << "    " << consMethod.order() << "-order remap analytical error : (L2) " 
                << remap_stat.errors[RemapStat::Errors::REMAP_L2] << " (Lmax) "
                << remap_stat.errors[RemapStat::Errors::REMAP_LINF] << "\n";
    Log::info() << "    " << consMethod.order() << "-order global mass conservation error : " 
                << remap_stat.errors[RemapStat::Errors::REMAP_CONS] << "\n";
    outfile << std::setw(10) 
            << remap_stat.errors[RemapStat::Errors::REMAP_L2] 
            << std::setw(10)
            << remap_stat.errors[RemapStat::Errors::REMAP_LINF] 
            << std::setw(10) 
            << remap_stat.errors[RemapStat::Errors::REMAP_CONS];
}

void do_remapping_test(Grid src_grid, Grid tgt_grid, double func(const PointLonLat&), std::ofstream& outfile) {
    util::Config gmsh_config;
    // Allow command-line argument to change coordinates output to lonlat; e.g.
    //    <program> --coordinates lonlat
    gmsh_config.set("coordinates", eckit::Resource<std::string>("--coordinates","lonlat"));
    gmsh_config.set("ghost", true);
    util::Config config;
    auto cell_data = [](const std::string& resource, const Grid& grid) -> bool {
        bool resource_default = grid.name()[0] == 'H' ? true : false;
        bool retval = eckit::Resource<bool>(resource,resource_default);
        return retval;
    };
    config.set("src_cell_data", cell_data("--src-cell-data",src_grid));
    config.set("tgt_cell_data", cell_data("--tgt-cell-data",tgt_grid));
    config.set("matrix_free", eckit::Resource<bool>("--matrix-free",false));
    config.set("normalise_intersections", eckit::Resource<bool>("--normalise",true));
    ConservativeMethod consMethod(config);

    outfile << std::setw(10) << src_grid.name() << std::setw(10) << tgt_grid.name();

    Log::info() << "REMAPPING: " << src_grid.name() << " --> " << tgt_grid.name() 
            << ", matrix-free: " << consMethod.matrix_free()
            << ", normalise: " << consMethod.normalise_intersections() << std::endl;
    Log::info().indent();
    auto start = std::chrono::system_clock::now();
    consMethod.do_setup(src_grid, tgt_grid);
    std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start;
    Log::info() << "Setup (computing supermesh) took " << elapsed_seconds.count() << " seconds.\n";
    outfile << std::setw(10) << elapsed_seconds.count();

    const auto& src_mesh = consMethod.src_mesh();
    const auto& tgt_mesh = consMethod.tgt_mesh();
    const auto& src_fs   = consMethod.source();
    const auto& tgt_fs   = consMethod.target();
    auto src_field       = src_fs.createField<double>();
    auto tgt_field       = tgt_fs.createField<double>();
    auto src_vals        = array::make_view<double, 1>(src_field);
    auto tgt_vals        = array::make_view<double, 1>(tgt_field);

    consMethod.setup_stat();
    output::Gmsh("cons-remap_srcmesh.msh", gmsh_config).write(consMethod.src_mesh());
    output::Gmsh("cons-remap_tgtmesh.msh", gmsh_config).write(consMethod.tgt_mesh());

    for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
        PointLonLat pll;
        eckit::geometry::Sphere::convertCartesianToSpherical(1., consMethod.src_points(spt), pll);
        src_vals(spt) = func(pll);
    }
    output::Gmsh("cons-remap_srcfield.msh", gmsh_config).write(src_field);

    consMethod.set_order(1);
    start = std::chrono::system_clock::now();
    consMethod.do_execute(src_field, tgt_field);
    elapsed_seconds = std::chrono::system_clock::now() - start;

    // compute difference field
    auto diff_field = src_fs.createField<double>();
    auto diff_vals  = array::make_view<double, 1>(diff_field);
    consMethod.remap_stat(src_vals, tgt_vals, &diff_vals, func);

    // remap statistics
    auto& remap_stat = consMethod.remap_stat();
    Log::info() << "Created " << remap_stat.counts[RemapStat::Counts::SRC_PLG] 
                << " (sub)polygons from "
                << src_mesh.cells().size() << " source mesh cells.\n";
    Log::info() << "    Total sum of subpolygon over/undershoots : "
                << remap_stat.errors[RemapStat::Errors::SRC_PLG_L1] << "\n";
    Log::info() << "    Max over/undershoots per cell : "
                << remap_stat.errors[RemapStat::Errors::SRC_PLG_LINF] << "\n";
    Log::info() << "Created " << remap_stat.counts[RemapStat::Counts::TGT_PLG] 
                << " (sub)polygons from "
                << tgt_mesh.cells().size() << " target mesh cells.\n";
    Log::info() << "    Total sum of subpolygon over/undershoots : "
                << remap_stat.errors[RemapStat::Errors::TGT_PLG_L1] << "\n";
    Log::info() << "    Max over/undershoots per cell : "
                << remap_stat.errors[RemapStat::Errors::TGT_PLG_LINF] << "\n";
    Log::info() << "Intersection polygons : "
                << remap_stat.counts[RemapStat::Counts::INT_PLG] << "\n";
    Log::info() << "    Total mismatch area in polygons intersections : "
                << remap_stat.errors[RemapStat::Errors::GEO_L1] << "\n";
    Log::info() << "    Source-cell maximal mismatch area in polygons intersections : "
                << remap_stat.errors[RemapStat::Errors::GEO_LINF] << "\n";
    Log::info() << "Non covered source polygons : "
                << remap_stat.counts[RemapStat::Counts::UNCVR_SRC] << "\n";
    Log::info() << "Diff in source mesh vs target mesh coverage with polygons : "
                << remap_stat.errors[RemapStat::Errors::GEO_DIFF] << "\n";
    Log::info() << "  1-order remap took " << elapsed_seconds.count() << " seconds.\n";
    outfile << std::setw(10) << elapsed_seconds.count();
    outfile << std::setw(10) << remap_stat.errors[RemapStat::Errors::GEO_DIFF];
    outfile << std::setw(10) << remap_stat.errors[RemapStat::Errors::REMAP_L2]
            << std::setw(10) << remap_stat.errors[RemapStat::Errors::REMAP_LINF];
    print_remap_errors(consMethod, outfile);
    output::Gmsh("cons-remap_tgtfield-1ord.msh", gmsh_config).write(tgt_field);
    output::Gmsh("cons-remap_difffield-1ord.msh", gmsh_config).write(diff_field);

    consMethod.set_order(2);
    start = std::chrono::system_clock::now();
    consMethod.do_execute(src_field, tgt_field);
    elapsed_seconds = std::chrono::system_clock::now() - start;

    // remap statistics
    consMethod.remap_stat(src_vals, tgt_vals, &diff_vals, func);
	remap_stat = consMethod.remap_stat();
    Log::info() << "  2-order remap took " << elapsed_seconds.count() << " seconds.\n";
    outfile << std::setw(10) << elapsed_seconds.count();
    print_remap_errors(consMethod, outfile);
    output::Gmsh("cons-remap_tgtfield-2ord.msh", gmsh_config).write(tgt_field);
    output::Gmsh("cons-remap_difffield-2ord.msh", gmsh_config).write(diff_field);

    (outfile << "\n").flush();
    Log::info().unindent();
}

CASE("test_interpolation_conservative") {
    std::stringstream ss;
    ss << "# (1) s-grid   (2) t-grid   (3) setup [s]   (4) err.polygon.create";
    ss << "   (5) err.polygon.intersecting.L1\n#(6) err.polygon.intersecting.Lmax   (7) Time 1st-rmp [sec]";
    ss << "# (8) err.1st.ana.L2   (9) err.1st.ana.Lmax   (10) err.1st.global.cons\n#";
    ss << " (11) Time 2nd-remap [sec]   (12) 2nd-ana-err.L2   (13) 2nd-ana-err.Lmax   (14) err.2nd.global.cons\n";
    for (int i = 1; i < 15; ++i) {
        ss << std::setw(10) << i;
    }
    ss << "\n";

    SECTION("analytic constfunc") {
        auto func = [](const PointLonLat& p) { return 1.; };
        std::ofstream outfile;
//        int func_id = eckit::Resource<std::string>("--test", "Jones_Y22"));
        outfile.open("cons-remap_constfunc.dat", std::ios_base::app);
        outfile << "# Test -- analytic function = 1\n";
        outfile << std::scientific << std::setprecision(1);
        outfile << ss.str();
        // Allow to override via command-line, e.g.
        //     <program> --src-grid O16 --tgt-grid O32
        auto src_grid = Grid{eckit::Resource<std::string>("--src-grid", "H16")};
        auto tgt_grid = Grid{eckit::Resource<std::string>("--tgt-grid", "H32")};
 //       do_remapping_test(src_grid, tgt_grid, func, outfile);
        return;


        const int start_res            = 32;
        const int end_res              = start_res << 0;
        std::vector<std::string> grids = {"F", "N", "O", "H"};
        for (int i = start_res; i <= end_res; i *= 2) {
            for (int gi = 0; gi < grids.size(); gi++) {
                auto gridA = Grid(grids[gi] + std::to_string(i));
                for (int gj = 0; gj < grids.size(); gj++) {
                    auto gridB = Grid(grids[gj] + std::to_string(i));
                    do_remapping_test(gridA, gridB, func, outfile);
                }
            }
        }

        outfile.close();
    }

    SECTION("analytic Y_2^2 as in Jones - scaling") {
        std::ofstream outfile;
        outfile.open("cons-remap_JonesY22_scaling.dat", std::ios_base::app);
        outfile << "# Test -- analytic Y_2^2 as in Jones\n";
        outfile << std::scientific << std::setprecision(1);
        outfile << ss.str();
        auto func = [](const PointLonLat& p) {
            double cos = std::cos(0.025 * p[0]);
            return 2. + cos * cos * std::cos(2 * 0.025 * p[1]);
        };

        // Allow to override via command-line, e.g.
        //     <program> --src-grid O16 --tgt-grid O32
        auto src_grid = Grid{eckit::Resource<std::string>("--src-grid", "H16")};
        auto tgt_grid = Grid{eckit::Resource<std::string>("--tgt-grid", "H32")};
        do_remapping_test(src_grid, tgt_grid, func, outfile);
        return;

        const int start_res            = 32;
        const int end_res              = start_res << 0;
        std::vector<std::string> grids = {"N", "O", "H", "F"};
        for (int gi = 0; gi < grids.size(); gi++) {
            for (int gj = 0; gj < grids.size(); gj++) {
                for (int i = start_res; i <= end_res; i *= 2) {
                    for (int j = end_res; j <= end_res; j *= 2) {
                        auto gridA = Grid(grids[gi] + std::to_string(i));
                        auto gridB = Grid(grids[gj] + std::to_string(j));
                        do_remapping_test(gridA, gridB, func, outfile);
                    }
                }
            }
        }

        outfile.close();
    }

    SECTION("analytic Y_2^2 as in Jones - all2all") {
        return;
        std::ofstream outfile;
        outfile.open("cons-remap_JonesY22_all2all.dat", std::ios_base::app);
        outfile << "# Test -- analytic Y_2^2 as in Jones\n";
        outfile << std::scientific << std::setprecision(1);
        outfile << ss.str();
        auto func = [](const PointLonLat& p) {
            double cos = std::cos(0.025 * p[0]);
            return 2. + cos * cos * std::cos(2 * 0.025 * p[1]);
        };

        const int start_res            = 32;
        const int end_res              = start_res << 0;
        std::vector<std::string> grids = {"N", "O", "H"};
        for (int i = start_res; i <= end_res; i *= 2) {
            for (int gi = 0; gi < grids.size(); gi++) {
                auto gridA = Grid(grids[gi] + std::to_string(i));
                for (int gj = 0; gj < grids.size(); gj++) {
                    for (int j = start_res; j <= end_res; j *= 2) {
                        auto gridB = Grid(grids[gj] + std::to_string(j));
                        do_remapping_test(gridA, gridB, func, outfile);
                    }
                }
            }
        }

        outfile.close();
    }

    SECTION("analytic Hill as in Jones") {
        return;
        std::ofstream outfile;
        outfile.open("cons-remap_JonesHill.dat", std::ios_base::app);
        outfile << "# Test -- analytic Hill as in Jones\n";
        outfile << std::scientific << std::setprecision(1);
        outfile << ss.str();
        auto func = [](const PointLonLat& p) {
            PointXYZ c = {1., 0., 0.};
            PointXYZ p_sph;
            eckit::geometry::Sphere::convertSphericalToCartesian(1., p, p_sph);
            double r = PointXYZ::norm(p_sph - c);
            return 2. + std::cos(M_PI * r / 10.);
        };
        outfile.close();
    }
}

}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
