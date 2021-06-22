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

Grid localgrid( int nx, int ny ) {
    util::Config gridspec;
    gridspec.set( "type", "regional" );
    gridspec.set( "nx", nx );
    gridspec.set( "ny", ny );
    gridspec.set( "north", 80 );
    gridspec.set( "south", 0 );
    gridspec.set( "west", 0 );
    gridspec.set( "east", 90 );
    return Grid{gridspec};
}

void compute_geom_errors( const array::ArrayView<double, 1>& src_vals, const array::ArrayView<double, 1>& tgt_vals,
                          ConservativeMethod& conservativeMethod, double func( const PointLonLat& ),
                          std::ofstream& outfile ) {
    double src_sum = 0.;
    double tgt_sum = 0.;
    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        src_sum += conservativeMethod.src_area( scell );
    }
    for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
        tgt_sum += conservativeMethod.tgt_area( tcell );
    }
    Log::info() << "    cons err in polygon create     : " << std::abs( src_sum - tgt_sum ) * 0.25 * M_1_PI << "\n";
    outfile << std::setw( 10 ) << std::abs( src_sum - tgt_sum ) * 0.25 * M_1_PI;

    double err_1   = 0.;
    double err_max = 0.;
    size_t no_iplg = 0;
    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        double diff_cell   = conservativeMethod.src_area( scell );
        const auto& iparam = conservativeMethod.iparam()[scell];
        for ( idx_t icell = 0; icell < iparam.weights.size(); ++icell ) {
            diff_cell -= iparam.weights[icell];
        }
        no_iplg += iparam.weights.size();
        err_1 += std::abs( diff_cell );
        err_max = std::max( err_max, std::abs( diff_cell ) );
    }
    err_1 *= 0.25 * M_1_PI;
    Log::info() << "    size of src_grid, tgt_grid, supergrid: " << src_vals.size() << " " << tgt_vals.size() << " "
                << no_iplg << "\n";
    Log::info() << "    cons err in polygon intersect  : (L1) " << err_1 << " (Lmax) " << err_max << "\n";
    outfile << std::setw( 10 ) << err_1 << std::setw( 10 ) << err_max;
}

void compute_field_errors( const array::ArrayView<double, 1>& src_vals, const array::ArrayView<double, 1>& tgt_vals,
                           array::ArrayView<double, 1>& diff_vals, ConservativeMethod& conservativeMethod,
                           double func( const PointLonLat& ), std::ofstream& outfile ) {
    double global_cons_err = 0;
    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        diff_vals( scell ) = src_vals( scell ) * conservativeMethod.src_area( scell );
        global_cons_err += diff_vals( scell );
        const auto& iparam = conservativeMethod.iparam()[scell];
        for ( idx_t icell = 0; icell < iparam.weights.size(); ++icell ) {
            diff_vals( scell ) -= tgt_vals( iparam.tcell_id[icell] ) * iparam.weights[icell];
        }
        diff_vals( scell ) = std::abs( diff_vals( scell ) );
    }

    double err_2   = 0.;
    double err_max = 0.;
    for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
        global_cons_err -= tgt_vals( tcell ) * conservativeMethod.tgt_area( tcell );
        auto p = conservativeMethod.tgt_centroid( tcell );
        PointLonLat pll;
        eckit::geometry::Sphere::convertCartesianToSpherical( 1., p, pll );
        double err_l = std::abs( tgt_vals( tcell ) - func( pll ) );
        err_2 += err_l * err_l * conservativeMethod.tgt_area( tcell );
        err_max = std::max( err_max, err_l );
    }
    err_2           = std::sqrt( err_2 * 0.25 * M_1_PI );
    global_cons_err = std::sqrt( std::abs(global_cons_err) * 0.25 * M_1_PI );
    Log::info() << "    " << conservativeMethod.order() << "-order remap analytical error : (L2) " << err_2
                << " (Lmax) " << err_max << "\n";
    Log::info() << "    " << conservativeMethod.order() << "-order global remap error : " << std::abs( global_cons_err )
                << "\n";
    outfile << std::setw( 10 ) << err_2 << std::setw( 10 ) << err_max << std::setw( 10 ) << std::abs( global_cons_err );
}

void do_remapping_test( Grid src_grid, Grid tgt_grid, double func( const PointLonLat& ), std::ofstream& outfile ) {
    util::Config config;
    config.set( "include_pole", true );
    config.set( "matrix_free", false );
    config.set( "normalise_intersections", 1 );
    config.set( "triangulate", false );

    outfile << std::setw( 10 ) << src_grid.name() << std::setw( 10 ) << tgt_grid.name();

    ConservativeMethod conservativeMethod( config );

    auto start = std::chrono::system_clock::now();
    conservativeMethod.do_setup( src_grid, tgt_grid );
    std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start;
    Log::info() << "REMAPPING: " << src_grid.name() << " --> " << tgt_grid.name() << "\n";
    Log::info() << "  Setup (computing supermesh) took " << elapsed_seconds.count() << " seconds.\n";
    outfile << std::setw( 10 ) << elapsed_seconds.count();

    const auto& src_fs = conservativeMethod.source();
    const auto& tgt_fs = conservativeMethod.target();
    auto src_field     = src_fs.createField<double>();
    auto tgt_field     = tgt_fs.createField<double>();
    auto src_vals      = array::make_view<double, 1>( src_field );
    auto tgt_vals      = array::make_view<double, 1>( tgt_field );

    compute_geom_errors( src_vals, tgt_vals, conservativeMethod, func, outfile );
    output::Gmsh( "cons-remap_smesh.msh", util::Config( "coordinates", "lonlat" ) )
        .write( conservativeMethod.src_mesh() );
    output::Gmsh( "cons-remap_tmesh.msh", util::Config( "coordinates", "lonlat" ) )
        .write( conservativeMethod.tgt_mesh() );

    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        auto p = conservativeMethod.src_centroid( scell );
        PointLonLat pll;
        eckit::geometry::Sphere::convertCartesianToSpherical( 1., p, pll );
        src_vals( scell ) = func( pll );
    }
    output::Gmsh( "cons-remap_sfield.msh", util::Config( "coordinates", "lonlat" ) ).write( src_field );

    conservativeMethod.set_order( 1 );
    start = std::chrono::system_clock::now();
    conservativeMethod.do_execute( src_field, tgt_field );
    elapsed_seconds = std::chrono::system_clock::now() - start;
    Log::info() << "  1-order remap took " << elapsed_seconds.count() << " seconds.\n";
    outfile << std::setw( 10 ) << elapsed_seconds.count();
    output::Gmsh( "cons-remap_tfield-1ord.msh", util::Config( "coordinates", "lonlat" ) ).write( tgt_field );

    auto diff_field = src_fs.createField<double>();
    auto diff_vals  = array::make_view<double, 1>( diff_field );
    compute_field_errors( src_vals, tgt_vals, diff_vals, conservativeMethod, func, outfile );
    output::Gmsh( "cons-remap_dfield-1ord.msh", util::Config( "coordinates", "lonlat" ) ).write( diff_field );

    conservativeMethod.set_order( 2 );
    start = std::chrono::system_clock::now();
    conservativeMethod.do_execute( src_field, tgt_field );
    elapsed_seconds = std::chrono::system_clock::now() - start;
    Log::info() << "  2-order remap took " << elapsed_seconds.count() << " seconds.\n";
    outfile << std::setw( 10 ) << elapsed_seconds.count();
    output::Gmsh( "cons-remap_tfield-2ord.msh", util::Config( "coordinates", "lonlat" ) ).write( tgt_field );

    compute_field_errors( src_vals, tgt_vals, diff_vals, conservativeMethod, func, outfile );
    output::Gmsh( "cons-remap_dfield-2ord.msh", util::Config( "coordinates", "lonlat" ) ).write( diff_field );

    ( outfile << "\n" ).flush();
}

CASE( "test_interpolation_conservative" ) {
    std::stringstream ss;
    ss << "# (1) s-grid   (2) t-grid   (3) setup [s]   (4) err.polygon.create";
    ss << "   (5) err.polygon.intersecting.L1\n#(6) err.polygon.intersecting.Lmax   (7) Time 1st-rmp [sec]";
    ss << "# (8) err.1st.ana.L2   (9) err.1st.ana.Lmax   (10) err.1st.global.cons\n#";
    ss << " (11) Time 2nd-remap [sec]   (12) 2nd-ana-err.L2   (13) 2nd-ana-err.Lmax   (14) err.2nd.global.cons\n";
    for ( int i = 1; i < 15; ++i ) {
        ss << std::setw( 10 ) << i;
    }
    ss << "\n";

    SECTION( "analytic constfunc" ) {
        auto func = []( const PointLonLat& p ) { return 1.; };
        std::ofstream outfile;
        outfile.open( "cons-remap_constfunc.dat", std::ios_base::app );
        outfile << "# Test -- analytic function = 1\n";
        outfile << std::scientific << std::setprecision( 1 );
        outfile << ss.str();

        const int start_res            = 32;
        const int end_res              = start_res << 0;
        std::vector<std::string> grids = {"F", "N", "O", "H"};
        for ( int i = start_res; i <= end_res; i *= 2 ) {
            for ( int gi = 0; gi < grids.size(); gi++ ) {
                auto gridA = Grid( grids[gi] + std::to_string( i ) );
                for ( int gj = 0; gj < grids.size(); gj++ ) {
                    auto gridB = Grid( grids[gj] + std::to_string( i ) );
                    do_remapping_test( gridA, gridB, func, outfile );
                }
            }
        }

        outfile.close();
    }

    SECTION( "analytic Y_2^2 as in Jones - scaling" ) {
        std::ofstream outfile;
        outfile.open( "cons-remap_JonesY22_scaling.dat", std::ios_base::app );
        outfile << "# Test -- analytic Y_2^2 as in Jones\n";
        outfile << std::scientific << std::setprecision( 1 );
        outfile << ss.str();
        auto func = []( const PointLonLat& p ) {
            double cos = std::cos( 0.025 * p[0] );
            return 2. + cos * cos * std::cos( 2 * 0.025 * p[1] );
        };

        const int start_res            = 32;
        const int end_res              = start_res << 0;
        std::vector<std::string> grids = {"N", "O", "H", "F"};
        for ( int gi = 0; gi < grids.size(); gi++ ) {
            for ( int gj = 0; gj < grids.size(); gj++ ) {
                for ( int i = start_res; i <= end_res; i *= 2 ) {
                    for ( int j = end_res; j <= end_res; j *= 2 ) {
                        auto gridA = Grid( grids[gi] + std::to_string( i ) );
                        auto gridB = Grid( grids[gj] + std::to_string( j ) );
                        do_remapping_test( gridA, gridB, func, outfile );
                    }
                }
            }
        }

        outfile.close();
    }

    SECTION( "analytic Y_2^2 as in Jones - all2all" ) {
        std::ofstream outfile;
        outfile.open( "cons-remap_JonesY22_all2all.dat", std::ios_base::app );
        outfile << "# Test -- analytic Y_2^2 as in Jones\n";
        outfile << std::scientific << std::setprecision( 1 );
        outfile << ss.str();
        auto func = []( const PointLonLat& p ) {
            double cos = std::cos( 0.025 * p[0] );
            return 2. + cos * cos * std::cos( 2 * 0.025 * p[1] );
        };

        const int start_res            = 32;
        const int end_res              = start_res << 0;
        std::vector<std::string> grids = {"N", "O", "H"};
        for ( int i = start_res; i <= end_res; i *= 2 ) {
            for ( int gi = 0; gi < grids.size(); gi++ ) {
                auto gridA = Grid( grids[gi] + std::to_string( i ) );
                for ( int gj = 0; gj < grids.size(); gj++ ) {
                    for ( int j = start_res; j <= end_res; j *= 2 ) {
                        auto gridB = Grid( grids[gj] + std::to_string( j ) );
                        do_remapping_test( gridA, gridB, func, outfile );
                    }
                }
            }
        }

        outfile.close();
    }

    SECTION( "analytic Hill as in Jones" ) {
        std::ofstream outfile;
        outfile.open( "cons-remap_JonesHill.dat", std::ios_base::app );
        outfile << "# Test -- analytic Hill as in Jones\n";
        outfile << std::scientific << std::setprecision( 1 );
        outfile << ss.str();
        auto func = []( const PointLonLat& p ) {
            PointXYZ c = {1., 0., 0.};
            PointXYZ p_sph;
            eckit::geometry::Sphere::convertSphericalToCartesian( 1., p, p_sph );
            double r = PointXYZ::norm( p_sph - c );
            return 2. + std::cos( M_PI * r / 10. );
        };
        outfile.close();
    }
}

}  // namespace test
}  // namespace atlas


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
