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
            no_iplg += iparam.weights.size();
        }
        if ( std::abs( diff_cell ) > 1e-8 ) {
            Log::info() << "diff_cell, scell: " << diff_cell << " " << scell << "\n";
            for ( idx_t i = 0; i < iparam.weights.size(); ++i ) {
                Log::info() << iparam.cell_id[i] << " ";
            }
            Log::info() << "\n";
        }
        err_1 += diff_cell;
        err_max = std::max( err_max, std::abs( diff_cell ) );
    }
    err_1 = std::abs( err_1 );
    Log::info() << "    Size of src_grid, tgt_grid, supergrid: " << src_vals.size() << " " << tgt_vals.size() << " "
                << no_iplg << "\n";
    Log::info() << "    cons err in polygon intersect  : (L2) " << err_1 << " (Lmax) " << err_max << "\n";
    outfile << std::setw( 10 ) << err_1 << std::setw( 10 ) << err_max;
}

void compute_field_errors( const array::ArrayView<double, 1>& src_vals, const array::ArrayView<double, 1>& tgt_vals,
                           array::ArrayView<double, 1>& diff_vals, ConservativeMethod& conservativeMethod,
                           double func( const PointLonLat& ), std::ofstream& outfile ) {
    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        diff_vals( scell ) = src_vals( scell ) * conservativeMethod.src_area( scell );
        const auto& iparam = conservativeMethod.iparam()[scell];
        for ( idx_t icell = 0; icell < iparam.weights.size(); ++icell ) {
            diff_vals( scell ) -= tgt_vals( iparam.cell_id[icell] ) * iparam.weights[icell];
        }
        diff_vals( scell ) = std::abs( diff_vals( scell ) );
        //diff_vals( scell ) = std::abs( diff_vals( scell ) ) / conservativeMethod.src_area( scell );
    }

    double err_2   = 0.;
    double err_max = 0.;
    for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
        auto p = conservativeMethod.tgt_centroid( tcell );
        PointLonLat pll;
        eckit::geometry::Sphere::convertCartesianToSpherical( 1., p, pll );
        double err_l = std::abs( tgt_vals( tcell ) - func( pll ) );
        err_2 += err_l * err_l * conservativeMethod.tgt_area( tcell );
        err_max = std::max( err_max, err_l );
    }
    err_2 = std::sqrt( err_2 * 0.25 * M_1_PI );
    Log::info() << "    " << conservativeMethod.order() << "-order remap analytical error : (L2) " << err_2
                << " (Lmax) " << err_max << "\n";
    outfile << std::setw( 10 ) << err_2 << std::setw( 10 ) << err_max;

    err_2   = 0.;
    err_max = 0.;
    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        const auto& iparam = conservativeMethod.iparam()[scell];
        for ( idx_t icell = 0; icell < iparam.weights.size(); ++icell ) {
            double err_l = std::abs( src_vals( scell ) - tgt_vals( iparam.cell_id[icell] ) );
            err_2 += err_l * err_l * iparam.weights[icell];
            err_max = std::max( err_max, err_l );
        }
    }
    err_2 = std::sqrt( err_2 * 0.25 * M_1_PI );
    Log::info() << "    " << conservativeMethod.order() << "-order remap mesh2mesh error  : (L2) " << err_2
                << " (Lmax) " << err_max << "\n";
    outfile << std::setw( 10 ) << err_2 << std::setw( 10 ) << err_max;
}

void do_remapping_test( Grid src_grid, Grid tgt_grid, double func( const PointLonLat& ), std::ofstream& outfile ) {
    util::Config config;
    config.set( "include_pole", true );
    //config.set( "triangulate", true );

    outfile << std::setw( 10 ) << src_grid.name() << std::setw( 10 ) << tgt_grid.name();
    bool src_healpix = ( src_grid.name()[0] == 'H' ) or ( src_grid.name()[0] == 'h' );
    bool tgt_healpix = ( tgt_grid.name()[0] == 'H' ) or ( tgt_grid.name()[0] == 'h' );
    auto src_meshgen = ( src_healpix ? MeshGenerator{"healpix"} : MeshGenerator{"structured", config} );
    auto tgt_meshgen = ( tgt_healpix ? MeshGenerator{"healpix"} : MeshGenerator{"structured", config} );
    Mesh src_mesh    = src_meshgen.generate( src_grid );
    Mesh tgt_mesh    = tgt_meshgen.generate( tgt_grid );

    ConservativeMethod conservativeMethod( config );

    functionspace::CellColumns src_fs( src_mesh );
    functionspace::CellColumns tgt_fs( tgt_mesh );
    auto src_field = src_fs.createField<double>();
    auto tgt_field = tgt_fs.createField<double>();
    auto src_vals  = array::make_view<double, 1>( src_field );
    auto tgt_vals  = array::make_view<double, 1>( tgt_field );

    auto start = std::chrono::system_clock::now();
    conservativeMethod.do_setup( src_mesh, tgt_mesh );
    std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start;
    Log::info() << "REMAPPING: " << src_grid.name() << " --> " << tgt_grid.name() << "\n";
    Log::info() << "  Setup (computing supermesh) took " << elapsed_seconds.count() << " seconds.\n";
    outfile << std::setw( 10 ) << elapsed_seconds.count();

    compute_geom_errors( src_vals, tgt_vals, conservativeMethod, func, outfile );
    output::Gmsh( "cons-remap_smesh.msh", util::Config( "coordinates", "lonlat" ) ).write( src_mesh );
    output::Gmsh( "cons-remap_tmesh.msh", util::Config( "coordinates", "lonlat" ) ).write( tgt_mesh );

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

    ( outfile << "\n\n" ).flush();
}

CASE( "test_interpolation_conservative" ) {
    std::stringstream ss;
    ss << "# (1) s-grid   (2) t-grid   (3) setup [s]   (4) err.polygon.create";
    ss << "   (5) err.polygon.intrsc.L1   (6) err.polygon.intrsc.Lmax   (7) 1st-rmp [s]\n";
    ss << "# (8) err.1st.ana.L2   (9) err.1st.ana.Lmax   (10) err.1st.mesh2mesh.err.L2";
    ss << "   (11) err.1st.mesh2mesh.err.Lmax   (12) 2nd-remap[s]   (13) 2nd-ana-err.L2\n";
    ss << "# (14) 2nd-ana-err.Lmax   (15) 2nd-mesh2mesh-err.L2   (16) 2nd-mesh2mesh-err.Lmax\n";
    for ( int i = 1; i < 17; ++i ) {
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
        do_remapping_test( Grid( "F32" ), Grid( "H32" ), func, outfile );
        do_remapping_test( Grid( "H32" ), Grid( "O32" ), func, outfile );
        do_remapping_test( Grid( "O32" ), Grid( "N32" ), func, outfile );
        do_remapping_test( Grid( "N32" ), Grid( "H32" ), func, outfile );
        outfile.close();
    }

    SECTION( "analytic Y_2^2 as in Jones" ) {
        std::ofstream outfile;
        outfile.open( "cons-remap_JonesY22.dat", std::ios_base::app );
        outfile << "# Test -- analytic Y_2^2 as in Jones\n";
        outfile << std::scientific << std::setprecision( 1 );
        outfile << ss.str();
        auto func = []( const PointLonLat& p ) {
            double cos = std::cos( 0.025 * p[0] );
            return 2. + cos * cos * std::cos( 2 * 0.025 * p[1] );
        };
        do_remapping_test( Grid( "O1" ), Grid( "H128" ), func, outfile );
        do_remapping_test( Grid( "O2" ), Grid( "H128" ), func, outfile );
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
