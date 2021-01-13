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

void compute_errors( array::ArrayView<double, 1>& src_vals, array::ArrayView<double, 1>& tgt_vals,
                     ConservativeMethod& conservativeMethod, double func( const PointLonLat& ) ) {
    double err_2   = 0.;
    double err_max = 0.;
    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        double tgt_sum_cell = 0.;
        const auto& iparam  = conservativeMethod.iparam()[scell];
        for ( idx_t icell = 0; icell < iparam.weights.size(); ++icell ) {
            tgt_sum_cell += iparam.weights[icell];
        }
        double err_l = conservativeMethod.src_area( scell ) - tgt_sum_cell;
        err_2 += err_l * err_l;
        err_max = std::max( err_max, std::abs( err_l ) );
    }
    err_2 = std::sqrt( err_2 * 0.25 * M_1_PI );
    err_max *= 0.25 * M_1_PI;
    Log::info() << "     local cons error : (L2) " << err_2 << " (Lmax) " << err_max << "\n";

    double src_sum = 0.;
    double tgt_sum = 0.;
    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        src_sum += conservativeMethod.src_area( scell );
    }
    for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
        tgt_sum += conservativeMethod.tgt_area( tcell );
    }
    Log::info() << "    global cons error : " << std::abs( src_sum - tgt_sum ) * 0.25 * M_1_PI << "\n";

    err_2   = 0.;
    err_max = 0.;
    for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
        auto p = conservativeMethod.tgt_centroid( tcell );
        PointLonLat pll;
        eckit::geometry::Sphere::convertCartesianToSpherical( 1., p, pll );
        double afunc = func( pll );
        double err_l = std::abs( tgt_vals( tcell ) - afunc ) * conservativeMethod.tgt_area( tcell );
        err_2 += err_l * std::abs( tgt_vals( tcell ) - afunc );
        err_max = std::max( err_max, err_l );
    }
    err_2 = std::sqrt( err_2 * 0.25 * M_1_PI );
    err_max *= 0.25 * M_1_PI;
    Log::info() << "	remap error : (L2) " << err_2 << " (Lmax) " << err_max << "\n";
}

void do_remapping_test( Grid src_grid, Grid tgt_grid, double func( const PointLonLat& ) ) {
    util::Config config;

    auto src_meshgen =
        ( src_grid.name() == "healpix" ? MeshGenerator{"healpix"}
                                       : MeshGenerator{"structured", util::Config( "include_pole", true )} );
    auto tgt_meshgen =
        ( tgt_grid.name() == "healpix" ? MeshGenerator{"healpix"}
                                       : MeshGenerator{"structured", util::Config( "include_pole", true )} );
    Mesh src_mesh = src_meshgen.generate( src_grid );
    Mesh tgt_mesh = tgt_meshgen.generate( tgt_grid );

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
    Log::info() << "ConservativeMethod::do_setup took " << elapsed_seconds.count() << " seconds.\n";

    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        auto p = conservativeMethod.src_centroid( scell );
        PointLonLat pll;
        eckit::geometry::Sphere::convertCartesianToSpherical( 1., p, pll );
        src_vals( scell ) = func( pll );
    }

    Log::info() << " ====== " << src_grid.name() << " --> " << tgt_grid.name() << " 1. order\n";
    conservativeMethod.set_order( 1 );
    start = std::chrono::system_clock::now();
    conservativeMethod.do_execute( src_field, tgt_field );
    elapsed_seconds = std::chrono::system_clock::now() - start;
    Log::info() << "ConservativeMethod::do_execute took " << elapsed_seconds.count() << " seconds.\n";

    compute_errors( src_vals, tgt_vals, conservativeMethod, func );

    Log::info() << " ====== " << src_grid.name() << " --> " << tgt_grid.name() << " 2. order\n";
    conservativeMethod.set_order( 2 );
    start = std::chrono::system_clock::now();
    conservativeMethod.do_execute( src_field, tgt_field );
    elapsed_seconds = std::chrono::system_clock::now() - start;
    Log::info() << "ConservativeMethod::do_execute took " << elapsed_seconds.count() << " seconds.\n";

    compute_errors( src_vals, tgt_vals, conservativeMethod, func );
}

CASE( "test_interpolation_conservative" ) {
    SECTION( "analytic function = 1" ) {
        auto func = []( const PointLonLat& p ) { return 1.; };
        do_remapping_test( Grid( "F8" ), Grid( "H13" ), func );
        do_remapping_test( Grid( "H2" ), Grid( "N32" ), func );
        do_remapping_test( Grid( "N32" ), Grid( "O7" ), func );
    }

    SECTION( "analytic Y_2^2 as in Jones" ) {
        auto func = []( const PointLonLat& p ) {
            double cos = std::cos( p[0] );
            return 2. + cos * cos * std::cos( 2 * p[1] );
        };
        do_remapping_test( Grid( "F8" ), Grid( "H13" ), func );
        do_remapping_test( Grid( "H2" ), Grid( "N32" ), func );
        do_remapping_test( Grid( "N32" ), Grid( "O7" ), func );
    }

    SECTION( "analytic Hill as in Jones" ) {
        auto func = []( const PointLonLat& p ) {
            PointXYZ c = {1., 0., 0.};
            PointXYZ p_sph;
            eckit::geometry::Sphere::convertSphericalToCartesian( 1., p, p_sph );
            double r = PointXYZ::norm( p_sph - c );
            return 2. + std::cos( M_PI * r / 0.2 );
        };
        do_remapping_test( Grid( "F8" ), Grid( "H13" ), func );
        do_remapping_test( Grid( "H2" ), Grid( "N32" ), func );
        do_remapping_test( Grid( "N32" ), Grid( "O7" ), func );
    }
}

}  // namespace test
}  // namespace atlas


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
