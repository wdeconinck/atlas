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
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"
#include "atlas/interpolation/method/knn/ConservativeMethod.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

using CSPolygon = util::ConvexSphericalPolygon;
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

double func( const double& lon, const double& lat ) {
    return lon + lat;
}

double func( const double& x, const double& y, const double& z ) {
    return 100 * x + 10 * y + z;
}


CASE( "test_interpolation_conservative" ) {
    Grid src_grid = localgrid( 3, 33 );
    Grid tgt_grid = localgrid( 33, 3 );
    MeshGenerator meshgen( "regular" );
    Mesh src_mesh = meshgen.generate( src_grid );
    Mesh tgt_mesh = meshgen.generate( tgt_grid );

	util::Config config;
	config.set( "order", 2 );
	ConservativeMethod conservativeMethod( config );

    functionspace::CellColumns src_fs( src_mesh );
    functionspace::CellColumns tgt_fs( tgt_mesh );
    auto src_field = src_fs.createField<double>();
    auto tgt_field = tgt_fs.createField<double>();
    auto src_vals  = array::make_view<double, 1>( src_field );
    auto tgt_vals  = array::make_view<double, 1>( tgt_field );

	conservativeMethod.do_setup( src_mesh, tgt_mesh );

    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        auto p            = conservativeMethod.src_centroid(scell);
        src_vals( scell ) = func( p[0], p[1], p[2] );
    }

	conservativeMethod.do_execute( src_field, tgt_field );

    // validate first order conservation property
    double err = 0.;
    for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
        double scell_area               = 0.;
        const auto& iparam = conservativeMethod.iparam();
        for ( idx_t icell = 0; icell < iparam[scell].weights.size(); ++icell ) {
            scell_area += iparam[scell].weights[icell];
        }
        Log::info() << " scell, scell_area, diff: " << conservativeMethod.src_area(scell) << " " << scell_area << " "
                    << conservativeMethod.src_area(scell) - scell_area << "\n";
        err += std::abs( conservativeMethod.src_area(scell) - scell_area );
    }

    // validate error
    err = 0.;
    for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
        auto p = conservativeMethod.tgt_centroid(tcell);
        err += std::abs( tgt_vals( tcell ) - func( p[0], p[1], p[2] ) );
    }
    err /= tgt_vals.size();
    Log::info() << "target field err: " << err << "\n";
}

}  // namespace test
}  // namespace atlas


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
