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
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"
#include "atlas/util/ConvexSphericalPolygon.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

using CSPolygon = util::ConvexSphericalPolygon;

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

struct InterpolationParameters {
	std::vector< int > cell_id;
	std::vector< PointXYZ > centroids;
	std::vector< double > weights;
	std::ostream& print( std::ostream& os) {
		os <<"centroids : " <<centroids <<"\n"
		   <<"weights   : " <<weights <<"\n";
	}
};


CASE( "test_interpolation_conservative" ) {
    Grid gridA = localgrid( 3, 3 );
    Grid gridB = localgrid( 4, 3 );
	MeshGenerator meshgen( "regular" );
	Mesh mA = meshgen.generate( gridA );
	Mesh mB = meshgen.generate( gridB );
	
	std::vector< PointLonLat > pts_ll;
	const idx_t nb_cells_a = mA.cells().size();
	const auto& node_connectivity_a = mA.cells().node_connectivity();
	std::vector< CSPolygon > cspA( nb_cells_a );
	const auto lonlat_a = array::make_view<double,2>( mA.nodes().lonlat() );
	for( idx_t jcell=0; jcell < nb_cells_a; ++jcell ) {
		const idx_t nb_nodes = node_connectivity_a.cols( jcell );
		pts_ll.clear();
		pts_ll.resize( nb_nodes );
		for( idx_t jnode = 0; jnode < nb_nodes; ++jnode ) {
			idx_t inode = node_connectivity_a( jcell, jnode );
			pts_ll[ nb_nodes - 1 - jnode ] = PointLonLat{ lonlat_a(inode,0), lonlat_a(inode,1) };
		}
		cspA[ jcell ] =  CSPolygon( pts_ll );
	}

	const idx_t nb_cells_b = mB.cells().size();
	const auto& node_connectivity_b = mB.cells().node_connectivity();
	std::vector< CSPolygon > cspB( nb_cells_b );
	const auto lonlat_b = array::make_view<double,2>( mB.nodes().lonlat() );
	for( idx_t jcell = 0; jcell < nb_cells_b; ++jcell ) {
		const idx_t nb_nodes = node_connectivity_b.cols( jcell );
		pts_ll.clear();
		pts_ll.resize( nb_nodes );
		for( idx_t jnode=0; jnode < nb_nodes; ++jnode ) {
			idx_t inode = node_connectivity_b( jcell, jnode );
			pts_ll[ nb_nodes - 1 - jnode ] = PointLonLat{ lonlat_b(inode,0), lonlat_b(inode,1) };
		}
		cspB[ jcell ] =  CSPolygon( pts_ll );
	}

	std::vector< InterpolationParameters > ip( nb_cells_b );
	// brute force (!) needs to be changed
	for( idx_t bcell = 0; bcell < nb_cells_b; ++bcell ) {
		for( idx_t acell=0; acell < nb_cells_a; ++acell ) {
			CSPolygon cspAB = cspA[ acell ].intersect( cspB[ bcell ] );
			if ( cspAB.area() > 0 ) {
				ip[ bcell ].cell_id.emplace_back( acell );
				ip[ bcell ].weights.emplace_back( cspAB.area() );
				ip[ bcell ].centroids.emplace_back( cspAB.centroid() );
			}
		}
	}

	for( idx_t bcell = 0; bcell < nb_cells_b; ++bcell ) {
		Log::info() <<"Grid A Polygon " <<cspB[ bcell ] <<" intersects Grid B polygons:\n";
		for( idx_t acell = 0; acell < ip[ bcell ].cell_id.size(); ++acell ) {
			Log::info() <<"   " <<cspA[ acell ] <<"\n";
		}
		ip[ bcell ].print( Log::info() );
	}
}


}  // namespace test
}  // namespace atlas


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
