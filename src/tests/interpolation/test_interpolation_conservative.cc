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

#include "eckit/types/FloatCompare.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildNode2CellConnectivity.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/mesh/detail/AccumulateFacets.h"
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


CASE( "test_interpolation_conservative" ) {
    Grid gridA = localgrid( 3, 3 );
    Grid gridB = localgrid( 2, 2 );
	std::cout <<"grid A:";
	for ( auto& p : gridA.lonlat() ) {
		std::cout <<p <<" ";
	}
	std::cout <<"\ngrid B:";
	for ( auto& p : gridB.lonlat() ) {
		std::cout <<p <<" ";
	}
	std::cout <<"\n";
	MeshGenerator meshgen( "regular" );
	Mesh mA = meshgen.generate( gridA );
	Mesh mB = meshgen.generate( gridB );
	auto node_glb_idx = array::make_view<gidx_t, 1>( mA.nodes().global_index() );
	auto cell_glb_idx = array::make_view<gidx_t, 1>( mA.cells().global_index() );
	std::cout <<" mA: " <<mA.nodes().size() <<" " <<mA.cells().size() <<" " <<mA.edges().size() <<"\n";
	std::cout <<" mB: " <<mB.nodes().size() <<"\n";
	
	std::vector< PointLonLat > pts_ll;
	const idx_t nb_cells_a = mA.cells().size();
	const auto& node_connectivity_a = mA.cells().node_connectivity();
	CSPolygon* ctpA = new CSPolygon[ nb_cells_a ];
	const auto lonlat_a = array::make_view<double,2>( mA.nodes().lonlat() );
	for( idx_t jcell=0; jcell < nb_cells_a; ++jcell ) {
		const idx_t nb_nodes = node_connectivity_a.cols( jcell );
		pts_ll.clear();
		pts_ll.reserve( nb_nodes );
		for( idx_t jnode=0; jnode < nb_nodes; ++jnode ) {
			idx_t inode = node_connectivity_a( jcell, jnode );
			pts_ll.insert( pts_ll.begin(), PointLonLat{ lonlat_a(inode,0), lonlat_a(inode,1) } );
		}
		std::cout <<" Forming A-CSPolygon from: " <<pts_ll <<"\n";
		ctpA[ jcell ] =  CSPolygon( pts_ll );
	}

	const idx_t nb_cells_b = mB.cells().size();
	const auto& node_connectivity_b = mB.cells().node_connectivity();
	CSPolygon* ctpB = new CSPolygon[ nb_cells_b ];
	const auto lonlat_b = array::make_view<double,2>( mB.nodes().lonlat() );
	for( idx_t jcell=0; jcell < nb_cells_b; ++jcell ) {
		const idx_t nb_nodes = node_connectivity_b.cols( jcell );
		pts_ll.clear();
		pts_ll.reserve( nb_nodes );
		for( idx_t jnode=0; jnode < nb_nodes; ++jnode ) {
			idx_t inode = node_connectivity_b( jcell, jnode );
			pts_ll.insert( pts_ll.begin(), PointLonLat{ lonlat_b(inode,0), lonlat_b(inode,1) } );
		}
		std::cout <<" Forming B-CSPolygon from: " <<pts_ll <<"\n";
		ctpB[ jcell ] =  CSPolygon( pts_ll );
	}

	CSPolygon* ctpAB = new CSPolygon[ nb_cells_a * nb_cells_b ];
	for( idx_t bcell=0; bcell < nb_cells_b; ++bcell ) {
		for( idx_t acell=0; acell < nb_cells_a; ++acell ) {
			std::cout <<" Intersecting polygon\n" <<ctpB[ bcell ] <<"\nand\n";
			std::cout <<" polygon\n" <<ctpA[ acell ] <<"\n";
			ctpAB[ acell * nb_cells_b + bcell ] = ctpA[ acell ].intersect( ctpB[ bcell ] );
			std::cout <<" and got polygon\n" <<ctpAB[ acell * nb_cells_b + bcell ] <<"\n";
			std::cout <<" 	of area " <<ctpAB[ acell * nb_cells_b + bcell ].area();
			std::cout <<" and centroid " <<ctpAB[ acell * nb_cells_b + bcell ].centroid() <<"\n\n";
		}
	}

	delete [] ctpA;
	delete [] ctpB;
}


}  // namespace test
}  // namespace atlas


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
