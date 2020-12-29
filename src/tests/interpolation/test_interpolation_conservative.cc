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
		os <<"centroids         : " <<centroids <<"\n"
		   <<"fractional area   : " <<weights <<"\n";
	}
};

double func( const double& lon, const double& lat ) {
	return 1.; //lon + lat;
}


CASE( "test_interpolation_conservative" ) {
    Grid src_grid = localgrid( 2, 2 );
    Grid tgt_grid = localgrid( 2, 3 );
	MeshGenerator meshgen( "regular" );
	Mesh src_mesh = meshgen.generate( src_grid );
	Mesh tgt_mesh = meshgen.generate( tgt_grid );
	
	std::vector< PointLonLat > pts_ll;
	const idx_t src_nb_cells = src_mesh.cells().size();
	const auto& src_node_connectivity = src_mesh.cells().node_connectivity();
	std::vector< CSPolygon > src_csp( src_nb_cells );
	const auto src_lonlat = array::make_view<double,2>( src_mesh.nodes().lonlat() );
	for( idx_t jcell=0; jcell < src_nb_cells; ++jcell ) {
		const idx_t nb_nodes = src_node_connectivity.cols( jcell );
		pts_ll.clear();
		pts_ll.resize( nb_nodes );
		for( idx_t jnode = 0; jnode < nb_nodes; ++jnode ) {
			idx_t inode = src_node_connectivity( jcell, jnode );
			pts_ll[ nb_nodes - 1 - jnode ] = PointLonLat{ src_lonlat(inode,0), src_lonlat(inode,1) };
		}
		src_csp[ jcell ] =  CSPolygon( pts_ll );
	}

	const idx_t tgt_nb_cells = tgt_mesh.cells().size();
	const auto& tgt_node_connectivity = tgt_mesh.cells().node_connectivity();
	std::vector< CSPolygon > tgt_csp( tgt_nb_cells );
	const auto tgt_lonlat = array::make_view<double,2>( tgt_mesh.nodes().lonlat() );
	for( idx_t jcell = 0; jcell < tgt_nb_cells; ++jcell ) {
		const idx_t nb_nodes = tgt_node_connectivity.cols( jcell );
		pts_ll.clear();
		pts_ll.resize( nb_nodes );
		for( idx_t jnode=0; jnode < nb_nodes; ++jnode ) {
			idx_t inode = tgt_node_connectivity( jcell, jnode );
			pts_ll[ nb_nodes - 1 - jnode ] = PointLonLat{ tgt_lonlat(inode,0), tgt_lonlat(inode,1) };
		}
		tgt_csp[ jcell ] =  CSPolygon( pts_ll );
	}

	std::vector< InterpolationParameters > ip( src_nb_cells );
	// brute force (!) needs to be changed
	for( idx_t scell = 0; scell < src_nb_cells; ++scell ) {
		for( idx_t tcell=0; tcell < tgt_nb_cells; ++tcell ) {
			CSPolygon csp_i = src_csp[ scell ].intersect( tgt_csp[ tcell ] );
			if ( csp_i.area() > 0 ) {
				ip[ scell ].cell_id.emplace_back( tcell );
				ip[ scell ].weights.emplace_back( csp_i.area()/src_csp[ scell ].area() );
				ip[ scell ].centroids.emplace_back( csp_i.centroid() );
			}
		}
	}

	for( idx_t scell = 0; scell < src_nb_cells; ++scell ) {
		Log::info() <<"Source-Polygon " <<src_csp[ scell ] <<" intersects Target-Polygons:\n";
		for( idx_t tcell = 0; tcell < ip[ scell ].cell_id.size(); ++tcell ) {
			Log::info() <<"   " <<tgt_csp[ tcell ] <<"\n";
		}
		ip[ scell ].print( Log::info() );
	}

	functionspace::CellColumns src_fs( src_mesh );
	functionspace::CellColumns tgt_fs( tgt_mesh );
	auto src_field = src_fs.createField< double >();
	auto tgt_field = tgt_fs.createField< double >();
	auto src_vals = array::make_view< double, 1 >( src_field );

	Log::info() <<"   src_fs.size: " <<src_fs.size() <<"\n";
	Log::info() <<"   tgt_fs.size: " <<tgt_fs.size() <<"\n";
	Log::info() <<"   src_vals.size: " <<src_vals.size() <<"\n";

	for ( idx_t jnode = 0; jnode < src_vals.size(); ++jnode ) {
		src_vals( jnode ) = func( src_lonlat(jnode,0), src_lonlat(jnode,1) );
	}
}

CASE( "test_interpolation_conservative" ) {
}


}  // namespace test
}  // namespace atlas


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
