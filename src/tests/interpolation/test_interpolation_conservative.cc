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
#include "atlas/mesh.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
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
	return lon+lat;
}

double func( const double& x, const double& y, const double& z ) {
	return 100*x + 10*y + z;
}


CASE( "test_interpolation_conservative" ) {
    Grid src_grid = localgrid( 9, 9 );
    Grid tgt_grid = localgrid( 33, 33 );
	MeshGenerator meshgen( "regular" );
	Mesh src_mesh = meshgen.generate( src_grid );
	Mesh tgt_mesh = meshgen.generate( tgt_grid );
	
	std::vector< PointLonLat > pts_ll;
	const idx_t src_nb_cells = src_mesh.cells().size();
	const auto& src_node_connectivity = src_mesh.cells().node_connectivity();
	std::vector< CSPolygon > src_csp;
	src_csp.resize( src_nb_cells );
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
	std::vector< CSPolygon > tgt_csp;
	tgt_csp.resize( tgt_nb_cells );
	const auto tgt_lonlat = array::make_view<double,2>( tgt_mesh.nodes().lonlat() );
	for( idx_t jcell = 0; jcell < tgt_nb_cells; ++jcell ) {
		const idx_t nb_nodes = tgt_node_connectivity.cols( jcell );
		pts_ll.clear();
		pts_ll.resize( nb_nodes );
		for( idx_t jnode = 0; jnode < nb_nodes; ++jnode ) {
			idx_t inode = tgt_node_connectivity( jcell, jnode );
			pts_ll[ nb_nodes - 1 - jnode ] = PointLonLat{ tgt_lonlat(inode,0), tgt_lonlat(inode,1) };
		}
		tgt_csp[ jcell ] =  CSPolygon( pts_ll );
	}

	std::vector< InterpolationParameters > interpolationParameters;
	interpolationParameters.resize( src_nb_cells );
	// brute force (!) needs to be changed
	for( idx_t scell = 0; scell < src_nb_cells; ++scell ) {
		for( idx_t tcell = 0; tcell < tgt_nb_cells; ++tcell ) {
			CSPolygon csp_i = src_csp[ scell ].intersect( tgt_csp[ tcell ] );
			if ( csp_i.area() > 0 ) {
				interpolationParameters[ scell ].cell_id.emplace_back( tcell );
				interpolationParameters[ scell ].weights.emplace_back( csp_i.area() );
				interpolationParameters[ scell ].centroids.emplace_back( csp_i.centroid() );
			}
		}
	}

	functionspace::CellColumns src_fs( src_mesh );
	functionspace::CellColumns tgt_fs( tgt_mesh );
	auto src_field = src_fs.createField< double >();
	auto tgt_field = tgt_fs.createField< double >();
	auto src_vals = array::make_view< double, 1 >( src_field );
	auto tgt_vals = array::make_view< double, 1 >( tgt_field );


	mesh::actions::build_edges( src_mesh );
	const auto& src_cell2edge = src_mesh.cells().edge_connectivity();
	const auto& src_edge2cell = src_mesh.edges().cell_connectivity();

	// assign field values on source mesh
	for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
		auto p = src_csp[ scell ].centroid();
		src_vals( scell ) = func( p[0], p[1], p[2] );
	}

	for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
		tgt_vals( tcell ) = 0.;
	}
	for ( idx_t scell = 0; scell < src_vals.size(); ++scell ) {
		auto p = src_csp[ scell ].centroid();
		// get cell neighbours
		idx_t src_nb_edges = src_cell2edge.cols( scell );
		std::vector< idx_t > src_neighbour_cells;
		src_neighbour_cells.reserve( src_nb_edges );
		for( idx_t sedge = 0; sedge < src_nb_edges; ++sedge ) {
			idx_t iedge = src_cell2edge( scell, sedge );
			idx_t cell0 = src_edge2cell( iedge, 0 );
			idx_t cell1 = src_edge2cell( iedge, 1 );
			if ( cell0 != src_cell2edge.missing_value() && cell0 != scell ) {
				src_neighbour_cells.emplace_back( cell0 );
			}
			else if ( cell1 != src_cell2edge.missing_value() ) {
				src_neighbour_cells.emplace_back( cell1 );
			}
			else {
				src_neighbour_cells.emplace_back( scell );
			}
		}
		// calculate gradient
		PointXYZ grad = { 0., 0., 0. };
		for ( idx_t nid = 0; nid < src_neighbour_cells.size(); ++nid ) {
			idx_t ncell = src_neighbour_cells[ nid ];
			idx_t nncell = src_neighbour_cells[ nid != src_neighbour_cells.size()-1 ? nid+1 : 0 ];
			if ( src_vals( ncell ) != scell && src_vals( nncell ) != scell ) {
				double coeff = src_vals( ncell );
				coeff += src_vals( nncell );
				coeff = 0.5*coeff - src_vals( scell );
				grad = PointXYZ::add( grad, PointXYZ::mul( PointXYZ::cross( src_csp[ ncell ].centroid(), src_csp[ nncell ].centroid() ), coeff ) );
			}
		}
		grad = PointXYZ::mul( grad, src_csp[ scell ].area() ); // this is WRONG we need area of Fig 2. in Kritsikis et al. (2017) -> overestimation of gradient but still conservative
		InterpolationParameters& iparam = interpolationParameters[ scell ];
		for( idx_t icell = 0; icell < iparam.centroids.size(); ++icell ) {
			// NOTE: check if this is correct barycenter of the source cell!!
			tgt_vals( iparam.cell_id[ icell ] ) += iparam.weights[ icell ] * ( src_vals( scell ) + PointXYZ::dot( grad, iparam.centroids[ icell ] - src_csp[ scell ].centroid() ) );
		}
	}
	for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
		tgt_vals( tcell ) /= tgt_csp[ tcell ].area();
	}

	// validate
	double err = 0;
	for ( idx_t tcell = 0; tcell < tgt_vals.size(); ++tcell ) {
		auto p = tgt_csp[ tcell ].centroid();
		err += std::abs( tgt_vals( tcell ) - func( p[0],p[1], p[2] ) );
	}
	err /= tgt_vals.size();
	Log::info() <<"target field err: " <<err <<"\n";
}

}  // namespace test
}  // namespace atlas


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
