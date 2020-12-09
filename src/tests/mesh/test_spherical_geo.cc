#include <iostream>
#include <cmath>
#include <vector>

#include "eckit/geometry/Sphere.h"
#include "atlas/util/Point.h"
#include "atlas/util/ConvexSphericalPolygon.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using ConvexSphericalPolygon = util::ConvexSphericalPolygon;

void overlappingPolygons( ConvexSphericalPolygon* plg1, ConvexSphericalPolygon* plg2 ) {
	std::vector< PointLonLat > p1;
	p1.emplace_back( PointLonLat( 0., 70. ) );
	p1.emplace_back( PointLonLat( 0., 60. ) );
	p1.emplace_back( PointLonLat( 10., 60. ) );
	p1.emplace_back( PointLonLat( 10, 70. ) );
	p1.emplace_back( PointLonLat( 0., 70. ) );
	std::vector< PointLonLat > p2;
	p2.emplace_back( PointLonLat( 0., 60. ) );
	p2.emplace_back( PointLonLat( 10., 60. ) );
	p2.emplace_back( PointLonLat( 10., 70. ) );
	p2.emplace_back( PointLonLat( 0., 60. ) );
	plg1 = new ConvexSphericalPolygon( p1 );
	plg2 = new ConvexSphericalPolygon( p2 );
}

CASE( "test_spherical_segment_intersection" ) {
	std::vector< PointLonLat > p1;
	p1.emplace_back( PointLonLat( 0., 70. ) );
	p1.emplace_back( PointLonLat( 0., 60. ) );
	p1.emplace_back( PointLonLat( 40., 60. ) );
	p1.emplace_back( PointLonLat( 40., 70. ) );
	p1.emplace_back( PointLonLat( 0., 70. ) );
	std::vector< PointLonLat > p2;
	p2.emplace_back( PointLonLat( 20., 70. ) );
	p2.emplace_back( PointLonLat( 0.1, 70. ) );
	p2.emplace_back( PointLonLat( 20., 50. ) );
	p2.emplace_back( PointLonLat( 40., 50. ) );
	p2.emplace_back( PointLonLat( 20., 70. ) );
	ConvexSphericalPolygon plg1( p1 );
	ConvexSphericalPolygon plg2( p2 );
	int i = 0;
	int j = 1;
	PointXYZ ip;
	bool inters = plg2.intersect( plg1[i], plg1[j], ip, 0 );
	if ( inters ) {
		PointLonLat ip_ll;
		eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ip_ll );
		std::cout <<"Intersection at " <<ip_ll <<"\n";
		std::cout <<"  with seg. [" <<plg1[i] <<", " <<plg1[j] <<"]\n";
	} else {
		std::cout <<"No intersection with seg. [" <<plg1[i] <<", " <<plg1[j] <<"]\n";
	}
}

CASE( "test_spherical_polygon_intersection" ) {
	std::vector< PointLonLat > p1;
	p1.emplace_back( PointLonLat( 0., 80. ) );
	p1.emplace_back( PointLonLat( 0., 60. ) );
	p1.emplace_back( PointLonLat( 40., 60. ) );
	p1.emplace_back( PointLonLat( 40., 80. ) );
	p1.emplace_back( PointLonLat( 0., 80. ) );
	std::vector< PointLonLat > p2;
	p2.emplace_back( PointLonLat( 20., 70. ) );
	p2.emplace_back( PointLonLat( 20., 50. ) );
	p2.emplace_back( PointLonLat( 60., 50. ) );
	p2.emplace_back( PointLonLat( 60., 70. ) );
	p2.emplace_back( PointLonLat( 20., 70. ) );
	ConvexSphericalPolygon plg1( p1 );
	ConvexSphericalPolygon plg2( p2 );
	//ConvexSphericalPolygon* ipol = plg1.intersect( plg2 );
	ConvexSphericalPolygon* ipol = plg2.intersect( plg1 );
}

//-----------------------------------------------------------------------------

} // end namespace test
} // end namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
