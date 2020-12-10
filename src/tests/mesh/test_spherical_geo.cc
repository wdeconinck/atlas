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

CASE( "test_spherical_segment_intersection" ) {
	std::vector< PointLonLat > p;
	p.emplace_back( PointLonLat( 0., 80. ) );
	p.emplace_back( PointLonLat( 0., 60. ) );
	p.emplace_back( PointLonLat( 40., 60. ) );
	p.emplace_back( PointLonLat( 40., 80. ) );
	p.emplace_back( PointLonLat( 0., 80. ) );
	ConvexSphericalPolygon plg( p );
	std::vector< PointLonLat > pp;
	pp.emplace_back( PointLonLat( 60., 50. ) );
	pp.emplace_back( PointLonLat( 60., 70. ) );
	pp.emplace_back( PointLonLat( 60., 50. ) );
	ConvexSphericalPolygon pplg( pp );
	PointXYZ ip;
	int inters = -1 + plg.intersect( pplg[0], pplg[1], ip, 0 );
	if ( inters != -1 ) {
		PointLonLat ip_ll;
		eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ip_ll );
		std::cout <<"Intersection at " <<ip_ll <<" on edge " <<inters <<"\n";
	} else {
		std::cout <<"No intersection with seg. [" <<pplg[0] <<", " <<pplg[1] <<"]\n";
	}
}

CASE( "test_spherical_polygon_intersection" ) {
	std::vector< PointLonLat > p1;
	p1.emplace_back( PointLonLat( 0., 70. ) );
	p1.emplace_back( PointLonLat( 0., 60. ) );
	p1.emplace_back( PointLonLat( 40., 60. ) );
	p1.emplace_back( PointLonLat( 40., 70. ) );
	p1.emplace_back( PointLonLat( 0., 70. ) );
	std::vector< PointLonLat > p1a;
	p1a.emplace_back( PointLonLat( 0., 90. ) );
	p1a.emplace_back( PointLonLat( 0., 0. ) );
	p1a.emplace_back( PointLonLat( 40., 0. ) );
	p1a.emplace_back( PointLonLat( 0., 90. ) );
	std::vector< PointLonLat > p2;
#if 0
	// see Test 1)
	p2.emplace_back( PointLonLat( 0., 60. ) );
	p2.emplace_back( PointLonLat( 0., 50. ) );
	p2.emplace_back( PointLonLat( 40., 60. ) );
	p2.emplace_back( PointLonLat( 0., 60. ) );
#elif 0
	// see Test 2)
	p2.emplace_back( PointLonLat( 0., 60. ) );
	p2.emplace_back( PointLonLat( 0., 50. ) );
	p2.emplace_back( PointLonLat( 20., 60. ) );
	p2.emplace_back( PointLonLat( 0., 60. ) );
#elif 0
	// see Test 3)
	p2.emplace_back( PointLonLat( 10., 60. ) );
	p2.emplace_back( PointLonLat( 10., 50. ) );
	p2.emplace_back( PointLonLat( 30., 60. ) );
	p2.emplace_back( PointLonLat( 10., 60. ) );
#elif 0
	// see Test 4)
	p2.emplace_back( PointLonLat( 40., 80. ) );
	p2.emplace_back( PointLonLat( 0., 60. ) );
	p2.emplace_back( PointLonLat( 40., 60. ) );
	p2.emplace_back( PointLonLat( 40., 80. ) );
#elif 0
	// see Test 5)
	p2.emplace_back( PointLonLat( 0., 80. ) );
	p2.emplace_back( PointLonLat( 0., 60. ) );
	p2.emplace_back( PointLonLat( 40., 60. ) );
	p2.emplace_back( PointLonLat( 0., 80. ) );
#elif 0
	// see Test 6)
	p2.emplace_back( PointLonLat( 20., 80. ) );
	p2.emplace_back( PointLonLat( 0., 60. ) );
	p2.emplace_back( PointLonLat( 40., 60. ) );
	p2.emplace_back( PointLonLat( 20., 80. ) );
#elif 0
	// see Test 7)
	p2.emplace_back( PointLonLat( 20., 70. ) );
	p2.emplace_back( PointLonLat( 0., 50. ) );
	p2.emplace_back( PointLonLat( 40., 50. ) );
	p2.emplace_back( PointLonLat( 20., 70. ) );
#elif 0
	// see Test 8)
	p2.emplace_back( PointLonLat( 0., 90. ) );
	p2.emplace_back( PointLonLat( 0., 60. ) );
	p2.emplace_back( PointLonLat( 40., 60. ) );
	p2.emplace_back( PointLonLat( 0., 90. ) );
#elif 0
	// see Test 9)
	p2.emplace_back( PointLonLat( -10., 80. ) );
	p2.emplace_back( PointLonLat( -10., 50. ) );
	p2.emplace_back( PointLonLat( 50., 50. ) );
	p2.emplace_back( PointLonLat( 50., 80. ) );
	p2.emplace_back( PointLonLat( -10., 80. ) );
#elif 0
	// see Test 10)
	p2.emplace_back( PointLonLat( 0., 80. ) );
	p2.emplace_back( PointLonLat( 0., 50. ) );
	p2.emplace_back( PointLonLat( 40., 50. ) );
	p2.emplace_back( PointLonLat( 40., 80. ) );
	p2.emplace_back( PointLonLat( 0., 80. ) );
#elif 1
	// see Test 11)
	p2.emplace_back( PointLonLat( 0., 65. ) );
	p2.emplace_back( PointLonLat( 20., 55. ) );
	p2.emplace_back( PointLonLat( 40., 60. ) );
	p2.emplace_back( PointLonLat( 20., 65. ) );
	p2.emplace_back( PointLonLat( 0., 65. ) );
#elif 1
	// see Test 12)
	p2.emplace_back( PointLonLat( 20., 65. ) );
	p2.emplace_back( PointLonLat( 0., 60. ) );
	p2.emplace_back( PointLonLat( 20., 55. ) );
	p2.emplace_back( PointLonLat( 40., 60. ) );
	p2.emplace_back( PointLonLat( 20., 65. ) );
#elif 1
	// see Test 13)
	p2.emplace_back( PointLonLat( 10., 63. ) );
	p2.emplace_back( PointLonLat( 20., 55. ) );
	p2.emplace_back( PointLonLat( 30., 63. ) );
	p2.emplace_back( PointLonLat( 20., 65. ) );
	p2.emplace_back( PointLonLat( 10., 63. ) );
#elif 1
	// see Test 14)
	p2.emplace_back( PointLonLat( 20., 75. ) );
	p2.emplace_back( PointLonLat( 0., 70. ) );
	p2.emplace_back( PointLonLat( 5., 65. ) );
	p2.emplace_back( PointLonLat( 10., 0. ) );
	p2.emplace_back( PointLonLat( 20., 0. ) );
	p2.emplace_back( PointLonLat( 40., 70. ) );
	p2.emplace_back( PointLonLat( 20., 75. ) );
#elif 1
	// see Test 15)
	p2.emplace_back( PointLonLat( 0., 50. ) );
	p2.emplace_back( PointLonLat( 0., 40. ) );
	p2.emplace_back( PointLonLat( 5., 45. ) );
	p2.emplace_back( PointLonLat( 0., 50. ) );
#elif 1
	// see Test 16)
	p2.emplace_back( PointLonLat( 0., 90. ) );
	p2.emplace_back( PointLonLat( 0., 80. ) );
	p2.emplace_back( PointLonLat( 20., 0. ) );
	p2.emplace_back( PointLonLat( 40., 80. ) );
	p2.emplace_back( PointLonLat( 0., 90. ) );
#endif
	ConvexSphericalPolygon plg1( p1 );
	ConvexSphericalPolygon plg1a( p1a );
	ConvexSphericalPolygon plg2( p2 );
	std::cout <<"\n ======== 0\n\n";
	ConvexSphericalPolygon* ipol1 = plg1.intersect( plg2 );
	std::cout <<"\n ======== 1\n\n";
	ConvexSphericalPolygon* ipol2 = plg2.intersect( plg1 );
	std::cout <<"\n ======== 2\n\n";
	ConvexSphericalPolygon* ipol1a = plg1a.intersect( plg2 );
	std::cout <<"\n ======== 3\n\n";
	ConvexSphericalPolygon* ipol2a = plg2.intersect( plg1a );
}

//-----------------------------------------------------------------------------

} // end namespace test
} // end namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
