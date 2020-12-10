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


ConvexSphericalPolygon& getCSPolygon( const int np, PointLonLat* p ) {
	std::vector< PointLonLat > pts;
	for( int i = 0; i < np; i++ ) {
		pts.emplace_back( p[i] );
	}
	pts.emplace_back( p[0] );
	return *(new ConvexSphericalPolygon( pts ));
}

CASE( "test_spherical_segment_intersection" ) {
	ConvexSphericalPolygon plg[2] = { 
		getCSPolygon( 4, new PointLonLat[4]{{0, 80}, {0, 60}, {40, 60}, {40,80}} ),
		getCSPolygon( 2, new PointLonLat[2]{{60, 50}, {60, 70}} )
	};
	PointXYZ ip;
	int inters = -1 + plg[0].intersect( plg[1][0], plg[1][1], ip, 0 );
	if ( inters != -1 ) {
		PointLonLat ip_ll;
		eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ip_ll );
		std::cout <<"Intersection at " <<ip_ll <<" on edge " <<inters <<"\n";
	} else {
		std::cout <<"No intersection with seg. [" <<plg[1][0] <<", " <<plg[1][1] <<"]\n";
	}
}

CASE( "test_spherical_polygon_intersection" ) {
	const int nplg_f = 2;
	const int nplg_g = 17;
	ConvexSphericalPolygon plg_f[ nplg_f ] = { 
		getCSPolygon( 4, new PointLonLat[4]{{0, 70}, {0, 60}, {40, 60}, {40,70}} ),
		getCSPolygon( 3, new PointLonLat[3]{{0, 90}, {0, 0}, {40, 0}} )
	};
	ConvexSphericalPolygon plg_g[ nplg_g ] = {
		getCSPolygon( 3, new PointLonLat[3]{{0, 60}, {0, 50}, {40, 60}} ), // 1
		getCSPolygon( 3, new PointLonLat[3]{{0, 60}, {0, 50}, {20, 60}} ), // 
		getCSPolygon( 3, new PointLonLat[3]{{10, 60}, {10, 50}, {30, 60}} ), // 3
		getCSPolygon( 3, new PointLonLat[3]{{40, 80}, {0, 60}, {40, 60}} ),
		getCSPolygon( 3, new PointLonLat[3]{{0, 80}, {0, 60}, {40, 60}} ), // 5
		getCSPolygon( 3, new PointLonLat[3]{{20, 80}, {0, 60}, {40, 60}} ),
		getCSPolygon( 3, new PointLonLat[3]{{20, 70}, {0, 50}, {40, 50}} ), // 7
		getCSPolygon( 3, new PointLonLat[3]{{0, 90}, {0, 60}, {40, 60}} ),
		getCSPolygon( 3, new PointLonLat[3]{{-10, 80}, {-10, 50}, {50, 80}} ), // 9
		getCSPolygon( 4, new PointLonLat[4]{{0, 80}, {0, 50}, {40, 50}, {40, 80}} ),
		getCSPolygon( 4, new PointLonLat[4]{{0, 65}, {20, 55}, {40, 60}, {20, 65}} ), // 11
		getCSPolygon( 4, new PointLonLat[4]{{20, 65}, {0, 60}, {20, 55}, {40, 60}} ),
		getCSPolygon( 4, new PointLonLat[4]{{10, 63}, {20, 55}, {30, 63}, {20, 65}} ), // 13
		getCSPolygon( 6, new PointLonLat[6]{{20, 75}, {0, 70}, {5, 65}, {10, 0}, {20, 0}, {40, 70}} ),
		getCSPolygon( 3, new PointLonLat[3]{{0, 50}, {0, 40}, {5, 45}} ), // 15
		getCSPolygon( 4, new PointLonLat[4]{{0, 90}, {0, 80}, {20, 0}, {40, 80}} ),
		getCSPolygon( 4, new PointLonLat[4]{{0, 65}, {0, 55}, {40, 65}, {40, 75}} ), // 17
	};
	for( int i = 0; i < nplg_f; i++ ) {
		for( int j = 0; j < nplg_g; j++ ) {
			std::cout <<"\n\n("<<i*nplg_g+j <<") Intersecting polygon\n    ";
			plg_f[i].print( std::cout );
			std::cout <<"\nwith polygon\n    ";
			plg_g[j].print( std::cout );
			ConvexSphericalPolygon* plg_fg = plg_f[i].intersect( plg_g[j] );
			ConvexSphericalPolygon* plg_gf = plg_g[j].intersect( plg_f[i] );
			std::cout <<"\ngot polygon\n    ";
			if ( plg_fg ) {
				plg_fg->print( std::cout );
			}
			else {
				std::cout <<"	empty";
			}
		}
	}
	std::cout <<"\n";
}

//-----------------------------------------------------------------------------

} // end namespace test
} // end namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
