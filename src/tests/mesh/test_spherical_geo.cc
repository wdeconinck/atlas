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

ConvexSphericalPolygon getCSPolygon( const int np, PointLonLat* p ) {
	if ( !p ) {
        return ConvexSphericalPolygon();
	}
	std::vector< PointLonLat > pts;
	pts.reserve( np + 1 );
	for( int i = 0; i < np; i++ ) {
		pts.emplace_back( p[i] );
	}
	pts.emplace_back( p[0] );
    return ConvexSphericalPolygon( pts );
}

CASE( "test default constructor" ) {
    ConvexSphericalPolygon p;
    EXPECT( bool(p) == false );
    EXPECT( p.validate() == false );
}

CASE( "test_spherical_polygon_area" ) {
    auto plg1 = getCSPolygon( 3, new PointLonLat[3]{{0, 90}, {0, 0}, {90, 0}});
    ATLAS_ASSERT( plg1.area() == M_PI_2 );
    std::cout <<"area: " <<plg1.area() <<"\n";
    auto plg2 = getCSPolygon( 4, new PointLonLat[4]{{0, 45}, {0, 0}, {90,0},{90, 45}});
    std::cout <<"area: " <<plg2.area() <<"\n";

    auto plg3 = getCSPolygon( 3, new PointLonLat[3]{{0, 90}, {0, 45}, {90,45}} );
    std::cout <<"area: " <<plg3.area() <<"\n";
    std::cout <<"area diff: " <<plg1.area() - plg2.area() - plg3.area() <<"\n";
	std::cout.flush();
    ATLAS_ASSERT( std::abs( plg1.area() - plg2.area() - plg3.area() ) < 1e-15 );
}

CASE( "test_spherical_polygon_intersection" ) {
    constexpr int nplg_f = 2;
    constexpr int nplg_g = 17;
    constexpr int nplg_i = nplg_f * nplg_g;
    std::array<ConvexSphericalPolygon,nplg_f> plg_f = {
		getCSPolygon( 4, new PointLonLat[4]{{0, 70}, {0, 60}, {40, 60}, {40,70}} ),
		getCSPolygon( 3, new PointLonLat[3]{{0, 90}, {0, 0}, {40, 0}} )
	};
    std::array<ConvexSphericalPolygon,nplg_g> plg_g = {
		getCSPolygon( 3, new PointLonLat[3]{{0, 60}, {0, 50}, {40, 60}} ),//0
		getCSPolygon( 3, new PointLonLat[3]{{0, 60}, {0, 50}, {20, 60}} ),
		getCSPolygon( 3, new PointLonLat[3]{{10, 60}, {10, 50}, {30, 60}} ),//3
		getCSPolygon( 3, new PointLonLat[3]{{40, 80}, {0, 60}, {40, 60}} ),
		getCSPolygon( 3, new PointLonLat[3]{{0, 80}, {0, 60}, {40, 60}} ),//5
		getCSPolygon( 3, new PointLonLat[3]{{20, 80}, {0, 60}, {40, 60}} ),
		getCSPolygon( 3, new PointLonLat[3]{{20, 70}, {0, 50}, {40, 50}} ),//7
		getCSPolygon( 3, new PointLonLat[3]{{0, 90}, {0, 60}, {40, 60}} ),
		getCSPolygon( 3, new PointLonLat[3]{{-10, 80}, {-10, 50}, {50, 80}} ),//9
		getCSPolygon( 4, new PointLonLat[4]{{0, 80}, {0, 50}, {40, 50}, {40, 80}} ),
		getCSPolygon( 4, new PointLonLat[4]{{0, 65}, {20, 55}, {40, 60}, {20, 65}} ),//11
		getCSPolygon( 4, new PointLonLat[4]{{20, 65}, {0, 60}, {20, 55}, {40, 60}} ),
		getCSPolygon( 4, new PointLonLat[4]{{10, 63}, {20, 55}, {30, 63}, {20, 65}} ),//13
		getCSPolygon( 6, new PointLonLat[6]{{20, 75}, {0, 70}, {5, 5}, {10, 0}, {20, 0}, {40, 70}} ),
		getCSPolygon( 3, new PointLonLat[3]{{0, 50}, {0, 40}, {5, 45}} ),//15
		getCSPolygon( 4, new PointLonLat[4]{{0, 90}, {0, 80}, {20, 0}, {40, 80}} ),
		getCSPolygon( 4, new PointLonLat[4]{{0, 65}, {0, 55}, {40, 65}, {40, 75}} ),//17
	};
    std::array<ConvexSphericalPolygon,nplg_i> plg_i = {
		getCSPolygon( 2, new PointLonLat[2]{{0, 60}, {40, 60}} ),//0
		getCSPolygon( 0, nullptr ),
		getCSPolygon( 0, nullptr ),//2
		getCSPolygon( 4, new PointLonLat[4]{{0, 60}, {40, 60}, {40, 70}, {10, 70.8}} ),
		getCSPolygon( 4, new PointLonLat[4]{{0, 70}, {0, 60}, {40, 60}, {30, 70.8}} ),//4
		getCSPolygon( 4, new PointLonLat[4]{{0, 60}, {40, 60}, {34.6, 70.5}, {5.3, 70.5}} ),
		getCSPolygon( 3, new PointLonLat[3]{{7.5, 60.9}, {32.5, 60.9}, {20, 70}} ),//6
		getCSPolygon( 4, new PointLonLat[4]{{0, 70}, {0, 60}, {40, 60}, {40, 70}} ),
		getCSPolygon( 3, new PointLonLat[3]{{0, 65.5}, {6.9, 70.6}, {0, 70}} ),//8
		getCSPolygon( 4, new PointLonLat[4]{{0, 70}, {0, 60}, {40, 60}, {40, 70}} ),
		getCSPolygon( 4, new PointLonLat[4]{{0, 65}, {10, 61}, {40, 60}, {20, 65}} ),//10
		getCSPolygon( 3, new PointLonLat[3]{{0, 60}, {40, 60}, {20, 65}} ),
		getCSPolygon( 5, new PointLonLat[5]{{12.6, 61}, {27.4, 61.3}, {30, 63}, {20, 65}, {10, 63}} ),//12
		getCSPolygon( 4, new PointLonLat[4]{{0, 70}, {1.9, 60.3}, {32.9, 60.9}, {40, 70}} ),
		getCSPolygon( 3, new PointLonLat[3]{{0, 50}, {0, 40}, {5, 45}} ),//14
		getCSPolygon( 4, new PointLonLat[4]{{13.7, 61.4}, {26.3, 61.4}, {30, 70.8}, {10, 70.8}} ),
		getCSPolygon( 6, new PointLonLat[6]{{0, 65}, {0, 60}, {16.8, 61.5}, {40, 65}, {40, 70}, {15, 71.1}} ),//16
		getCSPolygon( 3, new PointLonLat[3]{{0, 60}, {0, 50}, {40, 60}} ),
		getCSPolygon( 3, new PointLonLat[3]{{0, 60}, {0, 50}, {20, 60}} ),//18
		getCSPolygon( 3, new PointLonLat[3]{{10, 60}, {10, 50}, {30, 60}} ),
		getCSPolygon( 3, new PointLonLat[3]{{0, 60}, {40, 60}, {40, 80}} ),//20
		getCSPolygon( 3, new PointLonLat[3]{{0, 80}, {0, 60}, {40, 60}} ),
		getCSPolygon( 3, new PointLonLat[3]{{0, 60}, {40, 60}, {20, 80}} ),//22
		getCSPolygon( 3, new PointLonLat[3]{{0, 50}, {40, 50}, {20, 70}} ),
		getCSPolygon( 3, new PointLonLat[3]{{0, 90}, {0, 60}, {40, 60}} ),//24
		getCSPolygon( 4, new PointLonLat[4]{{0, 65.5}, {40, 79.2}, {40, 80.8}, {0, 80.8}} ),
		getCSPolygon( 4, new PointLonLat[4]{{0, 80}, {0, 50}, {40, 50}, {40, 80}} ),//26
		getCSPolygon( 4, new PointLonLat[4]{{0, 65}, {20, 55}, {40, 60}, {20, 65}} ),
		getCSPolygon( 4, new PointLonLat[4]{{0, 60}, {20, 55}, {40, 60}, {20, 65}} ),//28
		getCSPolygon( 4, new PointLonLat[4]{{10, 63}, {20, 55}, {30, 63}, {20, 65}} ),
		getCSPolygon( 6, new PointLonLat[6]{{0, 70}, {5, 5}, {10, 0}, {20, 0}, {40, 70}, {20,75}} ),//30
		getCSPolygon( 3, new PointLonLat[3]{{0, 50}, {0, 40}, {5, 45}} ),
		getCSPolygon( 4, new PointLonLat[4]{{0, 90}, {0, 80}, {20, 0}, {40, 80}} ),//32
		getCSPolygon( 4, new PointLonLat[4]{{0, 65}, {0, 55}, {40, 65}, {40, 75}} )
	};
	for( int i = 0; i < nplg_f; i++ ) {
		for( int j = 0; j < nplg_g; j++ ) {
			std::cout <<"\n\n("<<i*nplg_g+j <<") Intersecting polygon\n    ";
            plg_f[i].print( std::cout );
            std::cout <<"\nof area: " <<plg_f[i].area() <<", convex: " <<plg_f[i].validate();
			std::cout <<"\nwith polygon\n    ";
            plg_g[j].print( std::cout );
            std::cout <<"\nof area: " <<plg_g[j].area() <<", convex: " <<plg_g[j].validate();
            auto plg_fg = plg_f[i].intersect( plg_g[j] );
            auto plg_gf = plg_g[j].intersect( plg_f[i] );
			std::cout <<"\ngot polygon\n    ";
			if ( plg_fg ) {
                plg_fg.print( std::cout );
                std::cout <<"\nof area: " <<plg_fg.area()
                    <<", convex: " <<plg_g[j].validate();
                ATLAS_ASSERT( plg_fg.equals( plg_gf ) );
                ATLAS_ASSERT( plg_fg.equals( plg_i[ i*nplg_g + j ], 0.5 ) );
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
