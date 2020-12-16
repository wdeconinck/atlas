#include <cmath>
#include <iostream>
#include <vector>

#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/Point.h"
#include "eckit/geometry/Sphere.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using ConvexSphericalPolygon = util::ConvexSphericalPolygon;

ConvexSphericalPolygon getCSPolygon( std::initializer_list<PointLonLat> list ) {
    if ( list.size() == 0 ) {
        return ConvexSphericalPolygon();
    }
    std::vector<PointLonLat> pts;
    pts.reserve( list.size() + 1 );
    for ( auto& p : list ) {
        pts.emplace_back( p );
    }
    return ConvexSphericalPolygon( pts );
}


CASE( "test default constructor" ) {
    ConvexSphericalPolygon p;
    EXPECT( bool( p ) == false );
    EXPECT( p.validate() == false );
}

CASE( "test_spherical_polygon_area" ) {
    auto plg1 = getCSPolygon( {{0, 90}, {0, 0}, {90, 0}} );
    EXPECT_APPROX_EQ( plg1.area(), M_PI_2 );
    Log::info() << "area: " << plg1.area() << "\n";
    auto plg2 = getCSPolygon( {{0, 45}, {0, 0}, {90, 0}, {90, 45}} );
    Log::info() << "area: " << plg2.area() << "\n";

    auto plg3 = getCSPolygon( {{0, 90}, {0, 45}, {90, 45}} );
    Log::info() << "area: " << plg3.area() << "\n";
    Log::info() << "area diff: " << plg1.area() - plg2.area() - plg3.area() << std::endl;
    EXPECT_APPROX_EQ( std::abs( plg1.area() - plg2.area() - plg3.area() ), 0, 1e-15 );
}

CASE( "test_spherical_polygon_intersection" ) {
    constexpr int nplg_f                             = 2;
    constexpr int nplg_g                             = 17;
    constexpr int nplg_i                             = nplg_f * nplg_g;
    std::array<ConvexSphericalPolygon, nplg_f> plg_f = {getCSPolygon( {{0, 70}, {0, 60}, {40, 60}, {40, 70}} ),
                                                        getCSPolygon( {{0, 90}, {0, 0}, {40, 0}} )};
    std::array<ConvexSphericalPolygon, nplg_g> plg_g = {
        getCSPolygon( {{0, 60}, {0, 50}, {40, 60}} ),  //0
        getCSPolygon( {{0, 60}, {0, 50}, {20, 60}} ),
        getCSPolygon( {{10, 60}, {10, 50}, {30, 60}} ),  //3
        getCSPolygon( {{40, 80}, {0, 60}, {40, 60}} ),
        getCSPolygon( {{0, 80}, {0, 60}, {40, 60}} ),  //5
        getCSPolygon( {{20, 80}, {0, 60}, {40, 60}} ),
        getCSPolygon( {{20, 70}, {0, 50}, {40, 50}} ),  //7
        getCSPolygon( {{0, 90}, {0, 60}, {40, 60}} ),
        getCSPolygon( {{-10, 80}, {-10, 50}, {50, 80}} ),  //9
        getCSPolygon( {{0, 80}, {0, 50}, {40, 50}, {40, 80}} ),
        getCSPolygon( {{0, 65}, {20, 55}, {40, 60}, {20, 65}} ),  //11
        getCSPolygon( {{20, 65}, {0, 60}, {20, 55}, {40, 60}} ),
        getCSPolygon( {{10, 63}, {20, 55}, {30, 63}, {20, 65}} ),  //13
        getCSPolygon( {{20, 75}, {0, 70}, {5, 5}, {10, 0}, {20, 0}, {40, 70}} ),
        getCSPolygon( {{0, 50}, {0, 40}, {5, 45}} ),  //15
        getCSPolygon( {{0, 90}, {0, 80}, {20, 0}, {40, 80}} ),
        getCSPolygon( {{0, 65}, {0, 55}, {40, 65}, {40, 75}} ),  //17
    };
    std::array<ConvexSphericalPolygon, nplg_i> plg_i = {
        getCSPolygon( {{0, 60}, {40, 60}} ),  //0
        getCSPolygon( {} ),
        getCSPolygon( {} ),  //2
        getCSPolygon( {{0, 60}, {40, 60}, {40, 70}, {10, 70.8}} ),
        getCSPolygon( {{0, 70}, {0, 60}, {40, 60}, {30, 70.8}} ),  //4
        getCSPolygon( {{0, 60}, {40, 60}, {34.6, 70.5}, {5.3, 70.5}} ),
        getCSPolygon( {{7.5, 60.9}, {32.5, 60.9}, {20, 70}} ),  //6
        getCSPolygon( {{0, 70}, {0, 60}, {40, 60}, {40, 70}} ),
        getCSPolygon( {{0, 65.5}, {6.9, 70.6}, {0, 70}} ),  //8
        getCSPolygon( {{0, 70}, {0, 60}, {40, 60}, {40, 70}} ),
        getCSPolygon( {{0, 65}, {10, 61}, {40, 60}, {20, 65}} ),  //10
        getCSPolygon( {{0, 60}, {40, 60}, {20, 65}} ),
        getCSPolygon( {{12.6, 61}, {27.4, 61.3}, {30, 63}, {20, 65}, {10, 63}} ),  //12
        getCSPolygon( {{0, 70}, {1.9, 60.3}, {32.9, 60.9}, {40, 70}} ),
        getCSPolygon( {{0, 50}, {0, 40}, {5, 45}} ),  //14
        getCSPolygon( {{13.7, 61.4}, {26.3, 61.4}, {30, 70.8}, {10, 70.8}} ),
        getCSPolygon( {{0, 65}, {0, 60}, {16.8, 61.5}, {40, 65}, {40, 70}, {15, 71.1}} ),  //16
        getCSPolygon( {{0, 60}, {0, 50}, {40, 60}} ),
        getCSPolygon( {{0, 60}, {0, 50}, {20, 60}} ),  //18
        getCSPolygon( {{10, 60}, {10, 50}, {30, 60}} ),
        getCSPolygon( {{0, 60}, {40, 60}, {40, 80}} ),  //20
        getCSPolygon( {{0, 80}, {0, 60}, {40, 60}} ),
        getCSPolygon( {{0, 60}, {40, 60}, {20, 80}} ),  //22
        getCSPolygon( {{0, 50}, {40, 50}, {20, 70}} ),
        getCSPolygon( {{0, 90}, {0, 60}, {40, 60}} ),  //24
        getCSPolygon( {{0, 65.5}, {40, 79.2}, {40, 80.8}, {0, 80.8}} ),
        getCSPolygon( {{0, 80}, {0, 50}, {40, 50}, {40, 80}} ),  //26
        getCSPolygon( {{0, 65}, {20, 55}, {40, 60}, {20, 65}} ),
        getCSPolygon( {{0, 60}, {20, 55}, {40, 60}, {20, 65}} ),  //28
        getCSPolygon( {{10, 63}, {20, 55}, {30, 63}, {20, 65}} ),
        getCSPolygon( {{0, 70}, {5, 5}, {10, 0}, {20, 0}, {40, 70}, {20, 75}} ),  //30
        getCSPolygon( {{0, 50}, {0, 40}, {5, 45}} ),
        getCSPolygon( {{0, 90}, {0, 80}, {20, 0}, {40, 80}} ),  //32
        getCSPolygon( {{0, 65}, {0, 55}, {40, 65}, {40, 75}} )};
    for ( int i = 0; i < nplg_f; i++ ) {
        for ( int j = 0; j < nplg_g; j++ ) {
            Log::info() << "\n(" << i * nplg_g + j << ") Intersecting polygon\n    " << plg_f[i] << std::endl;
            Log::info() << "of area: " << plg_f[i].area() << ", convex: " << plg_f[i].validate() << std::endl;
            Log::info() << "with polygon\n    " << plg_g[j] << std::endl;
            Log::info() << "of area: " << plg_g[j].area() << ", convex: " << plg_g[j].validate() << std::endl;
            auto plg_fg = plg_f[i].intersect( plg_g[j] );
            auto plg_gf = plg_g[j].intersect( plg_f[i] );
            Log::info() << "got polygon\n    ";
            if ( plg_fg ) {
                Log::info() << plg_fg << std::endl;
                Log::info() << "of area: " << plg_fg.area() << ", convex: " << plg_g[j].validate() << std::endl;
                EXPECT( plg_fg.equals( plg_gf ) );
                EXPECT( plg_fg.equals( plg_i[i * nplg_g + j], 0.5 ) );
            }
            else {
                Log::info() << "	empty" << std::endl;
            }
        }
    }
}

CASE( "Size of ConvexSphericalPolygon" ) {
    // This test illustrates that ConvexSphericalPolygon is allocated on the stack completely,
    // as sizeof(ConvexSphericalPolygon) includes space for MAX_SIZE coordinates of type PointXYZ
    EXPECT( sizeof( PointXYZ ) == sizeof( double ) * 3 );
    size_t expected_size = 0;
    expected_size += ConvexSphericalPolygon::MAX_SIZE * sizeof( PointXYZ );
    expected_size += sizeof( size_t );
    expected_size += sizeof( bool );
    EXPECT( sizeof( ConvexSphericalPolygon ) >= expected_size );  // greater because compiler may add some padding
}


//-----------------------------------------------------------------------------

}  // end namespace test
}  // end namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
