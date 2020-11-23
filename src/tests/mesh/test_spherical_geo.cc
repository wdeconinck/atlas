#include <iostream>
#include <cmath>
#include <vector>

#include "eckit/geometry/Sphere.h"
#include "atlas/util/Point.h"
#include "atlas/util/ConvexSphericalPolygon.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

static double eps_point() { return 1e-32; }
static double deg2rad() { return M_PI / 180.; }
static double rad2deg() { return 180. * M_1_PI; }

/*
	// return the intersecting polygon of "this" and "p"
	Pol* intersect( const Pol& pol ) const {
		std::vector< Point > lv;
		// is there a vertex of "pol" inside "this"?
		int i1 = 0;
		int inside = 0;
		for ( ; i1 < pol.nv; i1++ ) {
			if ( this->inside( pol.v[i1] ) ) {
				inside = 1;
				std::cout << "\nPoint inside: " << i1 <<"\n";
				break;
			}
		}
		if ( inside ) {
			lv.emplace_back( pol.v[i1] );
			std::cout << "  inside lv emplace back: " <<pol.v[i1] <<"\n";
		}
		else {
			// no vertex of "p" inside "this", loop over edges of p 
			// and find one intersecting polygon "this"
			// if no intersection found return null-polygon
			Point* ip1;
			for( int i = 0; i < pol.nv; i++ ) {
				//ip1 = this->intersect( pol.v[i], pol.v[(i!=pol.nv-1 ? i+1 : 0)] );
				if ( ip1 ) {
					lv.emplace_back( *ip1 );
					std::cout << "  lv emplace back: " <<*ip1 <<"\n";
					i1 = i;
					break;
				}
			}
			if ( ! ip1 ) {
				cout <<"Polygons do not intersect\n";
				return nullptr;
			}
		}
		// find remaining vertices of the intersection polygon
		cout <<" Find remaining points...\n";
		for( Point& p = lv[0]; p != lv[0]; ) {
			//Point* ip = this->nextIntersection( p, pol );
			//p = *ip;
		}

		Pol* rp = new Pol( lv.size() );
		rp->set( &lv[0] );
		return rp;	
	}
*/


CASE( "test_spherical_polygon_intersection" ) {
	std::vector< PointLonLat > p1;
	p1.emplace_back( PointLonLat( 0., 80. ) );
	p1.emplace_back( PointLonLat( 0., 60. ) );
	p1.emplace_back( PointLonLat( 40., 60. ) );
	p1.emplace_back( PointLonLat( 20., 80. ) );
	p1.emplace_back( PointLonLat( 0., 80. ) );
	std::vector< PointLonLat > p2;
	p2.emplace_back( PointLonLat( 20., 70. ) );
	p2.emplace_back( PointLonLat( 0., 70. ) );
	p2.emplace_back( PointLonLat( 20., 50. ) );
	p2.emplace_back( PointLonLat( 40., 50. ) );
	p2.emplace_back( PointLonLat( 20., 70. ) );

	atlas::PointXYZ sp1[5];
	atlas::PointXYZ sp2[5];
	eckit::geometry::Sphere sphere;
	for( int i=0; i<5; i++ ) {
		sphere.convertSphericalToCartesian( 1., p1[i], sp1[i] );
		sphere.convertSphericalToCartesian( 1., p2[i], sp2[i] );
	}

	util::ConvexSphericalPolygon pol1( p1 );
	util::ConvexSphericalPolygon pol2( p2 );

	int i = 1;
	int j = 2;
	PointXYZ ip;
	bool inters = pol2.intersect( sp1[i], sp1[j], ip, 2 );
	PointLonLat ip_ll;
  	sphere.convertCartesianToSpherical( 1., ip, ip_ll );
	if ( inters ) {
		std::cout <<"Intersection at " <<ip_ll <<"\n";
		std::cout <<"  with seg. [" <<p1[i] <<", " <<p1[j] <<"]\n";
	} else {
		std::cout <<"No intersection" <<"\n";
		std::cout <<"  with seg. [" <<p1[i] <<", " <<p1[j] <<"]\n";
	}

	//Pol* pol_insc = pol1.intersect( pol2 );
	//pol_insc->info();
}

//-----------------------------------------------------------------------------

} // end namespace test
} // end namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
