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
	// return 0:outside, -1:on_edge, 1:strictly_inside
	inline int inside( const Point &p ) const {
		for ( int i = 0; i < nv; i++ ) {
			const Point& cp = v[i].cross( v[(i!=nv-1 ? i+1 : 0)] );
			double dp = cp.dot( p );
			if (  dp < -eps_point() ) {
				return 0;
			} 
			else if ( abs(dp) < eps_point() ) {
				return -1;
			}
		}
		return 1;
	}

	Point* nextIntersection() const {
		if ( !ip ) {
			lv.emplace_back( *ip );
		}
		Point* ip;
		do {
			// find next vertex
			for( int i = 0; i < nv; i++ ) {
				ip = this->intersect( p.v[i], p.v[(i!=p.nv-1 ? i+1 : 0)] );
				if ( ip ) {
					lv.emplace_back( *ip );
					cout << "  found lv emplace back: " <<*ip <<"\n";
					break;
				}
			}
		} while ( ip && (*ip != lv[0]) );
	}

	// intersect first? segment [p1,p2] with this polygon=[v0,v1,...,v_nv]
	// first? in the order [v0,v1], [v1,v2],...
	// return type 0:no_intersection, 1:stricly_inside, 2:on_edge
	int intersect( const Point& p1, const Point& p2, Point& p3 ) const {
		const Pol& tp = *this;
		Point& sp1 = p1.cross( p2 );
		for( int i = 0; i < tp.nv; i++ ) {
			int ip = ( i != tp.nv-1 ? i+1 : 0 );
			Point& tpi = tp.v[i].cross( tp.v[ip] );
			Point& p3 = sp1.cross( tpi );
			//p3 /= sqrt( v_insc.dot( p3 ) );
			cout <<" p3 = " <<p3 <<"\n";
			//type = p3.onSegment( p1, p2 )
			//if ( type ) {
		//		return type;
	//		}
		}
		return 0;
	}

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
	p1.emplace_back( PointLonLat( 10., 0. ) );
	p1.emplace_back( PointLonLat( 30., 0. ) );
	p1.emplace_back( PointLonLat( 30., 40. ) );
	p1.emplace_back( PointLonLat( 10., 20. ) );
	std::vector< PointLonLat > p2;
	p2.emplace_back( PointLonLat( 20., 20. ) );
	p2.emplace_back( PointLonLat( 20., 0. ) );
	p2.emplace_back( PointLonLat( 40., 20. ) );
	p2.emplace_back( PointLonLat( 40., 40. ) );

	atlas::PointXYZ sp1[4];
	atlas::PointXYZ sp2[4];
	eckit::geometry::Sphere sphere;
	for( int i=0; i<4; i++ ) {
		sphere.convertSphericalToCartesian( 1., p1[i], sp1[i] );
		sphere.convertSphericalToCartesian( 1., p2[i], sp2[i] );
	}

	util::ConvexSphericalPolygon pol1( p1 );
	util::ConvexSphericalPolygon pol2( p2 );

	//Pol* pol_insc = pol1.intersect( pol2 );
	//pol_insc->info();
}

//-----------------------------------------------------------------------------

} // end namespace test
} // end namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
