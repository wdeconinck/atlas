/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "eckit/geometry/Sphere.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/runtime/Exception.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/NormaliseLongitude.h"

#define DEBUG_OUTPUT 0
#define DEBUG_OUTPUT_DETAIL 0

namespace atlas {
namespace util {

bool approx_eq( const double& v1, const double& v2, const double& t ) {
    return eckit::types::is_approximately_equal( v1, v2, t );
}

bool approx_eq( const PointXYZ& v1, const PointXYZ& v2, const double& t ) {
    return approx_eq( v1[0], v2[0], t ) && approx_eq( v1[1], v2[1], t ) && approx_eq( v1[2], v2[2], t );
}

//------------------------------------------------------------------------------------------------------

//ConvexSphericalPolygon::ConvexSphericalPolygon( const PartitionPolygon& partition_polygon ) :
//    PolygonCoordinates( partition_polygon.xy(), false ) {}

/*
ConvexSphericalPolygon::ConvexSphericalPolygon() : 
	ConvexSphericalPolygon( std::vector<PointXYZ>{PointXYZ(0,0,0),PointXYZ(0,0,0),PointXYZ(0,0,0)} ) {
	valid_ = false;
}

// TODO: earth radius set to 1 !!
ConvexSphericalPolygon::ConvexSphericalPolygon( const std::vector<PointXYZ>& points ) : PolygonCoordinates( points ){
	sph_coords_.clear();
	sph_coords_.resize( points.size() );
	for ( size_t i = 0; i < points.size(); ++i ) {
		sph_coords_[i] = points[i];
		eckit::geometry::Sphere::convertCartesianToSpherical( 1., sph_coords_[i], coordinates_[i] );
	}
	valid_ = true; // assume all are convex
#ifndef NDEBUG
	validate();
#endif
}
*/

// TODO: earth radius set to 1 !!
ConvexSphericalPolygon::ConvexSphericalPolygon( const std::vector<PointLonLat>& points ) :
    size_(points.size()) {
	sph_coords_.clear();
	sph_coords_.resize( points.size() );
	for ( size_t i = 0; i < points.size(); ++i ) {
		eckit::geometry::Sphere::convertSphericalToCartesian( 1., points[i], sph_coords_[i] );
	}
	valid_ = true; // assume all are convex
#ifndef NDEBUG
	validate();
#endif
}

bool ConvexSphericalPolygon::validate() {
	valid_ = true;
	for( int i = 0; i < size()-1; i++ ) {
		int ni = ( i != size()-1 ? i+1 : 0 );
		int nni = ( ni != size()-1 ? ni+1 : 0 );
		if ( ! leftOf( sph_coords_[nni], sph_coords_[i], sph_coords_[ni] ) ) {
			valid_ = false;
			return false;
		}
	}
	return true;
}

bool ConvexSphericalPolygon::equals( const ConvexSphericalPolygon& plg, const double prec ) const {
	const int sz = size();
	if ( sz != plg.size() ) {
		return false;
	}
	int i = 0;
	for( ; i < sz; i++ ) {
        if ( approx_eq( plg.sph_coords_[i], this->sph_coords_[i], prec ) ) {
			break;
		}
	}
	for( int j = 0; j < sz; j++ ) {
		int idx = (i+j)%sz;
        if ( ! approx_eq( plg.sph_coords_[idx], this->sph_coords_[idx], prec ) ) {
			return false;
		}
	}
	return true;
}

// note: two diameterly opposite points forming segment are not allowed
// return tangential angle between [pl,p] and [p,pr]
inline double ConvexSphericalPolygon::angle( const PointXYZ& pl, const PointXYZ& p, const PointXYZ& pr ) const {
	const PointXYZ& plp = PointXYZ( PointXYZ::cross( pl, p ) );
	const PointXYZ& ppr = PointXYZ( PointXYZ::cross( p, pr ) );
	double s1p = PointXYZ::dot( plp, ppr ) / ( PointXYZ::norm( plp ) * PointXYZ::norm( ppr ) );
	s1p = std::acos( ( s1p < 0. ? -1 : 1 ) * std::min( 1., std::abs(s1p) ) );
	return M_PI - s1p;
}


// note: unit sphere!
// I. Todhunter (1886), Paragr. 99
double ConvexSphericalPolygon::area() const {
	const int sz = size()-1;
	double a = M_PI*(2-sz);
	for ( int i = 0; i < sz; i++ ) {
		int im1 = (i != 0) ? i-1 : sz-1;
		int ip1 = (i != sz-1) ? i+1 : 0;
		a += angle( sph_coords_[im1], sph_coords_[i], sph_coords_[ip1] );
	}
	ATLAS_ASSERT( a > -eps_ );
	return a;
}

// return 0:P_right_of_[p1,p2], -1:overlap_of_[P,p1]_and_[P,p2], 1:P_left_of_[p1,p2]
inline int ConvexSphericalPolygon::leftOf( const PointXYZ& P, const PointXYZ& p1, const PointXYZ& p2 ) const {
	const PointXYZ& cp = PointXYZ( PointXYZ::cross( p1, p2 ) );
	ATLAS_ASSERT( PointXYZ::norm(cp) > eps_ );
	const double dp = PointXYZ::dot( cp, P );
	return ( dp > eps_ ? 1 : ( dp < -eps_ ? 0 : -1 ) );
}

// return 0:outside, -1:on_edge, 1:strictly_inside
int ConvexSphericalPolygon::contains( const PointXYZ& P ) const {
    ATLAS_ASSERT( sph_coords_.size() >= 2 );
	const size_t ncoord = sph_coords_.size()-1;
	for ( size_t i = 0; i < ncoord; ++i ) {
		const PointXYZ sp1 = sph_coords_[i];
		const PointXYZ sp2 = sph_coords_[i!=ncoord-1 ? i+1 : 0];
		int dp = leftOf( P, sp1, sp2 );
		if ( dp != 1 ) {
			return dp;
		}
	}
	return 1;
}

// note: two diameterly opposite points forming segment are not allowed
// note: [s1,s2] is always the smaller part of THE great circle through s1 and s2.
bool ConvexSphericalPolygon::onSegment( const PointXYZ& P, const PointXYZ& s1, const PointXYZ& s2 ) const {
	ATLAS_ASSERT( s1 != PointXYZ::mul(s2,-1.) ); // should be done in "validate" because s1,s2,p1,p2 are polydon vertices
	ATLAS_ASSERT( std::abs(PointXYZ::dot(P,P) - 1) < eps_ );
	ATLAS_ASSERT( std::abs(PointXYZ::dot(s1,s1) - 1) < eps_ ); // should be done in "validate" because s1,s2,p1,p2 are polydon vertices
	ATLAS_ASSERT( std::abs(PointXYZ::dot(s2,s2) - 1) < eps_ ); // should be done in "validate" because s1,s2,p1,p2 are polydon vertices
	double s1p = PointXYZ::dot(s1,P);
	double s1s2 = PointXYZ::dot(s1,s2);
	double ps2 = PointXYZ::dot(P,s2);
	s1p = ( s1p < 0. ? -1 : 1 ) * std::min( 1., std::abs(s1p) );
	s1s2 = ( s1s2 < 0. ? -1 : 1 ) * std::min( 1., std::abs(s1s2) );
	ps2 = ( ps2 < 0. ? -1 : 1 ) * std::min( 1., std::abs(ps2) );
	double angl12 = acos( s1p );
	double angl13 = acos( s1s2 );
	double angl23 = acos( ps2 );
	return std::abs(angl12 + angl23 - angl13) < eps_;
}

// intersect segment [s1,s2] with this polygon
// @param[in] s1, s2 segment end points
// @param[out] ip intersection point if any
// @return -1:overlap, 0:no_intersect, 1+(id of this polygon's segment intersecting the given segment)
int ConvexSphericalPolygon::intersect( const PointXYZ& s1, const PointXYZ& s2, PointXYZ& ip, int start ) const {
	const PointXYZ& cp1 = static_cast< PointXYZ >( PointXYZ::cross( s1, s2 ) );
	ATLAS_ASSERT( PointXYZ::norm(cp1) > eps_ );
	int ncoord = sph_coords_.size() - 1;

#if DEBUG_OUTPUT_DETAIL
	PointLonLat s1ll, s2ll;
	eckit::geometry::Sphere::convertCartesianToSpherical( 1., s1, s1ll );
	eckit::geometry::Sphere::convertCartesianToSpherical( 1., s2, s2ll );
	std::cout <<"   doing intersection with [" <<s1ll <<", " <<s2ll <<"]\n";
	std::cout <<"   	ncoord: " <<ncoord <<"\n";
	std::cout.flush();
#endif
	for( int i = start; i < start+ncoord; i++ ) {
		const PointXYZ& sp1 = sph_coords_[ i%ncoord ];
		const PointXYZ& sp2 = sph_coords_[ (i+1)%ncoord ];
#if DEBUG_OUTPUT_DETAIL
		std::cout <<"     check edge " <<i%ncoord <<" :" <<coordinates_[i%ncoord] <<", " <<coordinates_[(i+1)%ncoord] <<"]\n";
		std::cout.flush();
#endif
		const PointXYZ& cp2 = static_cast< PointXYZ >( PointXYZ::cross( sp1, sp2 ) );
		ATLAS_ASSERT( PointXYZ::norm(cp2) > eps_ );
		ip = static_cast< PointXYZ >( PointXYZ::cross( cp1, cp2 ) );

		if ( PointXYZ::norm( ip ) > eps_ ) {
#if DEBUG_OUTPUT_DETAIL
			PointLonLat ipp;
			eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ipp );
			std::cout <<"       try intersect " <<ipp <<"\n";
			std::cout.flush();
#endif
			ip = PointXYZ::div( ip, PointXYZ::norm( ip ) );
			if ( onSegment( ip, s1, s2 ) && onSegment( ip, sp1, sp2 ) ) {
				return 1+i%ncoord;
			}
			ip = PointXYZ::mul( ip, -1. );
#if DEBUG_OUTPUT_DETAIL
			eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ipp );
			std::cout <<"       try intersect " <<ipp <<"\n";
			std::cout.flush();
#endif
			if ( onSegment( ip, s1, s2 ) && onSegment( ip, sp1, sp2 ) ) {
				return 1+i%ncoord;
			}
		}
		else {
			//overlap
			if ( onSegment( s1, sp1, sp2 ) && (onSegment( s2, s1, sp2 ) || onSegment( sp2, s1, s2 )) ) {
				ip = *(new PointXYZ( s1 ));
#if DEBUG_OUTPUT_DETAIL
				PointLonLat ipp;
				eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ipp );
				std::cout <<"       		got first point " <<ipp <<"\n";
				std::cout.flush();
#endif
			}
			else if ( onSegment( sp1, s1, s2 ) && (onSegment( sp2, sp1, s2 ) || onSegment( s2, sp1, sp2 )) ) {
				ip = *(new PointXYZ( sp1 ));
#if DEBUG_OUTPUT_DETAIL
				PointLonLat ipp;
				eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ipp );
				std::cout <<"       		got first point " <<ipp <<"\n";
				std::cout.flush();
#endif
			}
			else {
				ip = *(new PointXYZ( {0.,0.,0.} ));
#if DEBUG_OUTPUT_DETAIL
				std::cout <<"       		no intersection in overlap\n";
				std::cout.flush();
#endif
				return 0;
			}
			return 1+i%ncoord;
		}
	}
	return 0;
}

// intersect a polygon with this polygon
// @param[in] pol clipping polygon
// @param[out] intersecting polygon
ConvexSphericalPolygon* ConvexSphericalPolygon::intersect( const ConvexSphericalPolygon& plg ) const {
	std::vector< PointXYZ > iplg_p;
	int ii = 0; // "this" vertex counter
	int jj = 0; // "plg" vertex counter
	PointXYZ ip;
	const int n_plg = this->size() - 1;
	for( ; ii < n_plg; ii++ ) {
		jj = -1 + plg.intersect( sph_coords_[ii], sph_coords_[(ii+1)%n_plg], ip );
		if ( jj != -1 ) {
#if DEBUG_OUTPUT_DETAIL
			PointLonLat ipp;
			eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ipp );
			std::cout <<"  " <<jj <<"th edge intersects with [" <<coordinates_[ii] <<" "
<<coordinates_[(ii+1)%n_plg] <<"] at " <<ipp <<"\n";
			std::cout.flush();
#endif
			break;
		}
#if DEBUG_OUTPUT_DETAIL
		std::cout <<"  polygon does not intersects with [" <<coordinates_[ii] <<" " <<coordinates_[(ii+1)%n_plg] <<"]\n";
		std::cout.flush();
#endif
	}
	if ( jj != -1 ) {
		iplg_p.emplace_back( ip );
		int intersect = nextIntersect( iplg_p, plg, ii, jj, 0 );
		std::vector< PointLonLat > iplg_p_ll;
		iplg_p_ll.resize( iplg_p.size() );
		for( int i = 0; i < iplg_p.size(); i++ ) {
            eckit::geometry::Sphere::convertCartesianToSpherical( 1., iplg_p[i], iplg_p_ll[i] );
		}
		return ( intersect ? new ConvexSphericalPolygon( iplg_p_ll ) : nullptr );
	}
	else {
#if DEBUG_OUTPUT_DETAIL
		std::cout <<" polygons edges do not intersect with edges of the other polygon.\n";
	    std::cout.flush();
#endif
		if ( this->contains( plg.sph_coords_[0] ) == 1 ) {
#if DEBUG_OUTPUT_DETAIL
			std::cout <<" this contains " <<plg.coordinates_[0] <<" -> plg inside this.\n";
			std::cout.flush();
#endif
			return new ConvexSphericalPolygon( plg );
		}
		else if ( plg.contains( sph_coords_[0] ) == 1 ) {
#if DEBUG_OUTPUT_DETAIL
			std::cout <<" plg contains " <<coordinates_[0] <<" -> this inside plg.\n";
			std::cout.flush();
#endif
			return new ConvexSphericalPolygon( *this );
		}
		else {
#if DEBUG_OUTPUT_DETAIL
			std::cout <<" this contains " <<plg.coordinates_[0] <<" -> plg inside this.\n";
			std::cout <<" plg NOT inside this && this NOT inside plg\n";
			std::cout.flush();
#endif
			return nullptr;
		}
	}
}

int ConvexSphericalPolygon::nextIntersect( std::vector< PointXYZ >& iplg_p,
							 const ConvexSphericalPolygon& plg,
							 const int ii0,
          					 const int jj0,
							 const int inside
						   ) const {
	const int n_iplg = iplg_p.size() - 1;
    const int n_plg1 = size() - 1;
	const int n_plg2 = plg.size() - 1;
	const PointXYZ P = iplg_p[ n_iplg ];
	int ii = ii0;
	int jj = jj0;
	int nii = (ii0 != n_plg1-1) ? ii0+1 : 0;
	int njj = (jj0 != n_plg2-1) ? jj0+1 : 0;

	const PointXYZ* pnp1 = &sph_coords_[nii];
	const PointXYZ* pnp2 = &plg.sph_coords_[njj];
	if ( !inside && PointXYZ::norm(P - *pnp1) < eps_ ) {
		int nnii = (nii != n_plg1-1) ? nii+1 : 0;
		ii = nii;
		nii = nnii;
		pnp1 = &sph_coords_[nii];
	}
	if ( !inside && PointXYZ::norm(P - *pnp2) < eps_ ) {
		int nnjj = (njj != n_plg2-1) ? njj+1 : 0;
		jj = njj;
		njj = nnjj;
		pnp2 = &(plg.sph_coords_[njj]);
	}
	const PointXYZ& np1 = *pnp1;
	const PointXYZ& np2 = *pnp2;

#if DEBUG_OUTPUT_DETAIL
	std::cout <<"\n == doing nextIntersect(i0,ii,nii,j0,jj,njj,inside): "<<ii0<<", "<<ii<<", "<<nii<<", "<<jj0<<", "<<jj<<", "<<njj<<", " <<inside<<"\n";
	std::cout <<" iplg_p: ";
	for( int i = 0; i <= n_iplg; i++ ) {
		PointLonLat ip_ll;
		eckit::geometry::Sphere::convertCartesianToSpherical( 1., iplg_p[i], ip_ll );
		std::cout <<", " <<ip_ll;
	}
	std::cout <<"\n";
	PointLonLat ip_ll;
	eckit::geometry::Sphere::convertCartesianToSpherical( 1., P, ip_ll );
	std::cout <<"P = " <<ip_ll;
	std::cout.flush();
#endif
	if ( n_iplg > 1 && PointXYZ::norm( iplg_p[0] - iplg_p[n_iplg] ) < eps_ ) {
		return 1;
	}

	PointXYZ ip; //next vertex of the iplg
	int new_ii = ii;
	int new_jj = jj;
	const bool P_eq_plgi = (PointXYZ::norm(P - sph_coords_[ii]) < eps_);
	const bool P_eq_plgj = (PointXYZ::norm(P - plg.sph_coords_[jj]) < eps_);
#if DEBUG_OUTPUT_DETAIL
	std::cout <<"   P_eq_plgi, P_eq_plgj: "<<P_eq_plgi<<", "<<P_eq_plgj <<"\n";
	std::cout.flush();
#endif
	
	if ( inside ) {
		if ( P_eq_plgi ) {
			new_jj = -1 + plg.intersect( P, np1, ip, njj );
			bool ip_eq_P = ( PointXYZ::norm( ip - iplg_p[n_iplg] ) < eps_ );
			if ( new_jj == -1 || ip_eq_P ) { // no intersection of plg with [P,np1] means np1 is inside plg
#if DEBUG_OUTPUT_DETAIL
				std::cout <<"   inside of plg: no intersection with polygon.\n";
				std::cout.flush();
#endif
				if ( plg.contains( np1 ) == 1 ) {
					iplg_p.emplace_back( np1 );
					return nextIntersect( iplg_p, plg, nii, jj, 1 );
				}
				else {
#if DEBUG_OUTPUT_DETAIL
					std::cout <<"   inside of plg: np1 is outside no intersection with polygon.\n";
					std::cout.flush();
#endif
					return 0;
				}
			} else {
#if DEBUG_OUTPUT_DETAIL
				PointLonLat ip_ll;
				eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ip_ll );
				PointLonLat ip_ll2;
				eckit::geometry::Sphere::convertCartesianToSpherical( 1., iplg_p[n_iplg], ip_ll2 );
				std::cout <<" inside of plg: intersection at " <<ip_ll <<", last point: " <<ip_ll2 <<"\n";
				std::cout.flush();
#endif
				iplg_p.emplace_back( ip );
				return nextIntersect( iplg_p, plg, ii, new_jj, 0 );
			}
		}
		else {
			new_ii = -1 + this->intersect( P, np2, ip, nii );
			bool ip_eq_P = ( PointXYZ::norm( ip - iplg_p[n_iplg] ) < eps_ );
			if ( new_ii == -1 || ip_eq_P ) { // no intersection of this with [P,np2]
#if DEBUG_OUTPUT_DETAIL
				std::cout <<"   inside of this: no intersection with polygon.\n";
				std::cout.flush();
#endif
				if ( this->contains( np2 ) == 1 ) {
					iplg_p.emplace_back( np2 );
					return nextIntersect( iplg_p, plg, ii, njj, 1 );
				}
				else {
#if DEBUG_OUTPUT_DETAIL
					std::cout <<"   inside of this: np2 is outside no intersection with polygon.\n";
					std::cout.flush();
#endif
					return 0;
				}
			} else {
#if DEBUG_OUTPUT_DETAIL
				PointLonLat ip_ll;
				eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ip_ll );
				PointLonLat ip_ll2;
				eckit::geometry::Sphere::convertCartesianToSpherical( 1., iplg_p[n_iplg], ip_ll2 );
				std::cout <<" inside of this: intersection at " <<ip_ll <<", last point: " <<ip_ll2 <<"\n";
				std::cout.flush();
#endif
				iplg_p.emplace_back( ip );
				return nextIntersect( iplg_p, plg, new_ii, jj, 0 );
			}
		}
		ATLAS_ASSERT( false );
	}

	if ( P_eq_plgi || P_eq_plgj ) {
#if DEBUG_OUTPUT_DETAIL
		std::cout <<" edges could overlap\n";
		std::cout <<"  " <<onSegment( np1, P, np2 ) <<", " <<onSegment( np2, P, np1 ) <<", " <<onSegment( P, np1, np2 ) <<"\n";
		std::cout.flush();
#endif
		// polygon edges overlap
		if ( onSegment( np1, P, np2 ) ) {
			iplg_p.emplace_back( np1 );
			return nextIntersect( iplg_p, plg, ii, jj, 0 );
		}
		else if ( onSegment( np2, P, np1 ) ) {
			iplg_p.emplace_back( np2 );
			return nextIntersect( iplg_p, plg, ii, jj, 0 );
		}
		else if ( onSegment( P, np1, np2 ) ) {
			// polygons intersect only on segment --> no polygon intersection
			return 0;
		}
	}

	// not an inside point & no edges-overlap
#if DEBUG_OUTPUT_DETAIL
	PointLonLat Pll, np1ll, np2ll;
	eckit::geometry::Sphere::convertCartesianToSpherical( 1., P, Pll );
	eckit::geometry::Sphere::convertCartesianToSpherical( 1., np1, np1ll );
	eckit::geometry::Sphere::convertCartesianToSpherical( 1., np2, np2ll );
	std::cout <<"   "<<np2ll <<" leftOf ["<<Pll <<"," <<np1ll <<" -> " <<leftOf( np2, P, np1 ) <<"\n";
	std::cout <<"   "<<np1ll <<" leftOf ["<<Pll <<"," <<np2ll <<" -> " <<leftOf( np1, P, np2 ) <<"\n";
	std::cout.flush();
#endif
	if ( leftOf( np2, P, np1 ) ) { // 2a) 3b)
		new_ii = -1 + this->intersect( P, np2, ip, nii );
		bool ip_eq_P = ( PointXYZ::norm( ip - iplg_p[n_iplg] ) < eps_ );
		if ( new_ii == -1 || ip_eq_P ) { // no intersection of plg with [P,np2] means np2 is inside plg
#if DEBUG_OUTPUT_DETAIL
			std::cout <<"   no intersection with polygon.\n";
			std::cout.flush();
#endif
			iplg_p.emplace_back( np2 );
			return nextIntersect( iplg_p, plg, ii, njj, 1 );
		} else {
#if DEBUG_OUTPUT_DETAIL
			PointLonLat ip_ll;
			eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ip_ll );
			PointLonLat ip_ll2;
	        eckit::geometry::Sphere::convertCartesianToSpherical( 1., iplg_p[n_iplg], ip_ll2 );
			std::cout <<" leftOf( np2, P, np1 ): intersection at " <<ip_ll <<", last point: " <<ip_ll2 <<"\n";
			std::cout.flush();
#endif
			iplg_p.emplace_back( ip );
			return nextIntersect( iplg_p, plg, new_ii, new_jj, 0 );
		}
	}
	else { //if ( leftOf( np1, P, np2 ) ) { // 3a) 2b)
		new_jj = -1 + plg.intersect( P, np1, ip, njj ); 
		bool ip_eq_P = ( PointXYZ::norm( ip - iplg_p[n_iplg] ) < eps_ );
		if ( new_jj == -1 || ip_eq_P) {
#if DEBUG_OUTPUT_DETAIL
			std::cout <<"   no intersection with polygon.\n";
			std::cout.flush();
#endif
			// no intersection of plg2 with [P,np1] means np1 is inside plg
			iplg_p.emplace_back( np1 );
			return nextIntersect( iplg_p, plg, nii, jj, 1 );
		} else {
#if DEBUG_OUTPUT_DETAIL
			PointLonLat ip_ll;
			eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ip_ll );
			PointLonLat ip_ll2;
	        eckit::geometry::Sphere::convertCartesianToSpherical( 1., iplg_p[n_iplg], ip_ll2 );
			std::cout <<" leftOf( np1, P, np2 ): intersection at " <<ip_ll <<", last point: " <<ip_ll2 <<"\n";
			std::cout.flush();
#endif
			iplg_p.emplace_back( ip );
			return nextIntersect( iplg_p, plg, new_ii, new_jj, 0 );
		}
	}
}

void ConvexSphericalPolygon::print( std::ostream& out ) const {
    out << "[";
    for ( size_t i = 0; i < size(); ++i ) {
        if ( i > 0 ) {
            out << " ";
        }
        out << sph_coords_[i];
    }
    out << "]";
}


//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
