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

#define DEBUG_OUTPUT

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

//ConvexSphericalPolygon::ConvexSphericalPolygon( const PartitionPolygon& partition_polygon ) :
//    PolygonCoordinates( partition_polygon.xy(), false ) {}

// TODO: earth radius set to 1 !!
ConvexSphericalPolygon::ConvexSphericalPolygon( const std::vector<PointLonLat>& points ) : PolygonCoordinates( points ) {
	eckit::geometry::Sphere sphere;
	sph_coords_.clear();
	sph_coords_.resize( points.size() );
	for ( size_t i = 0; i < points.size(); ++i ) {
		sphere.convertSphericalToCartesian( 1., points[i], sph_coords_[i] );
	}	
}

// return 0:outside, -1:on_edge, 1:strictly_inside
inline int ConvexSphericalPolygon::leftOf( const PointXYZ& P, const PointXYZ& p1, const PointXYZ& p2 ) const {
	const PointXYZ cp = PointXYZ( PointXYZ::cross( p1, p2 ) );
	ATLAS_ASSERT( PointXYZ::norm(cp) > eps_ );
	const double dp = PointXYZ::dot( cp, P );
	if ( dp < eps_ ) {
		return 0;
	}
	else if ( fabs(dp) < eps_ ) {
		return -1;
	}
	return 1;
}

// return 0:outside, -1:on_edge, 1:strictly_inside
int ConvexSphericalPolygon::contains( const PointXYZ& P ) const {
    ATLAS_ASSERT( coordinates_.size() >= 2 );
	const size_t ncoord = sph_coords_.size();
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
	ATLAS_ASSERT( s1 != PointXYZ::mul(s2,-1.) );
	ATLAS_ASSERT( fabs(PointXYZ::dot(P,P) - 1) < eps_ );
	ATLAS_ASSERT( fabs(PointXYZ::dot(s1,s1) - 1) < eps_ );
	ATLAS_ASSERT( fabs(PointXYZ::dot(s2,s2) - 1) < eps_ );
	double angl12 = acos( std::min( 1., PointXYZ::dot(s1,P) ) );
	double angl13 = acos( std::min( 1., PointXYZ::dot(s1,s2) ) );
	double angl23 = acos( std::min( 1., PointXYZ::dot(P,s2) ) );
	return fabs(angl12 + angl23 - angl13) < eps_;
}

// intersect segment [s1,s2] with this polygon
// @param[in] s1, s2 segment end points
// @param[out] ip intersection point if any
// @return -1:overlap, 0:no_intersect, 1+(id of this polygon's segment intersecting the given segment)
int ConvexSphericalPolygon::intersect( const PointXYZ& s1, const PointXYZ& s2, PointXYZ& ip, int start ) const {
	const PointXYZ& cp1 = static_cast< PointXYZ >( PointXYZ::cross( s1, s2 ) );
	ATLAS_ASSERT( PointXYZ::norm(cp1) > eps_ );
	int ncoord = sph_coords_.size() - 1;

	for( int i = start; i < start+ncoord; i++ ) {
		std::cout <<"Check with edge " <<i%ncoord <<"\n";
		const PointXYZ& sp1 = sph_coords_[ i%ncoord ];
		const PointXYZ& sp2 = sph_coords_[ (i+1)%ncoord ];
#ifdef DEBUG_OUTPUT
		PointLonLat pp1, pp2, ipp;
		eckit::geometry::Sphere::convertCartesianToSpherical( 1., sp1, pp1 );
		eckit::geometry::Sphere::convertCartesianToSpherical( 1., sp2, pp2 );
		std::cout <<"  with seg. [" <<pp1 <<", " <<pp2 <<"]\n";
		std::cout.flush();
#endif
		const PointXYZ& cp2 = static_cast< PointXYZ >( PointXYZ::cross( sp1, sp2 ) );
		ATLAS_ASSERT( PointXYZ::norm(cp2) > eps_ );
		ip = static_cast< PointXYZ >( PointXYZ::cross( cp1, cp2 ) );
#ifdef DEBUG_OUTPUT
		eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ipp );
		std::cout <<"  try intersect " <<ipp <<"\n";
		std::cout.flush();
#endif

		if ( PointXYZ::norm( ip ) > eps_ ) {
			ip = PointXYZ::div( ip, PointXYZ::norm( ip ) );
			if ( onSegment( ip, s1, s2 ) && onSegment( ip, sp1, sp2 ) ) {
				return 1+i;
			}
			ip = PointXYZ::mul( ip, -1. );
#ifdef DEBUG_OUTPUT
			eckit::geometry::Sphere::convertCartesianToSpherical( 1., ip, ipp );
			std::cout <<"  try intersect " <<ipp <<"\n";
			std::cout.flush();
#endif
			if ( onSegment( ip, s1, s2 ) && onSegment( ip, sp1, sp2 ) ) {
				return 1+i;
			}
		}
		else {
			//overlap
			return -1;
		}
	}
	return 0;
}

// intersect a polygon with this polygon
// @param[in] pol clipping polygon
// @param[out] intersecting polygon
ConvexSphericalPolygon* ConvexSphericalPolygon::intersect( const ConvexSphericalPolygon& plg ) const {
		std::cout <<" Entered ConvexSphericalPolygon::intersect\n";
	    std::cout.flush();
	std::vector< PointXYZ > plg_points;
	int ii = 0;
	int jj = 0;
	PointXYZ ip;
	const int plgsize = plg.size() - 1;
	for( ; ii < plgsize; ii++ ) {
		jj = this->intersect( plg.sph_coords_[ii], plg.sph_coords_[(ii+1)%(plgsize-1)], ip );
		if ( jj ) {
			break;
		}
	}
	if ( jj ) {
		plg_points.emplace_back( ip );
		nextIntersect( plg_points, plg, 0, *this, 0 );
	}
	else {
		std::cout <<" plg edges not intersecting: " <<"\n";
	    std::cout.flush();
		ATLAS_ASSERT( false );
	}
	return nullptr;
	//return new ConvexSphericalPolygon( plg_points );
}

int ConvexSphericalPolygon::nextIntersect( std::vector< PointXYZ >& plg_points,
							 const ConvexSphericalPolygon& plg1,
							 int ii,
					         const ConvexSphericalPolygon& plg2, 
          					 int jj
						   ) const {
	const int nplgi = plg_points.size();
	if ( nplgi > 2 && plg_points[0] == plg_points[nplgi-1] ) {
		return 1;
	}
	const int nplg1 = plg1.size() - 1;
    const int nplg2 = plg2.size() - 1;
	const PointXYZ P = plg_points[nplgi-1];
	int nii = (ii != nplg1-1) ? ii+1 : 0;
	int njj = (jj != nplg2-1) ? jj+1 : 0;
	const PointXYZ& np1 = plg1.sph_coords_[nii];
	const PointXYZ& np2 = plg2.sph_coords_[njj];

#ifdef DEBUG_OUTPUT
	std::cout <<" == doing nextIntersect( "<<ii<<", "<<jj<<" )\n";
	std::cout <<" plg_points: ";
	for( int i = 0; i < nplgi; i++ ) {
		PointLonLat ip_ll;
		eckit::geometry::Sphere::convertCartesianToSpherical( 1., plg_points[i], ip_ll );
		std::cout <<", " <<ip_ll <<"\n";
	}
	std::cout.flush();
#endif


	if ( P == plg2.sph_coords_[jj] || P == plg1.sph_coords_[ii] ) {
#ifdef DEBUG_OUTPUT
		std::cout <<" edges overlap\n";
		std::cout.flush();
#endif
		// polygon edges overlap
		if ( onSegment( np1, P, np2 ) ) {
			plg_points.emplace_back( np1 );
			return nextIntersect( plg_points, plg1, nii, plg2, njj );
		}
		else if ( onSegment( np2, P, np1 ) ) {
			plg_points.emplace_back( np2 );
			return nextIntersect( plg_points, plg1, nii, plg2, njj );
		}
		else if ( onSegment( P, np1, np2 ) ) {
			// polygons intersect only on segment --> no polygon intersection
			return 0;
		}
	}

	return 0;

	// no edges-overlap
	PointXYZ ip; //next intersection-polygon vertex
	int new_ii;
	int new_jj;
	if ( leftOf( np2, P, np1 ) ) { // 2a) 3b)
		new_jj = jj + 1;
		new_ii = -1 + plg1.intersect( P, np2, ip, ii );
		if ( ! new_ii ) { // no intersection of plg1 with [P,np2] means np2 is inside plg1
			new_ii = ii + 1;
			plg_points.emplace_back( np2 );
			return nextIntersect( plg_points, plg1, new_ii, plg2, new_jj );
		} else {
			plg_points.emplace_back( ip );
			return nextIntersect( plg_points, plg1, new_ii, plg2, new_jj );
		}
	}
	else { //if ( leftOf( np1, P, np2 ) ) { // 3a) 2b)
		new_jj = jj + 1;
		new_jj = -1 + plg2.intersect( P, np1, ip, jj );
		if ( ! new_jj ) {
			new_jj = jj + 1;
			// no intersection of plg2 with [P,np1] means np1 is inside plg1
			plg_points.emplace_back( np1 );
			return nextIntersect( plg_points, plg1, new_ii, plg2, new_jj );
		} else {
			plg_points.emplace_back( ip );
			return nextIntersect( plg_points, plg1, new_ii, plg2, new_jj );
		}
	}
}

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
