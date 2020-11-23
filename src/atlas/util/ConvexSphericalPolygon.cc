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

#include "eckit/types/FloatCompare.h"

#include "atlas/runtime/Exception.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/NormaliseLongitude.h"

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

ConvexSphericalPolygon::ConvexSphericalPolygon( const PartitionPolygon& partition_polygon ) :
    PolygonCoordinates( partition_polygon.xy(), false ) {}

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
int ConvexSphericalPolygon::contains( const PointXYZ& P ) const {
    ATLAS_ASSERT( coordinates_.size() >= 2 );
	const size_t ncoord = sph_coords_.size();
	for ( size_t i = 0; i < ncoord; ++i ) {
		const PointXYZ cp = PointXYZ( PointXYZ::cross( sph_coords_[i], sph_coords_[ i!=ncoord-1 ? i+1 : 0 ] ) );
		const double dp = PointXYZ::dot( cp, P );
		if ( dp < eps_ ) {
			return 0;
		}
		else if ( abs(dp) < eps_ ) {
			return -1;
		}
	}
	return 1;
}

bool ConvexSphericalPolygon::onSegment( const PointXYZ& P, const PointXYZ& s1, const PointXYZ& s2 ) const {
	ASSERT( abs(PointXYZ::dot(P,P) - 1) < eps_ );
	ASSERT( abs(PointXYZ::dot(s1,s1) - 1) < eps_ );
	ASSERT( abs(PointXYZ::dot(s2,s2) - 1) < eps_ );
	double angl12 = acos( std::min( 1., PointXYZ::dot(s1,P) ) );
	double angl13 = acos( std::min( 1., PointXYZ::dot(s1,s2) ) );
	double angl23 = acos( std::min( 1., PointXYZ::dot(P,s2) ) );
	return angl12 + angl23 - angl13 < eps_;
}

// intersect segment [s1,s2] with this polygon
// @param[in] s1, s2 segment end points
// @param[out] ip intersection point if any
bool ConvexSphericalPolygon::intersect( const PointXYZ& s1, const PointXYZ& s2, PointXYZ& ip, int start ) const {
	const PointXYZ& cp1 = static_cast< PointXYZ >( PointXYZ::cross( s1, s2 ) );
	int ncoord = sph_coords_.size() - 1;
	for( int i = start; i < start+ncoord; i++ ) {
		std::cout <<"Check with edge " <<i <<"\n";
		const PointXYZ cp2 = static_cast< PointXYZ >( PointXYZ::cross( sph_coords_[ i%ncoord ], sph_coords_[ (i+1)%ncoord ] ) );
		ip = static_cast< PointXYZ >( PointXYZ::cross( cp1, cp2 ) );
		ip = PointXYZ::div( ip, PointXYZ::norm( ip ) );
		if ( onSegment( ip, s1, s2 ) ) {
			return true;
		}
		ip = PointXYZ::mul( ip, -1. );
		if ( onSegment( ip, s1, s2 ) ) {
			return true;
		}
	}
	return false;
}

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
