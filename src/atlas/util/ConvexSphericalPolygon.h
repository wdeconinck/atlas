/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <vector>

#include "atlas/util/detail/Debug.h"
#include "atlas/util/Point.h"
#include "atlas/util/Polygon.h"

namespace atlas {
namespace util {

class PartitionPolygon;
//------------------------------------------------------------------------------------------------------

class ConvexSphericalPolygon : public PolygonCoordinates {
public:
    ConvexSphericalPolygon( const PartitionPolygon& );

    ConvexSphericalPolygon( const std::vector<PointLonLat>& points );
    //ConvexSphericalPolygon( const std::vector<PointXYZ>& points );

    /*
   * Point-in-polygon test on sphere with spherical polygons
   * @param[in] P given point in (x,y,z) coordinates
   * @return 0:outside, -1:on_edge, 1:strictly_inside
   */
    int contains( const PointXYZ& P ) const;

    bool contains( const Point2& P ) const override {
		ATLAS_ASSERT( false );
	}

	static constexpr double eps_ = 1e-16;

//protected:
    /*
   * Point-on-segment test on great circle segments
   * @param[in] P given point in (x,y,z) coordinates
   * @return 
   */
    bool onSegment( const PointXYZ& P, const PointXYZ& s1, const PointXYZ& p2 ) const;

    /*
   * Segment-sph_polygon intersection
   * @param[in] s1, s2 segment endpoints in (x,y,z) coordinates
   * @param[in] start start with polygon segments [pol[start],pol[start+1]],...
   * @param[out] ip intersection point or nullptr
   * @return 0:no_intersection, 1:
   */
    bool intersect( const PointXYZ& s1, const PointXYZ& s2, PointXYZ& ip, int start = 0 ) const;


private:
	std::vector<PointXYZ> sph_coords_;
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
