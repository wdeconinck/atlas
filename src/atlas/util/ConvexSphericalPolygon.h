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

#include "atlas/util/Point.h"
#include "atlas/util/Polygon.h"
#include "atlas/util/detail/Debug.h"

namespace atlas {
namespace util {

class PartitionPolygon;
//------------------------------------------------------------------------------------------------------

class ConvexSphericalPolygon {
public:
    static constexpr int MAX_SIZE = 8;

    ConvexSphericalPolygon();
    ConvexSphericalPolygon( const std::vector<PointXYZ>& points );
    ConvexSphericalPolygon( const std::vector<PointLonLat>& points );
    //ConvexSphericalPolygon( const PartitionPolygon& );

    /*
   * Point-in-polygon test on sphere with spherical polygons
   * @param[in] P given point in (x,y,z) coordinates
   * @return 0:outside, -1:on_edge, 1:strictly_inside
   */
    int contains( const PointXYZ& P ) const;

    operator bool() const { return valid_; }


    /*
   * Point left of [p1,p2]
   * @param[in] P, p1, p2 given point in xyz-coordinates
   * @return 0:P_right_of_[p1,p2], -1:overlap_of_[P,p1]_and_[P,p2], 1:P_left_of_[p1,p2]
   */
    inline int leftOf( const PointXYZ& P, const PointXYZ& p1, const PointXYZ& p2 ) const;

    double area() const;

	const PointXYZ& centroid() const;

	// return tangential angle between [pl,p] and [p,pr]
	inline double angle( const PointXYZ& pl, const PointXYZ& p, const PointXYZ& pr ) const;

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
    int intersect( const PointXYZ& s1, const PointXYZ& s2, PointXYZ& ip, int start = 0 ) const;

    /*
   * intersect a polygon with this polygon
   * @param[in] pol clipping polygon
   * @param[out] intersecting polygon
   */
    ConvexSphericalPolygon intersect( const ConvexSphericalPolygon& pol ) const;

    /*
   * @param[in] P given point in (x,y,z) coordinates
   * @return true if equal vertices
   */
    bool equals( const ConvexSphericalPolygon& plg, const double prec = eps_ ) const;

    /*
   * @return true:polygon is convex
   */
    bool validate();

    size_t size() const { return size_; }

    void print( std::ostream& ) const;

    friend std::ostream& operator<<( std::ostream& out, const ConvexSphericalPolygon& p ) {
        p.print( out );
        return out;
    }

private:
    /*
   * find next vertex of intersection polygon of plg1 and plg2 polygons
   * @param[out] plg_points collected vertices of intersection polygon
   * @param[in] i starting edge of plg1
   * @param[in] j starting edge of plg2
   */
    int nextIntersect( std::vector<PointXYZ>& plg_points, const ConvexSphericalPolygon& plg1, int i, int j,
                       int inside = 0 ) const;


private:
    std::array<PointXYZ, MAX_SIZE> sph_coords_;
	PointXYZ centroid_;
    size_t size_;
    bool valid_;

    static constexpr double eps_ = 1e-7;  // 1e-8 did not work for i=1,j=13 of test_spherical_geo.cc !!
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
