/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iomanip>
#include <iostream>

#include "eckit/geometry/Sphere.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/runtime/Exception.h"
#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/NormaliseLongitude.h"

#define DEBUG_OUTPUT 1
#define DEBUG_OUTPUT_DETAIL 1

namespace atlas {
namespace util {

bool approx_eq( const double& v1, const double& v2, const double& tol ) {
    return eckit::types::is_approximately_equal( v1, v2, tol );
}

bool approx_eq( const PointXYZ& v1, const PointXYZ& v2, const double& tol ) {
    //return approx_eq( v1[0], v2[0], t ) && approx_eq( v1[1], v2[1], t ) && approx_eq( v1[2], v2[2], t );
    return PointXYZ::norm( v1 - v2 ) < tol;
}

bool approx_eq_null( const PointXYZ& v1, const double& tol ) {
    //return approx_eq( v1[0], 0., t ) && approx_eq( v1[1], 0., t ) && approx_eq( v1[2], 0., t );
    return PointXYZ::norm( v1 ) < tol;
}

double cart_diff( const PointXYZ& v1, const PointXYZ& v2 ) {
    return PointXYZ::norm( v1 - v2 );
}

PointLonLat sph_to_lonlat( const PointXYZ& p ) {
    PointLonLat pp;
    eckit::geometry::Sphere::convertCartesianToSpherical( 1., p, pp );
    return pp;
}

//------------------------------------------------------------------------------------------------------

ConvexSphericalPolygon::ConvexSphericalPolygon() : valid_( false ), size_( 0 ), area_( 0 ) {}

ConvexSphericalPolygon::ConvexSphericalPolygon( const std::vector<PointLonLat>& points ) :
    size_( points.size() ), area_( 0 ) {
    ATLAS_ASSERT( size_ < MAX_SIZE );
    eckit::geometry::Sphere::convertSphericalToCartesian( 1., points[0], sph_coords_[0] );
    centroid_  = sph_coords_[0];
    size_t isp = 1;
    for ( size_t i = 1; i < points.size() - 1; ++i ) {
        eckit::geometry::Sphere::convertSphericalToCartesian( 1., points[i], sph_coords_[isp] );
        if ( approx_eq( sph_coords_[isp], sph_coords_[isp - 1], 1e-10 ) ) {
            continue;
        }
        centroid_ = centroid_ + sph_coords_[isp];
        ++isp;
    }
    eckit::geometry::Sphere::convertSphericalToCartesian( 1., points[points.size() - 1], sph_coords_[isp] );
    if ( approx_eq( sph_coords_[isp], sph_coords_[0], 1e-10 ) or
         approx_eq( sph_coords_[isp], sph_coords_[isp - 1], 1e-10 ) ) {
    }
    else {
        centroid_ = centroid_ + sph_coords_[isp];
        ++isp;
    }
    size_  = isp;
    valid_ = size_ > 2;
    if ( valid_ ) {
        ATLAS_ASSERT( validate() );
        ATLAS_ASSERT( not approx_eq_null( centroid_, 1e-10 ) );
        centroid_ = PointXYZ::div( centroid_, PointXYZ::norm( centroid_ ) );
        compute_area();
        cell_radius_ = 0.;
        for ( size_t i = 0; i < size_; ++i ) {
            cell_radius_ = std::max( cell_radius_, cart_diff( sph_coords_[i], centroid_ ) );
        }
    }
}

ConvexSphericalPolygon::ConvexSphericalPolygon( const std::vector<PointXYZ>& points, const bool debug ) :
    size_( points.size() ), area_( 0 ) {
#if DEBUG_OUTPUT
    if ( debug ) {
        std::cout << " \n\nConvexSphericalPolygon got points: ";
        for ( int i = 0; i < points.size(); ++i ) {
            std::cout << std::setprecision( 10 ) << sph_to_lonlat( points[i] ) << ",";
        }
        ( std::cout << "\n" ).flush();
    }
#endif
    ATLAS_ASSERT( size_ < MAX_SIZE );
    if ( size_ < 3 ) {
        size_  = 0;
        valid_ = false;
        return;
    }
    idx_t isp = 0;
    idx_t i   = 0;
    idx_t j;
    idx_t k;
    bool search_3pts = true;
    for ( ; i < points.size() && search_3pts; ++i ) {
        const PointXYZ& P0 = points[i];
        for ( j = i + 1; j < points.size() && search_3pts; ++j ) {
            const PointXYZ& P1 = points[j];
            if ( approx_eq( P0, P1, 1e-10 ) ) {
                continue;
            }
            for ( k = j + 1; k < points.size() && search_3pts; ++k ) {
                const PointXYZ& P2 = points[k];
                if ( approx_eq( P1, P2, 1e-10 ) or approx_eq( P0, P2, 1e-10 ) ) {
                    continue;
                }
                if ( leftOf( P2, P0, P1, 1e-14, debug ) ) {
                    sph_coords_[isp++] = P0;
                    sph_coords_[isp++] = P1;
                    sph_coords_[isp++] = P2;
                    search_3pts        = false;
                }
            }
        }
    }
    if ( search_3pts ) {
        valid_ = false;
        size_  = 0;
        return;
    }
#if DEBUG_OUTPUT
    if ( debug ) {
        std::cout << " ConvexSphericalPolygon 3 FIRST: " << std::setprecision( 10 ) << sph_to_lonlat( sph_coords_[0] )
                  << " " << sph_to_lonlat( sph_coords_[1] ) << " " << sph_to_lonlat( sph_coords_[2] ) << " "
                  << "\n";
    }
#endif

    centroid_ = sph_coords_[i] + sph_coords_[j] + sph_coords_[k];
    for ( ; k < points.size() - 1; ++k ) {
        if ( approx_eq( points[k], sph_coords_[isp - 1], 1e-10 ) or
             ( not leftOf( points[k], sph_coords_[isp - 2], sph_coords_[isp - 1], 1e-14, debug ) ) ) {
            continue;
        }
        sph_coords_[isp] = points[k];
        centroid_        = centroid_ + sph_coords_[isp];
        isp++;
    }
    const PointXYZ& Pl2 = sph_coords_[isp - 2];
    const PointXYZ& Pl1 = sph_coords_[isp - 1];
    const PointXYZ& P0  = sph_coords_[0];
    const PointXYZ& P   = points[size_ - 1];
    if ( ( not approx_eq( P, P0, 1e-14 ) ) and ( not approx_eq( P, Pl1, 1e-14 ) ) and leftOf( P, Pl2, Pl1, 1e-14, debug ) and
         leftOf( P0, Pl1, P, 1e-14, debug ) ) {
        sph_coords_[isp] = P;
        centroid_        = centroid_ + P;
        ++isp;
    }
    size_  = isp;
    valid_ = size_ > 2;
    if ( valid_ ) {
        ATLAS_ASSERT( not approx_eq_null( centroid_, 1e-10 ) );
        centroid_ = PointXYZ::div( centroid_, PointXYZ::norm( centroid_ ) );
        compute_area();
    }
}

bool ConvexSphericalPolygon::validate() {
    if ( valid_ ) {
        for ( int i = 0; i < size(); i++ ) {
            int ni                = ( i != size() - 1 ? i + 1 : 0 );
            int nni               = ( ni != size() - 1 ? ni + 1 : 0 );
            const PointXYZ& P     = sph_coords_[i];
            const PointXYZ& nextP = sph_coords_[ni];
            ATLAS_ASSERT( std::abs( PointXYZ::dot( P, P ) - 1. ) < 1e-14 );
            ATLAS_ASSERT( not approx_eq( P, PointXYZ::mul( nextP, -1. ), 1e-10 ) );
            valid_ = valid_ && leftOf( sph_coords_[nni], P, nextP, 1e-14, 0 );
        }
    }
    return valid_;
}

bool ConvexSphericalPolygon::equals( const ConvexSphericalPolygon& plg, const double deg_prec ) const {
    const double le = 2. * std::sin( M_PI * deg_prec / 360. );
    if ( ( not plg.valid_ ) || ( not valid_ ) || size_ != plg.size() ) {
        ( Log::info() << " ConvexSphericalPolygon::equals == not compatible\n" ).flush();
        return false;
    }
    int offset = 0;
    for ( ; offset < size_; ++offset ) {
        if ( PointXYZ::norm( plg.sph_coords_[0] - sph_coords_[offset] ) < le ) {
            break;
        }
    }
    if ( offset == size_ ) {
        ( Log::info() << "ConvexSphericalPolygon::equals == no point equal"
                      << "\n" )
            .flush();
        return false;
    }

    for ( int j = 0; j < size_; j++ ) {
        int idx   = ( offset + j ) % size_;
        auto dist = PointXYZ::norm( plg.sph_coords_[j] - sph_coords_[idx] );
        if ( dist > le ) {
            ( Log::info() << "  ConvexSphericalPolygon::equals == point distance " << dist << "\n" ).flush();
            return false;
        }
    }
    return true;
}

// note: unit sphere!
// I. Todhunter (1886), Paragr. 99
void ConvexSphericalPolygon::compute_area() {
    area_ = 0.;
    if ( size_ < 3 ) {
        return;
    }
    int valid_angle = 0;
    for ( int i = 0; i < size_; i++ ) {
        int im1            = ( i != 0 ) ? i - 1 : size_ - 1;
        int ip1            = ( i != size_ - 1 ) ? i + 1 : 0;
        const PointXYZ& pl = sph_coords_[im1];
        const PointXYZ& p  = sph_coords_[i];
        const PointXYZ& pr = sph_coords_[ip1];
        if ( not leftOf( pr, pl, p, -1e-14, 0 ) ) {
            continue;
        }
        PointXYZ ppl          = PointXYZ( PointXYZ::cross( p, pl ) );
        PointXYZ ppr          = PointXYZ( PointXYZ::cross( p, pr ) );
        const double ppl_norm = PointXYZ::norm( ppl );
        const double ppr_norm = PointXYZ::norm( ppr );
        //ATLAS_ASSERT( ppl_norm > std::numeric_limits<double>::epsilon()
        //	&& ppr_norm > std::numeric_limits<double>::epsilon() );
        if ( ppl_norm < std::numeric_limits<double>::epsilon() or ppr_norm < std::numeric_limits<double>::epsilon() ) {
            continue;
            Log::info() << " p, pl, pr : " << p << " " << pl << " " << pr << "\n";
            Log::info() << " ppl, ppr : " << ppl << " " << ppr << "\n";
            Log::info() << " ppl_norm, ppr_norm : " << ppl_norm << " " << ppr_norm << "\n";
            Log::info() << " compute_area fails for plg: ";
            for ( int aa = 0; aa < size_; aa++ ) {
                Log::info() << " " << sph_to_lonlat( sph_coords_[aa] );
            }
            Log::info().flush();
            ATLAS_ASSERT( false );
        }
        valid_angle++;
        double s1p = PointXYZ::dot( ppl, ppr ) / ( ppl_norm * ppr_norm );
        s1p        = std::acos( ( s1p < 0. ? -1 : 1 ) * std::min( 1., std::abs( s1p ) ) );
        area_ += s1p;
    }
    area_ = ( valid_angle > 2 ? std::abs( area_ ) + M_PI * ( 2. - valid_angle ) : 0. );
}

// return 0:P_right_of_[p1,p2], -1:overlap_of_[P,p1]_and_[P,p2], 1:P_left_of_[p1,p2]
int ConvexSphericalPolygon::leftOf( const PointXYZ& P, const PointXYZ& p1, const PointXYZ& p2, const double tol,
                                    const int debug ) const {
    const PointXYZ& cp = PointXYZ( PointXYZ::cross( p1, p2 ) );
    //if ( approx_eq_null( cp, 1e-15 ) ) {
    //	return true;
    //}
    double cpP = PointXYZ::dot( cp, P );
#if DEBUG_OUTPUT
    if ( debug ) {
        ( Log::info() << " p, p1, p2: " << sph_to_lonlat( P ) << ", " << sph_to_lonlat( p1 ) << ", "
                      << sph_to_lonlat( p2 ) << "\n" )
            .flush();
        ( Log::info() << " p1 x p2, p * (p1 x p2): " << cp << ", " << cpP << "\n" ).flush();
    }
#endif
    return ( cpP > -tol );
}

double ConvexSphericalPolygon::norm_max( const PointXYZ& p, const PointXYZ& q ) {
    double n01 = std::max( std::abs( p[0] - q[0] ), std::abs( p[1] - q[1] ) );
    return std::max( n01, std::abs( p[2] - q[2] ) );
}

bool ConvexSphericalPolygon::between( const PointXYZ& p, const PointXYZ& p1, const PointXYZ& p2, const int debug ) {
    PointXYZ p12 = PointXYZ::cross( p1, p2 );
    double p12n  = PointXYZ::norm( p12 ) - std::numeric_limits<double>::epsilon();
    double pp1n  = PointXYZ::norm( p - p1 ) - 1e+7 * std::numeric_limits<double>::epsilon();
    double pp2n  = PointXYZ::norm( p - p2 ) - 1e+7 * std::numeric_limits<double>::epsilon();
    if ( p12n < 0. && pp1n < 0. ) {
        return true;
    }
    p12             = PointXYZ::div( p12, p12n );
    const double dp = 1e+3 * std::numeric_limits<double>::epsilon() - std::abs( PointXYZ::dot( p, p12 ) );
    if ( dp < 0. ) {
#if DEBUG_OUTPUT_DETAIL
        if ( debug ) {
            ( Log::info() << "NOT between: point not in the plane, " << dp << "\n" ).flush();
        }
#endif
        return false;
    }
    double pp = PointXYZ::norm( p1 - p2 );
    pp        = std::min( pp - pp1n, pp - pp2n );
#if DEBUG_OUTPUT_DETAIL
    if ( debug ) {
        ( Log::info() << "  between pp, pp2n = " << pp << ", " << pp2n << "\n" ).flush();
    }
#endif
    return ( pp > 0. && pp2n > 0. );
}

PointXYZ ConvexSphericalPolygon::common( const PointXYZ& s1, const PointXYZ& s2, const PointXYZ& p1, const PointXYZ& p2,
                                         const int debug ) {
    PointXYZ s  = static_cast<PointXYZ>( PointXYZ::cross( s1, s2 ) );
    PointXYZ p  = static_cast<PointXYZ>( PointXYZ::cross( p1, p2 ) );
    PointXYZ sp = static_cast<PointXYZ>( PointXYZ::cross( s, p ) );

#if 0
	Log::info() << "s1-s2 = " << s1 - s2 << ", |s1-s2|: " << PointXYZ::norm( s1 - s2 ) << "\n";
	Log::info() << "p1-p2 = " << p1 - p2 << ", |p1-p2|: " << PointXYZ::norm( p1 - p2 ) << "\n";
	Log::info() << " s = " << s << ", |s| = " << PointXYZ::norm( s ) << "\n";
	Log::info() << " s = " << sph_to_lonlat(s) << "\n";
	Log::info() << " p = " << p << ", |p| = " << PointXYZ::norm( p ) << "\n";
	Log::info() << " p = " << sph_to_lonlat(p) << "\n";
	Log::info() << " sp = " << sp << ", |sp| = " << PointXYZ::norm( sp ) << "\n";
	Log::info() << " sp = " << sph_to_lonlat(sp) << "\n\n";
#endif

    double sp_norm = PointXYZ::norm( sp ) - std::numeric_limits<double>::epsilon();
#if DEBUG_OUTPUT_DETAIL
    if ( debug ) {
        ( Log::info() << " Parallel: " << sp_norm << " < 0 ?\n" ).flush();
    }
#endif
    if ( sp_norm > 0. ) {
        sp = PointXYZ::div( sp, PointXYZ::norm( sp ) );
#if DEBUG_OUTPUT_DETAIL
        if ( debug ) {
            ( Log::info() << "Test sp = " << sph_to_lonlat( sp ) << "\n" ).flush();
            ( Log::info() << "	p1-sp-p2: " << between( sp, p1, p2 ) << "\n" ).flush();
        }
#endif
        if ( between( sp, p1, p2 ) ) {
#if DEBUG_OUTPUT_DETAIL
            if ( debug ) {
                ( Log::info() << ",	GOT  p: " << sph_to_lonlat( sp ) << "\n" ).flush();
            }
#endif
            return sp;
        }
        sp = PointXYZ::mul( sp, -1 );
#if DEBUG_OUTPUT_DETAIL
        if ( debug ) {
            ( Log::info() << "Test sp = " << sph_to_lonlat( sp ) << "\n" ).flush();
        }
#endif
        if ( between( sp, p1, p2 ) ) {
#if DEBUG_OUTPUT_DETAIL
            if ( debug ) {
                ( Log::info() << ",	GOT -p: " << sph_to_lonlat( sp ) << "\n" ).flush();
            }
#endif
            return sp;
        }
#if DEBUG_OUTPUT_DETAIL
        if ( debug ) {
            ( Log::info() << " Intersection not on [p1, p2).\n" ).flush();
        }
#endif
        return PointXYZ( {0, 0, 0} );
    }
    else {
#if DEBUG_OUTPUT_DETAIL
        if ( debug ) {
            ( Log::info() << " Overlap. \n" ).flush();
        }
#endif
        return PointXYZ( {1, 1, 1} );
    }
}

// @return -1: overlap with one of polygon edges,
//  		0: no_intersect,
//          1 + (id of this polygon's segment intersecting [s1,s2]): otherwise
int ConvexSphericalPolygon::intersect( const PointXYZ& s1, const PointXYZ& s2, PointXYZ& I, int start,
                                       const bool debug ) const {
    for ( int i = start; i < size_; i++ ) {
        const int id0      = i;
        const int id1      = ( id0 == size_ - 1 ) ? 0 : id0 + 1;
        const PointXYZ& p1 = sph_coords_[id0];
        const PointXYZ& p2 = sph_coords_[id1];
#if DEBUG_OUTPUT
        if ( debug ) {
            ( Log::info() << "	** edge: " << sph_to_lonlat( p1 ) << " " << sph_to_lonlat( p2 ) << "\n" ).flush();
        }
#endif
        I = common( s1, s2, p1, p2 );
        if ( I[0] == 0 && I[1] == 0 && I[2] == 0 ) {
            // intersection not on [p1,p2)
            continue;
        }
        if ( I[0] == 1 && I[1] == 1 ) {
            // overlap
            return -1;
        }
        return 1 + id0;
    }
    return 0;
}

void ConvexSphericalPolygon::clip( const PointLonLat& s1, const PointLonLat& s2, const int debug ) {
    PointXYZ p1;
    PointXYZ p2;
    eckit::geometry::Sphere::convertSphericalToCartesian( 1., s1, p1 );
    eckit::geometry::Sphere::convertSphericalToCartesian( 1., s2, p2 );
    clip( p1, p2, debug );
}

void ConvexSphericalPolygon::clip( const PointXYZ& s1, const PointXYZ& s2, const int debug ) {
    ConvexSphericalPolygon deleteme = *this;
    if ( PointXYZ::norm( s1 - s2 ) < 1e-14 ) {
        return;
    }
    PointXYZ i1;
    PointXYZ i2 = PointXYZ( {0, 0, 0} );
    int f1;
    int f2 = -1;
    f1     = -1 + intersect( s1, s2, i1, 0, debug );
    if ( f1 >= 0 ) {
        // intersection with no-overlap
        f2 = -1 + intersect( s1, s2, i2, f1 + 1, debug );
    }
#if DEBUG_OUTPUT
    if ( debug ) {
        ( Log::info() << "f1, f2: " << f1 << " " << f2 << "\n" ).flush();
    }
#endif
    bool no_intersection = ( f1 == -1 );
    bool touches         = ( f1 == -2 ) && ( f2 == -1 );

    bool point_intrs      = ( f1 < 0 or ( f1 >= 0 && f2 < 0 ) );
    const PointXYZ& s1xs2 = PointXYZ( PointXYZ::cross( s1, s2 ) );
    double triangle_left  = PointXYZ::dot( s1xs2, sph_coords_[0] + sph_coords_[1] + sph_coords_[2] );
    if ( point_intrs && triangle_left > 0. ) {
        // edge overlap, this polygon on the inside
#if DEBUG_OUTPUT
        if ( debug ) {
            ( Log::info() << "this polygon on the inside\n" ).flush();
        }
#endif
        return;
    }
    else if ( point_intrs ) {
        // edge overlap, this polygon on the outside
#if DEBUG_OUTPUT
        if ( debug ) {
            ( Log::info() << "this polygon on the outside\n" ).flush();
        }
#endif
        valid_ = false;
        size_  = 0;
        return;
    }
#if DEBUG_OUTPUT
    if ( debug ) {
        ( Log::info() << "f1, f2, i1, i2: " << f1 << " " << f2 << " " << sph_to_lonlat( i1 ) << ", "
                      << sph_to_lonlat( i2 ) << "\n" )
            .flush();
    }
#endif
    int f1n       = ( f1 != size_ - 1 ) ? f1 + 1 : 0;
    int f2n       = ( f2 != size_ - 1 ) ? f2 + 1 : 0;
    int i1_is_f1  = ( PointXYZ::norm( i1 - sph_coords_[f1] ) < 1e-8 );
    int i1_is_f1n = ( PointXYZ::norm( i1 - sph_coords_[f1n] ) < 1e-8 );
    int i2_is_f2  = ( PointXYZ::norm( i2 - sph_coords_[f2] ) < 1e-8 );
    int i2_is_f2n = ( PointXYZ::norm( i2 - sph_coords_[f2n] ) < 1e-8 );

    std::vector<PointXYZ> cp_sph_coords;
    for ( int i = 0; i <= f1; i++ ) {
        cp_sph_coords.emplace_back( sph_coords_[i] );
    }
    if ( not i1_is_f1 && not i1_is_f1n ) {
        cp_sph_coords.emplace_back( i1 );
    }
    for ( int i = f1 + 1; i <= f2; i++ ) {
        cp_sph_coords.emplace_back( sph_coords_[i] );
    }
    if ( not i2_is_f2 && not i2_is_f2n ) {
        cp_sph_coords.emplace_back( i2 );
    }
    for ( int i = f2 + 1; i < size_; i++ ) {
        cp_sph_coords.emplace_back( sph_coords_[i] );
    }

    int id = 0;
    for ( int i = 0; i < cp_sph_coords.size(); i++ ) {
        PointXYZ& cpi = cp_sph_coords[i];
        if ( leftOf( cpi, s1, s2, 1e-14, debug ) ) {
            sph_coords_[id++] = cp_sph_coords[i];
        }
    }
    size_ = id;

#if DEBUG_OUTPUT
    if ( debug ) {
        ( Log::info() << " new size " << size_ << "\n" ).flush();
        ( Log::info() << " New plg: " ).flush();
        for ( int i = 0; i < size_; i++ ) {
            Log::info() << " " << sph_to_lonlat( sph_coords_[i] );
        }
        ( Log::info() << "\n\n" ).flush();
    }
#endif
    // check
    for ( int i = 0; i < size_; i++ ) {
        if ( PointXYZ::norm( sph_coords_[i] ) < 1e-10 ) {
            Log::info() << " clip == size too large when clipping\n";
            for ( int i = 0; i < deleteme.size_; i++ ) {
                Log::info() << " " << sph_to_lonlat( deleteme.sph_coords_[i] );
            }
            Log::info() << "\n with " << sph_to_lonlat( s1 ) << " " << sph_to_lonlat( s2 ) << "\n";
            Log::info() << " clip got size = " << size_ << " from size = " << deleteme.size_ << "\n";
            for ( int i = 0; i < size_; i++ ) {
                Log::info() << " " << sph_to_lonlat( sph_coords_[i] );
            }
            ATLAS_ASSERT( false );
        }
    }
}

// intersect a polygon with this polygon
// @param[in] pol clipping polygon
// @param[out] intersecting polygon
ConvexSphericalPolygon ConvexSphericalPolygon::intersect( const ConvexSphericalPolygon& plg, const int debug ) const {
    ConvexSphericalPolygon obj = *this;
    for ( int i = 0; i < plg.size_ && bool( obj ); i++ ) {
        const PointXYZ& s1 = plg.sph_coords_[i];
        const PointXYZ& s2 = plg.sph_coords_[( i != plg.size_ - 1 ) ? i + 1 : 0];
#if DEBUG_OUTPUT
        if ( debug ) {
            ( Log::info() << std::setprecision( 8 ) << "\n 	plg " << i << ": " << obj << "\n" ).flush();
            ( Log::info() << " 	now clip with " << sph_to_lonlat( s1 ) << " " << sph_to_lonlat( s2 ) << "\n" ).flush();
        }
#endif
        obj.clip( s1, s2, debug );
    }
    obj.area_ = 0;
    if ( bool( obj ) ) {
        obj.centroid_ = PointXYZ( {0, 0, 0} );
        for ( int i = 0; i < size_; i++ ) {
            obj.centroid_ = obj.centroid_ + obj.sph_coords_[i];
        }
        obj.centroid_ = PointXYZ::div( obj.centroid_, PointXYZ::norm( obj.centroid_ ) );
        obj.compute_area();
    }
    return obj;
}

void ConvexSphericalPolygon::print( std::ostream& out ) const {
    out << "{";
    for ( size_t i = 0; i < size(); ++i ) {
        if ( i > 0 ) {
            out << ",";
        }
        PointLonLat ip_ll;
        eckit::geometry::Sphere::convertCartesianToSpherical( 1., sph_coords_[i], ip_ll );
        out << ip_ll;
    }
    out << "}";
}


//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
