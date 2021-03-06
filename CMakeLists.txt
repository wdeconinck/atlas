# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

############################################################################################
# Atlas

cmake_minimum_required( VERSION 3.6 FATAL_ERROR )
if( POLICY CMP0074 )
  cmake_policy( SET CMP0074 NEW )
  # This policy allows to search for packages with <package>_ROOT variables
  #                        (only supported with CMake 3.12 and above)
  # This policy can be removed once cmake_minimum_required( VERSION 3.12 ) is used
endif()

find_package( ecbuild 3.1.0 REQUIRED )

################################################################################
# Initialise project Atlas

project( atlas LANGUAGES CXX )

set( ATLAS_BITS_GLOBAL 64 ) # bits used to define a global index (atlas::gidx_t)
set( ATLAS_BITS_LOCAL 32 )  # bits used to define a local index  (atlas::idx_t)

################################################################################
# Required packages

ecbuild_find_package( NAME eckit VERSION 1.4.0 REQUIRED )
ecbuild_debug( "   eckit_FEATURES : [${eckit_FEATURES}]" )

################################################################################
# Features that can be enabled / disabled with -DENABLE_<FEATURE>

include( features/FORTRAN )
include( features/MPI )
include( features/OMP )
include( features/FFTW )
include( features/TRANS )
include( features/TESSELATION )
include( features/GRIDTOOLS_STORAGE )
include( features/ACC )
include( features/EIGEN )
include( features/PROJ )
include( features/SANDBOX )
include( features/CLANG_TIDY )
include( features/INCLUDE_WHAT_YOU_USE )
include( features/DOCS )

################################################################################
# sources

include( atlas_compile_flags )

add_subdirectory( src )

add_subdirectory( doc )

################################################################################
# Export and summarize

include( atlas_export )

ecbuild_print_summary()

