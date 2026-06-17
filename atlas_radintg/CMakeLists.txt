cmake_minimum_required(VERSION 3.20)

project(atlas_radintg LANGUAGES Fortran CXX VERSION 0.0.0 )

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/atlas_radintg)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#### Dependency

find_package(atlas REQUIRED)
find_package(ecrad REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -finit-real=snan -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-integer=-2147483648")

add_library(atlas_radintg STATIC 
  src/util/MultiFieldCreatorRad.h
  src/util/MultiFieldCreatorRad.cc
  src/radiation/atlas_ecrad_mod.F90
)
target_include_directories(atlas_radintg PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)
target_include_directories(atlas_radintg PUBLIC ${HDF5_C_INCLUDE_DIRS})
target_link_libraries(atlas_radintg PUBLIC atlas atlas_f ecrad.dp ifs.dp ${HDF5_C_LIBRARIES})

add_library(atlas_ifs STATIC
  src/ifs/type_geometry.F90
  src/ifs/yomdim.F90
  src/ifs/yommp.F90
  src/ifs/yomgem.F90
  src/ifs/yomprad.F90
)
target_link_libraries(atlas_ifs atlas_radintg)

#### Example

add_executable(atlas-radintg
  src/programs/atlas-radintg.cc
)
target_link_libraries(atlas-radintg PRIVATE atlas_radintg)
# target_include_directories(atlas-radintg PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)
# target_include_directories(atlas_radintg PRIVATE ${HDF5_C_INCLUDE_DIRS})


add_executable(atlas-radintg-from-ifs
  src/programs/atlas-radintg-from-ifs.F90
)
target_link_libraries(atlas-radintg-from-ifs PRIVATE atlas_radintg atlas_ifs pluto_f)
