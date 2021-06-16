/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/Tiles.h"
#include "atlas/grid/detail/partitioner/CubedSpherePartitioner.h"
#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using grid::detail::partitioner::CubedSpherePartitioner;

CASE( "cubedsphere_grid_mesh_field_test" ) {
    // THIS IS TEMPORARY!
    // I expect this will be replaced by some more aggressive tests.

    // Set grid.
    const auto grid = atlas::Grid( "CS-EA-24" );

    atlas::Log::info() << grid->type() << std::endl;
    atlas::Log::info() << grid.size() << std::endl;


    // Set mesh.
    auto meshGen = atlas::MeshGenerator( "cubedsphere" );
    auto mesh    = meshGen.generate( grid );

    // Set functionspace
    auto functionSpace = atlas::functionspace::NodeColumns(
        mesh, atlas::util::Config( "levels", 1 ) | atlas::util::Config( "periodic_points", true ) );

    // Set field
    auto field     = functionSpace.createField<idx_t>( atlas::option::name( "indices" ) );
    auto fieldView = atlas::array::make_view<idx_t, 2>( field );

    for ( idx_t i = 0; i < fieldView.shape()[0]; ++i ) {
        fieldView( i, 0 ) = i;
    }
}

CASE( "cubedsphere_partitioner_test" ) {
    int resolution( 4 );
    std::vector<std::string> grid_names{
        "CS-LFR-C-" + std::to_string( resolution ),
    };
    Grid grid{grid_names[0]};

    using grid::detail::partitioner::CubedSpherePartitioner;

    if ( mpi::size() == 1 ) {
        // factory based constructor
        {
            std::vector<int> globalProcStartPE{0, 0, 0, 1, 1, 1};
            std::vector<int> globalProcEndPE{0, 0, 0, 1, 1, 1};
            std::vector<int> nprocx{1, 1, 1, 1, 1, 1};
            std::vector<int> nprocy{1, 1, 1, 1, 1, 1};

            atlas::util::Config conf;
            conf.set( "starting rank on tile", globalProcStartPE );
            conf.set( "final rank on tile", globalProcEndPE );
            conf.set( "nprocx", nprocx );
            conf.set( "nprocy", nprocy );

            grid::Partitioner partitioner( "cubedsphere", conf );
            grid::Distribution d_cs = partitioner.partition( grid );

    // Set gmsh config.
    auto gmshConfigXy     = atlas::util::Config( "coordinates", "xy" );
    auto gmshConfigXyz    = atlas::util::Config( "coordinates", "xyz" );
    auto gmshConfigLonLat = atlas::util::Config( "coordinates", "lonlat" );

    // Set source gmsh object.
    const auto gmshXy     = atlas::output::Gmsh( "cs_xy_mesh.msh", gmshConfigXy );
    const auto gmshXyz    = atlas::output::Gmsh( "cs_xyz_mesh.msh", gmshConfigXyz );
    const auto gmshLonLat = atlas::output::Gmsh( "cs_lonlat_mesh.msh", gmshConfigLonLat );

    // Write gmsh.
    gmshXy.write( mesh );
    gmshXy.write( field );
    gmshXyz.write( mesh );
    gmshXyz.write( field );
    gmshLonLat.write( mesh );
    gmshLonLat.write( field );
}
}  // namespace test
}  // namespace atlas

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
