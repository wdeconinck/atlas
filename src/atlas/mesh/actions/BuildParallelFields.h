/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date   June 2014

#ifndef BuildParallelFields_h
#define BuildParallelFields_h

#include "atlas/internals/atlas_config.h"

namespace atlas {
namespace mesh {
    class Mesh;
    class Nodes;
} }

namespace atlas {
namespace mesh {
namespace actions {

/*
 * Build all parallel fields in the mesh
 *  - calls build_nodes_parallel_fields()
 */
void build_parallel_fields( Mesh& mesh );

/*
 * Build parallel fields for the "nodes" function space if they don't exist.
 * - glb_idx:    create unique indices for non-positive values
 * - partition:  set to parallel::mpi::comm().rank() for negative values
 * - remote_idx: rebuild from scratch
 */
void build_nodes_parallel_fields( mesh::Nodes& nodes );

/*
 * Build parallel fields for the "edges" function space if they don't exist.
 * - glb_idx:    create unique indices for non-positive values
 * - partition:  set to partition of node with lowest glb_idx
 * - remote_idx: rebuild from scratch
 *
 *  TODO: Make sure that the edge-partition is at least one of the partition numbers of the
 *        neighbouring elements.
 *        Because of this problem, the size of the halo should be set to 2 instead of 1!!!
 */
void build_edges_parallel_fields( Mesh& mesh );

void renumber_nodes_glb_idx (mesh::Nodes& nodes);

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
#define mesh_Nodes mesh::Nodes
extern "C"
{
  void atlas__build_parallel_fields (Mesh* mesh);
  void atlas__build_nodes_parallel_fields (mesh_Nodes* nodes);
  void atlas__build_edges_parallel_fields (Mesh* mesh);
  void atlas__renumber_nodes_glb_idx (mesh_Nodes* nodes);
}
#undef mesh_Nodes
// ------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas


#endif // BuildParallelFields_h
