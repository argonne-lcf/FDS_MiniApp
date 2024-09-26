# Loop3D

This repository contains several variants of the 'Loop3D' mini-app,
with the purpose being to explore GPU offload strategies for the Fire Dynamics Simulator (FDS). These codes can be broadly categorized as:

- Automatic data migration with map clauses.
- Unified shared memory.
- Explicit control of data, with C bindings.

Further categorizations:
- Single grid/single task/shared memory. These codes are basically the original Loop3D, but with some OpenMP pragmas thrown in. Run on a single compute node.
- Single grid/MPI. Decomposes the rectangular domain among an arbitrary number of MPI ranks, with GPU-aware message-passing for halo exchange. This is basically just a distributed version of the original Loop3D. Can run across multiple compute nodes.
- Single mesh/single task/shared memory. Like the original Loop3D, but with field variables encapsulated in a MESH data type as in FDS. Custom mappers for automatically handling data migration for this derived type.
- Multi-mesh/single task/shared memory. Like the above version, but with a list of meshes, and the "point-to-mesh" feature from FDS. Uses unified shared memory.

All programs were compiled and tested on Sunspot, a testbed for the Aurora supercomputer (Argonne National Laboratory).

## References
[The Fire Dynamics Simulator](https://pages.nist.gov/fds-smv/) 

## Contributors
- Aristotle Martin (Duke University) ([@aris134](https://github.com/aris134))
- Marcos Vanella (National Institutes of Standards and Technology) ([@marcosvanella](https://github.com/marcosvanella))
- Saumil S. Patel (Argonne National Laboratory) 
