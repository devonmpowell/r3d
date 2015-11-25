# `r3d`

---

## About:

Routines for fast, geometrically robust clipping operations and analytic volume/moment computations 
over polytopes in 2D and 3D (as well as experimental ND). This software forms the kernel for an exact 
general remeshing scheme. Also includes physically conservative voxelization 
(volume sampling) of 3D polyhedra to a Cartesian grid.

As described in 
[Powell & Abel (2015)](http://www.sciencedirect.com/science/article/pii/S0021999115003563) and
[LA-UR-15-26964](la-ur-15-26964.pdf). For information about the API itself, see
[LA-UR-15-26964](la-ur-15-26964.pdf).

---

### Features:

- Robustly clip polytopes against planes.

- Compute volumes and moments over polytopes using the optimal recursive method of
[Koehl (2012)](https://www.computer.org/csdl/trans/tp/2012/11/ttp2012112158.pdf).

- Voxelize 3D polyhedra onto a Cartesian grid by calculating the exact coordinate moments
  of the intersections between the tetrahedron and each underlying grid cell.

- Utility functions for orientation tests, box initialization, conversion between polyhedral
  representations, and more.

- A set of rigorous unit-tests, located in `tests`. These tests also serve as examples of how to
  use `r3d`. 

- All declarations and documentation are located in `r3d.h` and `r2d.h`.

---

### Usage:

- To build, type

`make`

- To compile into your code,

`#include <r3d.h>`
or
`#include <r2d.h>`

- To link,

`-lr3d`

---

### Licensing: 

`r3d.c`, `r3d.h`, `r2d.c`, `r2d.h`, `rNd.c`, `rNd.h`, and contents of `tests` 
Copyright (C) 2015, DOE and Los Alamos National Security, LLC.

Contents of `deprecated` Copyright (C) 2015, The Board of Trustees of the Leland Stanford Junior University, 
through SLAC National Accelerator Laboratory (subject to receipt of any required approvals 
from the U.S. Dept. of Energy). 

See source file headers for full license text. All code is open-source, subject to terms of the
respective license. We request that you cite 
[Powell & Abel (2015)](http://www.sciencedirect.com/science/article/pii/S0021999115003563) and
[LA-UR-15-26964](la-ur-15-26964.pdf) when using this code for research purposes.


