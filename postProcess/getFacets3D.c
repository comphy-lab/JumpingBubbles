/**
# Facet Output Visualization Tool

This module provides functionality for extracting and outputting interface
facets from volume fraction fields in computational fluid dynamics simulations.
It is designed to work with the Basilisk C framework for solving partial
differential equations on adaptive octree/quadtree grids.

## Overview

The tool reads volume fraction data from saved simulation files and extracts
the interface geometry using facet reconstruction algorithms. The extracted
facets can be visualized or processed for further analysis of multiphase
flow interfaces.

## Authors
- Vatsal Sanjay (vatsalsanjay@gmail.com)
- Physics of Fluids Group

## Dependencies
- Basilisk C framework
- Octree grid implementation
- Navier-Stokes centered solver
- Volume fraction utilities
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f[];
char filename[80];

/**
### output_facets_v2

Extracts and outputs the facets (interface segments/polygons) from a volume
fraction field. This function implements a marching squares/cubes-like
algorithm to reconstruct interfaces from volume fraction data.

The function iterates through all grid cells and identifies those containing
an interface (where the volume fraction is between 0 and 1). For each
interface cell, it:
1. Computes the interface normal using neighboring cell data
2. Calculates the plane constant α for the linear interface approximation
3. Extracts the intersection points between the interface and cell edges
4. Outputs the facet vertices to the specified file

#### Parameters
- `c`: Volume fraction scalar field containing the interface data
- `fp`: Output file pointer for writing facet coordinates (default: stdout)
- `s`: Face vector field for surface tension calculations (default: {{-1}})

#### Output Format
The output format depends on the spatial dimension:
- **1D**: Single x-coordinate of the interface position
- **2D**: Line segments as pairs of (x,y) coordinates, blank line separated
- **3D**: Polygon vertices as (x,y,z) coordinates, blank line separated

#### Implementation Notes
- Only processes cells with volume fraction in range (1e-6, 1-1e-6) to
  avoid numerical issues near pure phases
- Uses the PLIC (Piecewise Linear Interface Calculation) method for
  interface reconstruction
- Facet extraction tolerance is set to 1.1 for robust polygon generation
  in 3D cases
*/
trace
void output_facets_v2(scalar c, FILE * fp = stdout, 
                      face vector s = {{-1}})
{
  foreach()
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = facet_normal(point, c, s);
      double alpha = plane_alpha(c[], n);
#if dimension == 1
      fprintf(fp, "%g\n", x + Delta*alpha/n.x);
#elif dimension == 2
      coord segment[2];
      if (facets(n, alpha, segment) == 2)
        fprintf(fp, "%g %g\n%g %g\n\n", 
                x + segment[0].x*Delta, y + segment[0].y*Delta, 
                x + segment[1].x*Delta, y + segment[1].y*Delta);
#else // dimension == 3
      coord v[12];
      int m = facets(n, alpha, v, 1.1);
      for (int i = 0; i < m; i++)
        fprintf(fp, "%g %g %g\n",
                x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      if (m > 0)
        fputc('\n', fp);
#endif
    }

  fflush(fp);
}

/**
### main

Entry point for the facet extraction tool. Processes a saved simulation
file specified as a command-line argument and outputs the extracted facets
to standard error.

#### Command Line Usage
```
./program_name simulation_file.dump
```

#### Parameters
- `a`: Number of command line arguments
- `arguments`: Array of command line argument strings
  - `arguments[1]`: Path to the input simulation file

#### Processing Steps
1. Reads the filename from command line arguments
2. Sets boundary conditions for the volume fraction field
3. Restores the simulation state from the saved file
4. Ensures proper prolongation operators for adaptive mesh refinement
5. Calls output_facets_v2 to extract and output interface facets

#### Boundary Conditions
- Bottom boundary: Dirichlet condition with f = 1 (fully filled)

#### Tree Grid Handling
For adaptive (tree) grids, the function ensures that the volume fraction
field uses the appropriate `fraction_refine` prolongation operator. This
is crucial for maintaining conservation properties during mesh refinement
and coarsening operations.

#### Output
Interface facets are written to standard error (ferr) in a format suitable
for visualization with tools like gnuplot or paraview.
*/
int main(int a, char const *arguments[]) {
  sprintf(filename, "%s", arguments[1]);

  f[bottom] = dirichlet(1.);
  restore(file=filename);

  #if TREE
  // Ensure proper prolongation for volume fraction fields on adaptive grids
  void (* prolongation)(Point, scalar) = f.prolongation;
  if (prolongation != fraction_refine) {
    f.prolongation = fraction_refine;
    f.dirty = true;
  }
  #endif // TREE

  output_facets_v2(f, ferr);
}