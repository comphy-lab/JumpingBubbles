/**
# Boundary Visualization Tool

This program extracts and visualizes boundary data from Basilisk C simulation 
files, specifically designed for analyzing fluid dynamics simulations on 
adaptive octree grids. The tool outputs boundary cell coordinates in a format 
suitable for visualization software.

## Overview

The program reads a Basilisk simulation snapshot file and iterates through 
all cells at the bottom boundary of the computational domain. For each 
boundary cell, it outputs the four corner coordinates of the cell face, 
creating a quadrilateral representation suitable for 3D visualization.

## Authors
- Vatsal Sanjay (vatsalsanjay@gmail.com)
- Youssef
- Physics of Fluids Group

## Usage
```
./bview <simulation_file>
```

## Output Format
The program outputs to standard error (stderr) the coordinates of boundary 
cell faces as quadrilaterals. Each quadrilateral is defined by four corner 
points, with a blank line separating each cell.

## Dependencies
- Basilisk C framework
- Octree grid module for adaptive mesh refinement
- Navier-Stokes solver (centered scheme)
- Volume fraction tracking module
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

char filename[80];

/**
### Main Function

Entry point for the boundary visualization tool. Processes command-line 
arguments and extracts boundary data from the specified simulation file.

## Parameters
- `a`: Number of command-line arguments
- `arguments`: Array of command-line argument strings
  - `arguments[0]`: Program name (unused)
  - `arguments[1]`: Path to the Basilisk simulation file to process

## Process Flow
1. Extracts the simulation filename from command-line arguments
2. Restores the simulation state from the specified file
3. Iterates through all cells at the bottom boundary
4. Outputs quadrilateral coordinates for each boundary cell face

## Output Details
For each boundary cell, the program outputs four corner coordinates in the 
following order:
- Bottom-left corner: (x - Δ/2, y, z - Δ/2)
- Top-left corner: (x - Δ/2, y, z + Δ/2)
- Top-right corner: (x + Δ/2, y, z + Δ/2)
- Bottom-right corner: (x + Δ/2, y, z - Δ/2)

Where (x, y, z) represents the cell center and Δ (Delta) is the cell size.

## Notes
- Output is directed to stderr to separate from standard output
- Each quadrilateral is followed by two newlines for separation
- The coordinate ordering creates counter-clockwise vertices when viewed 
  from outside the domain (positive y direction)
*/
int main(int a, char const *arguments[]) {
  // Extract simulation filename from command-line arguments
  sprintf(filename, "%s", arguments[1]);
  
  // Restore simulation state from the specified file
  restore(file = filename);
  
  // Iterate through all cells at the bottom boundary
  foreach_boundary(bottom) {
    // Output quadrilateral coordinates for visualization
    // Format: four corner points of the boundary cell face
    fprintf(ferr, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n\n",
            x - Delta/2., y, z - Delta/2.,  // Bottom-left corner
            x - Delta/2., y, z + Delta/2.,  // Top-left corner
            x + Delta/2., y, z + Delta/2.,  // Top-right corner
            x + Delta/2., y, z - Delta/2.); // Bottom-right corner
  }
  
  return 0;
}