/**
 * @file JumpingBubbles.c
 * @brief Simulation of jumping bubbles sitting on a substrate using Basilisk C. 
 * @author Vatsal
 * @version 1.5 
 * @date Jan 5, 2025

# Changelog (v1.5) Jan 5, 2025
- Extended support for arbitary contact angle. 

 * This code simulates two bubbles coalescing and jumping off a substrate. It uses an adaptive octree grid for spatial discretization  and a two-phase flow model with surface tension. The code reads an STL file for the bubble geometry and captures interface evolution using the volume-of-fluid (VOF) method.
 *   Simulation Parameters:
 *   Oh: Ohnesorge number (ratio of viscous to inertial-capillary forces) 
 *   MAXlevel: Maximum refinement level, controlling the finest grid
 *   tmax: Maximum simulation time (default: 1e-2)
 *   
 *   Other Parameters kept fixed:
 *   Ldomain: Length of the simulation domain (default: 4)
 *   Rho21, Mu21: Density and viscosity ratios of the second phase with respect to the first (default: 1e-3)
 *   tsnap: Time interval between solution snapshots
 *   tsnap2: Time interval for log outputs
 *   FILTERED: Enable density and viscosity jump smoothing (default: false)
 *   fErr, KErr, VelErr: Adaptive refinement error tolerances for the interface, curvature, and velocity fields
 *   MINlevel: Minimum refinement level, controlling the coarsest grid
 *
 * Boundary Conditions:
 *   - On the bottom boundary, the tangential and radial velocity components are set to zero.
 *   - The interface fraction is set to 1 at the bottom, indicating liquid presence in cells touching that boundary.
 * The simulation proceeds via standard Basilisk events:
 *   - init: Restores from a dump file if available; otherwise constructs the initial interface from an STL file
 *   - adapt: Adaptive mesh refinement based on interface, curvature, and velocity field errors
 *   - writingFiles: Dumps solution snapshots at specified intervals
 *   - logWriting: Records kinetic energy to a log file at specified intervals
 *
 * Implementation details:
 *   - Basilisk C's navier-stokes/centered solver is used for momentum conservation
 *   - The two-phase flow model is combined with surface tension through the tension.h module
 *   - The distance() and fractions() functions handle geometric input and construct the volume fraction field
 *
 * Example:
 *   Running on a Linux system with OpenMP support:
 *   qcc -O2 -Wall -disable-dimensions -fopenmp JumpingBubbles.c -o JumpingBubbles -lm
 *   ./JumpingBubbles
 */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
//#define FILTERED
#include "contact-fixed.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"

#if !_MPI
#include "distance.h"
#endif

#include "reduced.h"

#define MINlevel 2                                              // maximum level

#define tsnap (1e-2)
#define tsnap2 (1e-4)
// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-4)                            // error tolerances in velocity

#define Mu21 (1.00e-3)
#define Rho21 (1.00e-3)

// domain
#define Ldomain 4                                // Dimension of the domain

// boundary conditions
u.t[bottom] = dirichlet(0.);
u.r[bottom] = dirichlet(0.);
uf.t[bottom] = dirichlet(0.);
uf.r[bottom] = dirichlet(0.);

double theta0, patchR;
vector h[];
h.t[bottom] = contact_angle (theta0*pi/180.);
h.r[bottom] = contact_angle (theta0*pi/180.);

double tmax, Oh;
int MAXlevel; // maximum level
char nameOut[80];

int main() {

#if !_MPI
  tmax = 1e-2;
#else
  tmax = 2e0;
#endif

  Oh = 0.0066; // <\eta_l/sqrt(\rho_l*\gamma*Requiv)>
  MAXlevel = 8;
  theta0 = 15; // contact angle in degrees
  patchR = 0.184; // Rcont/Requiv
  Bo = 0.016; // Bo = \rho_l*g*Requiv^2/\gamma

  init_grid (1 << MINlevel);
  L0=Ldomain;
  fprintf(ferr, "tmax = %g. Oh = %g\n",tmax, Oh);
  rho1 = 1.0; mu1 = Oh;
  rho2 = Rho21; mu2 = Mu21*Oh;
  
  f.height = h;
  f.sigma = 1.0;

  G.y = -Bo;

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  sprintf(dumpFile, "restartFile");

  run();

}

event init(t = 0){
#if _MPI // this is for supercomputers without OpenMP support
  if (!restore (file = dumpFile)){
    fprintf(ferr, "Cannot restored from a dump file!\n");
  }
#else // note that distance.h is incompatible with OpenMPI. So, the below code should not be used with MPI
  if(!restore (file = dumpFile)){
    char filename[60];
    sprintf(filename,"InitialCondition.stl");
    FILE * fp = fopen (filename, "r");
    if (fp == NULL){
      fprintf(ferr, "There is no file named %s\n", filename);
      return 1;
    }
    coord * p = input_stl (fp);
    fclose (fp);
    coord min, max;

    bounding_box (p, &min, &max);
    fprintf(ferr, "xmin %g xmax %g\nymin %g ymax %g\nzmin %g zmax %g\n", min.x, max.x, min.y, max.y, min.z, max.z);
    fprintf(ferr, "x0 = %g, y0 = %g, z0 = %g\n", 0., - 1.0, (min.z+max.z)/2.);
    origin (0., - 1.0 - 0.025, (min.z+max.z)/2.);

    scalar d[];
    distance (d, p);
    while (adapt_wavelet ((scalar *){f, d}, (double[]){1e-6, 1e-6*L0}, MAXlevel).nf);
    vertex scalar phi[];
    foreach_vertex(){
      phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1] +
  	     d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
    }
    fractions (phi, f);

    foreach () {
      foreach_dimension(){
        u.x[] = 0.0;
      }
    }

    dump (file = "dumpInit"); // to check the initial condition
    fprintf(ferr, "Done with initial condition!\n");
    // return 1;
  }
#endif
}

event adapt(i++) {
  adapt_wavelet_limited ((scalar *){f, u.x, u.y, u.z, h.x, h.y, h.z},
     (double[]){fErr, VelErr, VelErr, VelErr, hErr, hErr, hErr},
      MAXlevel, MINlevel);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax+tsnap) {
  dump (file = dumpFile);
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (t = 0; t += tsnap2; t <= tmax+tsnap) {

  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += 0.5*(sq(u.x[]) + sq(u.y[]) + sq(u.z[]))*clamp(1.-f[], 0., 1.)*cube(Delta);
  }

  if (pid() == 0){
    static FILE * fp;
    if (i == 0) {
      fprintf (ferr, "i dt t ke\n");
      fp = fopen ("log", "w");
      fprintf (fp, "i dt t ke\n");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    } else {
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
  }

}
