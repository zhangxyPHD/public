/**
 * Version 2.0
 * Author: Vatsal Sanjay
 * Last updated: Oct 11, 2024

# Introduction:
We investigate the classical problem of VP drop impacting a solid surface.
# Numerical code
Id 1 is for the Viscoplastic liquid drop, and Id 2 is Newtonian gas.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED // Smear density and viscosity jumps
/**
To model Viscoplastic liquids, we use a modified version of [two-phase.h](http://basilisk.fr/src/two-phase.h). [two-phaseVP.h](two-phaseVP.h) contains these modifications.
*/
#include "two-phaseVP-HB.h"
/**
 You can use: conserving.h as well. Even without it, I was still able to conserve the total energy (also momentum?) of the system if I adapt based on curvature and vorticity/deformation tensor norm (see the adapt even). I had to smear the density and viscosity anyhow because of the sharp ratios in liquid (Bingham) and the gas.
*/
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "distance.h"
#include "tag.h"
#include "adapt_wavelet_limited.h"

#define MINlevel 3

// Error tolerancs
#define fErr (1e-3)     // error tolerance in VOF
#define KErr (1e-2)     // error tolerance in KAPPA
#define VelErr (1e-2)   // error tolerances in velocity
#define DissErr (1e-2)  // error tolerances in dissipation
#define OmegaErr (1e-2) // error tolerances in vorticity

// gas properties!
#define RHO21 (1e-3)
#define MU21 (1e-2)

// Distance and radius of drop calculations
#define Xdist (1.02)
#define R2Drop(x, y) (sq(x - Xdist) + sq(y))

// boundary conditions
u.t[left] = dirichlet(0.0);
u.n[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

p[right] = dirichlet(0.);
u.n[right] = neumann(0.);

p[top] = dirichlet(0.);
u.n[top] = neumann(0.);

int MAXlevel, MAXlevel_f;
double We, Oh, J, Bo;
double tmax, Ldomain;
double tsnap = 0.01;
char nameOut[80], resultsName[80], dumpFile[80];
double ifmodel1 = 0.0;
int main(int argc, char const *argv[])
{
  origin(0., 0.);
  init_grid(1 << 5);
  MAXlevel = atoi(argv[1]);   // 10
  MAXlevel_f = atoi(argv[2]); // 10
  J = atof(argv[3]);          // 0 for Newtonian.
  We = atof(argv[4]);         // 10
  Oh = atof(argv[5]);         // 0.01
  Bo = atof(argv[6]);         // 0 without considering density
  epsilon = atof(argv[7]);    // 1e-2
  tmax = atof(argv[8]);       // 10
  Ldomain = atof(argv[9]);    // 8
  DT = atof(argv[10]);        // 1e-4
  CFL = atof(argv[11]);       // 0.1
  sprintf(resultsName, "%s", argv[12]);

  L0 = Ldomain;
  NITERMAX = 300;
  TOLERANCE = 1e-4;

  char comm[80];
  sprintf(comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1., rho2 = RHO21;
  mu1 = Oh / sqrt(We), mu2 = MU21 * Oh / sqrt(We);
  f.sigma = 1.0 / We;
  tauy = J / We;
  G.x = -Bo / We; // uncomment only if Gravity is needed!
  run();
}

event init(t = 0)
{
  if (!restore(file = "dump", list = all))
  {
    // log
    fprintf(ferr, "We,Oh,J,MAXlevel,epsilon,MU21,DT,Ldomain,tmax,Bo\n");
    fprintf(ferr, "%g,%g,%g,%d,%4.3e,%g,%g,%g,%g,%g\n", We, Oh, J, MAXlevel, epsilon, MU21, DT, Ldomain, tmax, Bo);
    refine(R2Drop(x, y) < 1 && (level < MAXlevel));
    refine(fabs(R2Drop(x, y) - 1) < 0.1 && (level < MAXlevel_f));
    fraction(f, 1. - R2Drop(x, y));
    foreach ()
    {
      u.x[] = -1.0 * f[];
      u.y[] = 0.0;
    }
    boundary((scalar *){f, u.x, u.y});
  }
}

int refRegion(double x, double y, double z)
{
  return (MAXlevel);
}

int refRegion_f(double x, double y, double z)
{
  return (MAXlevel_f);
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++)
{
  if (t > 1e-2)
  {
    scalar KAPPA[];
    curvature(f, KAPPA);
    scalar D2c[];
    foreach ()
    {
      double D11 = (u.y[0, 1] - u.y[0, -1]) / (2 * Delta);
      double D22 = (u.y[] / max(y, 1e-20));
      double D33 = (u.x[1, 0] - u.x[-1, 0]) / (2 * Delta);
      double D13 = 0.5 * ((u.y[1, 0] - u.y[-1, 0] + u.x[0, 1] - u.x[0, -1]) / (2 * Delta));
      double D2 = (sq(D11) + sq(D22) + sq(D33) + 2.0 * sq(D13));
      D2c[] = f[] * D2;
    }
    if (ifmodel1 < 0.5)
    {
      adapt_wavelet_limited((scalar *){u.x, u.y, D2c}, (scalar *){f, KAPPA}, (double[]){VelErr, VelErr, DissErr}, (double[]){fErr, KErr}, refRegion, refRegion, MINlevel);
    }
    else
    {
      adapt_wavelet_limited((scalar *){u.x, u.y, D2c}, (scalar *){f, KAPPA}, (double[]){VelErr, VelErr, DissErr}, (double[]){fErr, KErr}, refRegion, refRegion_f, MINlevel);
    }
    unrefine(x > 0.95 * Ldomain);
  }

  if (t > 5)
  {
    DT = 1e-3;
  }
  double time_now = perf.t / 60.0;
  fprintf(ferr, "%g,%d,%g\n", t, i, time_now);
}

double t_last = 0.0;
double DeltaT = 0.0;
int count_run = 0;
double ke_last = 100;
event postProcess(t += tsnap)
{
  double ke = 0., xMin = Ldomain;
  foreach (reduction(+ : ke) reduction(min : xMin))
  {
    ke += (2 * pi * y) * (sq(u.x[]) + sq(u.y[])) * rho(f[]) / 2. * sq(Delta);
    if (f[] > 1e-6)
    {
      if ((x < xMin))
      {
        xMin = x;
      }
    }
  }
  if (ke > ke_last)
  {
    ifmodel1 = 1.0;
  }
  // only run 5 step for the fine step
  if (DT < 5e-5 && count_run > 5)
  {
    return 1;
  }
  count_run = count_run + 1;

  /**
  ## Dumping snapshots, when ke > 0 && ke < 100, We store the dumpfile to ensure there is no problem for the current time step.
  */
  if (ke > 0 && ke < 2.1)
  {
    p.nodump = false;
    dump(file = "dump");
    sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
    dump(file = nameOut);
  }
  else
  {
    fprintf(ferr, "Ke_Error, Exit...\n");
    return 1;
  }
  // log
  DeltaT = perf.t / 60.0 - t_last;
  t_last = perf.t / 60.0;
  static FILE *fp1;
  if (pid() == 0)
  {
    if (i == 0)
    {
      fp1 = fopen("log_run", "w");
      fprintf(fp1, "t,i,Cell,Wallclocktime(min),CPUtime(min),ke\n");
      fflush(fp1);
    }
    fp1 = fopen("log_run", "a");
    fprintf(fp1, "%g,%d,%ld,%g,%g,%g\n", t, i, grid->tn, perf.t / 60.0, DeltaT, ke);
    fflush(fp1);
  }
  // stop condition
  if ((t > tmax - tsnap) || (t > 0.5 && (ke < 1e-5 || (xMin > 0.04))))
  {
    char comm[256];
    sprintf(comm, "cp log_run ../Results_Running/log_%s.csv", resultsName);
    system(comm);
    fprintf(ferr, "Kinetic energy is too small or droplet bounce off. Exiting...\n");
    return 1;
  }
  ke_last = ke;
}

event end(t = tmax)
{
}