/**
 * Version 3.0
 * Author: Vatsal Sanjay, Xiangyu Zhang
 * Last updated: April 8, 2025

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
 You can use: conserving.h as well. Even without it, It was still able to conserve the total energy (also momentum?) of the system if I adapt based on curvature and vorticity/deformation tensor norm (see the adapt even). I had to smear the density and viscosity anyhow because of the sharp ratios in liquid (Bingham) and the gas.
*/
#include "navier-stokes/conserving.h"
#include "tension.h"
// #include "reduced.h"
#include "distance.h"
#include "tag.h"

// gas properties!
#define RHO21 (1e-3)
#define MU21 (1e-2)
#define MINlevel 3

// Error tolerancs
#define fErr (1e-3)     // error tolerance in VOF
#define KErr (1e-2)     // error tolerance in KAPPA
#define VelErr (1e-2)   // error tolerances in velocity
#define DissErr (1e-2)  // error tolerances in dissipation
#define OmegaErr (1e-2) // error tolerances in vorticity

// Distance and radius of drop calculations
#define Xdist (1.1)
#define R2Drop(x, y) (sq(x - Xdist) + sq(y))

// boundary conditions
u.t[left] = dirichlet(0.0);
u.n[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

p[right] = dirichlet(0.);
u.n[right] = neumann(0.);

p[top] = dirichlet(0.);
u.n[top] = neumann(0.);

int MAXlevel;
double We, Oh, J, Bo, epsilon_1;
double tmax, Ldomain;
double tsnap = 0.01;
char nameOut[80], dumpFile[80];

int main(int argc, char const *argv[])
{
  origin(0., 0.);
  init_grid(1 << 5);
  MAXlevel = atoi(argv[1]);  // 10
  J = atof(argv[2]);         // 0 for Newtonian.
  We = atof(argv[3]);        // 10
  Oh = atof(argv[4]);        // 0.01
  Bo = atof(argv[5]);        // 0 without considering density
  epsilon_1 = atof(argv[6]); // 1e-2
  tmax = atof(argv[7]);      // 10
  Ldomain = atof(argv[8]);   // 8
  DT = atof(argv[9]);        // 1e-3
  CFL = atof(argv[10]);      // 1e-3

  L0 = Ldomain;
  NITERMAX = 1000;
  NITERMIN = 2;
  // TOLERANCE = 1e-4;

  char comm[80];
  sprintf(comm, "mkdir -p intermediate");
  system(comm);
  // To mimic a liquid-gas free surface, we set the Ohnesorge number of the surrounding gas Ohs 10−5 and the density ratio 10−3.
  rho1 = 1., rho2 = RHO21;
  mu1 = Oh / sqrt(We), mu2 = MU21 * Oh / sqrt(We);
  // mu1 = Oh / sqrt(We), mu2 = 1e-5;
  f.sigma = 1.0 / We;
  tauy = J / We;
  // G.x = -Bo / We; // uncomment only if Gravity is needed!
  // when t<tsnap, epsilon=0.1. and when t>tsnap epsilon=epsilon_1
  epsilon = 0.1;
  fprintf(ferr, "We,Oh,J,MAXlevel,epsilon,mu2,DT,Ldomain,tmax,Bo,CFL\n");
  fprintf(ferr, "%g,%g,%g,%d,%4.3e,%g,%g,%g,%g,%g,%g\n", We, Oh, J, MAXlevel, epsilon_1, mu2, DT, Ldomain, tmax, Bo, CFL);
  run();
}

event init(t = 0)
{
  if (!restore(file = "dump", list = all))
  {
    // log
    refine(R2Drop(x, y) < 1.05 && (level < MAXlevel));
    fraction(f, 1. - R2Drop(x, y));
    foreach ()
    {
      u.x[] = -1.0 * f[];
      u.y[] = 0.0;
    }
    boundary((scalar *){f, u.x, u.y});
  }
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++)
{
  if (t < 1e-2)
  {
    adapt_wavelet((scalar *){f, u.x, u.y},
                  (double[]){fErr, VelErr, VelErr},
                  MAXlevel, MINlevel);
  }
  else
  {
    epsilon = epsilon_1; // update the value of epsilon
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
    adapt_wavelet((scalar *){f, KAPPA, u.x, u.y, D2c},
                  (double[]){fErr, KErr, VelErr, VelErr, DissErr},
                  MAXlevel, MINlevel);
    unrefine(x > 0.95 * Ldomain);
  }
  double time_now = perf.t / 60.0;
  fprintf(ferr, "%g,%d,%g,%g\n", t, i, time_now, epsilon);
}

double t_last = 0.0;
double DeltaT = 0.0;
int count_run = 0;
event postProcess(t += 0.01)
{
  double ke = 0., xMin = Ldomain;
  face vector s[];
  s.x.i = -1;
  foreach (reduction(+ : ke) reduction(min : xMin))
  {
    ke += (2 * pi * y) * (sq(u.x[]) + sq(u.y[])) * rho(f[]) / 2. * sq(Delta);
    if (f[] > 1e-6 && f[] < 1. - 1e-6)
    {
      coord n1 = facet_normal(point, f, s);
      double alpha1 = plane_alpha(f[], n1);
      coord segment1[2];
      if (facets(n1, alpha1, segment1) == 2)
      {
        double x1 = x + (segment1[0].x + segment1[1].x) * Delta / 2.;
        // double y1 = y + (segment1[0].y + segment1[1].y) * Delta / 2.;       
        if ((x1 < xMin))
        {
          xMin = x1;
        }
      }
    }
  }

  if (isnan(ke)) {
    fprintf(ferr, "Ke is NaN\n");
    return 1;
  }
  else{
    // this is designed for the recalculate with finer time step, only run 5 step for the fine step

    // we check the value of kenetic energy.
    if (ke < 1e-6 || ke > 2.1)
    {
      fprintf(ferr, "Ke_Error, Exit...\n");
      return 1;
    }
    else
    {
      p.nodump = true;
      dump(file = "dump");
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
        fprintf(fp1, "t,i,Cell,Wallclocktime(min),CPUtime(min),ke,xMin\n");
        fflush(fp1);
      }
      fp1 = fopen("log_run", "a");
      fprintf(fp1, "%g,%d,%ld,%g,%g,%g,%g\n", t, i, grid->tn, perf.t / 60.0, DeltaT, ke, xMin);
      fflush(fp1);
    }
    // stop condition
    if ((t > tmax - tsnap) || (t > 0.5 && (ke < 1e-6 || (xMin > 0.04))))
    {
      fprintf(ferr, "Success, Reach Max time Or Kinetic energy is too small Or droplet bounce off. Exiting...\n");
      sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
      dump(file = nameOut);
      return 1;
    }
  }

  
}

event snapshot(t += 0.01; t <= 50)
{
  // when p.nodump=false the restore function will not work well. So when needed we use p.nodump=false, and the force can be calculated by d(mv)/dt.
  // p.nodump = false;
  p.nodump = true;
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file = nameOut);
}
