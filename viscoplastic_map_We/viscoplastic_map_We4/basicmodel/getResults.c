/* Title: Getting Energy
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last Update: Sep 06 2021
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "tag.h"
#include "adapt_wavelet_limited.h"
#include "heights.h"

scalar f[];
double Oh, mu1, J, tauy, epsOh, epsJ, We;
#define RHO21 (1e-3)
double rho1 = 1., rho2 = RHO21;
#ifndef rho
# define rho(f) (clamp(f,0,1)*(rho1 - rho2) + rho2)
#endif

int Level;
char filename[80], nameEnergy[80];

int main(int a, char const *arguments[])
{
  sprintf(filename, "%s", arguments[1]);
  Oh = atof(arguments[2]);
  We = atof(arguments[3]);
  J = atof(arguments[4]);
  Level = atof(arguments[5]);

  sprintf(filename, "%s", arguments[1]);
  restore(file = filename);

  // boundary conditions
  u.t[left] = dirichlet(0.0);
  u.n[left] = dirichlet(0.0);
  f[left] = dirichlet(0.0);

  p[right] = dirichlet(0.);
  u.n[right] = neumann(0.);

  p[top] = dirichlet(0.);
  u.n[top] = neumann(0.);
  mu1 = Oh / sqrt(We);
  tauy = J / We;
  f.prolongation = refine_bilinear;
  boundary((scalar *){f, u.x, u.y, p});

  // fprintf(ferr, "Oh %3.2e, We %g\n", Oh, We);
  // return 1;

  // tag all liquid parts starts
  scalar d[];
  double threshold = 1e-4;
  foreach ()
  {
    d[] = (f[] > threshold);
  }
  int n = tag(d), size[n];
  for (int i = 0; i < n; i++)
  {
    size[i] = 0;
  }
  foreach_leaf()
  {
    if (d[] > 0)
    {
      size[((int)d[]) - 1]++;
    }
  }
  int MaxSize = 0;
  int MainPhase = 0;
  for (int i = 0; i < n; i++)
  {
    // fprintf(ferr, "%d %d\n",i, size[i]);
    if (size[i] > MaxSize)
    {
      MaxSize = size[i];
      MainPhase = i + 1;
    }
  }
  // tag all liquid parts ends

  scalar sf[];
  foreach ()
    sf[] = (4. * f[] +
            2. * (f[0, 1] + f[0, -1] + f[1, 0] + f[-1, 0]) +
            f[-1, -1] + f[1, -1] + f[1, 1] + f[-1, 1]) /
           16.;
  sf.prolongation = refine_bilinear;
  boundary({sf});

  /*
  Do calculations start
  */
  double ke = 0., vcm = 0., vcm1 = 0., vcmR = 0., vcmR1 = 0., vc = 0., vc1 = 0., dm = 0., dm1 = 0.;
  face vector s[];
  s.x.i = -1;
  double yMax = 0;
  double xMax = 0;
  double xMin = HUGE;
  // double vR = 0., uH = 0.;
  double vR = 0., xTP = 0.;
  double uH = 0., yTP = 0.;

  foreach (){
    //vcm and vcm1 are used for the calculation of the normal force.
    vcm += (2 * pi * y) * (clamp(sf[], 0., 1.) * u.x[]) * sq(Delta);
    vcmR += (2 * pi * y) * (clamp(sf[], 0., 1.) * u.y[]) * sq(Delta);
    dm += (2 * pi * y) * (clamp(sf[], 0., 1.)) * sq(Delta);
    ke += sq(Delta) * (2 * pi * y) * (sq(u.x[]) + sq(u.y[])) * rho(sf[]) / 2.;
    if (d[] == MainPhase)
    {
      vcm1 += (2 * pi * y) * (clamp(sf[], 0., 1.) * u.x[]) * sq(Delta);      
      vcmR1 += (2 * pi * y) * (clamp(sf[], 0., 1.) * u.y[]) * sq(Delta);
      dm1 += (2 * pi * y) * (clamp(sf[], 0., 1.)) * sq(Delta);      
    }               
  }
  vc=vcm/dm;
  vc1=vcm1/dm1;

  foreach (){
    if (f[] > 1e-6 && f[] < 1. - 1e-6 && d[] == MainPhase)
    {
      coord n1 = facet_normal(point, f, s);
      double alpha1 = plane_alpha(f[], n1);
      coord segment1[2];
      if (facets(n1, alpha1, segment1) == 2)
      {
        double x1 = x + (segment1[0].x + segment1[1].x) * Delta / 2.;
        double y1 = y + (segment1[0].y + segment1[1].y) * Delta / 2.;
       
        if (y1 > yMax)
        {
          yMax = y1;
          xTP = x1;          
        }
        if (y < 0.01)
        {
          if (x1 > xMax)
          {
            xMax = x1;
            yTP = y1;            
          }
        }
        if (x1 < xMin)
        {
          xMin = x1;
        }
      }
    }    
  }
  vR = interpolate(u.y, xTP, yMax);
  uH = interpolate(u.x, xMax, yTP);

  double pdatum = 0, wt = 0;
  foreach_boundary(top, reduction(+:pdatum), reduction(+:wt)){
    pdatum += 2*pi*y*p[]*(Delta);
    wt += 2*pi*y*(Delta);
  }
  if (wt >0){
    pdatum /= wt;
  }

  double pforce = 0.,pforce1 = 0.;
  foreach_boundary(left, reduction(+:pforce), reduction(+:pforce1)){
    if (Delta<1.0/sq(Level-1)){
      pforce += 2*pi*y*(p[]-pdatum)*(Delta);
      pforce1 += 2*pi*y*(p[])*(Delta);
    }
  }

  /*
  Do calculations end
  */

  fprintf(ferr, "%6.5e %d %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e\n",
          t,n,ke,vcm,vc,vcm1,vc1,vcmR,vcmR1,yMax,xMax,xMin,vR,uH,pforce,pforce1);
}