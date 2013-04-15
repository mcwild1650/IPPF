#include "calc.h"
#include "tools.h"
#include "field.h"
#include "restart.h"
#include "config.h"
#include "grid.h"
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

static double field_max_difference(int Nx,int My,field A, field B)
{
  int i,j;
  double max_diff=0;
  for(i=0;i<My;++i)
  for(j=0;j<Nx;++j)
  {
    if(fabs(A[i][j]-B[i][j])>max_diff)
      max_diff=fabs(A[i][j]-B[i][j]);
  }
  return max_diff;
}

static int set_boundaries(state* s)
{
  int i,j;
  const double amewa=(JET(ia)*((CFG(Nx)+1)/2+1))*DRV(dx);
  const double f=sin(2*PI*JET(freq)*VOL(time));

  // Upper and Lower BCs
  for(j=0; j<CFG(Nx)+2; ++j) {
    FLD(Psi)[0][j]=0;
    if(j>=JET(ia)-1 && j<=JET(ib)-1)
      FLD(Psi)[0][j]=(-JET(c0)*(-0.5*amewa*sqrt(amewa*amewa+1)-0.5*sinh(amewa)
          +0.5*SPC(x)[j]*sqrt(SPC(x)[j]*SPC(x)[j]+1)+0.5*sinh(SPC(x)[j])))*f;
    if(j>JET(ib)-1)
      FLD(Psi)[0][j]=FLD(Psi)[0][JET(ib)-1];
    FLD(Omega)[0][j]=(7*FLD(Psi)[0][j]-8*FLD(Psi)[1][j]+
                     FLD(Psi)[2][j])/(2*DRV(dy2))/FLD(DMsq)[0][j];
    FLD(u)[0][j]=0;
    FLD(v)[0][j]=0;
    if(j>JET(ia)-1 && j<JET(ib)-1)
      FLD(v)[0][j]=JET(c0)*f;
    if(j>=JET(ia)-2 && j<=JET(ib))
      FLD(Omega)[0][j]+=(FLD(v)[0][j+1]*sqrt(SPC(x)[j+1]*SPC(x)[j+1]+1)
        -FLD(v)[0][j-1]*sqrt(SPC(x)[j-1]*SPC(x)[j-1]+1))/(2*DRV(dx))/FLD(DMsq)[0][j];
    FLD(Omega)[CFG(My)+1][j]=0;
    FLD(Psi)[CFG(My)+1][j]=(SPC(x)[j]+CFG(A))*(SPC(y)[CFG(My)+1]-1);
    FLD(u)[CFG(My)+1][j]=(SPC(x)[j]+CFG(A))/FLD(DM)[CFG(My)+1][j];
    FLD(v)[CFG(My)+1][j]=-(SPC(y)[CFG(My)+1]-1)/FLD(DM)[CFG(My)+1][j];
  }

  // Side BCs
  for(i=1; i<CFG(My)+1; ++i) {
    if(i<CFG(IBL)) {
      FLD(Omega)[i][0]=FLD(Omega)[i][1];
      FLD(Psi)[i][0]=FLD(Psi)[i][1];
      FLD(u)[i][0]=FLD(u)[i][1];
      FLD(v)[i][0]=FLD(v)[i][1];
      FLD(Omega)[i][CFG(Nx)+1]=FLD(Omega)[i][CFG(Nx)];
      FLD(Psi)[i][CFG(Nx)+1]=FLD(Psi)[i][CFG(Nx)];
      FLD(u)[i][CFG(Nx)+1]=FLD(u)[i][CFG(Nx)];
      FLD(v)[i][CFG(Nx)+1]=FLD(v)[i][CFG(Nx)];
    }
    else {
      FLD(Omega)[i][0]=0;
      FLD(Psi)[i][0]=(SPC(x)[0]+CFG(A))*(SPC(y)[i]-1);
      FLD(u)[i][0]=(SPC(x)[0]+CFG(A))/FLD(DM)[i][0];
      FLD(v)[i][0]=-(SPC(y)[i]-1)/FLD(DM)[i][0];
      FLD(Omega)[i][CFG(Nx)+1]=0;
      FLD(Psi)[i][CFG(Nx)+1]=(SPC(x)[CFG(Nx)+1]+CFG(A))*(SPC(y)[i]-1);
      FLD(u)[i][CFG(Nx)+1]=(SPC(x)[CFG(Nx)+1]+CFG(A))/FLD(DM)[i][CFG(Nx)+1];
      FLD(v)[i][CFG(Nx)+1]=-(SPC(y)[i]-1)/FLD(DM)[i][CFG(Nx)+1];
    }
  }
  return 0;
}

static int velocity_calc(state* s)
{
  // Calculate velocities
  int i,j;
  for(i=1; i<CFG(My)+1; ++i)
  for(j=1; j<CFG(Nx)+1; ++j)
  {
    FLD(u)[i][j]=(FLD(Psi)[i+1][j]-FLD(Psi)[i-1][j])/DRV(dx2)/FLD(DM)[i][j];
    FLD(v)[i][j]=-(FLD(Psi)[i][j+1]-FLD(Psi)[i][j-1])/DRV(dx2)/FLD(DM)[i][j];
  }
  return 0;
}

static int psi_calc(state* s)
{
  VOL(Psi_k)=0;
  field temp;
  int i,j;
  while(VOL(Psi_tol)>CFG(Tol)) {
   ++VOL(Psi_k);
   VOL(Psi_tol)=0;

   temp=FLD(Psi0i);
   FLD(Psi0i)=FLD(Psi);
   FLD(Psi)=temp;

   for(i=1; i<CFG(My)+1; ++i)
     for(j=1; j<CFG(Nx)+1; ++j) {
       FLD(Psi)[i][j]=DRV(KappaA)*(DRV(dxsq)*
                      FLD(Omega)[i][j]*FLD(DMsq)[i][j]+
                      FLD(Psi0i)[i][j+1]+
                      FLD(Psi0i)[i][j-1]+
                      DRV(Kappasq)*(FLD(Psi0i)[i+1][j]+
                      FLD(Psi0i)[i-1][j]));
       if(fabs(FLD(Psi)[i][j]-FLD(Psi0i)[i][j])>VOL(Psi_tol))
         VOL(Psi_tol)=fabs(FLD(Psi)[i][j]-FLD(Psi0i)[i][j]);
     }
  }
  return 0;
}

static int omega_calc(state* s)
{
  int i,j;
  for(i=1;i<CFG(My)+1;++i)
  for(j=1;j<CFG(Nx)+1;++j)
  {
    FLD(Omega)[i][j]=FLD(Omega0)[i][j]*(1-DRV(alpha)/FLD(DMsq)[i][j])+
                     FLD(Omega0)[i][j+1]*
                     (-DRV(Cxd2)*FLD(u)[i][j+1]*FLD(DM)[i][j+1]+
                     DRV(alphaX))/FLD(DMsq)[i][j]+
                     FLD(Omega0)[i][j-1]*(DRV(Cxd2)*
                     FLD(u)[i][j-1]*FLD(DM)[i][j-1]+
                     DRV(alphaX))/FLD(DMsq)[i][j]+
                     FLD(Omega0)[i+1][j]*
                     (-DRV(Cyd2)*FLD(v)[i+1][j]*FLD(DM)[i+1][j]+
                     DRV(alphaY))/FLD(DMsq)[i][j]+
                     FLD(Omega0)[i-1][j]*
                     (DRV(Cyd2)*FLD(v)[i-1][j]*FLD(DM)[i-1][j]+
                     DRV(alphaY))/FLD(DMsq)[i][j];
  }
  return 0;
}

int initialize_all(state* s)
{
  return 0;
}

int run_time_step(state* s)
{
  return 0;
}
