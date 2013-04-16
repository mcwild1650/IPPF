#include "calc.h"
#include "tools.h"
#include "field.h"
#include "restart.h"
#include "config.h"
#include "grid.h"
#include "driver.h"
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
  const double amewa=(JET(ia)*((CFG(Nx)+1)/2+1))*SPC(dx);
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
        -FLD(v)[0][j-1]*sqrt(SPC(x)[j-1]*SPC(x)[j-1]+1))/(2*SPC(dx))/FLD(DMsq)[0][j];
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
  int Psi_k=0;
  int My=CFG(My);
  int Nx=CFG(Nx);
  double Psi_tol=1;
  double KappaA=DRV(KappaA);
  double Kappasq=DRV(Kappasq);
  double dxsq=DRV(dxsq);

  field Psi=FLD(Psi);
  field Psi0i=FLD(Psi0i);
  field Omega=FLD(Omega);
  field DMsq=FLD(DMsq);
  
  field temp;
  int i,j;
  while((Psi_tol)>CFG(Tol)) {
   ++(Psi_k);
   (Psi_tol)=0;

   temp=(Psi0i);
   (Psi0i)=(Psi);
   (Psi)=temp;

   for(i=1; i<(My)+1; ++i)
     for(j=1; j<(Nx)+1; ++j) {
       (Psi)[i][j]=(KappaA)*((dxsq)*
                      (Omega)[i][j]*(DMsq)[i][j]+
                      (Psi0i)[i][j+1]+
                      (Psi0i)[i][j-1]+
                      (Kappasq)*((Psi0i)[i+1][j]+
                      (Psi0i)[i-1][j]));
       if(fabs((Psi)[i][j]-(Psi0i)[i][j])>(Psi_tol))
         (Psi_tol)=fabs((Psi)[i][j]-(Psi0i)[i][j]);
     }
  }
  VOL(Psi_k)=Psi_k;
  VOL(Psi_tol)=Psi_tol;
  return 0;
}

static int omega_calc(state* s)
{
  int i,j;
  field Omega=FLD(Omega);
  field Omega0=FLD(Omega0);
  field DMsq=FLD(DMsq);
  field u=FLD(u);
  field v=FLD(v);
  field DM=FLD(DM);
  double alpha=DRV(alpha);
  double alphaX=DRV(alphaX);
  double alphaY=DRV(alphaY);
  double Cxd2=DRV(Cxd2);
  double Cyd2=DRV(Cyd2);
  int My=CFG(My);
  int Nx=CFG(Nx);
  for(i=1;i<My+1;++i)
  for(j=1;j<Nx+1;++j)
  {
    Omega[i][j]=Omega0[i][j]*(1-alpha/(DMsq)[i][j])+
                     Omega0[i][j+1]*
                     (-(Cxd2)*(u)[i][j+1]*(DM)[i][j+1]+
                     alphaX)/(DMsq)[i][j]+
                     Omega0[i][j-1]*((Cxd2)*
                     (u)[i][j-1]*(DM)[i][j-1]+
                     alphaX)/(DMsq)[i][j]+
                     Omega0[i+1][j]*
                     (-(Cyd2)*(v)[i+1][j]*(DM)[i+1][j]+
                     alphaY)/(DMsq)[i][j]+
                     Omega0[i-1][j]*
                     ((Cyd2)*(v)[i-1][j]*(DM)[i-1][j]+
                     alphaY)/(DMsq)[i][j];
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

void init_derived(config* c, space* s, derived* d)
{
  double dx = s->dx;
  double dy = s->dy;
  double dxx = square(dx);
  double dyy = square(dy);
  d->dxsq = dxx;
  d->dysq = dyy;
  double dx2 = 2*dx;
  double dy2 = 2*dy;
  d->dx2 = dx2;
  d->dy2 = dy2;
  double dt = c->dt;
  double Kappa2 = square(dx/dy);
  d->Kappasq = Kappa2;
  d->KappaA = 1/(2*(1.0+Kappa2));
  double Cx = dt/dx;
  double Cy = dt/dy;
  d->Cx = Cx;
  d->Cy = Cy;
  double C = max(Cx,Cy);
  d->C = C;
  double Cx2 = 0.5*Cx;
  double Cy2 = 0.5*Cy;
  d->Cxd2 = Cx2;
  d->Cyd2 = Cy2;
  double Re = c->Re;
  double alphaX = dt/(dxx*Re);
  double alphaY = dt/(dyy*Re);
  double alpha = 2*alphaX+2*alphaY;
  d->alpha = alpha;
  d->alphaX = alphaX;
  d->alphaY = alphaY;
}

void init_fields(config* c, space* s, derived* d, fields* f)
{
  int i,j;
  int My = f->y;
  int Nx = f->x;
  field u = f->u;
  field v = f->u;
  field DM = f->DM;
  field DM2 = f->DMsq;
  field Psi = f->Psi;
  field Omega = f->Omega;
  double* x = s->x;
  double* y = s->y;
  double A = c->A;
  double dyy = d->dysq;
  for(i=0; i<My+2; ++i)
    for(j=0; j<Nx+2; ++j) {
      DM2[i][j]=square(x[j])+square(y[i]);
      DM[i][j]=sqrt(DM2[i][j]);
    }
  for(i=1; i<My+2; ++i)
    for(j=0; j<Nx+2; ++j) {
      v[i][j]=-(y[i]-1)/DM[i][j];
      u[i][j]=(x[j]+A)/DM[i][j];
      Psi[i][j]=(x[j]+A)*(y[i]-1);
      Omega[i][j]=0;
    }
  for(j=0; j<Nx+2; ++j) {
    v[0][j]=-(y[0]-1)/DM[0][j];
    u[0][j]=0;
    Psi[0][j]=(x[j]+A)*(y[0]-1);
    Omega[0][j]=((7*Psi[0][j]-8*Psi[1][j]+Psi[2][j])/(2*dyy))/DM2[0][j];
  }
}
