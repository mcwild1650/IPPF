#include "calc.h"
#include "tools.h"
#include "field.h"
#include "restart.h"
#include "config.h"
#include "grid.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static void copy_boundaries(int Nx, int My, field a, field b)
{
  int i,j;
  for (i=1;i<=My;++i) {
    b[i][0]=a[i][0];
    b[i][Nx+1]=a[i][Nx+1];
  }
  for (j=1;j<=Nx;++j) {
    b[0][j]=a[0][j];
    b[My+1][j]=a[My+1][j];
  }
}

static int setBoundaries(
    config* c,
    space* s,
    double t,
    fields* fld,
    derived* d)
{
  int i,j;
  jet_config* jet = &(c->jet);
  int ia = jet->ia;
  int ib = jet->ib;
  double c0 = jet->c0;
  double freq = jet->ia;
  int Nx = c->Nx;
  int My = c->My;
  double A = c->A;
  int IBL = c->IBL;
  double dx = s->dx;
  double* x = s->x;
  double* y = s->y;
  field Psi = fld->Psi;
  field Omega = fld->Omega;
  field u = fld->u;
  field v = fld->v;
  field DMsq = fld->DMsq;
  field DM = fld->DM;
  const double amewa=(ia*((Nx+1)/2+1))*dx;
  const double f=sin(2*PI*freq*t);
  double dy2 = d->dy2;
  // Upper and Lower BCs
  for(j=0; j<Nx+2; ++j) {
    Psi[0][j]=0;
    if(j>=ia-1 && j<=ib-1)
      Psi[0][j]=(-c0*(-0.5*amewa*sqrt(amewa*amewa+1)-0.5*sinh(amewa)
          +0.5*x[j]*sqrt(x[j]*x[j]+1)+0.5*sinh(x[j])))*f;
    if(j>ib-1)
      Psi[0][j]=Psi[0][ib-1];
    Omega[0][j]=(7*Psi[0][j]-8*Psi[1][j]+Psi[2][j])/(2*dy2)/DMsq[0][j];
    u[0][j]=0;
    v[0][j]=0;
    if(j>ia-1 && j<ib-1)
      v[0][j]=c0*f;
    if(j>=ia-2 && j<=ib)
      Omega[0][j]+=(v[0][j+1]*sqrt(x[j+1]*x[j+1]+1)
        -v[0][j-1]*sqrt(x[j-1]*x[j-1]+1))/(2*dx)/DMsq[0][j];
    Omega[My+1][j]=0;
    Psi[My+1][j]=(x[j]+A)*(y[My+1]-1);
    u[My+1][j]=(x[j]+A)/DM[My+1][j];
    v[My+1][j]=-(y[My+1]-1)/DM[My+1][j];
  }

  // Side BCs
  for(i=1; i<My+1; ++i) {
    if(i<IBL) {
      Omega[i][0]=Omega[i][1];
      Psi[i][0]=Psi[i][1];
      u[i][0]=u[i][1];
      v[i][0]=v[i][1];
      Omega[i][Nx+1]=Omega[i][Nx];
      Psi[i][Nx+1]=Psi[i][Nx];
      u[i][Nx+1]=u[i][Nx];
      v[i][Nx+1]=v[i][Nx];
    }
    else {
      Omega[i][0]=0;
      Psi[i][0]=(x[0]+A)*(y[i]-1);
      u[i][0]=(x[0]+A)/DM[i][0];
      v[i][0]=-(y[i]-1)/DM[i][0];
      Omega[i][Nx+1]=0;
      Psi[i][Nx+1]=(x[Nx+1]+A)*(y[i]-1);
      u[i][Nx+1]=(x[Nx+1]+A)/DM[i][Nx+1];
      v[i][Nx+1]=-(y[i]-1)/DM[i][Nx+1];
    }
  }
  return 0;
}

static void velocityCalc(
    config* c,
    fields* f,
    derived* d)
{
  int i,j;
  field u = f->u;
  field v = f->v;
  field Psi = f->Psi;
  field DM = f->DM;
  const int Nx = c->Nx;
  const int My = c->My;
  const double dx2 = d->dx2;
  const double dy2 = d->dy2;
  for(i=1; i<My+1; ++i)
  for(j=1; j<Nx+1; ++j)
  {
    u[i][j]= (Psi[i+1][j]-Psi[i-1][j])/dy2/DM[i][j];
    v[i][j]=-(Psi[i][j+1]-Psi[i][j-1])/dx2/DM[i][j];
  }
}

static double onePsiCalc(
    config* c,
    fields* f,
    field Psi,
    field Psi0,
    derived* d)
{
  int Nx = c->Nx;
  int My = c->My;
  field DMsq = f->DMsq;
  field Omega = f->Omega;
  double KappaA = d->KappaA;
  double Kappasq = d->Kappasq;
  double dxsq = d->dxsq;
  double Psi_tol = 0;
  int i,j;
  for(i=1; i<My+1; ++i)
    for(j=1; j<Nx+1; ++j) {
      Psi[i][j]=
        KappaA*
//todo: build a matrix of these:
         (dxsq*Omega[i][j]*DMsq[i][j]+
          Psi0[i][j+1]+
          Psi0[i][j-1]+
          Kappasq*(Psi0[i+1][j]+
                   Psi0[i-1][j]));
      if(fabs(Psi[i][j]-Psi0[i][j])>Psi_tol)
        Psi_tol=fabs(Psi[i][j]-Psi0[i][j]);
    }
  copy_boundaries(Nx,My,Psi0,Psi);
  return Psi_tol;
}

static void psiCalc(
    config* c,
    fields* f,
    derived* d,
    vol* vl)
{
  int Psi_k;
  double Tol = c->Tol;
  double Psi_tol = onePsiCalc(c,f,f->Psi,f->Psi0,d);
  Psi_k = 1;
  while (Psi_tol > Tol) {
     swapFields(&(f->Psi),&(f->Psi0i));
     Psi_tol = onePsiCalc(c,f,f->Psi,f->Psi0i,d);
     ++(Psi_k);
  }
  vl->Psi_k=Psi_k;
  vl->Psi_tol=Psi_tol;
}

static void omegaCalc(
    config* c,
    fields* f,
    derived* d)
{
  int i,j;
  field Omega=f->Omega;
  field Omega0=f->Omega0;
  field DMsq=f->DMsq;
  field u=f->u;
  field v=f->v;
  field DM=f->DM;
  double alpha=d->alpha;
  double alphaX=d->alphaX;
  double alphaY=d->alphaY;
  double Cxd2=d->Cxd2;
  double Cyd2=d->Cyd2;
  int My=c->My;
  int Nx=c->Nx;

  for(i=1;i<My+1;++i)
  for(j=1;j<Nx+1;++j)
  {
    Omega[i][j]=Omega0[i][j]*(1-alpha/DMsq[i][j])+
                Omega0[i][j+1]*(-Cxd2*u[i][j+1]*DM[i][j+1]+alphaX)/DMsq[i][j]+
                Omega0[i][j-1]*( Cxd2*u[i][j-1]*DM[i][j-1]+alphaX)/DMsq[i][j]+
                Omega0[i+1][j]*(-Cyd2*v[i+1][j]*DM[i+1][j]+alphaY)/DMsq[i][j]+
                Omega0[i-1][j]*( Cyd2*v[i-1][j]*DM[i-1][j]+alphaY)/DMsq[i][j];
  }
  copy_boundaries(Nx,My,Omega0,Omega);
}

void initDerived(config* c, space* s, derived* d)
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

void initVolatile(
    config* c,
    int k,
    vol* v)
{
  v->step = k;
  v->time = k*(c->dt);
}

void initFields(config* c, space* s, derived* d, fields* f)
{
  int i,j;
  int Nx = c->Nx;
  int My = c->My;
  field u = f->u;
  field v = f->v;
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

static double fieldMaxDifference(int Nx,int My,field A, field B)
{
  int i,j;
  double max_diff=0;
  for(i=0;i<My+2;++i)
  for(j=0;j<Nx+2;++j)
  {
    if(fabs(A[i][j]-B[i][j])>max_diff)
      max_diff=fabs(A[i][j]-B[i][j]);
  }
  return max_diff;
}

static void maxDiffCalc(
    config* c,
    fields* f,
    vol* vl)
{
  int Nx = c->Nx;
  int My = c->My;
  vl->Omega_tol = fieldMaxDifference(Nx,My,f->Omega,f->Omega0);
  vl->Psi_tol = fieldMaxDifference(Nx,My,f->Psi,f->Psi0);
}

void oneTimeStep(
    config* c,
    space* s,
    fields* f,
    derived* d,
    vol* v)
{
  omegaCalc(c,f,d);
  psiCalc(c,f,d,v);
  setBoundaries(c,s,v->time,f,d);
  velocityCalc(c,f,d);
  maxDiffCalc(c,f,v);
}

void calculate(
    config* c,
    space* s,
    fields* f,
    derived* d,
    vol* v)
{
  for (; v->step < c->Ot; ++(v->step))
  {
    swapFields(&(f->Psi),&(f->Psi0));
    swapFields(&(f->Omega),&(f->Omega0));
    oneTimeStep(c,s,f,d,v);
    v->time += c->dt;
  }
}
