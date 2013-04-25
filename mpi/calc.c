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
#include <mpi.h>

int pxpy2rank(int px, int py, grid* g)
{
  int n= g->n, m= g->m;
  if(px>n-1 || py >m-1 || px<0 || py<0) return -1;
  else return px+py*n;
}

void mpi_copy_boundaries(int Nx, int My, field a, field b, grid* g)
{
  int px,py;  //position of the cell
  int i;
  int dest_rank;
  double *send_column1,*send_column_1,*rec_column1,*rec_column_1;
  MPI_Request request[8];
  px= g->px;
  py= g->py;


  ALLOCATE(send_column1,My);
  ALLOCATE(send_column_1,My);
  ALLOCATE(rec_column1,My);
  ALLOCATE(rec_column_1,My);

  dest_rank=pxpy2rank(px+1,py,g);
  if(dest_rank!=-1)
  { //receive right boundary from right processor or previous
    for(i=1;i<=My;i++)
    { //send right column to right processor 
      send_column1[i-1]=b[i][Nx];
    }
    MPI_Isend(send_column1,My,MPI_DOUBLE,dest_rank,0,MPI_COMM_WORLD,&request[0]);
    MPI_Irecv(rec_column1 ,My,MPI_DOUBLE,dest_rank,0,MPI_COMM_WORLD,&request[1]);
  }
  else
  {
    for(i=0;i<My+2;i++) //boundaries are copied entirely to maintain consistent corners
      b[i][Nx+1]=a[i][Nx+1];
    request[0]=request[1]=MPI_REQUEST_NULL;
  }

  dest_rank=pxpy2rank(px-1,py,g);
  if(dest_rank!=-1)
  { //receive left boundary from left processor or previous
    for(i=1;i<=My;i++)
    { //send left column to left processor
      send_column_1[i-1]=b[i][1];
    }
    MPI_Isend(send_column_1,My,MPI_DOUBLE,dest_rank,1,MPI_COMM_WORLD,&request[2]);
    MPI_Irecv(rec_column_1 ,My,MPI_DOUBLE,dest_rank,1,MPI_COMM_WORLD,&request[3]);
  }
  else
  {
    for(i=0;i<My+2;i++)
      b[i][0]=a[i][0];
    request[2]=request[3]=MPI_REQUEST_NULL;
  }

  dest_rank=pxpy2rank(px,py-1,g);
  if(dest_rank!=-1)
  {
    //send bottom row to below processor
    MPI_Isend(&b[1][1],Nx,MPI_DOUBLE,dest_rank,2,MPI_COMM_WORLD,&request[6]);
    //receive bottom boundary from below processor or previous
    MPI_Irecv(&b[0][1],Nx,MPI_DOUBLE,dest_rank,2,MPI_COMM_WORLD,&request[5]);
  }
  else
  {
    printf("Copying field0[0][0]=%lf to field[0][0]\n",a[0][0]);
    for(i=0;i<Nx+2;i++)
      b[0][i]=a[0][i];
    printf("field[0][0]=%lf\n",b[0][0]);
    request[4]=request[5]=MPI_REQUEST_NULL;
  }

  dest_rank=pxpy2rank(px,py+1,g);
  if(dest_rank!=-1)
  {
    //send top row to above processor
    MPI_Isend(&b[My][1]  ,Nx,MPI_DOUBLE,dest_rank,3,MPI_COMM_WORLD,&request[4]);
    //receive top boundary from above processor or previous
    MPI_Irecv(&b[My+1][1],Nx,MPI_DOUBLE,dest_rank,3,MPI_COMM_WORLD,&request[7]);
  }
  else
  {
    for(i=0;i<Nx+2;i++)
      b[My+1][i]=a[My+1][i];
    request[6]=request[7]=MPI_REQUEST_NULL;
  }

  MPI_Waitall(8,request,MPI_STATUSES_IGNORE);

  //copy the columns where they belong
  for(i=1;i<=My;i++)
    b[i][Nx+1] = rec_column1[i-1];
  for(i=1;i<=My;i++)
    b[i][0] = rec_column_1[i-1];

  deallocate(send_column_1);
  deallocate(rec_column_1);
  deallocate(send_column1);
  deallocate(rec_column1);
}

static int setBoundaries(
    config* c,
    space* s,
    double t,
    fields* fld,
    derived* d,
    grid* g)
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
  int px=g->px;
  int py=g->py;
  int x=g->x;
  int y=g->y;
  int m=g->m;
  int n=g->n;
  int lastx=g->lastx;
  int lasty=g->lasty;
  int start_x=g->x*g->dx;
  int len_x=(px!=n-1)?g->dx:g->lastx;
  int len_y=(py!=m-1)?g->dy:g->lasty;
    

  const double amewa=(ia-((x+1)/2+1))*dx;
  double jet_val;
  //const double amewa=(ia-((Nx+1)/2+1))*dx;
  const double f=sin(2*PI*freq*t);
  double dy2 = d->dy2;
  // Upper and Lower BCs
  if(px==0)
  {
    for(j=0; j<len_x+2; ++j) {
      Psi[0][j]=0;
      if(start_x+j>=ia-1 && start_x+j<=ib-1)
        Psi[0][j]=(-c0*(-0.5*amewa*sqrt(amewa*amewa+1)-0.5*sinh(amewa)
            +0.5*x[j]*sqrt(x[j]*x[j]+1)+0.5*sinh(x[j])))*f;
      if(start_x+j>ib-1)
        Psi[0][j]=Psi[0][ib-1];
      Omega[0][j]=(7*Psi[0][j]-8*Psi[1][j]+Psi[2][j])/(2*dy2)/DMsq[0][j];
      u[0][j]=0;
      v[0][j]=0;
      if(j>ia-1 && j<ib-1)
        v[0][j]=c0*f;
      if(j>=ia-2 && j<=ib)
        Omega[0][j]+=(v[0][j+1]*sqrt(x[j+1]*x[j+1]+1)
          -v[0][j-1]*sqrt(x[j-1]*x[j-1]+1))/(2*dx)/DMsq[0][j];
    }
  }
  if(px==n)
  {
    for(j=0; j<len_x+2; ++j) {
      Omega[My+1][j]=0;
      Psi[My+1][j]=(x[j]+A)*(y[My+1]-1);
      u[My+1][j]=(x[j]+A)/DM[My+1][j];
      v[My+1][j]=-(y[My+1]-1)/DM[My+1][j];
    }
  }
  for(j=0; j<len_x+2; ++j) {
  //for(j=0; j<Nx+2; ++j) {
    Psi[0][j]=0;
    if(start_x+j>=ia-1 && start_x+j<=ib-1)
      Psi[0][j]=(-c0*(-0.5*amewa*sqrt(amewa*amewa+1)-0.5*sinh(amewa)
          +0.5*x[j]*sqrt(x[j]*x[j]+1)+0.5*sinh(x[j])))*f;
    if(start_x+j>ib-1)
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
    derived* d,
    grid* g)
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
  mpi_copy_boundaries(Nx,My,Psi0,Psi,g);
  return Psi_tol;
}

static void psiCalc(
    config* c,
    fields* f,
    derived* d,
    vol* vl,
    grid* g)
{
  int Psi_k;
  double Tol = c->Tol;
  double Psi_tol = onePsiCalc(c,f,f->Psi,f->Psi0,d,g);
  Psi_k = 1;
  while (Psi_tol > Tol) {
     swapFields(&(f->Psi),&(f->Psi0i));
     Psi_tol = onePsiCalc(c,f,f->Psi,f->Psi0i,d,g);
     ++(Psi_k);
  }
  vl->Psi_k=Psi_k;
  vl->Psi_tol=Psi_tol;
}

static void omegaCalc(
    config* c,
    fields* f,
    derived* d,
    grid* g)
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
  fprintf(stderr,"copying Omega\n");
  mpi_copy_boundaries(Nx,My,Omega0,Omega,g);
  fprintf(stderr,"done copying Omega\n");
  fprintf(stderr,"Omega[0][0]=%lf\n",Omega[0][0]);
}


void initDerived(config* c, space* s, derived* d)
{
  double dx = s->dx;
  double dy = s->dy;
  double dxx = SQUARE(dx);
  double dyy = SQUARE(dy);
  d->dxsq = dxx;
  d->dysq = dyy;
  double dx2 = 2*dx;
  double dy2 = 2*dy;
  d->dx2 = dx2;
  d->dy2 = dy2;
  double dt = c->dt;
  double Kappa2 = SQUARE(dx/dy);
  d->Kappasq = Kappa2;
  d->KappaA = 1/(2*(1.0+Kappa2));
  double Cx = dt/dx;
  double Cy = dt/dy;
  d->Cx = Cx;
  d->Cy = Cy;
  double C = MAX(Cx,Cy);
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
      DM2[i][j]=SQUARE(x[j])+SQUARE(y[i]);
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
    grid* g,
    fields* f,
    derived* d,
    vol* v)
{
  omegaCalc(c,f,d,g);
  psiCalc(c,f,d,v,g);
  printf("before BCs Omega[0][0]=%lf\n",f->Omega[0][0]);
  setBoundaries(c,s,v->time,f,d);
  printf("after BCs Omega[0][0]=%lf\n",f->Omega[0][0]);
  velocityCalc(c,f,d);
  maxDiffCalc(c,f,v);
}

void calculate(
    config* c,
    space* s,
    grid* g,
    fields* f,
    derived* d,
    vol* v)
{
  for (; v->step < c->Ot; ++(v->step))
  {
    swapFields(&(f->Psi),&(f->Psi0));
    swapFields(&(f->Omega),&(f->Omega0));
    oneTimeStep(c,s,g,f,d,v);
    v->time += c->dt;
  }
}
