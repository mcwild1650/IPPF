#include "field.h"
#include "tools.h"
#include <math.h>

field make_field(int x, int y, int size)
{
  int i;
  field f;
  ALLOCATE(f,y);
  ALLOCATE(f[0],size);
  for (i=1; i < y; ++i)
    f[i] = f[i-1] + x;
  return f;
}

void swap_fields(field* a, field* b)
{
  field temp = *a;
  *a = *b;
  *b = temp;
}

void free_field(field f)
{
  deallocate(f[0]);
  deallocate(f);
}

void make_fields(fields* f, int x, int y)
{
  int x2 = x+2;
  int y2 = y+2;
  f->x=x;
  f->y=y;
  f->size=x*y;
  f->Omega = make_field(x2,y2,f->size);
  f->Omega0 = make_field(x2,y2,f->size);
  f->Psi = make_field(x2,y2,f->size);
  f->Psi0 = make_field(x2,y2,f->size);
  f->Psi0i = make_field(x2,y2,f->size);
  f->u = make_field(x2,y2,f->size);
  f->v = make_field(x2,y2,f->size);
  f->DM = make_field(x2,y2,f->size);
  f->DMsq = make_field(x2,y2,f->size);
}

void free_fields(fields* f)
{
  free_field(f->Omega);
  free_field(f->Omega0);
  free_field(f->Psi);
  free_field(f->Psi0);
  free_field(f->Psi0i);
  free_field(f->u);
  free_field(f->v);
  free_field(f->DM);
  free_field(f->DMsq);
}

void init_fields(config* c, space* s, fields* f)
{
  int i,j;
  int My = f->y;
  int Nx = f->x;
  field u = f->u;
  field v = f->u;
  field DM = s->f->DM;
  field DM2 = s->f->DMsq;
  field Psi = s->f->Psi;
  field Omega = s->f->Omega;
  double* x = s->sp->x;
  double* y = s->sp->y;
  double A = s->c->A;
  double dyy = s->dr->dysq;
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

