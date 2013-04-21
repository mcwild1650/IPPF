#include "field.h"
#include "tools.h"
#include <math.h>
#include <assert.h>

field makeField(int x, int y, int size)
{
  int i;
  field f;
  ALLOCATE(f,y);
  ALLOCATE(f[0],size);
  for (i=1; i < y; ++i)
    f[i] = f[i-1] + x;
  return f;
}

void swapFields(field* a, field* b)
{
  field temp = *a;
  *a = *b;
  *b = temp;
}

void freeField(field f)
{
  deallocate(f[0]);
  deallocate(f);
}

void testFieldSanity(field f, int My,int Nx)
{
  int i,j;
  for(i=0;i<My+2;++i)
  for(j=0;j<Nx+2;++j)
  {
     assert((f[i][j]>0 || f[i][j]<0 || f[i][j]==0)&&
    f[i][j]!=INFINITY && f[i][j]!=-INFINITY);
  }
}

void makeFields(fields* f, int x, int y)
{
  int x2 = x+2;
  int y2 = y+2;
  f->x=x;
  f->y=y;
  f->size=x2*y2;
  f->Omega = makeField(x2,y2,f->size);
  f->Omega0 = makeField(x2,y2,f->size);
  f->Psi = makeField(x2,y2,f->size);
  f->Psi0 = makeField(x2,y2,f->size);
  f->Psi0i = makeField(x2,y2,f->size);
  f->u = makeField(x2,y2,f->size);
  f->v = makeField(x2,y2,f->size);
  f->DM = makeField(x2,y2,f->size);
  f->DMsq = makeField(x2,y2,f->size);
}

void freeFields(fields* f)
{
  freeField(f->Omega);
  freeField(f->Omega0);
  freeField(f->Psi);
  freeField(f->Psi0);
  freeField(f->Psi0i);
  freeField(f->u);
  freeField(f->v);
  freeField(f->DM);
  freeField(f->DMsq);
}

