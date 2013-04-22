#include "field.h"
#include "tools.h"
#include <math.h>
#include <assert.h>

field makeField(int x, int y)
{
  int i;
  field f;
  ALLOCATE(f,y);
  ALLOCATE(f[0],x*y);
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
  if (f) deallocate(f[0]);
  deallocate(f);
}

void testFieldSanity(field f, int My,int Nx)
{
  int i,j;
  for(i=0;i<My+2;++i)
  for(j=0;j<Nx+2;++j)
  {
    testDoubleSanity(f[i][j]);
  }
}

void makeFields(fields* f, int x, int y)
{
  int x2 = x+2;
  int y2 = y+2;
  f->Omega = makeField(x2,y2);
  f->Omega0 = makeField(x2,y2);
  f->Psi = makeField(x2,y2);
  f->Psi0 = makeField(x2,y2);
  f->Psi0i = makeField(x2,y2);
  f->u = makeField(x2,y2);
  f->v = makeField(x2,y2);
  f->DM = makeField(x2,y2);
  f->DMsq = makeField(x2,y2);
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

