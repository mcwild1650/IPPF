#include "field.h"
#include "tools.h"
#include <math.h>

field make_field(int x, int y)
{
  int i;
  field f;
  ALLOCATE(f,y);
  ALLOCATE(f[0],x*y);
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
  if (f) deallocate(f[0]);
  deallocate(f);
}

void make_fields(fields* f, int x, int y)
{
  int x2 = x+2;
  int y2 = y+2;
  f->Omega = make_field(x2,y2);
  f->Omega0 = make_field(x2,y2);
  f->Psi = make_field(x2,y2);
  f->Psi0 = make_field(x2,y2);
  f->Psi0i = make_field(x2,y2);
  f->u = make_field(x2,y2);
  f->v = make_field(x2,y2);
  f->DM = make_field(x2,y2);
  f->DMsq = make_field(x2,y2);
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

