#include "field.h"
#include "tools.h"

field make_field(int x, int y)
{
  int i;
  field f;
  ALLOCATE(f,y);
  ALLOCATE(f[0],y*x);
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
  int x2 = x;
  int y2 = y;
  f->Omega = make_field(x2,y2);
  f->Psi = make_field(x2,y2);
  f->u = make_field(x2,y2);
  f->v = make_field(x2,y2);
}

void free_fields(fields* f)
{
  free_field(f->Omega);
  free_field(f->Psi);
  free_field(f->u);
  free_field(f->v);
}
