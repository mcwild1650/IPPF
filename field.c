#include "field.h"
#include "tools.h"

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
  int x2 = x;
  int y2 = y;
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
