#include "field.h"

field make_field(int x, int y)
{
  field f;
  ALLOCATE(f,y);
  ALLOCATE(f[0],y*x);
  for (int i=1; i < y; ++i)
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
