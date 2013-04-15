//#include "calc.h"
#include "tools.h"
#include "space.h"

void make_space(space* s, int x, double dx, int y,double dy)
{
  int i;
  ALLOCATE(s->x,x+y);
  s->y=s->x+x;
  s->x[0]=0;
  s->y[0]=0;
  for(i=1;i<x;++i)
    s->x[i]=s->x[i-1]+dx;
  for(i=1;i<y;++i)
    s->y[i]=s->y[i-1]+dy;
}
  
void free_space(space* s)
{
  deallocate(s->x);
}
