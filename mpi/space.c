#include "tools.h"
#include "space.h"

static void linspace(const double min, const double max, const int N,
        int start, int len, double* x, double* dx) {
  int i;
  *dx=(max-min)/(N-1);
  if(start==0)
    x[0]=min;
  else
    x[0]=min+(*dx)*start;
  if(start+len==N)
    x[N-1]=max;
  else
    x[N-1]=min+(*dx)*(start+len-1);
  for(i=1; i<len-1; ++i)
    x[i]=x[i-1]+*dx;
}

void initSpace(config* c, grid* g, space* s)
{
  int Xmin = c->Xmin;
  int Xmax = c->Xmax;
  int Ymin = c->Ymin;
  int Ymax = c->Ymax;
  int Nx=c->Nx;
  int My=c->My;
  int start_x=g->px*g->dx;
  int start_y=g->py*g->dy;
  double* x;
  double* y;
  ALLOCATE(x,Nx+My+4);
  y=x+Nx+2;
  double dx,dy;
  linspace(Xmin,Xmax,g->x+2,start_x,Nx,x,&dx);
  linspace(Ymin,Ymax,g->y+2,start_y,My,y,&dy);
  s->x = x;
  s->y = y;
  s->dx = dx;
  s->dy = dy;
}

void freeSpace(space* s)
{
  deallocate(s->x);
}
