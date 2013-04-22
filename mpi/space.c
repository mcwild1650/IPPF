#include "tools.h"
#include "space.h"

static void linspace(const double min, const double max, const int N,
        double* x, double* dx) {
  int i;
  *dx=(max-min)/(N-1);
  x[0]=min;
  x[N-1]=max;
  for(i=1; i<N-1; ++i)
    x[i]=min+(*dx)*i;
}

void initSpace(config* c, space* s)
{
  int Xmin = c->Xmin;
  int Xmax = c->Xmax;
  int Ymin = c->Ymin;
  int Ymax = c->Ymax;
  int Nx = c->Nx; 
  int My = c->My; 
  double* x;
  double* y;
  ALLOCATE(x,Nx+My+4);
  y=x+Nx+2;
  double dx,dy;
  linspace(Xmin,Xmax,Nx+2,x,&dx);
  linspace(Ymin,Ymax,My+2,y,&dy);
  s->x = x;
  s->y = y;
  s->dx = dx;
  s->dy = dy;
}

void freeSpace(space* s)
{
  deallocate(s->x);
}
