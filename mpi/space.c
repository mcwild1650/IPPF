#include "tools.h"
#include "space.h"
#include <stdio.h>

static void linspace(const double min, const double max,
        int px, int dxdp, int tx, int Nx, double* x, double* dx)
{
  int start = dxdp*px;
  int i;
  *dx = (max-min)/(Nx+1);
  for (i=0; i < Nx+2; ++i)
    x[i] = (i+start)*(*dx) + min;
}

void initSpace(config* c, grid* g, space* s)
{
  int Xmin = c->Xmin;
  int Xmax = c->Xmax;
  int Ymin = c->Ymin;
  int Ymax = c->Ymax;
  int Nx=c->Nx;
  int My=c->My;
  double* x;
  double* y;
  ALLOCATE(x,Nx+My+4);
  y=x+Nx+2;
  double dx,dy;
  linspace(Xmin,Xmax,g->px,g->dx,g->x,Nx,x,&dx);
  linspace(Ymin,Ymax,g->py,g->dy,g->y,My,y,&dy);
  s->x = x;
  s->y = y;
  s->dx = dx;
  s->dy = dy;
}

void freeSpace(space* s)
{
  deallocate(s->x);
}
