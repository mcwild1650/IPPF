#include "tools.h"
#include "space.h"

static void linspace(const double min, const double max, const int N,
        double* x, double* dx,int x_start,int x_end) {
  int i;
  *dx=(max-min)/(N-1);
  x[0]=min;
  //x[N-1]=max;
  for(i=1; i<=N-1; ++i)
  {
  	if(i>=x_start && i <= x_end)x[i-x_start]=min+(*dx)*i;
  }
    
}

void initSpace(config* c, space* s)
{
  int Xmin = c->Xmin;
  int Xmax = c->Xmax;
  int Ymin = c->Ymin;
  int Ymax = c->Ymax;
  int Nx = c->tot_Nx; 
  int My = c->tot_My; 
  double* x;
  double* y;
  ALLOCATE(x,Nx+My+4);
  y=x+Nx+2;
  double dx,dy;
  linspace(Xmin,Xmax,Nx+2,x,&dx,c-> mpi_para.my_x_start,c-> mpi_para.my_x_end);
  linspace(Ymin,Ymax,My+2,y,&dy,c-> mpi_para.my_y_start,c-> mpi_para.my_y_end);
  s->x = x;
  s->y = y;
  s->dx = dx;
  s->dy = dy;
}

void freeSpace(space* s)
{
  deallocate(s->x);
}
