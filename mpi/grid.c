#include "grid.h"
#include "tools.h"
#include <limits.h>
#include <stdio.h>

#define MIN(a,b) (((b)<(a))?(b):(a))

static int ceildiv(int a, int b)
{
  if (a%b)
    return a/b+1;
  else
    return a/b;
}

void partition(config* c, grid* g)
{
  int x = c->Nx;
  int y = c->My;
  int good_m = 0;
  int good_n = 0;
  int min_comm = INT_MAX;
  int min_comp = INT_MAX;
  int p = parallelSize();
  int maxx = MIN(p,x);
  int maxy = MIN(p,y);
  int m,n;
  int rank = parallelRank();
  for (m=1; m <= maxy; ++m)
  for (n=1; n <= maxx; ++n)
    if (m*n == p)
    {
      int dx = ceildiv(x,n);
      if ((n*dx-x) > dx)
        continue;
      int dy = ceildiv(y,m);
      if ((m*dy-y) > dy)
        continue;
      int comm=(m-1)*x + (n-1)*y;
      int comp=dx*dy;
      if ((comm <= min_comm)&&
          (comp <= min_comp))
      {
        min_comm = comm;
        min_comp = comp;
        good_n = n;
        good_m = m;
      }
    }
  if ((!good_m)||(!good_n))
    die("failed to partition\n");
  g->x = x;
  g->y = y;
  g->n = good_n;
  g->m = good_m;
  g->dx = ceildiv(x,good_n);
  g->dy = ceildiv(y,good_m);
  printf("rank %d dx=%d dy=%d\n",rank,g->dx,g->dy);
  g->lastx = x -((g->n-1) * g->dx);
  g->lasty = y -((g->m-1) * g->dy);
  g->py = rank/(g->n);
  g->px = rank%(g->n);
  if (g->px == g->n-1)
    c->Nx = g->lastx;
  else
    c->Nx = g->dx;
  if (g->py == g->m-1)
    c->My = g->lasty;
  else
    c->My = g->dy;
  printf("rank %d Nx=%d My=%d\n",rank,c->Nx,c->My);
  g->x0 = (g->px == 0)?(0):(1);
  g->x1 = (g->px == g->n-1)?(c->Nx+2):(c->Nx+1);
  g->y0 = (g->py == 0)?(0):(1);
  g->y1 = (g->py == g->m-1)?(c->My+2):(c->My+1);
  printf("rank %d is [%d,%d]X[%d,%d]\n",rank,g->x0,g->x1,g->y0,g->y1);
}
