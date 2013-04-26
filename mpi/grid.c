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
  for (m=1; m <= maxy; ++m)
  for (n=1; n <= maxx; ++n)
    if (m*n == p)
    {
      int dx = ceildiv(x,n);
      int dy = ceildiv(y,m);
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
  g->lastx = x -((g->n-1) * g->dx);
  g->lasty = y -((g->m-1) * g->dy);
  int rank = parallelRank();
  g->py = rank/(g->n);
  g->px = rank%(g->n);
  g->len_x=(g->px!=n-1)?g->dx:g->lastx;
  g->len_y=(g->py!=m-1)?g->dy:g->lasty;
  if (g->py == g->m-1)
    c->Nx = g->lastx;
  else
    c->Nx = g->dx;
  if (g->px == g->n-1)
    c->My = g->lasty;
  else
    c->My = g->dy;
}
