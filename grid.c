#include "grid.h"
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

void computeGrid(int x, int y, int p, grid* g)
{
  int good_px = 0;
  int good_py = 0;
  int min_comm = INT_MAX;
  int min_comp = INT_MAX;
  int maxx = MIN(p,x);
  int maxy = MIN(p,y);
  for (int px=1; px <= maxx; ++px)
  for (int py=1; py <= maxy; ++py)
    if (px*py == p)
    {
      int dx = ceildiv(x,px);
      int dy = ceildiv(y,py);
      int comm=0;
      int comp=dx*dy;
      for (int c=dx; c < x; c += dx)
        comm += y;
      for (int r=dy; r < y; r += dy)
        comm += x;
      printf("try %d %d: %d, %d\n",px,py,comm,comp);
      if ((comm <= min_comm)&&
          (comp <= min_comp))
      {
        min_comm = comm;
        min_comp = comp;
        good_px = px;
        good_py = py;
      }
    }
  if ((!good_px)&&(!good_py))
  {
    fprintf(stderr,"fail\n");
    return;
  }
  printf("comm %d\n",min_comm);
  printf("comp %d\n",min_comp);
  g->px = good_px;
  g->py = good_py;
  g->dx = ceildiv(x,good_px);
  g->dy = ceildiv(x,good_py);
  g->lastx = x - (g->px-1) * g->dx;
  g->lasty = y - (g->py-1) * g->dy;
}
