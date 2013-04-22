#include <stdio.h>
#include <limits.h>

int ceildiv(int a, int b)
{
  if (a % b)
    return a/b+1;
  return a/b;
}

int main()
{
  int x,y,p;
  scanf("%d%d%d",&x,&y,&p);
  int good_px;
  int good_py;
  int min_comm = INT_MAX;
  int min_comp = INT_MAX;
  for (int px=1; px <= p; ++px)
  for (int py=1; py <= p; ++py)
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
      if ((comm <= min_comm)&&
          (comp <= min_comp))
      {
        printf("%d %d\n",px,py);
        printf("%d %d\n",dx,dy);
        printf("comm %d\n",comm);
        printf("comp %d\n",comp);
        min_comm = comm;
        min_comp = comp;
        good_px = px;
        good_py = py;
      }
    }
  return 0;
}
