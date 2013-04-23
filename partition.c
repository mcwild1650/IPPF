#include <stdio.h>
#include "grid.h"

int main()
{
  int x,y,p;
  scanf("%d%d%d",&x,&y,&p);
  grid g;
  computeGrid(x,y,p,&g);
  printf("%d %d\n",g.px,g.py);
  return 0;
}
