#ifndef GRID_H
#define GRID_H

#include "config.h"

typedef struct {
  int x,y;
  int x0,y0;
  int x1,y1;
  int n,m;
  int dx,dy;
  int lastx,lasty;
  int px,py;
} grid;

void partition(config* c, grid* g);

#endif
