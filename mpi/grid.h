#ifndef GRID_H
#define GRID_H

#include "config.h"

typedef struct {
  int x,y;
  int n,m;
  int dx,dy;
  int lastx,lasty;
  int px,py;
  int len_x,len_y;
} grid;

void partition(config* c, grid* g);

#endif
