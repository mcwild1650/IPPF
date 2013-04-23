#ifndef GRID_H
#define GRID_H

#include "config.h"

typedef struct {
  int x,y;
  int m,n;
  int dx,dy;
  int lastx,lasty;
  int px,py;
} grid;

void partition(config* c, grid* g);

#endif
