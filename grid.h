#ifndef GRID_H
#define GRID_H

typedef struct {
  int px;
  int py;
  int dx;
  int dy;
  int lastx;
  int lasty;
} grid;

void compute_grid(int x, int y, int p, grid* g);

#endif
