#ifndef SPACE_H
#define SPACE_H

#include "config.h"
#include "grid.h"

typedef struct{
  double dx,dy;
  double* x;
  double* y;
} space;

void initSpace(config* c, grid* g, space* s);
void freeSpace(space* s);

#endif
