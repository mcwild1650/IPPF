#ifndef SPACE_H
#define SPACE_H

#include "config.h"

typedef struct{
  double dx,dy;
  double* x;
  double* y;
} space;

void initSpace(config* c, space* s);
void freeSpace(space* s);

#endif
