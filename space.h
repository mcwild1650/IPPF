#ifndef SPACE_H
#define SPACE_H

#include "config.h"

typedef struct{
  double dx,dy;
  double* x;
  double* y;
} space;

void init_space(config* c, space* s);
void free_space(space* s);

#endif
