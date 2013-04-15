#ifndef SPACE_H
#define SPACE_H

typedef struct{
  double* x;
  double* y;
} space;

void make_space(space* s, int x, double dx, int y, double dy);
void free_space(space* s);
#endif
