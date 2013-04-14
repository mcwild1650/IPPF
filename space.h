#ifndef SPACE_H
#define SPACE_H

typedef struct{
  double* x;
  double* y;
} space;

make_space(space* s, int x, int y);
free_space(space* s);
#endif
