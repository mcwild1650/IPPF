#ifndef FIELD_H
#define FIELD_H

#include "config.h"
#include "space.h"

typedef double** field;

field make_field(int x, int y, int total);
void swap_fields(field* a, field* b);
void free_field(field f);

typedef struct {
  int x;
  int y;
  int size;
  field Omega;
  field Omega0;
  field Psi;
  field Psi0;
  field Psi0i;
  field u;
  field v;
  field DM;
  field DMsq;
} fields;

void make_fields(fields* f, int x, int y);
void free_fields(fields* f);

#endif
