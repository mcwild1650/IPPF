#ifndef FIELD_H
#define FIELD_H

#include "config.h"

typedef double** field;

field make_field(int x, int y);
void swap_fields(field* a, field* b);
void free_field(field f);

typedef struct {
  field Omega;
  field Psi;
  field u;
  field v;
} fields;

void make_fields(fields* f, config* c);
void free_fields(fields* f, config* c);

#endif
