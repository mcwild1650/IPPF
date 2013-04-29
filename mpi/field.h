#ifndef FIELD_H
#define FIELD_H

#include "config.h"
#include "space.h"

typedef double** field;

field makeField(int x, int y);
void swapFields(field* a, field* b);
void freeField(field f);

typedef struct {
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

void makeFields(fields* f, int x, int y);
void freeFields(fields* f);
void testFieldSanity(field f, int My,int Nx);

#endif
