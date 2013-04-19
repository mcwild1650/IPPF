#ifndef CALC_H
#define CALC_H

#include "config.h"
#include "space.h"
#include "field.h"

#define PI 3.14159265358979

//derived constants
typedef struct {
  double dxsq,dysq;
  double dx2,dy2;
  double Kappasq,KappaA;
  double C;
  double Cx,Cy;
  double Cxd2,Cyd2;
  double alpha,alphaX,alphaY;
} derived;

//volatile variables, changing from iteration to iteration
typedef struct {
  int step;
  double time;
  int Psi_k;
  double Omega_tol;
  double Psi_tol;
} vol;

void init_derived(config* c, space* s, derived* d);
void init_fields(config* c, space* s, derived* d, fields* f);
void one_time_step(
    config* c,
    space* s,
    fields* f,
    derived* d,
    vol* v);
void calculate(
    config* c,
    space* s,
    fields* f,
    derived* d,
    vol* v);

#endif
