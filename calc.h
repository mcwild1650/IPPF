#ifndef CALC_H
#define CALC_H

#include "space.h"
#include "field.h"

//a series of macros to make data access a little cleaner
//assumes that the state is named s

#define PI 3.14159265358979
#define FLD(x) (s->f->x)
#define DRV(x) (s->dr->x)
#define JET(x) (s->c->jet->x)
#define CFG(x) (s->c->x)
#define SPC(x) (s->sp->x)
#define VOL(x) (s->x->v)

//derived constants
typedef struct {
  int step;
  double time;
  double Kappasq,KappaA;
  double C,Cx,Cy,Cxd2,Cyd2;
  double dx,dy,dx2,dxsq,dysq,dy2;
  double alpha,alphaX,alphaY;
} derived;

//volatile variables, changing from iteration to iteration
typdef struct {
  int Psi_k,tic1,tic2,ct;
  double Omega_tol;
  double Psi_tol;
} vol;

//struct containing heirarchy of data for calculation
typedef struct {
  fields* f;
  space* sp;
  config* c;
  vol* v;
} state;

int initialize_all(state* state);
int run_time_step(state* state);

#endif
