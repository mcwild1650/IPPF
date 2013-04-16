#ifndef CALC_H
#define CALC_H

#include "config.h"
#include "space.h"
#include "field.h"

//a series of macros to make data access a little cleaner
//assumes that the state is named s

#define PI 3.14159265358979
#define FLD(x) (s->f->x)
#define DRV(x) (s->dr->x)
#define JET(x) (s->c->jet.x)
#define CFG(x) (s->c->x)
#define SPC(x) (s->sp->x)
#define VOL(x) (s->v->x)

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
  int Psi_k,tic1,tic2,ct;
  double Omega_tol;
  double Psi_tol;
} vol;

void init_derived(config* c, space* s, derived* d);
void init_fields(config* c, space* s, derived* d, fields* f);

#endif
