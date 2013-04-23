#ifndef DRIVER_H
#define DRIVER_H

#include "field.h"
#include "space.h"
#include "config.h"
#include "calc.h"
#include "grid.h"

//struct containing heirarchy of data for calculation
typedef struct {
  config* c;
  grid* g;
  space* sp;
  vol* v;
  derived* dr;
  fields* f;
} state;

typedef enum {
  config_mode,
  run_mode,
  restart_mode
} mode;

typedef struct {
  mode m;
  char configfile[FILENAME_LENGTH];
} args;

void makeState(state* s);
void freeState(state* s);
int getArgs(int argc, char** argv, args* a);

#endif
