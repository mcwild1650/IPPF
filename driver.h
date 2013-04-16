#ifndef DRIVER_H
#define DRIVER_H

#include "field.h"
#include "space.h"
#include "config.h"
#include "calc.h"

//struct containing heirarchy of data for calculation
typedef struct {
  fields* f;
  space* sp;
  config* c;
  derived* dr;
  vol* v;
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

void make_state(state* s);
void free_state(state* s);
int get_args(int argc, char** argv, args* a);

#endif
