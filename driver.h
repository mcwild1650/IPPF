#ifndef DRIVER_H
#define DRIVER_H

#include "config.h"

typedef enum {
  config_mode,
  run_mode,
  restart_mode
} mode;

typedef struct {
  mode m;
  char configfile[FILENAME_LENGTH];
} args;

int get_args(int argc, char** argv, args* a);

#endif
