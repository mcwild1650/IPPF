#ifndef RESTART_H
#define RESTART_H

#include "field.h"

void read_restart(const char* filename, config* c, fields* f);
void write_restart(const char* filename, config* c, fields* f);

#endif
