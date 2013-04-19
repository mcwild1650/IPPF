#ifndef RESTART_H
#define RESTART_H

#include "field.h"
#include "calc.h"
#include "space.h"

void read_restart(const char* filename, config* c, fields* f);
void write_old_restart(
    config* c,
    space* s,
    fields* f,
    derived* d,
    vol* vl);

#endif
