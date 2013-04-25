#ifndef RESTART_H
#define RESTART_H

#include "field.h"
#include "calc.h"
#include "space.h"

void readRestart(const char* filename, config* c, fields* f);
void writeOldRestart(
    config* c,
    space* s,
    fields* f,
    derived* d,
    vol* vl);
void writeOutput(
    config* c,
    space* s,
    fields* f,
    vol* vl,
    grid* g);

#endif
