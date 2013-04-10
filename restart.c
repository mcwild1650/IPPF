#include "restart.h"

void read_restart(const char* filename, config* c, fields* f)
{
  int Nx = c->Nx;
  int My = c->My;
  make_fields(f,Nx+2,My+2);
  //...
}

void write_restart(const char* filename, config* c, fields* f)
{
}
