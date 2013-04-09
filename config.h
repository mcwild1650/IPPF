#ifndef CONFIG_H
#define CONFIG_H

typedef struct {
  int exists;
  int ia;
  int ib;
  double c0;
  double freq;
} jet_config;

#define FILENAME_LENGTH 64

typedef struct {
  int Nx, Ny;
  double Re;
  double A;
  int Ot;
  int report;
  double dt;
  double Tol;
  int IBL;
  jet_config jet;
  const char filename[FILENAME_LENGTH];
  int psave;
} config;

void ask_config(config* c);
void read_config(const char* filename, config* c);
void broadcast_config(config* c);

#endif
