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
  int Nx, My;
  double Re;
  double A;
  int Ot;
  int report;
  double dt;
  double Tol;
  int IBL;
  jet_config jet;
  char filename[FILENAME_LENGTH];
  int psave;
  double Xmin,Xmax;
  double Ymin,Ymax;
} config;

void ask_config(config* c);
void write_config(const char* filename, const config* c);
void read_config(const char* filename, config* c);
void broadcast_config(config* c);

#endif
