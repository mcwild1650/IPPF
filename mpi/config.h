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
 int myrank;
 int mpi_size;
 int my_x_start,my_x_end;
 int my_y_start,my_y_end;
 int m,n;  // totoal partition number in x and y direction
 int px,py;   //position of this cell 
}  mpipara;

typedef struct {
  int Nx, My;    // Nx,My in a rank
  int tot_Nx,tot_My;   /////// the same with the original Nx, My
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
  mpipara mpi_para;
} config;

void askConfig(config* c);
void writeConfig(const char* filename, const config* c);
void readConfig(const char* filename, config* c);
void broadcastConfig(config* c);

#endif
