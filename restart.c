#include "restart.h"
#include "tools.h"
#include <stdio.h>

// Program output to file
static int writeFile(
    double** Omega,
    double** Psi,
    double** u,
    double** v,
    int* OutIN,
    double* OutDP,
    const int Nx,
    const int My,
    char* outfile)
{
  int i,j;

  FILE* f_out=fopen(outfile,"w");

  if(f_out==NULL)
    return -1;

  for(i=0; i<My+2; ++i) {
    for(j=0; j<Nx+2; ++j)
      fprintf(f_out,"%30.13le",Omega[i][j]);
    fprintf(f_out,"\n");
  }

  for(i=0; i<My+2; ++i) {
    for(j=0; j<Nx+2; ++j)
      fprintf(f_out,"%30.13le",Psi[i][j]);
    fprintf(f_out,"\n");
  }

  for(i=0; i<My+2; ++i) {
    for(j=0; j<Nx+2; ++j)
      fprintf(f_out,"%30.13le",u[i][j]);
    fprintf(f_out,"\n");
  }

  for(i=0; i<My+2; ++i) {
    for(j=0; j<Nx+2; ++j)
      fprintf(f_out,"%30.13le",v[i][j]);
    fprintf(f_out,"\n");
  }

  for(j=0; j<Nx+2; ++j)
    fprintf(f_out,"%30d",OutIN[j]);
  fprintf(f_out,"\n");

  for(j=0; j<Nx+2; ++j)
    fprintf(f_out,"%30.13le",OutDP[j]);

  fclose(f_out);
  return 0;
}


void writeOldRestart(
    config* c,
    space* s,
    fields* f,
    derived* d,
    vol* vl)
{
  int* OutIN;
  double* OutDP;
  int j;
  ALLOCATE(OutIN,c->Nx+2);
  ALLOCATE(OutDP,c->Nx+2);
  for(j=0; j<c->Nx+2; ++j) {
    OutIN[j]=0;
    OutDP[j]=0;
  }
  OutIN[0]=c->Nx+2;
  OutIN[1]=c->My+2;
  OutIN[2]=vl->time;

  OutIN[4]=c->IBL;
  OutIN[5]=c->jet.ia;
  OutIN[6]=c->jet.ib;

  OutDP[0]=c->Re;
  OutDP[1]=vl->Omega_tol;
  OutDP[2]=vl->Psi_tol;
  OutDP[3]=s->dx;
  OutDP[4]=s->dy;
  OutDP[5]=c->dt;
  OutDP[6]=c->Tol;
  OutDP[7]=c->Xmin;
  OutDP[8]=c->Xmax;
  OutDP[9]=c->Ymin;
  OutDP[10]=c->Ymax;
  OutDP[11]=c->A;
  OutDP[12]=c->jet.freq;
  OutDP[13]=c->jet.c0;

  writeFile(f->Omega,f->Psi,f->u,f->v,
      OutIN,OutDP,c->Nx,c->My,c->filename);
  deallocate(OutIN);
  deallocate(OutDP);
}
