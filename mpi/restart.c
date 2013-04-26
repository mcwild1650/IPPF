#include "restart.h"
#include "tools.h"
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#define VALUE_CHARS 30
#define NUM_INTEGERS 7
#define NUM_DOUBLES 14

static void writeField(
    MPI_File file,
    int x0,
    int x1,
    int y0,
    int y1,
    MPI_Offset start,
    MPI_Offset rc,
    int eol,
    field f)
{
  int i,j;
  MPI_Offset ro,o;
  char v[VALUE_CHARS+1];
  for (i = y0; i < y1; ++i)
  {
    ro = start + rc*(i-y0);
    o = ro;
    for (j = x0; j < x1; ++j)
    {
      snprintf(v,VALUE_CHARS+1,"%30.13le",f[i][j]);
      MPI_File_write_at(file,o,v,VALUE_CHARS,MPI_BYTE,MPI_STATUS_IGNORE);
      o += VALUE_CHARS;
    }
    if (eol)
    {
      snprintf(v,2,"\n");
      MPI_File_write_at(file,o,v,1,MPI_BYTE,MPI_STATUS_IGNORE);
    }
  }
}

static void makeIntegers(
    config* c,
    vol* v,
    grid* g,
    int* ints)
{
  ints[0] = g->x+2;
  ints[1] = g->y+2;
  ints[2] = v->step;
  ints[3] = 0; //???
  ints[4] = c->IBL;
  ints[5] = c->jet.ia;
  ints[6] = c->jet.ib;
}

static void writeIntegers(
    MPI_File file,
    MPI_Offset o,
    int x,
    int* ints)
{
  int i;
  char v[VALUE_CHARS+1];
  for (i=0; i < NUM_INTEGERS; ++i)
  {
    snprintf(v,VALUE_CHARS+1,"%30d",ints[i]);
    MPI_File_write_at(file,o,v,VALUE_CHARS,MPI_BYTE,MPI_STATUS_IGNORE);
    o += VALUE_CHARS;
  }
  snprintf(v,VALUE_CHARS+1,"%30d",0);
  for (; i < x+2; ++i)
  {
    MPI_File_write_at(file,o,v,VALUE_CHARS,MPI_BYTE,MPI_STATUS_IGNORE);
    o += VALUE_CHARS;
  }
  snprintf(v,2,"\n");
  MPI_File_write_at(file,o,v,1,MPI_BYTE,MPI_STATUS_IGNORE);
}

static void makeDoubles(
    config* c,
    vol* v,
    space* s,
    double* doubles)
{
  doubles[0] = c->Re;
  doubles[1] = v->Omega_tol;
  doubles[2] = v->Psi_tol;
  doubles[3] = s->dx;
  doubles[4] = s->dy;
  doubles[5] = c->dt;
  doubles[6] = c->Tol;
  doubles[7] = c->Xmin;
  doubles[8] = c->Xmax;
  doubles[9] = c->Ymin;
  doubles[10] = c->Ymax;
  doubles[11] = c->A;
  doubles[12] = c->jet.freq;
  doubles[13] = c->jet.c0;
}

static void writeDoubles(
    MPI_File file,
    MPI_Offset o,
    int x,
    double* doubles)
{
  int i;
  char v[VALUE_CHARS+1];
  for (i=0; i < NUM_DOUBLES; ++i)
  {
    snprintf(v,VALUE_CHARS+1,"%30.13le",doubles[i]);
    MPI_File_write_at(file,o,v,VALUE_CHARS,MPI_BYTE,MPI_STATUS_IGNORE);
    o += VALUE_CHARS;
  }
  snprintf(v,VALUE_CHARS+1,"%30.13le",0.0);
  for (; i < x+2; ++i)
  {
    MPI_File_write_at(file,o,v,VALUE_CHARS,MPI_BYTE,MPI_STATUS_IGNORE);
    o += VALUE_CHARS;
  }
  snprintf(v,2,"\n");
  MPI_File_write_at(file,o,v,1,MPI_BYTE,MPI_STATUS_IGNORE);
}

void writeOutput(
    config* c,
    space* s,
    fields* f,
    vol* vl,
    grid* g)
{
  int dx = g->dx;
  int dy = g->dy;
  int n = g->n;
  int px = g->px;
  int py = g->py;
  int x = g->x;
  int y = g->y;
  int x0 = g->x0;
  int x1 = g->x1;
  int y0 = g->y0;
  int y1 = g->y1;
  MPI_Offset rc,fc,fs,o;
  MPI_File file;
  int amode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
  char* fn = c->filename;
  int eol;
  MPI_File_open(MPI_COMM_WORLD,fn,amode,MPI_INFO_NULL,&file);
  rc = VALUE_CHARS*(x+2) + 1;
  fc = rc*(y+2);
  fs = fc*4 + rc*2;
  o = rc*(py*dy+y0) + VALUE_CHARS*(px*dx+x0);
  MPI_File_preallocate(file,fs);
  eol = (px == n-1);
  writeField(file,x0,x1,y0,y1,o,rc,eol,f->Omega);
  o += fc;
  writeField(file,x0,x1,y0,y1,o,rc,eol,f->Psi);
  o += fc;
  writeField(file,x0,x1,y0,y1,o,rc,eol,f->u);
  o += fc;
  writeField(file,x0,x1,y0,y1,o,rc,eol,f->v);
  if ( ! parallelRank())
  {
    int ints[NUM_INTEGERS];
    double doubles[NUM_DOUBLES];
    makeIntegers(c,vl,g,ints);
    makeDoubles(c,vl,s,doubles);
    o = fc*4;
    writeIntegers(file,o,x,ints);
    o += rc;
    writeDoubles(file,o,x,doubles);
  }
  MPI_File_close(&file);
}

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

  printf("in restart Omega[0][0]=%lf\n",Omega[0][0]);
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
  fprintf(f_out,"\n");

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
  OutIN[2]=vl->step;

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
