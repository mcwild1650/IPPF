#include "config.h"
#include "tools.h"
#include <stdio.h>
#include <mpi.h>

void askConfig(config* c)
{
  printf("Nx= ");
  scanf("%d",&(c->Nx));
  printf("My= ");
  scanf("%d",&(c->My));
  printf("What Reynolds Number would you like to run?\n");
  scanf("%lf",&(c->Re));
  printf("What value of circulation parameter (A-Tilde),\n");
  printf("would you like to use?\n");
  scanf("%lf",&(c->A));
  printf("How many time steps would you like to run?\n");
  scanf("%d",&(c->Ot));
  printf("How many time steps between reports?\n");
  scanf("%d",&(c->report));
  printf("dt= ");
  scanf("%lf",&(c->dt));
  printf("To what tolerance level would you like to iterate?\n");
  scanf("%lf",&(c->Tol));
  printf("There are %12d grid lines in the vertical direction.\n",c->My);
  printf("How many do you want to be in the boundary layer?\n");
  scanf("%d",&(c->IBL));
  printf("Would you like to include a jet in the simulation? (y/n): ");
  char yn;
  scanf(" %c",&yn);
  printf("read in '%c'\n",yn);
  if (yn == 'y')
  {
    c->jet.exists = 1;
    printf("What is the start location of the jet?\n");
    printf("(Range of 1 to %d)\n",c->Nx);
    scanf("%d",&(c->jet.ia));
    printf("What is the end location?\n");
    scanf("%d",&(c->jet.ib));
    printf("What is the amplitude of the jet?\n");
    scanf("%lf",&(c->jet.c0));
    printf("What is the frequency of the jet?\n");
    scanf("%lf",&(c->jet.freq));
  }
  else
  {
    c->jet.exists = 0;
    c->jet.ia = c->Nx+4;
    c->jet.ib = c->Nx+4;
    c->jet.c0 = 0;
    c->jet.freq = 0;
  }
  printf("What would you like to call the output file?\n");
  scanf("%s",c->filename);
  printf("How many iterations between incremental file writes?\n");
  printf("(use -1 to turn off incremental save)\n");
  scanf("%d",&(c->psave));
  c->Xmin = -20;
  c->Xmax =  20;
  c->Ymin = 1;
  c->Ymax = 11;
}

void writeConfig(const char* filename, const config* c)
{
  FILE* file = fopen(filename,"w");
  fprintf(file,"Nx %d\n",c->Nx);
  fprintf(file,"My %d\n",c->My);
  fprintf(file,"Reynolds Number %lf\n",c->Re);
  fprintf(file,"Circulation %lf\n",c->A);
  fprintf(file,"Timesteps %d\n",c->Ot);
  fprintf(file,"Report Steps %d\n",c->report);
  fprintf(file,"dt %lf\n",c->dt);
  fprintf(file,"Tolerance %lf\n",c->Tol);
  fprintf(file,"Boundary Layer %d\n",c->IBL);
  fprintf(file,"Jet %d\n",c->jet.exists);
  fprintf(file,"Jet Start %d\n",c->jet.ia);
  fprintf(file,"Jet End %d\n",c->jet.ib);
  fprintf(file,"Jet Amplitude %lf\n",c->jet.c0);
  fprintf(file,"Jet Frequency %lf\n",c->jet.freq);
  fprintf(file,"Filename %s\n",c->filename);
  fprintf(file,"Save Steps %d\n",c->psave);
  fclose(file);
}

void readConfig(const char* filename, config* c)
{
  if (parallelRank()==0)
  {
    FILE* file = fopen(filename,"r");
    fscanf(file,"Nx %d\n",&(c->Nx));
    fscanf(file,"My %d\n",&(c->My));
    fscanf(file,"Reynolds Number %lf\n",&(c->Re));
    fscanf(file,"Circulation %lf\n",&(c->A));
    fscanf(file,"Timesteps %d\n",&(c->Ot));
    fscanf(file,"Report Steps %d\n",&(c->report));
    fscanf(file,"dt %lf\n",&(c->dt));
    fscanf(file,"Tolerance %lf\n",&(c->Tol));
    fscanf(file,"Boundary Layer %d\n",&(c->IBL));
    fscanf(file,"Jet %d\n",&(c->jet.exists));
    fscanf(file,"Jet Start %d\n",&(c->jet.ia));
    fscanf(file,"Jet End %d\n",&(c->jet.ib));
    fscanf(file,"Jet Amplitude %lf\n",&(c->jet.c0));
    fscanf(file,"Jet Frequency %lf\n",&(c->jet.freq));
    fscanf(file,"Filename %s\n",c->filename);
    fscanf(file,"Save Steps %d\n",&(c->psave));
    fclose(file);
    c->Xmin = -20;
    c->Xmax =  20;
    c->Ymin = 1;
    c->Ymax = 11;
  }
  broadcastConfig(c);
}

void broadcastConfig(config* c)
{
  MPI_Bcast(c,sizeof(*c),MPI_PACKED,0,MPI_COMM_WORLD);
}
