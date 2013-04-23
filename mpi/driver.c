#include "driver.h"
#include "tools.h"
#include "restart.h"
#include <string.h>
#include <stdio.h>

int getArgs(int argc, char** argv, args* a)
{
  if (argc == 2)
  {
    if (!strcmp(argv[1],"config"))
    {
      a->m = config_mode;
      a->configfile[0] = '\0';
      return 1;
    }
  }
  if (argc == 3)
  {
    if (!strcmp(argv[1],"run"))
    {
      a->m = run_mode;
      strncpy(a->configfile,argv[2],FILENAME_LENGTH);
      return 1;
    }
    if (!strcmp(argv[1],"restart"))
    {
      a->m = restart_mode;
      strncpy(a->configfile,argv[2],FILENAME_LENGTH);
      return 1;
    }
  }
  printf("usage:\n");
  printf("`%s config` start here\n",argv[0]);
  printf("`%s run     <configfile>` run from a config file\n",argv[0]);
  printf("`%s restart <configfile>` to restart from saved data\n",argv[0]);
  return 0;
}

void makeState(state* s)
{
  ALLOCATE(s->f,1);
  ZERO_OUT(*(s->f));
  ALLOCATE(s->sp,1);
  ZERO_OUT(*(s->sp));
  ALLOCATE(s->c,1);
  ZERO_OUT(*(s->c));
  ALLOCATE(s->dr,1);
  ZERO_OUT(*(s->dr));
  ALLOCATE(s->v,1);
  ZERO_OUT(*(s->v));
}

void freeState(state* s)
{
  freeFields(s->f);
  deallocate(s->f);
  freeSpace(s->sp);
  deallocate(s->sp);
  deallocate(s->c);
  deallocate(s->dr);
  deallocate(s->v);
}

void initState(state* s)
{
  makeFields(s->f,s->c->tot_Nx,s->c->tot_My);
  initSpace(s->c,s->sp);
  initVolatile(s->c,0,s->v);
  initDerived(s->c,s->sp,s->dr);
  initFields(s->c,s->sp,s->dr,s->f);
}

int main(int argc, char** argv)
{
  args a;
  if (!getArgs(argc,argv,&a))
    return 0;
  state s;
  makeState(&s);
  config* c = s.c;
  if (a.m == config_mode)
  {
    askConfig(c);
    printf("What would you like to call the configuration file?\n");
    scanf("%s",a.configfile);
    writeConfig(a.configfile,c);
  }
  else if (a.m == run_mode)
  {
    readConfig(a.configfile,c); 
    initState(&s); 
    initMpi(&s);
    calculate(s.c,s.sp,s.f,s.dr,s.v);
    writeOldRestart(s.c,s.sp,s.f,s.dr,s.v);
    freeState(&s);
  }
  else if (a.m == restart_mode)
  {
    readConfig(a.configfile,c);
  }
  return 0;
}


void initMpi(state* s)
{
 
  startParallel();
   //mpipara *mpi_para =*(s->c->mpi_para);
  s->c->mpi_para.myrank=parallelRank( );
  s->c->mpi_para.mpi_size=parallelSize( );

  int n=2,m=2; ////change later 
  
   int myrank=s->c->mpi_para.myrank;
 int mpi_size=s->c->mpi_para.mpi_size;

  int i,total=0;
  int total_Nx=s->c->tot_Nx;
  int total_My=s->c->tot_My;
 
  printf("\n my rank:%d",myrank);
  for (i=0;i<n;i++)
  	{
  		if(i==myrank%n)  s->c->mpi_para.my_x_start=total;
  		
  		if(i<(total_Nx+2)%n )
  			total+=(total_Nx+2)/n+1;
  		else
  			total+=(total_Nx+2)/n;
  	 	if(i==myrank%n) 
  	 	{
  	 		s->c->mpi_para.my_x_end=total-1;
  	 		s->c->Nx=s->c->mpi_para.my_x_end-s->c->mpi_para.my_x_start+1;
  	 		break;
  	 	}
  	}
 
    printf ("Characters: %d %d \n", s->c->mpi_para.my_x_end ,total);
   total=0;
  for (i=0;i<m;i++)
  	{
  		if(i==myrank/n)  s->c->mpi_para.my_y_start=total;
  		if(i<(total_My+2)%m )
  			total+=(total_My+2)/m+1;
  		else
  			total+=(total_My+2)/m;
  			
  		if(i==myrank/n) 
  		{
  			s->c->mpi_para.my_y_end=total-1;
  			s->c->My=s->c->mpi_para.my_y_end-s->c->mpi_para.my_y_start+1;
  			break;
  		}
  	 	 
  	}
  	s->c->mpi_para.m=m;
  	s->c->mpi_para.n=n;
  	s->c->mpi_para.px=myrank%n;
  	s->c->mpi_para.py=myrank/n;
 printf ("Cha : %d %d \n", s->c->mpi_para.my_x_end ,total);
 printf("\nmy rank:%d,x: %d - %d ;  y: %d - %d",myrank,s->c->mpi_para.my_x_start,s->c->mpi_para.my_x_end,s->c->mpi_para.my_y_start,s->c->mpi_para.my_y_end);
 
 
}
