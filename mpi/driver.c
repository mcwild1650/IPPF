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
  MAKE(s->c);
  MAKE(s->g);
  MAKE(s->sp);
  MAKE(s->v);
  MAKE(s->dr);
  MAKE(s->f);
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
  partition(s->c,s->g);
  initSpace(s->c,s->g,s->sp);
  initVolatile(s->c,0,s->v);
  initDerived(s->c,s->sp,s->dr);
  makeFields(s->f,s->g->len_x,s->g->len_y);
  initFields(s->c,s->sp,s->dr,s->f,s->g);
}

int main(int argc, char** argv)
{
  args a;
  state s;
  startParallel();
  if (!getArgs(argc,argv,&a))
    return 0;
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
    calculate(s.c,s.sp,s.g,s.f,s.dr,s.v);
    writeOutput(s.c,s.sp,s.f,s.v,s.g);
    freeState(&s);
  }
  else if (a.m == restart_mode)
  {
    readConfig(a.configfile,c);
  }
  stopParallel();
  return 0;
}

