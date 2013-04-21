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
  ALLOCATE(s->sp,1);
  ALLOCATE(s->c,1);
  ALLOCATE(s->dr,1);
  ALLOCATE(s->v,1);
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
  makeFields(s->f,s->c->Nx,s->c->My);
  initSpace(s->c,s->sp);
  initFields(s->c,s->sp,s->dr,s->f);
  initDerived(s->c,s->sp,s->dr);
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
