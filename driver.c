#include "driver.h"
#include <string.h>
#include <stdio.h>

int get_args(int argc, char** argv, args* a)
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

void make_state(state* s)
{
  ALLOCATE(s->f,1);
  ALLOCATE(s->sp,1);
  ALLOCATE(s->c,1);
  ALLOCATE(s->dr,1);
  ALLOCATE(s->v,1);
}

void free_state(state* s)
{
  free_fields(s->f);
  deallocate(s->f);
  free_space(s->sp);
  deallocate(s->sp);
  deallocate(s->c);
  deallocate(s->dr);
  deallocate(s->v);
}

void init_state(state* s)
{
  make_fields(s->f,s->c->Nx,c->My);
  init_space(c,s->sp);
  init_fields(c,s->sp,s->f);
  init_derived(s->dr);
}

int main(int argc, char** argv)
{
  args a;
  if (!get_args(argc,argv,&a))
    return 0;
  state s;
  make_state(&s);
  config* c = s->c;
  if (a.m == config_mode)
  {
    ask_config(c);
    printf("What would you like to call the configuration file?\n");
    scanf("%s",a.configfile);
    write_config(a.configfile,c);
  }
  else if (a.m == run_mode)
  {
    read_config(a.configfile,c);
    init_state(&s);
  }
  else if (a.m == restart_mode)
  {
    read_config(a.configfile,c);
  }
  free_state(&s);
  return 0;
}
