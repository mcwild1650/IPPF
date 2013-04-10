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

int main(int argc, char** argv)
{
  args a;
  if (!get_args(argc,argv,&a))
    return 0;
  return 0;
}
