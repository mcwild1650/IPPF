#include "tools.h"
#ifdef PARALLEL
#include <mpi.h>
#endif
#include <stdlib.h>
#include <stdio.h>

void die(const char* message)
{
  fprintf(stderr,"%s\n",message);
  abort();
}

void* allocate(size_t bytes)
{
  if (bytes == 0)
    return NULL;
  void* m = malloc(bytes);
  if (m == NULL)
    die("malloc failed");
  return m;
}

void deallocate(void* memory)
{
  if (memory != NULL)
    free(memory);
}

double square(double x)
{
  return x*x;
}

double max(double a, double b)
{
  if (b>a) return b;
  return a;
}

#ifdef PARALLEL
void start_parallel(void)
{
  MPI_Init(NULL,NULL);
}

int parallel_rank(void)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  return rank;
}

int parallel_size(void)
{
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  return size;
}

void stop_parallel(void)
{
  MPI_Finalize();
}

void broadcast(void* p, size_t s)
{
  MPI_Bcast(p,s,MPI_BYTE,0,MPI_COMM_WORLD);
}
#endif


