#include "tools.h"
#include <mpi.h>
#include <stdlib.h>

void die(const char* message) __attribute__((noreturn))
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
