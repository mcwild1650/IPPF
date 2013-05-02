#include "tools.h"
#ifdef PARALLEL
#include <mpi.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

void testDoubleSanity(double n)
{
  assert((n>0 || n<0 || n==0) && n!=INFINITY && n!=-INFINITY);
}

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

void zero_out(void* p, size_t s)
{
  memset(p,0,s);
}

#ifdef PARALLEL
void startParallel(void)
{
  MPI_Init(NULL,NULL);
}

int parallelRank(void)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  return rank;
}

int parallelSize(void)
{
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  return size;
}

void stopParallel(void)
{
  MPI_Finalize();
}
#endif

#if defined(__i386__)

static unsigned long long rdtsc(void)
{
  unsigned long long int x;
  __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
  return x;
}

#elif defined(__x86_64__)

static unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

#elif defined(__powerpc__)

static unsigned long long rdtsc(void)
{
  unsigned long long int result=0;
  unsigned long int upper, lower,tmp;
  __asm__ volatile(
      "0:                  \n"
      "\tmftbu   %0           \n"
      "\tmftb    %1           \n"
      "\tmftbu   %2           \n"
      "\tcmpw    %2,%0        \n"
      "\tbne     0b         \n"
      : "=r"(upper),"=r"(lower),"=r"(tmp)
      );
  result = upper;
  result = result<<32;
  result = result|lower;
  return(result);
}

#endif

static unsigned long long global_timer;

void startTimer(void)
{
  global_timer = rdtsc();
}

#ifndef CLOCK_RATE
#error "define CLOCK_RATE to 1.6e9 or something when compiling"
#endif

double stopTimer(void)
{
  unsigned long long end = rdtsc();
  return ((double)(end-global_timer))/CLOCK_RATE;
}

