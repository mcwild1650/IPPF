#ifndef TOOLS_H
#define TOOLS_H

#define PARALLEL

#include <stddef.h>

void testDoubleSanity(double n);
void die(const char* message) __attribute__((noreturn));
void* allocate(size_t bytes);
#define ALLOCATE(p,n) ((p)=allocate((n)*sizeof(*(p))))
void deallocate(void* memory);

void startParallel(void);
int parallelRank(void);
int parallelSize(void);
void stopParallel(void);
void zero_out(void* p, size_t s);
#define ZERO_OUT(o) zero_out(&(o),sizeof(o))
#define MAKE(p) do{ALLOCATE(p,1);ZERO_OUT(*(p));}while(0)

#define SQUARE(x) ((x)*(x))
#define MAX(a,b) (((b)>(a))?(b):(a))
#define MIN(a,b) (((b)<(a))?(b):(a))

#ifdef DEBUG
#define debug(x) {fprintf(stderr,"FILE: %s LINE: %d RANK: %d STATE: %s\n",__FILE__, __LINE__,parallelRank(),x);}
#else
#define debug(x) ((void)0)
#endif

void startTimer(void);
double stopTimer(void);

#endif
