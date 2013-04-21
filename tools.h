#ifndef TOOLS_H
#define TOOLS_H

#include <stddef.h>

void die(const char* message) __attribute__((noreturn));
void* allocate(size_t bytes);
#define ALLOCATE(p,n) ((p)=allocate((n)*sizeof(*(p))))
void deallocate(void* memory);

void start_parallel(void);
int parallel_rank(void);
int parallel_size(void);
void stop_parallel(void);
void broadcast(void* p, size_t s);
#define BROADCAST(o) broadcast(&(o),sizeof(o))
void zero_out(void* p, size_t s);
#define ZERO_OUT(o) zero_out(&(o),sizeof(o))

double square(double x);
double max(double a, double b);

#endif
