#ifndef TOOLS_H
#define TOOLS_H

#include <stddef.h>

void die(const char* message) __attribute__((noreturn));
void* allocate(size_t bytes);
#define ALLOCATE(p,n) ((p)=allocate((n)*sizeof(*(p))))
void deallocate(void* memory);

void startParallel(void);
int parallelRan(void);
int parallelSize(void);
void stopParallel(void);
void broadcast(void* p, size_t s);
#define BROADCAST(o) broadcast(&(o),sizeof(o))

double square(double x);
double max(double a, double b);

#endif
