#ifndef TOOLS_H
#define TOOLS_H

void die(const char* message) __attribute__((noreturn));
void* allocate(size_t bytes);
#define ALLOCATE(p,n) ((p)=allocate((n)*sizeof(*(p))))
void deallocate(void* memory);

void start_parallel(void);
int parallel_rank(void);
int parallel_size(void);
void stop_parallel(void);

#endif
