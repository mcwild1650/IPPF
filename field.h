#ifndef FIELD_H
#define FIELD_H

typedef double** field;

field make_field(int x, int y);
void swap_fields(field* a, field* b);
void free_field(field f);

#endif
