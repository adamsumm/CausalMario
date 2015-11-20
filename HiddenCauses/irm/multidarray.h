#ifndef MULTIDARRAY_H 
#define MULTIDARRAY_H

#include "irmutils.h"
typedef struct multidarray_str *multidarray;

multidarray multidarray_create(int ndimensions, int *maxsizes);
datael multidarray_get(multidarray mdarray, int *indices);
void multidarray_set(multidarray mdarray, datael val, int *indices);
/* increment mdarray */
void multidarray_inc(multidarray mdarray, datael val, int *indices);
/* decrement mdarray */
void multidarray_dec(multidarray mdarray, datael val, int *indices);

/* increment index into mdarray*/
void multidarray_incrementindex(int ndim, int *index, int *classsizes);
void multidarray_copy(multidarray mdas, multidarray mdat);
void multidarray_free(multidarray mda);

#endif /* MULTIDARRAY_H*/
