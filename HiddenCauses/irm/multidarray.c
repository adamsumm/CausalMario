#include <stdlib.h>
#include "config.h"
#include "parameters.h"
#include "irmutils.h"
#include "multidarray.h"

/*****************************************************************************
Multidimensional arrays. We need something general here so we can deal with
relations of arbitrary arity.
*****************************************************************************/

struct multidarray_str
{ int ndimensions; 
  int *dimsizesmax; /* max size of each dimension */
  int *elements;
} multidarray_str;

int *multidarray_init_elements(int ndimensions, int *maxsizes);
void multidarray_free_elements(int *array, int ndimensions, int *maxsizes);

multidarray multidarray_create(int ndimensions, int *maxsizes) {
  int i;
  multidarray mda;
  mda = (multidarray) my_malloc(sizeof(struct multidarray_str));
  mda->ndimensions = ndimensions;
  mda->dimsizesmax = my_malloc(ndimensions*sizeof(int));
  for (i = 0; i < ndimensions; i++) {
    mda->dimsizesmax[i] = maxsizes[i];
  }
  mda->elements = multidarray_init_elements( ndimensions, maxsizes);
  return mda;
}

datael multidarray_get(multidarray mda, int *indices) 
{ int i;
  int * elements_ptr;
  elements_ptr = mda ->elements;
  for (i = 0; i < mda->ndimensions-1; i++) {
    elements_ptr = (int *) elements_ptr[indices[i]]; 
  }
  /* at the end of the string of pointers we should have a datael */
  return((datael) elements_ptr[indices[i]]);
}

void multidarray_set(multidarray mda, datael val, int *indices) 
{ int i;
  int * elements_ptr;
  elements_ptr = mda ->elements;
  for (i = 0; i < mda->ndimensions - 1; i++) {
    elements_ptr = (int *) elements_ptr[indices[i]]; 
  }
  datael_copy(val, (datael) elements_ptr[indices[i]]);
}

void multidarray_inc(multidarray mda, datael val, int *indices)
{
  datael_add(multidarray_get(mda, indices), val);
}

void multidarray_dec(multidarray mda, datael val, int *indices)
{
  datael_minus(multidarray_get(mda, indices),  val);
}


/* allocate memory for elements and initialize to zero */
int *multidarray_init_elements(int ndimensions, int *maxsizes) 
{ int i;
  int *array_ptr;
  int *newsizes;
  datael zeroel;

  newsizes=  (int *) my_malloc(ps.maxdim*sizeof(int));
  if (ndimensions == 1) {
    array_ptr = (int *) my_malloc(maxsizes[0] * sizeof(datael)); 
    for (i = 0; i < maxsizes[0]; i++) {
      zeroel = datael_create();
      array_ptr[i] = (int) zeroel;
    }
  } else {
    array_ptr = (int *) my_malloc(maxsizes[0] * sizeof(int)); 
    for (i = 1; i < ndimensions; i++) {
      newsizes[i-1] = maxsizes[i];
    }
    for (i = 0; i < maxsizes[0]; i++) {
      array_ptr[i] = 
	    (int) multidarray_init_elements(ndimensions-1, newsizes);
    }
  }
  free(newsizes);
  return array_ptr;
}

void multidarray_free_elements(int *array_ptr, int ndimensions, int *maxsizes) 
{ int i;
  int *newsizes;
  newsizes=  (int *) my_malloc(ps.maxdim*sizeof(int));
  if (ndimensions == 1) {
    for (i = 0; i < maxsizes[0]; i++) {
      datael_free((datael) array_ptr[i]);
    }
  } else {
    for (i = 1; i < ndimensions; i++) {
      newsizes[i-1] = maxsizes[i];
    }
    for (i = 0; i < maxsizes[0]; i++) {
      multidarray_free_elements((int *) array_ptr[i], ndimensions-1, newsizes); 
    }
  }
  free(newsizes);
  free(array_ptr);
}

void multidarray_copy(multidarray mdas, multidarray mdat) {
  int i, ndim, cellnum;
  int *absindex; 
  datael countsel;

  absindex=  (int *) my_malloc(ps.maxdim*sizeof(int));
  ndim = mdas->ndimensions;
  mdat->ndimensions = mdas->ndimensions;
  cellnum = 1;
  for (i = 0; i < ndim; i++) {
    absindex[i] = 0;
    cellnum *= mdas->dimsizesmax[i];
  }
  for (i = 0; i < cellnum; i++) {
    countsel = multidarray_get(mdas, absindex);
    multidarray_set(mdat, countsel, absindex);
    multidarray_incrementindex(ndim, absindex, mdas->dimsizesmax);
  }
  free(absindex);
}

void multidarray_free(multidarray mda) {
  multidarray_free_elements(mda->elements, mda->ndimensions, mda->dimsizesmax);
  free(mda->dimsizesmax);
  free(mda);
}

void multidarray_incrementindex(int ndim, int *index, int *classsizes) {
  int i = ndim;
  for (i = ndim-1; i >= 0; i--) {
    index[i]++;
    if (index[i] == classsizes[i]) {
      index[i] = 0;
    } else {
      break;
    }
  }
}

