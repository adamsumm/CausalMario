#include <stdio.h>
#include <stdlib.h>
#include "item.h"
#include "config.h"
#include "parameters.h"
#include "irmutils.h"

/*****************************************************************************
An item is an object (member of a domain). It has a class (absclass), and
participates in some relations.

ISSUE: the edges involving an item are represented by every chain in a chain
collection. This is unnecessary, and could be bad for very large graphs if one
copy of the graph barely fits into memory.
*****************************************************************************/

struct item_str 
{ int absclass; 
  int **nedges; /* for each instance of the item's parent domain, number of 
		   1-edges involving this item */
  int ****objtuples; /* links involving this item for each instance of its 
			 parent domain */
  double ***edgeweights; /* values of edges */
} item_str;

/* NB: copies of items can share everything except absclass */

void item_setclass(item itm, int absc) {
  itm->absclass = absc;
}

int item_getclass(item itm) {
   void* ptr = itm;
   
  return(itm->absclass);
}

item item_create(void) {
  item itm; 
  int relind, dimind, tupleind;
  itm = (item) my_malloc(sizeof(struct item_str));
  itm->absclass   = -1;
  itm->nedges     = (int **) my_malloc(ps.maxrel*sizeof(int *));
  itm->objtuples  = (int ****) my_malloc(ps.maxrel*sizeof(int ***));
  itm->edgeweights= (double ***) my_malloc(ps.maxrel*sizeof(int **));

  for (relind = 0; relind < ps.maxrel; relind++) {
    itm->nedges[relind]      = (int *) my_malloc(ps.maxdim* sizeof(int));
    itm->objtuples[relind]   = (int ***) my_malloc(ps.maxdim* sizeof(int *));
    itm->edgeweights[relind] = (double **) my_malloc(ps.maxdim* sizeof(int *));

    for (dimind = 0; dimind < ps.maxdim; dimind++) {
      itm->objtuples[relind][dimind] = (int **) 
	    my_malloc(ps.maxobjtuples* sizeof(int *));
      itm->edgeweights[relind][dimind] = (double *) 
	    my_malloc(ps.maxobjtuples* sizeof(double));
      for (tupleind = 0; tupleind < ps.maxobjtuples; tupleind++) {
	itm->objtuples[relind][dimind][tupleind] = (int *) 
		my_malloc(ps.maxdim* sizeof(int *));
      }
    }
  }

  for (relind = 0; relind < ps.maxrel; relind++) {
    for (dimind = 0; dimind < ps.maxdim; dimind++) {
       itm->nedges[relind][dimind] = 0;
    }
  }
  return itm;
}

/* add the NDIMensional edge in PARTICIPANTS, which includes ITS in dimensions
 * DIM of relation REL and has weight VAL */

void item_addedge(item itm, int rel, int dim, double val, int ndim, 
	int *participants) {
  int i, newedgeind;
  int **edgelist;

  edgelist = (int **) itm->objtuples[rel][dim];
  newedgeind = itm->nedges[rel][dim]++;

  if (newedgeind > ps.maxobjtuples) {
    fprintf(stderr, "ERROR: more edges than expected");
  }
  for (i = 0; i<ndim; i++) {
    edgelist[newedgeind][i] = participants[i];
  }
  itm->edgeweights[rel][dim][newedgeind]=val;
}

/* get nedges involving ITS in dimension DIM of relation REL */
int item_getnedges(item itm, int rel, int dim) {
  return itm->nedges[rel][dim];
}

/* get NDIMensional edge IND involving ITS in dimension DIM of relation REL */

void item_getparticipants(item itm, int rel, int dim, 
	int ndim, int ind, int *participants) {
  int i, *edgemembers;
  edgemembers= (int *) itm->objtuples[rel][dim][ind];

  for (i = 0; i < ndim; i++) {
    participants[i] = edgemembers[i];
  }
}

void item_copy(item is, item it) {
  it->absclass = is->absclass;
  /* don't need to copy the rest*/
}

void item_free(item i) {
  int relind, dimind, tupleind;

  for (relind = 0; relind < ps.maxrel; relind++) {
    for (dimind = 0; dimind < ps.maxdim; dimind++) {
      for (tupleind = 0; tupleind < ps.maxobjtuples; tupleind++) {
        free(i->objtuples[relind][dimind][tupleind]);
      }
      free(i->objtuples[relind][dimind]);
      free(i->edgeweights[relind][dimind]);
    }
    free(i->nedges[relind]);
    free(i->objtuples[relind]);
    free(i->edgeweights[relind]);
  }
  free(i->nedges);
  free(i->objtuples);
  free(i->edgeweights);
  free(i);
}


void item_print(item itm, int rel, int dim, int ndim) {
  int i, j, nedges;
  int *participants;

  participants=  (int *) my_malloc(ps.maxdim*sizeof(int));

  nedges = item_getnedges(itm, rel, dim);
  for (i = 0; i < nedges; i++) {
    item_getparticipants(itm, rel, dim, ndim, i, participants);      
    fprintf(stderr, "edge %d: ", i);
    for (j = 0; j < ndim; j++) {
      fprintf(stderr, "%d ", participants[j]);
    }
    fprintf(stderr, "\n");
  }
  free(participants);
}


double item_getedgeweight(item itm, int rel, int dim, int ind) {
  return itm->edgeweights[rel][dim][ind];
}


