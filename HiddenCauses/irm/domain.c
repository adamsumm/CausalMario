#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "domain.h"
#include "frontlist.h"
#include "class.h"
#include "item.h"
#include "config.h"
#include "parameters.h"
#include "irmutils.h"

/*****************************************************************************
A domain (or type) is a set of items (objects), which we usually want to
cluster.  An instance of a domain is a place where it appears in a relation (eg
dimension 1 of relation 2).
*****************************************************************************/

struct domain_str
{ int nitem; 
  int clusterflag;	    /* cluster this domain? */
  int nclassmax;
  int ninstances;
  int **instances;	    /* (reln, dim) pairs where domain appears */
  frontlist classorder;
  itemclass *classes;
  item *items;
  double alpha;		    /* concentration hyperparameter for CRP */
  double alphahyp;	    /* prior on concentration hyperparameter */
  double prior;
  int label;
  int *init;
} domain_str;


domain domain_create(int nitem, int nclassmax, int clusterflag, double alpha,
	double alphahyp, int label, int *initclasses) {
  int i;
  domain dmn;
  dmn = (domain) my_malloc(sizeof(struct domain_str));
  dmn->nitem = nitem;
  dmn->nclassmax = nclassmax;
  dmn->clusterflag= clusterflag;
  dmn->alpha = alpha;
  dmn->alphahyp = alphahyp;
  dmn->ninstances = 0;

  dmn->instances= (int **) my_malloc(ps.maxdim*ps.maxrel*sizeof(int *));
  for (i = 0; i < ps.maxdim*ps.maxrel; i++)  {
    dmn->instances[i] = (int *) my_malloc(2*sizeof(int));
  }

  dmn->classorder = frontlist_create(nclassmax);
  dmn->classes = (itemclass *) my_malloc(nclassmax*sizeof(itemclass));
  for (i = 0; i < nclassmax; i++) {
    dmn->classes[i] = itemclass_create(nitem);
  }
  dmn->items= (item *) my_malloc(nitem*sizeof(item));
  for (i = 0; i < nitem; i++) {
    dmn->items[i] = item_create();
  }
  dmn->label = label;
  dmn->init = (int *) my_malloc(nitem*sizeof(int));
  if (initclasses != NULL) {
    for (i = 0; i < nitem; i++) {
      dmn->init[i] = initclasses[i];
    }
  }
  dmn->prior = 0.0;
  return dmn;
}

int domain_getlabel(domain dmn) {
  return dmn->label;
}

double domain_getalpha(domain dmn) {
  return dmn->alpha;
}

void domain_setalpha(domain dmn, double alpha) {
  dmn->alpha = alpha;
}

void domain_setalphahyp(domain dmn, double alphahyp) {
  dmn->alphahyp = alphahyp;
}

int domain_getninstance(domain dmn) {
  return dmn->ninstances;
}

void domain_addinstance(domain dmn, int relind, int dim) {
  dmn->instances[dmn->ninstances][0] = relind;
  dmn->instances[dmn->ninstances][1] = dim;
  dmn->ninstances++;
}

void domain_getinstance(domain dmn, int instind, int *relind, int*dim) {
  *relind = dmn->instances[instind][0];
  *dim = dmn->instances[instind][1];
}

int domain_getnclasses(domain dmn) {
  return frontlist_getnclasses(dmn->classorder);
}

int domain_getnclassmax(domain dmn) {
  return dmn->nclassmax;
}

int  domain_getclasssize(domain dmn, int relc) {
  int absc = domain_absclasslabel(dmn, relc);
  return itemclass_getnmembers(dmn->classes[absc]);
}

int domain_getclassmember(domain dmn, int relc, int membind) {
 int absc = domain_absclasslabel(dmn, relc);
 return itemclass_getmember(dmn->classes[absc], membind); 
}

void domain_addemptyclass(domain dmn) {
  frontlist_add_rellabel(dmn->classorder);
}

void domain_removeemptyclass(domain dmn) {
  frontlist_remove_rellabel(dmn->classorder, domain_getnclasses(dmn)-1);
}

void domain_additemtoclass(domain dmn, int itemind, int relc) {
  int absc;
  /* if class relc is empty */
  if (domain_getclasssize(dmn, relc) == 0) {
    if (relc != domain_getnclasses(dmn) ) {
      /* if we're adding an item to an empty class it should have relative
       * class label equal to the current class count */
      fprintf(stderr, "unexpected class addition\n"); exit(1); 
    }
    frontlist_add_rellabel(dmn->classorder);
  }
  absc = domain_absclasslabel(dmn, relc);
  itemclass_addmember(dmn->classes[absc], itemind); 
  item_setclass(dmn->items[itemind], absc);
}

void domain_removeitemfromclass(domain dmn, int itemind) {
  int absc;
  absc = item_getclass(dmn->items[itemind]);
  itemclass_removemember(dmn->classes[absc], itemind); 
  item_setclass(dmn->items[itemind], -1);
  /* if class absc is now empty */
  if (itemclass_getnmembers(dmn->classes[absc]) == 0) {
    frontlist_remove_abslabel(dmn->classorder, absc);
  }
}

int domain_getitemabsclass(domain dmn, int itemind) {
  return( item_getclass(dmn->items[itemind]));
}

int domain_getiteminitclass(domain dmn, int itemind) {
  return dmn->init[itemind];

}

int domain_getitemrelclass(domain dmn, int itemind) {
  int absc;
  absc= item_getclass(dmn->items[itemind]);
  return domain_relclasslabel(dmn, absc);
}


int domain_absclasslabel(domain dmn, int relc) {
  return(frontlist_abslabel(dmn->classorder, relc));
}

int domain_relclasslabel(domain dmn, int absc) {
  return(frontlist_rellabel(dmn->classorder, absc));
}

int domain_getnitems(domain dmn) {
  return dmn->nitem;
}

int domain_getnedges(domain dmn, int relind, int dim,  int itemind) {
  return( item_getnedges(dmn->items[itemind], relind, dim) ); 
}

double domain_getedge(domain dmn, int relind, int dim, int itemind, 
	int ndim, int edgeind, int *participants) {
  item_getparticipants(dmn->items[itemind], relind, dim, ndim, edgeind, 
	participants);
  return( item_getedgeweight(dmn->items[itemind], relind, dim, edgeind) );
}

void domain_addedge(domain dmn, int relind, int dim, double val, int itemind, 
	int ndim, int *participants) {
  item_addedge(dmn->items[itemind], relind, dim, val, ndim, participants);
}

void domain_updateprior(domain dmn) {
  double logprior;
  double alpha; 
  int i, nclasses, nitems;

  alpha = dmn->alpha;
  nclasses = domain_getnclasses(dmn);
  nitems = domain_getnitems(dmn);
  logprior = nclasses*log(alpha) + lgamma(alpha) - lgamma(alpha + nitems);
  for (i = 0; i < nclasses; i++) {
    logprior += lgamma((double) domain_getclasssize(dmn, i));
  }
  /* exponential prior on alpha */
  logprior += -alpha*1.0/ps.alphahyp;
  dmn->prior = logprior;
}

double domain_getprior(domain dmn) {
  return dmn->prior;
}

int domain_getclusterflag(domain dmn) {
  return dmn->clusterflag;
}

void domain_copy(domain ds, domain dt) {
  int i;
  /* don't copy things that will be unchanged */
  frontlist_copy(ds->classorder, dt->classorder);
  for (i = 0; i < ds->nclassmax; i++) {
    itemclass_copy(ds->classes[i], dt->classes[i]);   
  }
  for (i = 0; i < ds->nitem; i++) {
    item_copy(ds->items[i], dt->items[i]);   
  }
  dt->alpha = ds->alpha;
  dt->prior = ds->prior;
}

void domain_free(domain d) {
  int i;
  for (i = 0; i < ps.maxdim*ps.maxrel; i++)  {
    free(d->instances[i]);
  }
  free(d->instances);
  frontlist_free(d->classorder);
  for (i = 0; i < d->nclassmax; i++) {
    itemclass_free(d->classes[i]);
  }
  free(d->classes);
  for (i = 0; i < d->nitem; i++) {
    item_free(d->items[i]);
  }
  free(d->items);
  free(d->init);
  free(d);
}



void domain_printclasses(domain dmn) {
  int i, nclasses, abslabel;
  nclasses = domain_getnclasses(dmn);
  for (i = 0; i < nclasses; i++) {
    fprintf(stdout, "Class %d ", i);
    abslabel = domain_absclasslabel(dmn, i);
    itemclass_print(dmn->classes[abslabel]);  
  }
}

void domain_printassignments(FILE *fileptr, domain dmn) {
  int i;
  int abslabel, rellabel;
  for (i = 0; i < domain_getnitems(dmn); i++) {
    abslabel = domain_getitemabsclass(dmn, i);
    rellabel = domain_relclasslabel(dmn, abslabel);
    fprintf(fileptr, "%d ", rellabel);
  }
  fprintf(fileptr, "\n");
}
