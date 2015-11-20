#include <stdio.h>
#include <stdlib.h>
#include "frontlist.h"
#include "irmutils.h"

/*****************************************************************************
A frontlist is a data structure for representing a set of numbers between 0 and
nelementsmax. Each item in the set has a relative label (its position in the
set: takes value between 0 and k if the set has n members), and an absolute
label (which can be seen as its true identity). The absolute labels of the
elements currently in the set are maintained at the front of absorder.
invabsorder represents where each abslabel appears in absorder.

This business is a way to represent sets that grow and shrink in size without
worrying about allocating and freeing memory. We therefore require that the
maximum size of each set is specified in advance.

*****************************************************************************/

struct frontlist_str 
{ int nelementsmax;
  int nelements;
  int *absorder;	/* map from rel label to abs label*/
  int *invabsorder;	/* map from abs label to rel label*/
} frontlist_str;

frontlist frontlist_create(int nelementsmax)
{ 
  int i;
  frontlist fl;
  fl = (frontlist) my_malloc(sizeof(struct frontlist_str));
  fl->nelementsmax = nelementsmax;
  fl->nelements = 0;
  fl->absorder	  = (int *) my_malloc(nelementsmax*sizeof(int));
  fl->invabsorder = (int *) my_malloc(nelementsmax*sizeof(int));
  for (i=0; i < nelementsmax; i++) {
    fl->absorder[i] = i;
    fl->invabsorder[i] = i;
  }
  return fl;
}

int frontlist_getnclasses(frontlist fl) {
  return (fl->nelements);
}

/* find abslabel corresponding to relind*/
int frontlist_abslabel(frontlist fl, int relind) {
  return fl->absorder[relind];
}

/* find rellabel corresponding to absind*/
int frontlist_rellabel(frontlist fl, int absind) {
  return fl->invabsorder[absind];
}

/* add abslabel to the frontlist */
void frontlist_add_abslabel(frontlist fl, int abslabel) {
  int lastrel, lastabs, rellabel;

  rellabel = fl->invabsorder[abslabel];
  lastrel= fl->nelements;
  lastabs = fl->absorder[lastrel];

  if (rellabel < fl->nelements) {
    fprintf(stderr, "abslabel already  in frontlist\n"); exit(1); 
  }
  if ( abslabel >= fl->nelementsmax) {
    fprintf(stderr, "absind is too big for frontlist\n"); exit(1); 
  }


  fl->absorder[rellabel] = lastabs;
  fl->invabsorder[lastabs] = rellabel;
  fl->absorder[lastrel] = abslabel;
  fl->invabsorder[abslabel] = lastrel;
  fl->nelements++;
  return; 
}

/* remove abslabel from the frontlist */
void frontlist_remove_abslabel(frontlist fl, int absind) {
  int rellabel;
  if (fl->invabsorder[absind] >= fl->nelements) {
    fprintf(stderr, "absind isn't in frontlist\n"); exit(1); 
  }
  if ( absind >= fl->nelementsmax) {
    fprintf(stderr, "absind is too big for frontlist\n"); exit(1); 
  }

  rellabel = fl->invabsorder[absind];
  frontlist_remove_rellabel(fl, rellabel);
  return; 
}


/* increase the size of the list by one. Choose a random abslabel for the new
 * element */
int frontlist_add_rellabel(frontlist fl) {
  int newabs;
  if (fl->nelements == fl->nelementsmax) {
    fprintf(stderr, "frontlist already full\n"); exit(1); 
  }
  newabs = fl->absorder[fl->nelements];
  fl->invabsorder[newabs] = fl->nelements;
  fl->nelements++;
  return( fl->nelements - 1 ); 
}

/* remove rellabel from frontlist*/
void frontlist_remove_rellabel(frontlist fl, int rellabel) {
  int lastrel, lastabs, abslabel;
  lastrel = fl->nelements - 1;
  lastabs = fl->absorder[lastrel];
  abslabel = fl->absorder[rellabel];
  fl->absorder[rellabel] = lastabs;
  fl->invabsorder[lastabs] = rellabel;
  fl->absorder[lastrel] = abslabel;
  fl->invabsorder[abslabel] = lastrel;
  fl->nelements--;
}

void frontlist_copy(frontlist fs, frontlist ft) {
  int i;
  ft->nelementsmax = fs->nelementsmax;
  ft->nelements = fs->nelements;
  for (i = 0; i < fs->nelementsmax; i++) {
    ft->absorder[i] = fs->absorder[i];
    ft->invabsorder[i] = fs->invabsorder[i];
  }
}

void frontlist_free(frontlist f) {
  free(f->absorder);
  free(f->invabsorder);
  free(f);
}


void frontlist_print(frontlist fl) {
  int i;
  fprintf(stdout, "nelements: %d \n", fl->nelements);
  for (i = 0; i < fl->nelementsmax; i++) {
    fprintf(stdout, "%d ", fl->absorder[i]);
  }
  fprintf(stdout, "\n");
  for (i = 0; i < fl->nelementsmax; i++) {
    fprintf(stdout, "%d ", fl->invabsorder[i]);
  }
  fprintf(stdout, "\n");
}



