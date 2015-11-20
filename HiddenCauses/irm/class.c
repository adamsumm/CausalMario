#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "frontlist.h"
#include "irmutils.h"
#include "class.h"


/*****************************************************************************
an itemclass is a latent class or cluster
*****************************************************************************/

struct itemclass_str
{ int nmembers; 
  frontlist members;
} itemclass_str;

int itemclass_getnmembers(itemclass icl) {
  return icl->nmembers;
}

itemclass itemclass_create(int nelementsmax) {
  itemclass icl;
  icl = (itemclass) my_malloc(sizeof(struct itemclass_str));
  icl->nmembers = 0;
  icl->members = frontlist_create(nelementsmax);
  return icl;
}

int itemclass_getmember(itemclass icl, int relind) {
 return frontlist_abslabel(icl->members, relind);
}

void itemclass_addmember(itemclass icl, int absind) {
 frontlist_add_abslabel(icl->members, absind);
 icl->nmembers++;
}

void itemclass_removemember(itemclass icl, int absind) {
 frontlist_remove_abslabel(icl->members, absind);
 icl->nmembers--;
}

void itemclass_copy(itemclass is, itemclass it) {
  it->nmembers = is->nmembers;
  frontlist_copy(is->members, it->members);
}

void itemclass_free(itemclass i) {
  frontlist_free(i->members);
  free(i);
}

void itemclass_print(itemclass icl) {
  int i;
  fprintf(stdout, "nelements: %d \n", icl->nmembers);
  for (i = 0; i < icl->nmembers; i++) {
    fprintf(stdout, "%d ", itemclass_getmember(icl, i));
  }
  fprintf(stdout, "\n");
}









