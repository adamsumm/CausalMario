#ifndef ITEM_H 
#define ITEM_H 

typedef struct item_str* item;

item item_create(void);
void item_setclass(item itm, int absc);
int item_getclass(item itm);

void item_addedge(item itm, int rel, int dim, double val, int ndim, 
	int *participants);
int item_getnedges(item itm, int rel, int dim);
void item_getparticipants(item itm, int rel, int dim, 
	int ndim, int ind, int *participants);

void item_print(item itm, int rel, int dim, int ndim);
void item_copy(item is, item it);
void item_free(item is);

/* For frequency data */
double item_getedgeweight(item itm, int rel, int dim, int ind);

#endif /* ITEM_H*/
