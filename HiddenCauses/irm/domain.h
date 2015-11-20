#ifndef DOMAIN_H 
#define DOMAIN_H 
#include <stdio.h>
#include "types.h"


domain domain_create(int nitem, int nclassmax, int clusterflag, double alpha,
		     double alphahyp, int label, int *initclasses);

double domain_getalpha(domain dmn);
void domain_setalpha(domain dmn, double alpha);
void domain_setalphahyp(domain dmn, double alphahyp);
int domain_getlabel(domain dmn);

/* an instance is a (relation, dimension) pair */
int  domain_getninstance(domain dmn);
void domain_addinstance(domain dmn, int relind, int dimind);
void domain_getinstance(domain dmn, int instind, int *relind, int *dimind);

/* we use relative and absolute class labels: see frontlist.h */
int domain_getnclasses(domain dmn);
int  domain_getnclassmax(domain dmn);
int  domain_getclasssize(domain dmn, int relc);
int  domain_getclassmember(domain dmn, int relc, int membind);
void domain_additemtoclass(domain dmn, int itemind, int relc); 
void domain_removeitemfromclass(domain dmn, int itemind); 
int  domain_getitemabsclass(domain dmn, int itemind);
int  domain_getitemrelclass(domain dmn, int itemind);
int  domain_getiteminitclass(domain dmn, int itemind);
int  domain_absclasslabel(domain dmn, int relc);
int  domain_relclasslabel(domain dmn, int absc);
int  domain_getnitems(domain dmn);

/* for chain_classlike */
void domain_addemptyclass(domain dmn); 
void domain_removeemptyclass(domain dmn); 



/* edges have weights */
int  domain_getnedges(domain dmn, int relind, int dim, int itemind);
double domain_getedge(domain dmn, int relind, int dim, int itemind, 
	int ndim, int edgeind, int *participants);
void domain_addedge(domain dmn, int relind, int dim, double val, int itemind, 
	int ndim, int *participants);

void domain_updateprior(domain dmn);
double domain_getprior(domain dmn);

int domain_getclusterflag(domain dmn);

void domain_copy(domain ds, domain dt);
void domain_free(domain d);

void domain_printclasses(domain dmn);
void domain_printassignments(FILE *fileptr, domain dmn);

#endif /* DOMAIN_H */
