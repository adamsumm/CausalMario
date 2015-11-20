#ifndef RELATION_H 
#define RELATION_H 
#include "types.h"


relation relation_create(int ndim, double betaprob, double betamag, 
			 double *nig, domain *doms, int t);
int relation_getdim(relation r);
void relation_getbeta(relation r, double *b1, double *b2);
void relation_setbeta(relation r, double b1, double b2);
domain * relation_getdoms(relation r);
void relation_addedges(relation r, datael weight, int *absclasses);
void relation_removeedges(relation r, datael weight, int *absclasses);
datael relation_edgeweight(relation r, int *absclasses); 
void relation_updatelikelihood(relation r);
double relation_likelihood(relation r);
int relation_getmissing(relation r);
void relation_setmissing(relation r, int val);
int relation_getdtype(relation r);
void relation_setdtype(relation r, int val);
void relation_printcounts(relation r);
void relation_copy(relation rs, relation rt);
void relation_free(relation r);
void relation_hypupdate(relation r, double temp, 
	int (*proposal_choose) (double, double, double), int *changeflag ) ;

double relation_getbetaprop(relation r);
double relation_getbetamag(relation r);
double relation_getnig(relation r, int i);


#endif /* RELATION_H*/
