#ifndef CHAIN_H 
#define CHAIN_H 
#include "types.h"


chain chain_create(int ndomains, int nrelns, double temp);
int chain_getndomains(chain chn);
int chain_getnrelations(chain chn);
double chain_gettemp(chain chn);
void chain_settemp(chain chn, double t);

void chain_adddomain(chain chn, int domind, int nitem, int nclassmax, 
	int clusterflag, double alpha, double alphahyp, int *initclasses); 
void chain_addrelation(chain chn, int relind, int ndim, double betaprop, double
	betamag, double *nig, int *domind);  
domain chain_getdomain(chain chn, int domind);
relation chain_getrelation(chain chn, int relind);


void chain_initclasscountsOLD(chain chn, int domind, int itemind); 
void chain_initclasscounts(chain chn, int domind, int itemind); 
double chain_classlikeOLD(chain chn, int domind, int itemind, int relclass); 
double chain_classlike(chain chn, int domind, int itemind, int relclass); 


void chain_additemtoclass(chain chn, int domind, int itemind, int relc);
void chain_removeitemfromclass(chain chn, int domind, int itemind);

/* add an edge to relation RELIND 
 * (updates counts for that relation, and edges for items in the associated 
 * domains)*/
void chain_addedge(chain chn, int relind, double val, int *participants);


/* probabilities */
void chain_updateprob(chain chn);
void chain_updatedomprobs(chain chn, int domind);
void chain_updaterelationprob(chain chn, int relind);
double chain_getprob(chain chn);

void chain_itergibbs(chain chn);
void chain_itersm(chain chn, chain sparechain);
int chain_climbsplit(chain chn, chain sparechn);
int chain_climbsplitfast(chain chn, chain sparechn);
int chain_climbmerge(chain chn, chain sparechn);
int chain_climbmergefast(chain chn, chain sparechn);
int chain_climbscan(chain chn);
void chain_randomize(chain chn);

void chain_hypupdate(chain chn, 
		     int (*proposal_choose) (double, double, double), 
		     int *changeflag ) ;
void chain_copy(chain c1, chain c2);
void chain_free(chain chn);

/* output */
void chain_printhyps(chain chn, FILE *fileptr);
void chain_printdetailedstatus(chain chn, FILE *fileptr);
void chain_printbriefstatus(chain chn, FILE *fileptr);
void chain_print(chain chn);
void chain_printitemclasscounts(chain chn, int domind);

#endif /* CHAIN_H */
