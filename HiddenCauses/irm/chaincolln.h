#ifndef CHAINCOLLN_H 
#define CHAINCOLLN_H 
#include "chain.h"
#include "types.h"


chaincolln chaincolln_create(int nchains, int ndomains, int nrelns, 
	char *prefix);

chaincolln chaincolln_readdata(void);
void chaincolln_mcmc(chaincolln cc, int nloops);
void chaincolln_climb(chaincolln cc, int nloops);

char * chaincolln_getprefix(chaincolln cc);  

int chaincolln_getnchain(chaincolln cc);
chain chaincolln_getchain(chaincolln cc, int chainind);

int chaincolln_getitercount(chaincolln cc);

void chaincolln_print(chaincolln cc);
void chaincolln_sample(chaincolln cc, int niter);
void chaincolln_itergibbsm(chaincolln cc);
int chaincolln_hypupdates(chaincolln cc,
	int (*proposal_choose) (double, double, double) ) ;

void chaincolln_tryswaps(chaincolln cc);
void chaincolln_swapchains(chaincolln cc, int sw1, int sw2);
void chaincolln_free(chaincolln cc);

void chaincolln_printassignments(chaincolln cc);
void chaincolln_printstatus(chaincolln cc, FILE *fptr);

#endif /* CHAINCOLLN_H */
