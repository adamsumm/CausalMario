#ifndef BLOCKUTILS_H
#define BLOCKUTILS_H 
#include "types.h"

double betaln(double r, double s);
int *my_malloc(int n);
int sample_multinomial(double *lprobs, int nelements, double temp);
int sample_best(double *score, int nelements, double empty);
int randomitem(int nitems);
double gaussrand(double mean, double sd);
int proposal_mcmc(double newscore, double oldscore, double temperature);
int proposal_best(double newscore, double oldscore, double temperature);
void logstoprobs(double *lprobs, int nelements);

void datael_assign(datael s, double val);
void datael_copy(datael s, datael t);
void datael_add(datael a, datael b);
void datael_minus(datael a, datael b);
datael datael_create(void);
void datael_free(datael d);
int datael_getint(datael d);
double datael_getdouble(datael d);
double datael_getcount(datael d); 
void datael_setcount(datael d, double c); 
double datael_getval(datael d, int i);
void datael_setval(datael d, int i, double v);

double locall(datael ss, dist prior, double nobs);
double globall(double ll, datael ss1, dist prior, double nobstot); 
int doubleeq(double a, double b);

/* this isn't in <math.h>? */
double lgamma(double x);

void dist_free(dist n);
void dist_copy(dist s, dist t);
double dist_ll(dist n1, dist n0, datael d);
double dist_logZ(dist n);
int dist_gettype(dist n);
void dist_setval(dist d, int i, double val);
double dist_getval(dist d, int i);


dist dist_create(relation r);
void dist_init(dist d, relation r);
void dist_update(dist n, datael d);


#endif /* BLOCKUTILS_H */
