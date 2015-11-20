#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include "config.h"
#include "cokus.h"
#include "irmutils.h"
#include "relation.h"

double normrand2( void ); 

struct datael_str
{ 
  short tflag; 
  double count;
  struct value {
    double dn[DATAELSIZE]; /* double number*/
  } v;
} datael_str;

/* distribution */
struct dist_str
{ 
  int dtype;
  double th[DISTSIZE];
} dist_str;

dist dist_createnorn(void);


/******************************************************************************/
/* betaln(): log beta function 
*/
double betaln(double r, double s)
{
  double result;

  result = lgamma(r)+lgamma(s)-lgamma(r+s);
  return result;
}


int *my_malloc(int chunksize)
{ int *ptr;
  ptr = malloc(chunksize);
  if (ptr == NULL) {
    fprintf(stderr, "ERROR: malloc failed");
  }
  return ptr;
}

void logstoprobs(double *lprobs, int nelements) {
  int i;
  double logmax, totprob;

  logmax = lprobs[0];
  for (i = 1; i < nelements ; i++) {
    if (lprobs[i] > logmax) { 
      logmax = lprobs[i];
    }
  }
  totprob = 0;
  for (i = 0; i < nelements; i++) {
    lprobs[i] = exp(lprobs[i] - logmax); 
    totprob += lprobs[i];
  }

  /* XXX: necessary? */
  for (i = 0; i < nelements; i++) {
    lprobs[i] = lprobs[i]/totprob; 
  }
}



/* sample from multinomial according to log probabilities in LPROBS. Convert
 * LPROBS to probabilities and raise them to the power TEMP in the process */

int sample_multinomial(double *lprobs, int nelements, double temp) {
  int i, chooseind;
  double logmax, totprob, rn, sofar;

  logmax = lprobs[0];
  for (i = 1; i < nelements ; i++) {
    if (lprobs[i] > logmax) { 
      logmax = lprobs[i];
    }
  }
  totprob = 0;
  for (i = 0; i < nelements; i++) {
    lprobs[i] = pow( exp(lprobs[i] - logmax), temp ); 
    totprob += lprobs[i];
  }

  /* XXX: necessary? */
  for (i = 0; i < nelements; i++) {
    lprobs[i] = lprobs[i]/totprob; 
  }

  rn = myrand();
  sofar = lprobs[0];
  chooseind = 0;
  while (sofar < rn) {
    chooseind++;
    sofar += lprobs[chooseind];
  }
  /*fprintf(stderr, "rand: %f\n", rn);*/
  return chooseind;
}

int sample_best(double *score, int nscores, double empty) {
  int i, bestind;
  double bestscore;

  bestind = 0; bestscore = score[0];
  for (i = 1; i < nscores; i++) {
    /* add noise to deal with ties*/
    score[i] = score[i]-0.0001*myrand();
    if (score[i] > bestscore) { 
      bestscore = score[i];
      bestind = i;
    }
  }
  return bestind;
}

int randomitem(int nitems) {
  return (int) (myrand()*nitems) % nitems;
}

double gaussrand( double mean, double std) {
  double rn, proposal;
  
  rn = normrand2();
  proposal = mean + std * rn;
  /* if (proposal <= 0)  {
    proposal = -proposal;
  } */

  return(proposal);
}


double normrand2( void ) {
  double x1,x2,n1,n2,rad;

  /* generate two normal variates */
  rad = 2;
  while (rad >= 1) {
    x1 = myrand()*2-1;
    x2 = myrand()*2-1;
    rad = x1*x1+x2*x2;
  }
  n1 = x1*sqrt(-2*log(rad)/rad);
  n2 = x2*sqrt(-2*log(rad)/rad);

  return(n1);
}

int proposal_mcmc(double newscore, double oldscore, double temp) {
  double randsamp;

  randsamp = myrand();
  if (randsamp < pow( exp(newscore-oldscore), temp) ) {
    return 1;
  } else {
    return 0;
  }
}

int proposal_best(double newscore, double oldscore, double temp) {
  if (newscore > oldscore) {
    return 1;
  } else {
    return 0;
 }
}

datael datael_create(void) {
  datael d;
  int i;
  d = (datael) my_malloc(sizeof(struct datael_str));
  d->count = 0;
  for (i = 0; i < DATAELSIZE; i++) {
    d->v.dn[i] = 0.0;
  }
  return(d);
}

void datael_free(datael d) {
  free(d);
}


/* t = s */
void datael_copy(datael s, datael t)  {
  int i;
  t->count = s->count;
  for (i = 0; i < DATAELSIZE; i++) {
    t->v.dn[i] = s->v.dn[i];
  }
}

/* a = a + b */
void datael_add(datael a, datael b)  {
  int i;
  a->count += b->count;
  for (i = 0; i < DATAELSIZE; i++) {
    a->v.dn[i] += b->v.dn[i];
  }
}

/* a = a - b */
void datael_minus(datael a, datael b)  {
  int i;
  a->count -= b->count;
  for (i = 0; i < DATAELSIZE; i++) {
    a->v.dn[i] -= b->v.dn[i];
  }
}

void datael_setcount(datael d, double c) {
  d->count = c;
}


double datael_getcount(datael d) {
  return d->count;
}

double datael_getval(datael d, int i) {
  if (i >= DATAELSIZE || i < 0) {
    fprintf(stderr, "unexpected index\n"); exit(1); 
  }
  return (d->v.dn[i]);
}

void datael_setval(datael d, int i, double val) {
  if (i >= DATAELSIZE || i < 0) {
    fprintf(stderr, "unexpected index\n"); exit(1); 
  }
  d->v.dn[i] = val;
}



double locall(datael ss, dist prior, double nobs) {
  int dtype;
  double ll;
  double onecounts, zerocounts, beta1, beta2; /* for binary data */
  double counts, pseudocounts;		      /* for frequency data */
  dist posttmp;

  dtype = dist_gettype(prior);
  if (dtype == BIN) {
    onecounts  = ss->v.dn[0]; 
    zerocounts = nobs - onecounts;
    beta1 = prior->th[0]; beta2 = prior->th[1];
    ll = betaln(onecounts + beta1, zerocounts + beta2) 
	    - betaln(beta1, beta2);
  } else if (dtype == FREQ) {
    counts = ss->v.dn[0]; pseudocounts = nobs*prior->th[0];
    if (doubleeq(nobs, 0)) {
      ll = 0.0;
    } else {
      ll = lgamma(counts + pseudocounts) - lgamma(pseudocounts) -
	   counts*log(nobs); 
    }
  } else if (dtype == CONT) {
    posttmp = dist_createnorn();
    dist_copy(prior, posttmp);
    dist_update(posttmp, ss); 
    ll = dist_ll(posttmp, prior, ss);
    dist_free(posttmp);
  } else {
    /* unknown type*/
    fprintf(stderr, "unexpected type\n"); exit(1); 
  }
  return ll;
}

double globall(double ll, datael ss1, dist prior, double nobstot) {
  double countstot, pcsingle;		      /* for frequency data */
  int dtype;

  dtype = dist_gettype(prior);
  if (dtype == BIN) {
  } else if (dtype == FREQ) {
    countstot = ss1->v.dn[0]; pcsingle = prior->th[0];
    ll += lgamma(pcsingle * nobstot) - lgamma(countstot + nobstot*pcsingle);
  } else if (dtype == CONT) {
  } else {
    /* unknown type*/
    fprintf(stderr, "unexpected type\n"); exit(1); 
  }

  return ll;
}

int doubleeq(double a, double b) {
  if (fabs(a-b) < 0.000000001) {
    return 1;
  } else {
    return 0;
  }
}

void dist_init(dist d, relation rn)
{
  double betamag, betaprop, beta1, beta2;
  if (d->dtype == BIN) {
    betaprop = relation_getbetaprop(rn); betamag = relation_getbetamag(rn);
    beta2 = betamag / (betaprop+1);
    beta1 = betaprop * beta2;
    dist_setval(d, 0, beta1); dist_setval(d, 1, beta2); 
  } else if (d->dtype == FREQ) {
    betaprop = relation_getbetaprop(rn); betamag = relation_getbetamag(rn);
    beta2 = betamag / (betaprop+1);
    beta1 = betaprop * beta2;
    dist_setval(d, 0, beta1); 
  } else if (d->dtype == CONT) {
    /* we use the same parameterization as Denison (via Murphy's slides) */
    d->th[0] = relation_getnig(rn,0);   /*m: mean of mu is m */
    d->th[1] = relation_getnig(rn,1);   /*v: sd of mu is v*sigma^2 */
    d->th[2] = relation_getnig(rn,2);   /*a: sigma^2 distributed according to IG(a,b) */
    d->th[3] = relation_getnig(rn,3)  ; /*b*/
		    /* expected value of sigma^2 is b/(a-1) */
  } else {
    fprintf(stderr, "unexpected type\n"); exit(1); 
  }
}

dist dist_create(relation rn)
{ dist n;
  int i, type;
  type = relation_getdtype(rn);
  n = (dist) my_malloc(sizeof(dist_str));
  for (i = 0; i < DISTSIZE; i++) {
    n->th[i] = 0;
  }
  n->dtype = type;
  dist_init(n, rn);
  return n;
}

dist dist_createnorn(void)
{ dist n;
  n = (dist) my_malloc(sizeof(dist_str));
  return n;
}


void dist_free(dist n)
{
  free(n);
}

void dist_copy(dist s, dist t)
{ 
  int i = 0;
  t->dtype = s->dtype;
  for (i = 0; i < DISTSIZE; i++) {
    t->th[i] = s->th[i];
  }
}

void dist_update(dist n, datael d) 
{  dist ntmp;
   double sx, sx2, count;
   ntmp = dist_createnorn();
   dist_copy(n, ntmp);
   count = datael_getcount(d);

   if (n->dtype == BIN) {
     n->th[0] += datael_getval(d,0);    
     n->th[1] += count - datael_getval(d,0);
   } else if (n->dtype == CONT) {
     sx = datael_getval(d, 0);
     sx2 = datael_getval(d, 1);

     n->th[0] = (1/(1/ntmp->th[1] + count))*(ntmp->th[0]/ntmp->th[1] + sx);
     n->th[1] = 1/(1/ntmp->th[1] + count);
     n->th[2] = ntmp->th[2] + count/2;
     n->th[3] = ntmp->th[3] + 0.5 *(pow(ntmp->th[0],2)/ntmp->th[1] + sx2 
				- pow(n->th[0],2)/n->th[1]);
  }

  dist_free(ntmp);
}

double dist_logZ(dist n)
{ double lZ;
  lZ = 0.0;

  if (n->dtype == CONT) {
    lZ = 0.5*log(2*M_PI) + 0.5*n->th[1] 
       + lgamma(n->th[2]) - n->th[2]*log(n->th[3]);
  }
  return (lZ);
}

double dist_ll(dist n1, dist n0, datael d) 
{
  double n, ll;
  n = datael_getcount(d);
  ll = 0.0;

  if (n0->dtype == CONT) {
    ll = 0.5*log(n1->th[1]) + n0->th[2]*log(n0->th[3]) +lgamma(n1->th[2]) 
	- 0.5*log(n0->th[1]) - n1->th[2]*log(n1->th[3]) - lgamma(n0->th[2]) 
	- n/2*log(M_PI);
  }

  return( ll );
}

int dist_gettype(dist n) {
  return n->dtype;
}

void dist_setval(dist d, int i, double val) {
  if (i >= DISTSIZE || i < 0) {
    fprintf(stderr, "unexpected index\n"); exit(1); 
  }
  d->th[i] = val;
}

double dist_getval(dist d, int i) {
  if (i >= DATAELSIZE || i < 0) {
    fprintf(stderr, "unexpected index\n"); exit(1); 
  }
  return (d->th[i]);
}


