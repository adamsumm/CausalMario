#include "relation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "parameters.h"
#include "domain.h"
#include "irmutils.h"
#include "multidarray.h"

/*****************************************************************************
A relation involves one or more domains. CLASSCOUNTS represents how this
relation links classes in the participating domains relate to each other.
*****************************************************************************/

struct relation_str
{ int ndimensions; 
  domain *domains;
  int dtype;		    /* binary, frequency or continuous data*/
  multidarray classcounts;  /* 1-edges between the classes of the participating
			       domains */ 
  int missingdataflag;	    /* are any edges unobserved? */
  double likelihood;
			    /* beta1 and beta2 are hyperparameters for the beta
			     * prior on entries in the eta matrix */
  double betaprop;	    /* proportion: beta1/beta2 */
  double betamag;	    /* magnitude : beta1 + beta2*/
  double nig[DISTSIZE];     /* for continuous data*/
} relation_str;


relation relation_create(int ndimensions, double betaprop, double betamag, 
	double *nig, domain *doms, int t) {
  int i;
  int *maxsizes;
  relation rn;

  maxsizes=  (int *) my_malloc(ps.maxdim*sizeof(int));

  for (i = 0; i < ndimensions; i++) {
    maxsizes[i] = domain_getnclassmax(doms[i]);
  }
  rn = (relation) my_malloc(sizeof(struct relation_str));
  rn->dtype = t;
  rn->ndimensions = ndimensions;
  rn->domains = (domain *) my_malloc(ps.maxdim*sizeof(domain));
  rn->classcounts = multidarray_create(ndimensions, maxsizes);
  rn->missingdataflag = 0;
  rn->betaprop = betaprop;
  rn->betamag = betamag;
  for (i = 0; i < DISTSIZE; i++) {
    rn->nig[i] = nig[i];
  }
  for (i = 0; i < ndimensions; i++) {
    rn->domains[i] = doms[i];
  }

  free(maxsizes);
  return rn;
}

int relation_getdim(relation rn) {
  return rn->ndimensions;
}

void relation_getbeta(relation rn, double *b1, double *b2) {
  *b1 = rn->betaprop;
  *b2 = rn->betamag;
}

void relation_setbeta(relation rn, double b1, double b2) {
  rn->betaprop = b1;
  rn->betamag = b2;
}

domain * relation_getdoms(relation rn) {
  return rn->domains;
}

void relation_addedges(relation rn, datael weight, int *absclasses) {
  multidarray_inc(rn->classcounts, weight, absclasses);
}

void relation_removeedges(relation rn, datael weight, int *absclasses) {
  multidarray_dec(rn->classcounts, weight, absclasses);
}

datael relation_edgeweight(relation rn, int *absclasses) {
  datael w;
  w = multidarray_get(rn->classcounts, absclasses);
  return w;
}

void relation_updatelikelihood(relation rn) {
  int ndim, dtype;
  double nobs, nobstot, betaprop, betamag;
  int *classsizes, *relindex, *absindex;
  int i, j, cellnum;
  datael tallyd, celld; 
  dist priord;
  double ll;


  classsizes=  (int *) my_malloc(ps.maxdim*sizeof(int));
  relindex=    (int *) my_malloc(ps.maxdim*sizeof(int));
  absindex=    (int *) my_malloc(ps.maxdim*sizeof(int));

  ndim =  rn->ndimensions;
  dtype = rn->dtype;

  /* ckemp debug */
  if (0) {
    relation_printcounts(rn); 
  }

  cellnum = 1;
  for (i = 0; i < ndim; i++) {
    classsizes[i] = domain_getnclasses(rn->domains[i]);
    cellnum *= classsizes[i];
    relindex[i] = 0;
  }
  ll = 0; nobstot = 0; 
  tallyd = datael_create(); priord = dist_create(rn);
  dist_init(priord, rn);
  for (i = 0; i < cellnum; i++) {
    /* relindex contains relative labels */
    for (j = 0; j < ndim; j++) {
      absindex[j] = domain_absclasslabel(rn->domains[j], relindex[j]);
    }

    celld = multidarray_get(rn->classcounts, absindex);
    datael_add(tallyd, celld);

    if (rn->missingdataflag) {
      nobs = datael_getcount(celld);
    } else {		   /* stored count < nobs because of sparsity */
      nobs = 1;
      for (j = 0; j < ndim; j++) {
        /* XXX: not quite right if we rule out self links */
        nobs *= domain_getclasssize(rn->domains[j], relindex[j]);
      }
    }
    nobstot += nobs;
    ll += locall(celld, priord, nobs);
    multidarray_incrementindex(ndim, relindex, classsizes);
  }

  ll = globall(ll, tallyd, priord, nobstot);
  relation_getbeta(rn, &betaprop, &betamag);

  if (relation_getdtype(rn) == BIN) {
    ll += -2.5*log(betamag); /* prior on hyperparameters */
  } 

  rn->likelihood = ll;

  free(classsizes); free(relindex); free(absindex); free(tallyd); free(priord);
}

double relation_likelihood(relation rn) {
  return( rn->likelihood );
}

void relation_printcounts(relation rn) {
  int ndim;
  int *classsizes, *relindex, *absindex;
  int i, j, cellnum;
  double counts, v1, v2;
  datael ocd;


  classsizes =  (int *) my_malloc(ps.maxdim*sizeof(int));
  relindex   =  (int *) my_malloc(ps.maxdim*sizeof(int));
  absindex   =  (int *) my_malloc(ps.maxdim*sizeof(int));

  ndim = rn->ndimensions;
  cellnum = 1;
  for (i = 0; i < ndim; i++) {
    classsizes[i] = domain_getnclasses(rn->domains[i]);
    /* to see the entire classcounts matrix */ 
    /* classsizes[i] = domain_getnclassmax(rn->domains[i]); */
    cellnum *= classsizes[i];
    relindex[i] = 0;
  }
  
  /* relindex contains relative labels */

  for (i = 0; i < cellnum; i++) {
    for (j = 0; j < ndim; j++) {
      absindex[j] = domain_absclasslabel(rn->domains[j], relindex[j]);
    }
    ocd = multidarray_get(rn->classcounts, absindex);
    v1 = datael_getval(ocd,0);
    counts = 1;
    for (j = 0; j < ndim; j++) {
      /* XXX: not quite right if we rule out self links */
      counts *= domain_getclasssize(rn->domains[j], relindex[j]);
    }
    if (relation_getmissing(rn)) {
      counts = datael_getcount(ocd);
    } else {
    }
    if (relation_getdtype(rn) == BIN) {
      v2 = counts- v1;
    } else {
      v2 = datael_getval(ocd,1);
    }

    for (j = 0; j < ndim; j++) {
      fprintf(stdout, "%d ", relindex[j]); 
    }
    fprintf(stdout, ":  %.1f %.1f %.0f\n", v1, v2, counts); 
    multidarray_incrementindex(ndim, relindex, classsizes);
  }

  free(classsizes); free(relindex); free(absindex);
}

void relation_hypupdate(relation r, double temp,
	int (*proposal_choose) (double, double, double), int *changeflag) {
  double oldbetamag, oldbetaprop, currentscore, newscore, proposal1, proposal2; 

  currentscore = relation_likelihood(r);
  /* oldbeta1 is proportion; oldbeta2 is magnitude */
  relation_getbeta(r, &oldbetaprop, &oldbetamag);

  if (ps.betapropupdate) {
    proposal1 = exp(gaussrand(log(oldbetaprop), BETASD));
    if (proposal1 < MINBETAPROP) {
      proposal1 = MINBETAPROP + (MINBETAPROP- proposal1); 
    }
  } else {
    proposal1 = oldbetaprop;
  }

  if (ps.betamagupdate) {
    proposal2 = exp(gaussrand(log(oldbetamag), BETASD));
  } else {
    proposal2 = oldbetamag;
  }

  relation_setbeta(r, proposal1, proposal2);
  relation_updatelikelihood(r);
  newscore = relation_likelihood(r);
  if (!proposal_choose(newscore, currentscore, temp)) {
    /* proposal rejected */
    relation_setbeta(r, oldbetaprop, oldbetamag);
    relation_updatelikelihood(r);
  } else {
    (*changeflag)++;
  }
}


void relation_copy(relation rs, relation rt) {
  int i;
  rt->ndimensions = rs->ndimensions;
  multidarray_copy(rs->classcounts, rt->classcounts);
  rt->likelihood= rs->likelihood;
  rt->betaprop = rs->betaprop;
  rt->betamag = rs->betamag;
  rt->dtype= rs->dtype;
  rt->missingdataflag= rs->missingdataflag;
  for (i = 0; i < DISTSIZE; i++) {
    rt->nig[i] = rs->nig[i];
  }
}

void relation_free(relation r) {
  free(r->domains);
  multidarray_free(r->classcounts);
  free(r);
}

int relation_getmissing(relation r) {
  return r->missingdataflag;
}

void relation_setmissing(relation r, int val) {
  r->missingdataflag = val;
}

int relation_getdtype(relation r) {
  return r->dtype;
}

void relation_setdtype(relation r, int val) {
  r->dtype= val;
}

double relation_getbetaprop(relation r) {
  return r->betaprop; 
}

double relation_getbetamag(relation r) {
  return r->betamag; 
}
double relation_getnig(relation r, int i) {
  return r->nig[i];
}

