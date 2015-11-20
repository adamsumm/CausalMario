#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cokus.h"
#include "chain.h"
#include "domain.h"
#include "relation.h"
#include "irmutils.h"
#include "config.h"
#include "parameters.h"
#include "multidarray.h"

/*****************************************************************************
A chain is a collection of domains and relations. Each domain includes a set of
class assignments for its objects.
*****************************************************************************/


int chain_gibbswrap(chain chn, int (*samplefunc) (double *, int, double)); 

double chain_iteritemclasscounts(chain chn, int domind, int itemind, 
	int relclassind, 
	double (*innerfunc) 
	    (chain, const int *, int *, int, int, int, int, int) );

double chain_gibbsscan(chain chn, int domind, int *items, int nitem, 
	int *classes, int nclass, int *rigged, int restrictedflag, 
	int (*samplefunc) (double *, int, double), int *changeflag );

double chain_initinnerfunc(chain chn, const int *mdarelindex, int *mdaabsindex, 
	int domind, int relind, int dimind, int itemind, int relclassind); 

double chain_classlikeinner_addtorelcounts(chain chn, const int *mdarelindex, 
    int *mdaabsindex, int domind, int relind, int dimind, int itemind, 
    int relclassind); 

double chain_classlikeinner_takefromrelcounts(chain chn, 
    const int *mdarelindex, int *mdaabsindex, int domind, int relind, 
    int dimind, int itemind, int relclassind); 

void chain_seqscan(chain chn, int domind, int *items, int nitem, int *classes, 
		   int nclass,  int restrictedflag);  

void chain_splitclass(chain chn, int domind, int i);
void chain_mergeclass(chain chn, int domind, int relclass1, int relclass2);

void chain_randomscan(chain chn, int domind, int *items, int nitem, 
	int *classes, int nclass, int *rigged, int restrictedflag);  

int find_ritems(domain d, int i, int j, int ci, int cj, int *ritems);

void make_rigged(domain origd, int *ritems, int ritemcount, int ci, int cj, 
	int *rigged) ; 

void chain_chainswap(chain c1, chain c2, double logqratio); 

void chain_addedgedom(chain chn, domain dom, int relind, int dim, double val, 
	int itemind, int ndim, int *participants); 

void chain_adjustrelations(chain chn, domain d, int itemind, 
    void (*pf) (relation, datael, int *) ); 

/* find ITEMCLASSES corresponding to PARTICIPANTS */
int chain_itemclasses(domain *doms, int *participants, int *itemclasses, 
		int ndim);  

void chain_resethyperparameters(chain chn);

struct chain_str
{ int ndomains;
  domain *domains;
  int nrelations;
  relation *relations;
  double temperature;		    /* for MC^3 */ 
  double chainprob;		    /* overall prob of chain state */
  multidarray **itemclasscounts;    /* temporary storage for the links between
				       an item and the currently existing 
				       classes */ 

} chain_str;

chain chain_create(int ndomains, int nrelns, double temp) {
  int i;
  chain chn;

  chn  = (chain) my_malloc(sizeof(struct chain_str));
  chn->domains   = (domain *) my_malloc(ndomains*sizeof( domain));
  chn->relations = (relation*) my_malloc(nrelns*sizeof( relation));
  chn->ndomains = ndomains;
  chn->nrelations = nrelns;
  chn->temperature = temp;

  chn->itemclasscounts = (multidarray **) 
	my_malloc(nrelns*sizeof(multidarray *));
  
  for (i = 0; i < nrelns; i++) {
    chn->itemclasscounts[i]  = (multidarray *) 
	    my_malloc(ps.maxdim*sizeof(multidarray));
  }

  return chn;
}

int chain_getndomains(chain chn) {
  return chn->ndomains;
}

int chain_getnrelations(chain chn) {
  return chn->nrelations;
}

double chain_gettemp(chain chn) {
  return chn->temperature;
}

void chain_settemp(chain chn, double t) {
  chn->temperature = t;
}

/* add domain to chain */
void chain_adddomain(chain chn, int domind, int nitem, int nclassmax, 
	int clusterflag, double alpha, double alphahyp, int *initclasses) { 
  chn->domains[domind] = domain_create(nitem, nclassmax, clusterflag, alpha,
    alphahyp, domind, initclasses) ;
}

void chain_addrelation(chain chn, int relind, int ndim, double betaprop, 
	double betamag, double *nig, int *domind) { 
  int i, j, *maxsizes;
  domain *doms;

  maxsizes=  (int *) my_malloc(ps.maxdim*sizeof(int));
  doms =  (domain *) my_malloc(ps.maxdim*sizeof(domain));

  for (i= 0 ; i < ndim; i++) {
    doms[i] = chn->domains[domind[i]];  
    domain_addinstance(doms[i], relind, i);
    maxsizes[i] = domain_getnclassmax(doms[i]);
  }
  chn->relations[relind] = 
	relation_create(ndim, betaprop, betamag, nig, doms, ps.datatype);
	/* XXX: need empirical sigma here */

  for (j = 0; j < ndim; j++) {
    /* initialize and fill with zeros*/
    chn->itemclasscounts[relind][j]  = multidarray_create(ndim, maxsizes); 
  }
  free(doms); free(maxsizes);
}

domain chain_getdomain(chain chn, int domind) {
  return chn->domains[domind];
}

relation chain_getrelation(chain chn, int relind) {
  return chn->relations[relind];
}

/* careful: we don't immediately update the prior for domain DOMIND */
void chain_additemtoclass(chain chn, int domind, int itemind, int relc) {
  domain d;

  d = chn->domains[domind];
  /* adjust domain domind */
  domain_additemtoclass(d, itemind, relc);
  /* adjust relations */
  chain_adjustrelations(chn, d, itemind, &relation_addedges);
}

void chain_removeitemfromclass(chain chn, int domind, int itemind) {
  domain d;
  
  d = chn->domains[domind];
  /* adjust relations */
  chain_adjustrelations(chn, d, itemind, &relation_removeedges);

  /* adjust domain domind */
  domain_removeitemfromclass(chn->domains[domind], itemind); 
}


/* increment or decrement classcounts involving domain D */

void chain_adjustrelations(chain chn, domain d, int itemind, 
	void (*adjustfunc) (relation, datael, int *) ) {
  int i;
  int relind, dim, ndim, nedges, edgeind, ninstances, *participants, 
      *itemclasses;
  int allpresentflag;
  datael del1;
  double eweight;
  relation r;

  participants=  (int *) my_malloc(ps.maxdim*sizeof(int));
  itemclasses =  (int *) my_malloc(ps.maxdim*sizeof(int));
  del1 = datael_create(); datael_setcount(del1, 1);

  ninstances = domain_getninstance(d);
  for (i = 0; i < ninstances; i++) {
    domain_getinstance(d, i, &relind, &dim); 
    r = chn->relations[relind];
    ndim = relation_getdim(r); 
    nedges = domain_getnedges(d, relind, dim, itemind);
    for (edgeind = 0; edgeind < nedges; edgeind++) {
      eweight = domain_getedge(d, relind, dim, itemind, ndim, edgeind,
          participants); 
      datael_setval(del1, 0, eweight); datael_setval(del1, 1, pow(eweight,2));
      allpresentflag = 
           chain_itemclasses(relation_getdoms(r), participants, itemclasses,
      		       ndim); 
      if (allpresentflag == 1) { /* all participants are live */
        adjustfunc(r, del1, itemclasses);
      }
    }
  }
  free(participants); free(itemclasses); datael_free(del1);
}


/* find ITEMCLASSES (absolute) corresponding to PARTICIPANTS. Return 1 if all
** ITEMCLASSES are positive (ie the edge doesn't involve an item that is 
** temporarily out of play */

int chain_itemclasses(domain *doms, int *participants, int *itemclasses, 
		int ndim) {
  int i, ppresentflag;
  ppresentflag = 1;
  for (i = 0; i < ndim; i++) {
     itemclasses[i] = domain_getitemabsclass(doms[i], participants[i]); 
     if (itemclasses[i] < 0) {ppresentflag = 0;}
  }
  return ppresentflag;
}


/* add an edge to relation RELIND */
void chain_addedge(chain chn, int relind, double val, int *participants)
{  int i, ndim, ppresentflag;
   int *absclasses;
   domain *memberdomains;
   datael del1;

   del1 = datael_create();
   datael_setcount(del1, 1);
   datael_setval(del1, 0, val); datael_setval(del1, 1, pow(val, 2));

   absclasses =  (int *) my_malloc(ps.maxdim*sizeof(int));
   memberdomains = relation_getdoms(chn->relations[relind]);
   ndim = relation_getdim(chn->relations[relind]);

   ppresentflag = chain_itemclasses(memberdomains, participants, absclasses, 
	ndim) ;
   if ( ppresentflag != 0 ) {
     /* add edges to items in all member domains */
     for (i = 0; i < ndim; i++) {
       chain_addedgedom(chn, memberdomains[i], relind, i, val, participants[i], 
  	    ndim, participants); 
     }
     /* increment counts for relation relind */ /*XXX: change for continuous*/
     relation_addedges(chn->relations[relind], del1, absclasses);
   }

   free(absclasses); free(del1);
}

void chain_addedgedom(chain chn, domain dom, int relind, int dim, double val, 
	int itemind, int ndim, int *participants)  {
  domain_addedge(dom, relind, dim, val, itemind, ndim, participants);  
}

/* collect probs already computed for each domain and relation */
void chain_updateprob(chain chn) {
  double ll = 0;
  int dind, rind;
  for (dind = 0; dind < chn->ndomains; dind++) {
    ll += domain_getprior(chn->domains[dind]);
  }

  for (rind= 0; rind < chn->nrelations; rind++) {
    ll += relation_likelihood(chn->relations[rind]);
  }
  chn->chainprob= ll;
}

/* iterate through absindices for all instances of domain DOMIND.
   INNERFUNC says what to do with each absindex */ 

double chain_iteritemclasscounts(chain chn, int domind, int itemind, 
	int relclassind, 
	double (*innerfunc) (chain, const int *, int *, int, int, int, int, int)) {
  int i,j, ninstance, relind, dimind, rdim, celltot, ccind;
  int *domclasses, *mdarelindex, *mdaabsindex;
  double result;
  multidarray icmda;
  domain d, *rdoms;
  relation r;

  domclasses =  (int *) my_malloc(ps.maxdim*sizeof(int));
  mdarelindex=  (int *) my_malloc(ps.maxdim*sizeof(int));
  mdaabsindex=  (int *) my_malloc(ps.maxdim*sizeof(int));

  d = chain_getdomain(chn, domind);
  ninstance = domain_getninstance(d);
  result = 0;
  for (i = 0; i < ninstance; i++) {
    domain_getinstance(d, i, &relind, &dimind);
    icmda = chn->itemclasscounts[relind][dimind];
    r = chain_getrelation(chn, relind);
    rdim = relation_getdim(r);
    rdoms = relation_getdoms(r);
    for (j = 0; j < rdim; j++) {
      domclasses[j] = domain_getnclasses(rdoms[j]);     
      mdarelindex[j] = 0;
    }
    domclasses[dimind] = 1;
    celltot= 1;
    for (j = 0; j < rdim; j++) { celltot*= domclasses[j]; }
    for (ccind = 0; ccind < celltot; ccind++) {
      for (j = 0; j < rdim; j++) { /* create absolute indices */
        mdaabsindex[j] = domain_absclasslabel(rdoms[j], mdarelindex[j]);
      }
      mdaabsindex[dimind] = 0;
      /* NB: innerfunc mustn't permanently alter mdarelindex */
      result += innerfunc(chn, mdarelindex, mdaabsindex, domind, relind,
			 dimind, itemind, relclassind);
      multidarray_incrementindex(rdim, mdarelindex, domclasses);
    }
  }
  free(domclasses); free(mdarelindex); free(mdaabsindex);
  return result;
}

/* initialize itemclasscounts (number of links to the current classes) 
** for ITEMIND */
void chain_initclasscounts(chain chn, int domind, int itemind) {
  int i, j, k, nedges, relind, dimind, rdim, ninstance;
  int *participants, *pabsclasses;
  int missedgeflag;
  multidarray icmda;
  domain d, *rdoms;
  relation r;
  datael del1;
  double eweight;

  del1 = datael_create();
  datael_setcount(del1, 1);

  participants=  (int *) my_malloc(ps.maxdim*sizeof(int));
  pabsclasses =  (int *) my_malloc(ps.maxdim*sizeof(int));

  /* initialize icmda */
  chain_iteritemclasscounts(chn, domind, itemind, -1, &chain_initinnerfunc); 
  d = chain_getdomain(chn, domind);
  ninstance = domain_getninstance(d);
  for (i = 0; i < ninstance; i++) {
    domain_getinstance(d, i, &relind, &dimind);
    icmda = chn->itemclasscounts[relind][dimind];
    r = chain_getrelation(chn, relind);
    rdim = relation_getdim(r);
    rdoms = relation_getdoms(r);
    /* add to itemclasscounts matrix */
    nedges = domain_getnedges(d, relind, dimind, itemind);
    for (j = 0; j < nedges; j++) {
      eweight = domain_getedge(d, relind, dimind, itemind, rdim, j, 
      	       participants);
      datael_setval(del1, 0, eweight); datael_setval(del1, 1, pow(eweight,2));
      missedgeflag = 0;
      chain_itemclasses(rdoms, participants, pabsclasses, rdim);
      /* convention for storing current item counts
         also important for next test, since we've removed itemind already*/
      pabsclasses[dimind] = 0;
      for (k = 0; k < rdim; k++) {
        if (pabsclasses[k] == -1) {
          missedgeflag = 1; break;
        }
      }        
      if (missedgeflag == 0) {
        multidarray_inc(icmda, del1, pabsclasses);
      }
    }
  }   
  free(del1);  free(participants); free(pabsclasses);
}

double chain_initinnerfunc(chain chn, const int *mdarelindex, int *mdaabsindex, 
	int domind, int relind, int dimind, int itemind, int relclassind) {
  multidarray icmda;
  relation r;
  datael d;

  d = datael_create();

  icmda = chn->itemclasscounts[relind][dimind];
  r = chain_getrelation(chn, relind);
  multidarray_set(icmda, d, mdaabsindex);
  datael_free(d);
  return(1.0);
}


/* compute log probability that ITEMIND belongs in RELCLASS of domain DOMIND 
** We use all edges that *don't* involve ITEMIND to effectively create a new
** beta prior on each eta matrix. Then compute the probability of the ITEMIND
** edges (stored in itemclasscounts) using this "new" prior.
*/

double chain_classlike(chain chn, int domind, int itemind, int relclass) {
  int relclasssize, absclass;
  double lprior, ll;
  domain d;
  relation r;

  d = chain_getdomain(chn, domind);
  r = chain_getrelation(chn, 0);
  relclasssize = domain_getclasssize(d, relclass); 
  absclass = domain_absclasslabel(d, relclass);
  if (relclasssize == 0) {
    lprior = log(domain_getalpha(d));
    /* add phantom class so chain_iteritemclasscounts can iterate through it*/
    domain_addemptyclass(d); 
  } else {
    lprior = log(( double) relclasssize);
  }


  /* ckemp debug */
  if (0) {
  relation_printcounts(chn->relations[0]);
  }
  
  /* compute likelihood of data for new assignment. update relcounts along the
   * way */
  ll = chain_iteritemclasscounts(chn, domind, itemind, relclass, 
	&chain_classlikeinner_addtorelcounts); 
  
  /* ckemp debug */
  if (0) {
    fprintf(stderr, "----------------\n");
    relation_printcounts(chn->relations[0]);
  }

  /* undo changes to relcounts */
  chain_iteritemclasscounts(chn, domind, itemind, relclass, 
	&chain_classlikeinner_takefromrelcounts); 

  /* ckemp debug */
  if (0) {
    fprintf(stderr, "----------------\n");
    relation_printcounts(chn->relations[0]);
  }

  if (relclasssize == 0) {
    /* remove phantom class */ 
    domain_removeemptyclass(d); 
  }

  return(lprior + ll);
}


/* compute likelihood of data for new assignment. update relcounts along the
   * way */

double chain_classlikeinner_addtorelcounts(chain chn, const int *mdarelindex, 
    int *mdaabsindex, int domind, int relind, int dimind, int itemind, 
    int relclass) {

  int j, rdim, absclass, relclasssize;
  double nobsold, nobsnew, cellcountold, cellcountnew;
  datael obsoldtmp, obsoldd, obsnewd;
  dist priord;
  int *classsizes, *mdarelindexcopy;
  int origmdareldim;
  domain d, *rdoms;
  relation r;
  multidarray icmda;
  double ll;

  classsizes = (int *) my_malloc(ps.maxdim*sizeof(int));
  mdarelindexcopy =  (int *) my_malloc(ps.maxdim*sizeof(int));
  for (j = 0; j < ps.maxdim; j++) {
    mdarelindexcopy[j] = mdarelindex[j];
  }
  
  /* we don't allow self links of any sort (no edge tuple can involve the same
   * object twice). This holds for 1 tuples (if no data are missing) and for
   * 1-tuples and 0-tuples (if we have missing data) */ 

  d = chain_getdomain(chn, domind);
  absclass = domain_absclasslabel(d, relclass);
  relclasssize = domain_getclasssize(d, relclass);
  r = chain_getrelation(chn, relind);
  rdim = relation_getdim(r);
  rdoms = relation_getdoms(r);

  icmda = chn->itemclasscounts[relind][dimind];

  origmdareldim = mdarelindexcopy[dimind];
  /* dimind set to 0 in both mdarelindex, mdaabsindex */
  mdaabsindex[dimind] = absclass;
  mdarelindexcopy[dimind] = relclass;
  priord = dist_create(r); 
  obsoldtmp= datael_create(); 

  cellcountold = 1;
  for (j = 0; j < rdim; j++) { /* classsizes*/
    classsizes[j] = domain_getclasssize(rdoms[j], mdarelindexcopy[j]);
    /* don't want to double count cells involving a new class, but don't want
     * to ignore them altogether. Our solution is perhaps not elegant : for all 
     * but one instance of each of these cells, we set new1c == new0c == 0 */
    if (classsizes[j] == 0 && j <= dimind) {
      classsizes[j] = 1;
    /* need to adjust classsizes in case we've already added itemind to
     * some dimensions */
    } else if (domain_getlabel(rdoms[j]) == domind && 
        mdarelindexcopy[j] == relclass && j < dimind ) {
      classsizes[j]++;
    }
    cellcountold *= classsizes[j]; /* will be zero for all but one instance of 
				      each cell involving a new class */
  }

  /* existing counts for this cell */
  obsoldd = relation_edgeweight(r, mdaabsindex);
  datael_copy(obsoldd, obsoldtmp);
  if (relation_getmissing(r)) {
    nobsold = datael_getcount(obsoldd);
  } else {
    if ( relclasssize == 0 ) {
      nobsold = 0;
    } else {
      nobsold = cellcountold ; 
    }
  }
  datael_setcount(obsoldtmp, nobsold);
  cellcountnew = cellcountold / classsizes[dimind]; 

  /* new counts for this cell */
  mdaabsindex[dimind] = 0;
  obsnewd = multidarray_get(icmda, mdaabsindex);
  /* ckemp debug*/
  if (0) {
  chain_printitemclasscounts(chn, domind);
  }

  if (relation_getmissing(r)) {
    nobsnew = datael_getcount(obsnewd);
  } else {
    nobsnew = cellcountnew;
  }

  dist_init(priord, r);  /* restore original prior */
  dist_update(priord, obsoldtmp);  
  ll = locall(obsnewd, priord, nobsnew); 

  /* increment oldcounts for this cell */
  mdaabsindex[dimind] = absclass;

  
  if (0) { /* ckemp debug */
  fprintf(stderr, "----------------\n");
  relation_printcounts(chn->relations[0]);
  }

  relation_addedges(r, obsnewd, mdaabsindex); 

  if (0) {
  fprintf(stderr, "----------------\n");
  relation_printcounts(chn->relations[0]);
  }

  free(classsizes); free(mdarelindexcopy);  free(priord); free(obsoldtmp);
  return ll;
}

double chain_classlikeinner_takefromrelcounts(chain chn, const int *mdarelindex, 
    int *mdaabsindex, int domind, int relind, int dimind, int itemind, 
    int relclass) {

  int absclass;
  datael new1cd;
  multidarray icmda;
  relation r;
  domain d;

  icmda = chn->itemclasscounts[relind][dimind];
  r = chain_getrelation(chn, relind);
  d = chain_getdomain(chn, domind);
  absclass = domain_absclasslabel(d, relclass);

  /* new 1 counts for this cell */
  mdaabsindex[dimind] = 0;
  new1cd = multidarray_get(icmda, mdaabsindex);

  mdaabsindex[dimind] = absclass;
  relation_removeedges(r, new1cd, mdaabsindex); 
  /** test will be violated for continuous data with negative edge weights
  if (datael_getval(relation_edgeweight(r, mdaabsindex),0) < 0) {
    fprintf(stderr, "negative edges"); exit(1);
  }
  **/
  return (1.0);
}

/* update probs involving domain dom */

void chain_updatedomprobs(chain chn, int domind) {
  domain dom;
  int *relflag;
  int i, ninstance, relind, dimind;

  relflag=  (int *) my_malloc(ps.maxrel*sizeof(int));
  dom = chn->domains[domind];
  domain_updateprior(dom);
  for (i = 0; i < chn->nrelations; i++) {
    relflag[i] = 0;
  }
  ninstance = domain_getninstance(dom);
  for (i = 0; i < ninstance; i++) {
    domain_getinstance(dom, i, &relind, &dimind);
    if (relflag[relind] == 0) {
      relation_updatelikelihood(chn->relations[relind]);
      relflag[relind] = 1;
    }
  }
  chain_updateprob(chn);
  free(relflag);
}

/* get chain prob */

double chain_getprob(chain chn) {
  return chn->chainprob;
}


/* find members of classes ci and cj (other than i and j) */
int find_ritems(domain d, int i, int j, int ci, int cj, int *ritems) {
  int cicount, cjcount, ritemind, k, cmember;
  cicount = domain_getclasssize(d, ci);
  ritemind = 0;
  for (k = 0; k < cicount; k++) {
    cmember = domain_getclassmember(d, ci, k);  
    if (cmember != i) {
      ritems[ritemind] = cmember; ritemind++;
    }
  }	
  cjcount = domain_getclasssize(d, cj);
  for (k = 0; k < cjcount; k++) {
    cmember = domain_getclassmember(d, cj, k);  
    if (cmember != j) {
      ritems[ritemind] = cmember; ritemind++;
    }
  }
  return ritemind; 
}


void chain_chainswap(chain schn, chain chn, double logqratio) {
  double aswap, rn;

  aswap = exp( logqratio ) * 
	pow(exp( chain_getprob( schn ) - chain_getprob( chn ) ), 
		chain_gettemp(schn) ) ;
  if (aswap > 1) { aswap = 1; } 
  rn = myrand();
  if (rn < aswap) {
    chain_copy(schn, chn);
    fprintf(stderr, "****SM: %f\n", aswap);
  }
}


/* find assignments in ORIGD for items in RITEMS*/
void make_rigged(domain origd, int *ritems, int ritemcount, int ci, int cj, 
	int *rigged) {
  int k;
  int origclass;

  /* set up rigged array */
  for (k = 0; k < ritemcount; k++) {
    origclass= domain_getitemrelclass(origd, ritems[k]); 
    if (origclass== ci) {
      rigged[k] = 1;
    } else {
      rigged[k] = 0;
      if (origclass != cj) {fprintf(stderr, "rigged problem"); exit(1); }
    }
  }
}

void chain_hypupdate(chain chn, 
	int (*proposal_choose) (double, double, double), int *changeflag) {
  int ucount, domind, relind;
  double oldalpha;
  double currentscore, newscore; 
  double proposal;
  domain d;

  for (ucount = 0; ucount < ps.hypupdates; ucount++) {
    if (ps.alphaupdate) { /* update alphas? */
    for (domind = 0; domind < chn->ndomains; domind++) {
      d = chain_getdomain(chn, domind);
      if ( domain_getclusterflag(d)  ) {
        currentscore = chain_getprob(chn);
        oldalpha = domain_getalpha(d); 
        proposal = exp(gaussrand(log(oldalpha), ALPHASD));
	if (proposal < MINALPHA) { /* minimum value for alpha is active */
          proposal = MINALPHA + (MINALPHA - proposal); 
	}
        domain_setalpha(d, proposal);
        chain_updatedomprobs(chn, domind);
        newscore= chain_getprob(chn);
        if (!proposal_choose(newscore, currentscore, chn->temperature) ) {
          /* proposal rejected */
          domain_setalpha(d, oldalpha);
          chain_updatedomprobs(chn, domind);
        } else {
          (*changeflag)++;
        }
      }
    }    
    }

    for (relind = 0; relind < chn->nrelations; relind++) {
      relation_hypupdate(chn->relations[relind], chn->temperature, 
			 proposal_choose, changeflag); 
    }
    chain_updateprob(chn);
  }
}

/* split-merge proposals (see Jain and Neal) */
void chain_itersm(chain chn, chain schn) {
  int domind, nitem, nclasses, i, j, ci, cj, smtype, ritemcount, 
	rclasscount, rscanind, ndomains, changeflag;
  int *ritems, *rigged;
  int rclasses[2];
  double scanprob, logqratio;

  domain d, origd;

  ritems  =  my_malloc(ps.maxitem*sizeof(int));
  rigged  =  my_malloc(ps.maxitem*sizeof(int));

  changeflag = 0;
  ndomains = chn->ndomains;
  for (domind = 0 ; domind < ndomains; domind++) {
    chain_copy(chn, schn);
    d = chain_getdomain(schn, domind);
    if (domain_getclusterflag(d)) {
      nitem = domain_getnitems(d);
      nclasses = domain_getnclasses(d);
      if (nitem == 1) {continue;}
      i = randomitem(nitem); j = i;
      while (j == i) {j = randomitem(nitem); }
      ci = domain_getitemrelclass(d, i); 
      cj = domain_getitemrelclass(d, j); 
      smtype = MERGE;
      if (ci == cj) {
        smtype = SPLIT;
	chain_removeitemfromclass(schn, domind, i);
	/* add i to a class of its own */
	chain_additemtoclass(schn, domind, i, nclasses); 
	ci = nclasses;
      }
      ritemcount= find_ritems(d, i, j, ci, cj, ritems);
      rclasses[0] = cj; rclasses[1] = ci; rclasscount = 2;
      /* initialize randomly */
      chain_randomscan(schn, domind, ritems, ritemcount, rclasses, rclasscount,
	    NULL, 1);

      for (rscanind = 0; rscanind < RSCAN; rscanind++) {
        scanprob = chain_gibbsscan(schn, domind, ritems, ritemcount, rclasses, 
		rclasscount, NULL, 1, &sample_multinomial, &changeflag); 
      }

      if (smtype == SPLIT) {
        logqratio = - chain_gibbsscan(schn, domind, ritems, ritemcount, 
		rclasses, rclasscount, NULL, 1, &sample_multinomial,
		&changeflag); 
      } else {
        origd = chn->domains[domind];
        make_rigged(origd, ritems, ritemcount, ci, cj, rigged);
        logqratio = chain_gibbsscan(schn, domind, ritems, ritemcount, 
		rclasses, rclasscount, rigged, 1, &sample_multinomial,
		&changeflag); 
	/* assign ritems to cj (the first entry in rclasses */
        chain_randomscan(schn, domind, ritems, ritemcount, rclasses, 1,
	    NULL, 1);
        chain_removeitemfromclass(schn, domind, i);
	/* cj might have changed */
        cj = domain_getitemrelclass(d, j); 
        chain_additemtoclass(schn, domind, i, cj);      
        chain_updatedomprobs(schn, domind);
      }
      chain_chainswap(schn, chn, logqratio);
    }
  }
  free(ritems); free(rigged);
}

/* hill-climbing pass: for each item, try splitting the class it belongs to.
 * This will be slow for big datasets */

int chain_climbsplit(chain chn, chain schn) {
  int domind, i, ndomains, nitems, changeflag;
  domain d;

  changeflag = 0;
  ndomains = chn->ndomains;
  for (domind = 0 ; domind < ndomains; domind++) {
    d = chain_getdomain(schn, domind);
    nitems = domain_getnitems(d);
    if (domain_getclusterflag(d)) {
      for (i = 0; i < nitems; i++) {
        chain_copy(chn, schn);
        chain_splitclass(schn, domind, i);
	if ( chain_getprob(schn) > chain_getprob(chn) ) {
          chain_copy(schn, chn);
  	  changeflag++;
        }
      }
    }
  }
  return changeflag;
}

/* hill-climbing pass: try splitting each class (choose bigger classes with
   higher probability. This scales better than the previous function for big
   datasets  */

int chain_climbsplitfast(chain chn, chain schn) {
  int domind, i,ci,itemind, itemindind, ndomains, changeflag, nclasses;
  double *logclasssizes;
  domain d;

  logclasssizes=  (double *) my_malloc(ps.maxclass*sizeof(double));
  changeflag = 0;
  ndomains = chn->ndomains;
  for (domind = 0 ; domind < ndomains; domind++) {
    d = chain_getdomain(chn, domind);
    if (domain_getclusterflag(d)) {
      nclasses = domain_getnclasses(d);
      for (i = 0; i < nclasses; i++) {
        logclasssizes[i] = log( (double) domain_getclasssize(d, i));
      }
      for (i = 0; i < nclasses; i++) {
        chain_copy(chn, schn);
        ci = sample_multinomial(logclasssizes, nclasses, 1.0);
        itemindind = randomitem(domain_getclasssize(d, ci));
        itemind = domain_getclassmember(d, ci, itemindind);
        chain_splitclass(schn, domind, itemind);
	if ( chain_getprob(schn) > chain_getprob(chn) ) {
          chain_copy(schn, chn);
  	  changeflag++;
        }
      }
    }
  }
  free(logclasssizes);
  return changeflag;
}


/* split the class that I belongs to */
void chain_splitclass(chain chn, int domind, int i) {
  int relci, relcisize, j, jind, relcj, nclasses, ritemcount,
    rclasscount;
  int *ritems;
  int rclasses[2];
  domain d;

  ritems= (int *) my_malloc(ps.maxitem * sizeof(int));

  d = chain_getdomain(chn, domind);
  nclasses = domain_getnclasses(d);
  relci = domain_getitemrelclass(d, i);
  relcisize = domain_getclasssize(d, relci);
  if (relcisize == 1) { free(ritems); return; }
  j = i;
  while (j == i) { /* find another member of ci */
    jind = randomitem(relcisize); 
    j = domain_getclassmember(d, relci, jind);
  }
  chain_removeitemfromclass(chn, domind, j);
  chain_additemtoclass(chn, domind, j, nclasses); 
  relcj = nclasses; 
  ritemcount= find_ritems(d, i, j, relci, relcj, ritems);
  rclasses[0] = relcj; rclasses[1] = relci; rclasscount = 2;
  chain_seqscan(chn, domind, ritems, ritemcount, rclasses, rclasscount,
  	      1);
  free(ritems);
}

/* hill-climbing pass: try merging all pairs of classes in each domain. Doesn't
 * scale well if there are many classes */

int chain_climbmerge(chain chn, chain schn) {
  int domind, nclasses, ci, cj, ndomains, nitems, changeflag;
  domain d;

  changeflag = 0;
  ndomains = chn->ndomains;
  for (domind = 0 ; domind < ndomains; domind++) {
    d = chain_getdomain(schn, domind);
    nitems = domain_getnitems(d);
    if (domain_getclusterflag(d)) {
      nclasses = domain_getnclasses(d);
      for (ci = 0; ci < nclasses; ci++) {
        for (cj = 0; cj < ci; cj++) {
          chain_copy(chn, schn);
          chain_mergeclass(schn, domind, ci, cj);
          if ( chain_getprob(schn) > chain_getprob(chn) ) {
            chain_copy(schn, chn);
            d = chain_getdomain(schn, domind);
	    nclasses = domain_getnclasses(d);
	    changeflag++;
	  }
	}
      }
    }
  }  
  return changeflag;
} 

/* hill-climbing pass: try merging each class with one other. Scales 0K */

int chain_climbmergefast(chain chn, chain schn) {
  int domind, nclasses, ci, cj, ndomains, changeflag;
  domain d;

  changeflag = 0;
  ndomains = chn->ndomains;
  for (domind = 0 ; domind < ndomains; domind++) {
    d = chain_getdomain(chn, domind);
    if (domain_getclusterflag(d)) {
      nclasses = domain_getnclasses(d);
      if (nclasses == 1) { continue; }
      for (ci = 0; ci < nclasses; ci++) {
        cj = ci;
	while (cj == ci) {cj = randomitem(nclasses); }
        chain_copy(chn, schn);
        chain_mergeclass(schn, domind, ci, cj);
        if ( chain_getprob(schn) > chain_getprob(chn) ) {
          chain_copy(schn, chn);
          d = chain_getdomain(chn, domind);
	  nclasses = domain_getnclasses(d);
	  if (nclasses == 1) {break;}
	  changeflag++;
	}
      }
    }
  }  
  return changeflag;
}


/* put all the members of RELC1 into RELC2 */ 

void chain_mergeclass(chain chn, int domind, int relc1, int relc2) {
  int relc1size, itemind, i; 
  domain d;
  
  d = chain_getdomain(chn, domind);
  relc1size = domain_getclasssize(d, relc1);
  for (i = 0; i < relc1size; i++) {
    itemind = domain_getclassmember(d, relc1, i);
    chain_removeitemfromclass(chn, domind, itemind);
    chain_additemtoclass(chn, domind, itemind, relc2) ;
  }
  chain_updatedomprobs(chn, domind);
}


void chain_itergibbs(chain chn) {
  chain_gibbswrap(chn, &sample_multinomial); 
}

int chain_climbscan(chain chn) {
  int changeflag;
  changeflag = chain_gibbswrap(chn, &sample_best); 
  return changeflag;
}

/* randomize class assignments if the hill-climbing search has reached a local
 * maximum  */
void chain_randomize(chain chn) {
  int i, ndomains, domind, nitem, nclass, itemind, initclass;
  int *items, *classes;

  items= (int *) my_malloc(ps.maxitem * sizeof(int));
  classes= (int *) my_malloc(ps.maxclass * sizeof(int));

  ndomains = chn->ndomains;
  if (ps.outsideinit) {
    for (domind = 0 ; domind < ndomains; domind++) {
      nitem = domain_getnitems(chn->domains[domind]);
      for (itemind = 0; itemind < nitem; itemind++) {
        chain_removeitemfromclass(chn, domind, itemind);
      }
      for (itemind = 0; itemind < nitem; itemind++) {
        initclass = domain_getiteminitclass(chn->domains[domind], itemind);
        chain_additemtoclass(chn, domind, itemind, initclass);
      }
      chain_updatedomprobs(chn, domind);
    }
  } else {
    /* classes[] indicates classes to consider */
    for (i = 0; i < ps.maxclass; i++) {
      classes[i] = i;      
    }
    /* items[] indicates items to consider */
    for (i = 0; i < ps.maxitem; i++) {
      items[i] = i;      
    }
    for (domind = 0 ; domind < ndomains; domind++) {
      if (domain_getclusterflag(chn->domains[domind])) {
        nclass = domain_getnclasses(chn->domains[domind]);
        nitem = domain_getnitems(chn->domains[domind]);
	/* reinitialize with one class per domain */
        chain_randomscan(chn, domind, items, nitem, classes, 1, NULL, 1); 
      }
      chain_updatedomprobs(chn, domind);
    }
  }
  chain_resethyperparameters(chn);

  free(items); free(classes);
}

/* set hyperparameters to original values */

void chain_resethyperparameters(chain chn) {
  int i;

  for (i = 0; i < chn->ndomains; i++) {
    domain_setalpha(chn->domains[i], ps.alpha);
  }

  for (i = 0; i < chn->nrelations; i++) {
    relation_setbeta(chn->relations[i], ps.betaprop, ps.betamag);
  } 
}

/* Update class assignments using a gibbs-style pass. SAMPLEFUNC is function
 * for choosing new class assignments: can choose probabilistically
 * (sample_multinomial) or maximize (sample_best) */

int chain_gibbswrap(chain chn, int (*samplefunc) (double *, int, double)) {
  int i, ndomains, domind, nitem, changeflag; 
  int *items, *classes;

  items= (int *) my_malloc(ps.maxitem * sizeof(int));
  classes= (int *) my_malloc(ps.maxclass * sizeof(int));

  changeflag = 0;
  ndomains = chn->ndomains;
  /* classes[] indicates classes to consider */
  for (i = 0; i < ps.maxclass; i++) {
    classes[i] = i;      
  }
  /* items[] indicates items to consider */
  for (i = 0; i < ps.maxitem; i++) {
    items[i] = i;      
  }

  for (domind = 0 ; domind < ndomains; domind++) {
    if (domain_getclusterflag(chn->domains[domind])) {
      nitem = domain_getnitems(chn->domains[domind]);
      chain_gibbsscan(chn, domind, items, nitem, classes, -1, NULL, 0,
		samplefunc, &changeflag); 
    }
  }

  free(items); free(classes);
  return changeflag;
}

/* consider putting ITEMS into CLASSES. If RIGGED is not NULL it specifies
 * which class assignments we'll choose for each item */

/* RIGGED specifies indices into CLASSES for each itemind in ITEMS */
/* CLASSES contains relative class labels */

double chain_gibbsscan(chain chn, int domind, int *items, int nitem, 
	int *classes, int nclass, int *rigged, int restrictedflag, 
	int (*samplefunc) (double *, int, double), int *changeflag) {

  double scanprob, *score, *checkscore;
  int itemind, itemindind, oldclass, oldclasssize,  relclass, relclassind, 
      samplerelclassind, newclass, newclasssize;
  domain d;

  score= (double *) my_malloc(ps.maxclass* sizeof(double));
  checkscore = (double *) my_malloc(ps.maxclass * sizeof(double));

  scanprob = 0;
  d = chain_getdomain(chn, domind);

  for (itemindind = 0; itemindind < nitem; itemindind++) {
    itemind = items[itemindind];
    oldclass = domain_getitemrelclass(d, itemind);
    oldclasssize = domain_getclasssize(d, oldclass);
    /* ckemp debug */
    if (0) {
      relation_printcounts(chain_getrelation(chn, 0));
      fprintf(stderr, "---------------\n");
    }

    chain_removeitemfromclass(chn, domind, itemind);
    if (restrictedflag == 0) { /* nclass might now be smaller than argument */
      /* +1 so we consider starting new class */
      nclass = domain_getnclasses(chn->domains[domind])+1;
    }

    /*add empty class so we initialize as much of itemclasscounts as we'lluse */
    domain_addemptyclass(d); 
    chain_initclasscounts(chn, domind, itemind);
    domain_removeemptyclass(d); 

    for (relclassind = 0; relclassind < nclass; relclassind++) {
      relclass = classes[relclassind];
      if (0) {
        checkscore[relclassind] = chain_classlike(chn, domind, itemind, relclass);
      }
      /* It's conceptually simpler but slower to add the item to relclass and
         compute the probability of the resulting configuration directly. This
	 can be used to check chain_classlike() */
      if (1) {
        chain_additemtoclass(chn, domind, itemind, relclass);      
        chain_updatedomprobs(chn, domind);
        score[relclassind] = chain_getprob(chn);
        chain_removeitemfromclass(chn, domind, itemind);
        chain_updatedomprobs(chn, domind);
      }
    }

    samplerelclassind= samplefunc(score, nclass, chn->temperature); 
    if (rigged != NULL) {
      samplerelclassind = rigged[itemindind]; 
    }
    newclass = classes[samplerelclassind];
    scanprob += log(score[samplerelclassind]);

    /* ckemp debug */
    if (0) {
      relation_printcounts(chain_getrelation(chn, 0));
      fprintf(stderr, "---------------\n");
    }

    chain_additemtoclass(chn, domind, itemind, newclass);      
    newclasssize = domain_getclasssize(d, newclass);

    /* ckemp debug */
    if (0) {
      fprintf(stderr, "---------------\n");
      relation_printcounts(chain_getrelation(chn, 0));
    }

    /* has state changed? */
    if (newclass != oldclass && !(oldclasssize == 1 && newclasssize == 1) ) {
      *changeflag = *changeflag + 1;
    }
  }
  chain_updatedomprobs(chn, domind);
  free(score); free(checkscore);
  
  return scanprob;
}


/* random scan involving ITEMS: choose new classes randomly from CLASSES*/
void chain_randomscan(chain chn, int domind, int *items, int nitem, 
	int *classes, int nclass, int *rigged, int restrictedflag)  {

  int itemind, itemindind,  newrelclass, oldrelclass, oldrelclasssize, 
	newrelclassind;
  domain d;

  d = chn->domains[domind];

  for (itemindind = 0; itemindind < nitem; itemindind++) {
    itemind = items[itemindind];
    oldrelclass= domain_getitemrelclass(d, itemind);
    oldrelclasssize = domain_getclasssize(d, oldrelclass);
    if (restrictedflag == 0) { /* nclass might now be smaller than argument */
      if (oldrelclasssize == 1) {
        nclass = domain_getnclasses(chn->domains[domind]);
      } else {
	    /* +1 if we want to consider starting new class */
        nclass = domain_getnclasses(chn->domains[domind]);
      }
    }
    newrelclassind = randomitem(nclass) ;
    newrelclass = classes[newrelclassind];
    if (newrelclass != oldrelclass) {
      chain_removeitemfromclass(chn, domind, itemind);
      chain_additemtoclass(chn, domind, itemind, newrelclass);      
    }
  }
  chain_updatedomprobs(chn, domind);
}

/* sequential scan: remove ITEMS from current classes. Make one pass through,
   choosing new classes from NCLASSES */
void chain_seqscan(chain chn, int domind,int *items,int nitem, int *classes, 
		   int nclass,  int restrictedflag)  {
  int itemind, itemindind,  relclassind, relclass, newrelclassind, newrelclass;
  double *score;
  domain d;

  score = (double *) my_malloc(ps.maxclass* sizeof(double));

  d = chain_getdomain(chn, domind);

  /* remove ITEMS from classes */
  for (itemindind = 0; itemindind < nitem; itemindind++) {
    itemind = items[itemindind];
    chain_removeitemfromclass(chn, domind, itemind);
  }

  for (itemindind = 0; itemindind < nitem; itemindind++) {
    itemind = items[itemindind];
    chain_initclasscounts(chn, domind, itemind);

    for (relclassind = 0; relclassind < nclass; relclassind++) {
      relclass = classes[relclassind];
      score[relclassind] = chain_classlike(chn, domind, itemind, relclass);
    }
    newrelclassind = sample_best(score, nclass, 1.0); 
    newrelclass = classes[newrelclassind];
    chain_additemtoclass(chn, domind, itemind, newrelclass);      
  }
  chain_updatedomprobs(chn, domind);

  free(score);
}

/* policy: only copy things that need to be copied (ie that have some chance of
 * being different between cs and ct */

void chain_copy(chain cs, chain ct) {
  int i;

  ct->ndomains = cs->ndomains;  
  ct->nrelations= cs->nrelations;  
  ct->temperature= cs->temperature;  
  ct->chainprob= cs->chainprob;  

  for (i = 0; i < ct->ndomains; i++) {
    domain_copy(cs->domains[i], ct->domains[i]);
  }

  for (i = 0; i < ct->nrelations; i++) {
    relation_copy(cs->relations[i], ct->relations[i]);
  } 
}


void chain_free(chain chn) {
  int i, j;
 
  for (i = 0; i < chn->nrelations; i++) {
    for (j = 0; j < relation_getdim(chn->relations[i]); j++) {
      multidarray_free(chn->itemclasscounts[i][j]);
    }
    free(chn->itemclasscounts[i]);
  }
  free(chn->itemclasscounts);
 
  for (i = 0; i < chn->ndomains; i++) {
    domain_free(chn->domains[i]);
  }
  free(chn->domains); 
 
  for (i = 0; i < chn->nrelations; i++) {
    relation_free(chn->relations[i]);
  }
  free(chn->relations); 
  
  free(chn);
}

void chain_printhyps(chain chn, FILE *fileptr) {
  int i;
  relation r;
  domain d;
  double beta1, beta2;

   /* rel hyps*/
  for (i = 0; i < chain_getnrelations(chn); i++) {
    r = chain_getrelation(chn, i); 
    relation_getbeta(r, &beta1, &beta2); 
    fprintf(fileptr, "%7.3f %7.3f ", beta1, beta2); 
  }

  /* dom hyps*/
  for (i = 0; i < chain_getndomains(chn); i++) {
    d = chain_getdomain(chn, i); 
    fprintf(fileptr, "%7.3f ", domain_getalpha(d)); 
  }
}

void chain_printdetailedstatus(chain chn, FILE *fileptr) {
  int i; 
  relation r;
  domain d;

  /* relscores */
  for (i = 0; i < chain_getnrelations(chn); i++) {
    r = chain_getrelation(chn, i); 
    fprintf(fileptr, "%7.3f ", relation_likelihood(r)); 
  }

  /* domscores */
  for (i = 0; i < chain_getndomains(chn); i++) {
    d = chain_getdomain(chn, i); 
    fprintf(fileptr, "%7.3f ", domain_getprior(d)); 
  }

  chain_printbriefstatus(chn, fileptr);
}

void chain_printbriefstatus(chain chn, FILE *fileptr) {
  int i;
  domain d;

  for (i = 0; i < chain_getndomains(chn); i++) {
    d = chain_getdomain(chn, i); 
    fprintf(fileptr, "%4d ", domain_getnclasses(d)); 
  }
 fprintf(fileptr, "%7.3f ", chain_getprob(chn)); 
}


void chain_print(chain chn) {
  int i;
  for (i = 0; i < chn->ndomains; i++) {
    fprintf(stdout, "Domain %d\n", i);
    domain_printclasses(chn->domains[i]);
  }
  for (i = 0; i < chn->nrelations ; i++) {
    fprintf(stdout, "Relation %d\n", i);
    relation_printcounts(chn->relations[i]);
  }
}

void chain_printitemclasscounts(chain chn, int domind) {
  int i,j, ninstance, relind, dimind, rdim, celltot, ccind;
  double cnt;
  datael cntd;
  int *domclasses, *mdarelindex, *mdaabsindex;
  multidarray icmda;
  domain d, *rdoms;
  relation r;

  domclasses=  (int *) my_malloc(ps.maxdim*sizeof(int));
  mdarelindex=  (int *) my_malloc(ps.maxdim*sizeof(int));
  mdaabsindex=  (int *) my_malloc(ps.maxdim*sizeof(int));

  d = chain_getdomain(chn, domind);
  ninstance = domain_getninstance(d);
  for (i = 0; i < ninstance; i++) {
    domain_getinstance(d, i, &relind, &dimind);
    fprintf(stdout, "domain %d: relation %d: dimension %d\n", domind, relind,
		dimind); 
    icmda = chn->itemclasscounts[relind][dimind];
    r = chain_getrelation(chn, relind);
    rdim = relation_getdim(r);
    rdoms = relation_getdoms(r);
    for (j = 0; j < rdim; j++) {
      domclasses[j] = domain_getnclasses(rdoms[j]);     
      mdarelindex[j] = 0;
    }
    domclasses[dimind] = 1;
    celltot= 1;
    for (j = 0; j < rdim; j++) { celltot*= domclasses[j]; }
    /* initialize itemclasscounts matrix to 0*/
    for (ccind = 0; ccind < celltot; ccind++) {
      for (j = 0; j < rdim; j++) { /* create absolute indices */
        mdaabsindex[j] = domain_absclasslabel(rdoms[j], mdarelindex[j]);
	fprintf(stdout, "%d ", mdarelindex[j]);
      }
      mdaabsindex[dimind] = 0;
      cntd = multidarray_get(icmda, mdaabsindex);
      cnt = datael_getcount(cntd);
      fprintf(stdout, ": %0.f\n", cnt);
      multidarray_incrementindex(rdim, mdarelindex, domclasses);
    }
  }

  free(domclasses); free(mdarelindex); free(mdaabsindex); 
}
