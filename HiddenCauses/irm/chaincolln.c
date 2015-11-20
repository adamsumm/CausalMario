#include <stdio.h>
#include <string.h>
#include <math.h>
#include "chaincolln.h"
#include "cokus.h"
#include "chain.h"
#include "config.h"
#include "parameters.h"
#include "irmutils.h"
#include "relation.h"
#include "domain.h"

#ifdef GSL
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>
#endif


/*****************************************************************************
A collection of chains. We maintain several chains in case we want to run
Metropolis-coupled MCMC (MC^3). We also keep one spare chain so that
split-merge proposals can be tried without losing the current state. 
*****************************************************************************/

struct chaincolln_str
{ int nchains;
  chain *chains;
  chain sparechain;
  char prefix[MAXSTRING];	/* unused right now */
  int itercount;		/* iteration count */
} chaincolln_str;

/* nchains already includes one extra for sparechain */
chaincolln chaincolln_create(int nchains, int ndomains, int nrelns, char *pref)
{
  int i;
  chaincolln cc;
  double temp;
  cc  = (chaincolln) my_malloc(sizeof(struct chaincolln_str));
  cc->chains = (chain *) my_malloc(nchains*sizeof(chain ));
  /* include spare chain */
  for (i = 0; i < nchains; i++) {
    temp = 1.0/(1+ps.temp*i);
    cc->chains[i] = chain_create(ndomains, nrelns, temp);
  }
  cc->nchains = nchains-1;
  cc->sparechain = cc->chains[nchains-1];
  strcpy(cc->prefix, pref);
  cc->itercount = 0;
  return cc;
}

/* read the configuration file and the graph */
chaincolln chaincolln_readdata(void) {
  FILE *fileptr, *initzsfile;
  int i, j, k, ndom, nreln, d, r, nitem, dim, maxclass, initclass, relcl, ndim, 
	domlabel, clusterflag, itemind, nchains, cind, zind;
  int *domlabels, *participants, participant;
  double val;
  double nig[DISTSIZE];
  domain *doms;
  relation rn;
  int *initclasses, ***edgecounts, *relsizes;
  char prefix[MAXSTRING];

  chaincolln cc;
  chain c, c0;
#ifdef GSL
  gsl_rng *rng;
  const gsl_rng_type *T;
  gsl_permutation *perm ;
  size_t N;

  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);
#endif 

  nchains = ps.nchains+1;
  nig[0] = ps.m; nig[1] = ps.v; nig[2] = ps.a; nig[3] = ps.b; 
  
  fileptr = fopen(ps.configfile,"r");
  if (fileptr == NULL) {
    fprintf(stderr, "couldn't read config file\n"); exit(1); 
  }

  /* initial read of ps.configfile to get ps.maxdim, ps.maxrel, ps.maxitem, 
     ps.maxclass */
  fscanf(fileptr, "%s", prefix);
  fscanf(fileptr, "%d %d", &ndom, &nreln);
  relsizes=  (int *) my_malloc(nreln*sizeof(int));
  ps.maxrel = nreln;
  ps.maxitem = 0; ps.maxclass = 0;
  for (d = 0; d < ndom; d++) {
    fscanf(fileptr, "%d %d %d %d", &nitem, &maxclass, &initclass, &clusterflag);
    if (nitem > ps.maxitem) {
      ps.maxitem = nitem;
    }
    if (maxclass > ps.maxclass) {
      ps.maxclass= maxclass;
    }
  }
  ps.maxdim = 0;
  for (r = 0; r < nreln; r++) {
    fscanf(fileptr, "%d", &ndim);
    relsizes[r] = ndim;
    if (ndim > ps.maxdim) {
      ps.maxdim = ndim;
    }
    for (dim=0; dim < ndim; dim++) {
      fscanf(fileptr, "%d", &domlabel);
    }
  }
  fclose(fileptr);

  domlabels=	 (int *) my_malloc(ps.maxdim*sizeof(int));
  participants=  (int *) my_malloc(ps.maxdim*sizeof(int));
  initclasses =  (int *) my_malloc(ps.maxitem*sizeof(int));

  /* initial read of ps.graphname to get ps.maxobjtuples */
  edgecounts =  (int ***) my_malloc(ps.maxrel*sizeof(int **));
  for (i = 0; i < ps.maxrel; i++) {
    edgecounts[i] =  (int **) my_malloc(ps.maxdim*sizeof(int *));
    for (j = 0; j < ps.maxdim; j++) {
      edgecounts[i][j] =  (int *) my_malloc(ps.maxitem*sizeof(int));
      for (k = 0; k < ps.maxitem; k++) {
        edgecounts[i][j][k] = 0;
      }
    }
  }
  ps.maxobjtuples = 0;

  fileptr = fopen(ps.graphname,"r");
  if (fileptr == NULL) {
    fprintf(stderr, "couldn't read graph\n"); exit(1); 
  }
  while( fscanf( fileptr, " %d", &r)!=EOF ) {
    ndim = relsizes[r];
    for (dim = 0; dim < ndim; dim++) {
      fscanf(fileptr, "%d", &participant);
      participants[dim] = participant;
    }
    fscanf(fileptr, "%lf", &val); 

    for (dim = 0; dim < ndim; dim++) {
        edgecounts[r][dim][participants[dim]]++;
    }
  }
  fclose(fileptr);
  for (i = 0; i < ps.maxrel; i++) {
    for (j = 0; j < ps.maxdim; j++) {
      for (k = 0; k < ps.maxitem; k++) {
        if (edgecounts[i][j][k] > ps.maxobjtuples) {
          ps.maxobjtuples = edgecounts[i][j][k];
        }
        edgecounts[i][j][k]= 0;
      }
    }
  }

  free(relsizes); 
  for (i = 0; i < ps.maxrel; i++) {
    for (j = 0; j < ps.maxdim; j++) {
      free(edgecounts[i][j]);
    }
    free(edgecounts[i]);
  }
  free(edgecounts);


  /* second read of ps.configfile where we set up datastructures */

  fileptr = fopen(ps.configfile,"r");
  if (ps.outsideinit) {
    initzsfile= fopen(ps.initfile,"r");
    if (initzsfile == NULL) {
      fprintf(stderr, "couldn't read initzsfile\n"); exit(1); 
    }
  } else {
    initzsfile = NULL;
  }

  fscanf(fileptr, "%s", prefix);
  fscanf(fileptr, "%d %d", &ndom, &nreln);

  cc = chaincolln_create(nchains, ndom, nreln, prefix);
  c0 = chaincolln_getchain(cc, 0);

  /* read domains */
  /* input file: nitem maxclass initclass clusterflag*/
  for (d = 0; d < ndom; d++) {
    fscanf(fileptr, "%d %d %d %d", &nitem, &maxclass, &initclass, &clusterflag);
#ifdef GSL
    N = nitem; 
#endif
    if (ps.outsideinit) {
      for (zind = 0; zind < nitem; zind++) {
        fscanf(initzsfile, "%d", &initclasses[zind]);
      }
    }

    /* add domains and items to chains */
    for (cind = 0; cind < nchains; cind++) {
      c = chaincolln_getchain(cc, cind);
      chain_adddomain(c, d, nitem, maxclass, clusterflag, ps.alpha,
		      ps.alphahyp, initclasses);
#ifdef GSL
      perm =  gsl_permutation_alloc(N);
      gsl_permutation_init(perm);
      gsl_ran_shuffle(rng, perm->data, N, sizeof(size_t)); 
#endif
      /* assign items to classes */
      relcl = 0;
      for (i = 0; i < nitem; i++) {
        if (ps.outsideinit) {
	  chain_additemtoclass(c, d, i, initclasses[i]);
	} else { 
          if (relcl == initclass) relcl = 0; 

	  /* without the GNUSL, each chain gets initialized the same way. This
	   * is suboptimal */
	  itemind = i;
#ifdef GSL
          itemind = gsl_permutation_get(perm, i);
#endif
          chain_additemtoclass(c, d, itemind, relcl);
          relcl++;
        }
      }
#ifdef GSL
      gsl_permutation_free(perm);
#endif
    }
  }
#ifdef GSL
  gsl_rng_free(rng);
#endif
  
  /* read relations*/
  /* input file: ndim d0 ... dn */

  for (r = 0; r < nreln; r++) {
    fscanf(fileptr, "%d", &ndim);
    for (dim=0; dim < ndim; dim++) {
      fscanf(fileptr, "%d", &domlabel);
      domlabels[dim] = domlabel;
    }
    for (cind = 0; cind < nchains; cind++) {
      c = chaincolln_getchain(cc, cind);
      chain_addrelation(c, r, ndim, ps.betaprop, ps.betamag, nig, domlabels);
    }
  }
  if (ps.outsideinit) {
    fclose(initzsfile);    
  }

  fclose(fileptr);
  /* second read of ps.graphname: store edges*/
  fileptr = fopen(ps.graphname,"r");
  /* input file: relind p0 p1 p2 .. pn val */
  while( fscanf( fileptr, " %d", &r)!= EOF ) {
    ndim = relation_getdim( chain_getrelation(c0, r) );
    doms = relation_getdoms( chain_getrelation(c0, r) ); 
    for (dim = 0; dim < ndim; dim++) {
      fscanf(fileptr, "%d", &participant);
      participants[dim] = participant;
      domlabels[dim] = domain_getlabel(doms[dim]); 
    }
    
    for (i = 0; i < ndim; i++) {
      for (j = 0; j < i; j++) {
        if (participants[i] == participants[j] && 
	    domlabels[i] == domlabels[j]) {
	  fprintf(stderr, "Self links not allowed.\n"); exit(1);  
	}
      } 
    } 

    fscanf(fileptr, "%lf", &val);
    for (cind = 0; cind < nchains; cind++) {
      c = chaincolln_getchain(cc, cind);
      chain_addedge(c, r, val, participants); 
      rn = chain_getrelation(c, r);
      if (doubleeq(val, 0)) {
	relation_setmissing(rn, 1);	
      }
      if (val > 1.5 && relation_getdtype(rn) != CONT) {
	relation_setdtype(rn, FREQ);	
      }
      if (!doubleeq(val, (int) val)) {
	relation_setdtype(rn, CONT);	
	relation_setmissing(rn, 1); /* XXX: no sparse continuous matrices */	
      }	
    }
  }

  fclose(fileptr);

  for (cind = 0; cind < nchains; cind++) {
    c = chaincolln_getchain(cc, cind);
    for (i = 0; i < chain_getndomains(c); i++) {
      chain_updatedomprobs(c, i);
    }
  }

  free(domlabels); free(participants); free(initclasses);

  return cc;
}

int chaincolln_getnchain(chaincolln cc) {
  return cc->nchains;
}

chain chaincolln_getchain(chaincolln cc, int chainind) {
  return cc->chains[chainind];
}

char *chaincolln_getprefix(chaincolln cc) {
  return cc->prefix;
}

int chaincolln_getitercount(chaincolln cc) {
  return cc->itercount;
}

void chaincolln_print(chaincolln cc) {
  chaincolln_printassignments(cc);
  chaincolln_printstatus(cc, NULL);
  chaincolln_printstatus(cc, stdout);
}

/* run hillclimbing with random restarts */
void chaincolln_climb(chaincolln cc, int maxloops) {
  int i, j, changeflag, nrepeats; 
  chain chn, schn;
  double oldscore, newscore, currscore;

  chn = chaincolln_getchain(cc, 0);
  schn = cc->sparechain;
  changeflag = 1; 
  nrepeats = 0;
  oldscore = -MAXDOUBLE;
  currscore= -MAXDOUBLE;

  chaincolln_print(cc);
  for (i = 0; i < maxloops; i++) {
    if (changeflag > 0) { 
      /* maximizing scan: move each object to the best class for it */
      changeflag = chain_climbscan(chn);
    }
    /*chain_print(chn); fprintf(stdout, "splitting\n");*/
    changeflag = changeflag + chain_climbsplitfast(chn, schn);

    /*chain_print(chn); fprintf(stdout, "merging\n");*/
    changeflag = changeflag + chain_climbmergefast(chn, schn);

    for (j = 0; j < ps.hypupdates; j++) {
      changeflag = changeflag + chaincolln_hypupdates(cc, &proposal_best);
    }

    cc->itercount++;
    newscore = chain_getprob(chn);
    if ( fabs(newscore - oldscore) < 0.00001 * fabs(newscore) ) {
      nrepeats++;
      if (nrepeats == 8) {
        chaincolln_print(cc);
        fprintf(stderr, "*************** randomize\n");
        nrepeats = 0;
	chain_randomize(cc->chains[0]);
        currscore = chain_getprob(chn);
        newscore = chain_getprob(chn);
      }
    } else {
      nrepeats = 0;
    }
    oldscore = newscore;
    chaincolln_printstatus(cc, stdout);
  }
}

/* run gibbs sampler with Metropolis updates (MC^3, split-merge) */
void chaincolln_mcmc(chaincolln cc, int maxloops) {
  int i, j;
  chaincolln_print(cc);
  for (i = 0; i < maxloops; i++) {
    /* gibbs scans and split merge updates */
    chaincolln_sample(cc, 1);
    /* MC^3: try swaps between chains at different temperatures */
    chaincolln_tryswaps(cc);
    /* hyperparameter updates? */
    for (j = 0; j < ps.hypupdates; j++) {
      chaincolln_hypupdates(cc, &proposal_mcmc);
    }
    chaincolln_print(cc);
  }
}


/* update hyperparameters */
int chaincolln_hypupdates(chaincolln cc, 
	int (*proposal_choose) (double, double, double)) {
  int i, changeflag;
  changeflag = 0;
  for (i = 0; i < cc->nchains; i++) {
    chain_hypupdate(chaincolln_getchain(cc, i), proposal_choose, &changeflag);
  }
  return changeflag;
}


void chaincolln_itergibbsm(chaincolln cc) {
  int i;
  if (1) { /* regular gibbs scans */
    for (i = 0; i < cc->nchains; i++) {
      chain_itergibbs(chaincolln_getchain(cc, i));
    }
  }
  if (1) { /* split-merge updates */
    for (i = 0; i < cc->nchains; i++) {
      chain_itersm(chaincolln_getchain(cc, i), cc->sparechain);
    }
  }
  cc->itercount++;
}

void chaincolln_sample(chaincolln cc, int niter) {
  int i;
  for (i = 0; i < niter; i++) {
    chaincolln_itergibbsm(cc);
  }
}

/* try swapping chains at different temperatures */
void chaincolln_tryswaps(chaincolln cc) {
  int nchains;
  int sw1, sw2, swcount;
  chain c1, c2; 
  double t1, t2, aswap, rn;

  nchains = cc->nchains;
  if (nchains > 1) {
    for (swcount = 0; swcount < 1; swcount++) {
      sw1 = randomitem(nchains);
      sw2 = sw1;
      while (sw2 == sw1) { sw2 = randomitem(nchains); }

      c1 = cc->chains[sw1]; c2 = cc->chains[sw2]; 
      t1 = chain_gettemp(c1); t2 = chain_gettemp(c2); 
      aswap = exp( (t1 - t2)* (chain_getprob(c2) - chain_getprob(c1)) );
      if (aswap > 1) { aswap=1; }
      rn = myrand();
      if (rn < aswap) {
        chaincolln_swapchains(cc, sw1, sw2);
        fprintf(stdout, "****CS: %f\n", aswap);
      }
    }
  }
}

/* swap chains SW1 and SW2 */
void chaincolln_swapchains(chaincolln cc, int sw1, int sw2) {
  chain c1, c2;
  double t1, t2;
  c1 = chaincolln_getchain(cc, sw1);
  c2 = chaincolln_getchain(cc, sw2);
  t1 = chain_gettemp(c1); t2 = chain_gettemp(c2); 
  chain_settemp(c1, t2); chain_settemp(c2, t1); 
  cc->chains[sw1] = c2;
  cc->chains[sw2] = c1;
}


void chaincolln_printassignments(chaincolln cc) {
  int i;
  chain chn; 
  char filename[MAXSTRING];
  FILE * domfileptr;

  chn = chaincolln_getchain(cc, 0);
  for (i = 0; i < chain_getndomains(chn); i++) {
    sprintf(filename, "%s_dom%d", ps.outroot, i);
    domfileptr = fopen(filename, "a");  
    if (domfileptr == NULL) {
      fprintf(stderr, "couldn't read domain file\n"); 
      exit(1); 
    }
    domain_printassignments(domfileptr, chain_getdomain(chn, i));
    fclose(domfileptr);
  }
}


void chaincolln_printstatus(chaincolln cc, FILE *fptr) {
  int i, nchains;
  chain chn; 
  char filename[MAXSTRING];
  FILE * domfileptr;

  nchains = chaincolln_getnchain(cc);
  if (fptr == NULL) {
    sprintf(filename, "%s_status", ps.outroot);
    domfileptr = fopen(filename, "a");  
    if (domfileptr == NULL) {
      fprintf(stderr, "couldn't read statusfile\n"); 
      exit(1); 
    }
  } else { 
    domfileptr = fptr;
  }
  fprintf(domfileptr, "%4d: ", chaincolln_getitercount(cc));
  chn = chaincolln_getchain(cc, 0);
  chain_printdetailedstatus(chn, domfileptr);
  for (i = 1; i < nchains; i++) {
    chn = chaincolln_getchain(cc, i);
    chain_printbriefstatus(chn, domfileptr);
  }
  chn = chaincolln_getchain(cc,0);
  chain_printhyps(chn, domfileptr);
  fprintf(domfileptr, "\n");


  if (fptr == NULL) {
    fclose(domfileptr);
  }
}

void chaincolln_free(chaincolln cc)
{
  int i;
  /* include spare chain */
  for (i = 0; i < cc->nchains+1; i++) {
    chain_free(cc->chains[i]);
  }
  free(cc->chains);
  free(cc);
}

