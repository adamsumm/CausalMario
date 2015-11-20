#ifndef PARAMETERS_H 
#define PARAMETERS_H 
#include "config.h"

struct parameter_str
{ int loops;
  int nchains;
  double temp;
  int mcmcflag;
  int hypupdates;   /* number of times to update hyperparameters per
		       iteration*/
  double alpha;	    /* concentration parameter for CRP */
  double alphahyp;   /* parameter for exponential prior on alpha */
  int alphaupdate;
  double betaprop;  /* beta1/beta2, where beta1 and beta2 determine the prior
		       on entries in each eta matrix   */

  int betapropupdate;
  double betamag;   /* beta1 + beta2 */
  int betamagupdate;
  double m;	    /* parameters for Normal Inverse Gamma */
  double v;
  double a; 
  double b;
  int datatype;
  char outroot[MAXSTRING];
  char configfile[MAXSTRING];
  char graphname[MAXSTRING];
  char initfile[MAXSTRING];
  int outsideinit;
  int maxitem;    /* max number of items in any domain */
  int maxdim;	  /* max number of dimensions in any relation*/
  int maxrel;	  /* number of relations */
  int maxobjtuples; /* max number of objtuples for any (reln, dimension) pair*/
  int maxclass;   /* max number of classes in any domain*/
} parameter_str;

/* global variable */
extern struct parameter_str ps;

void parameters_getps(int *argcp, char ***argvp); 

#endif /* PARAMETERS_H */
