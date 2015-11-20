#include "parameters.h"
#include "config.h"
#include "string.h"
#include "opt.h"

/*****************************************************************************
The global variable ps specifies parameters for the current run.
*****************************************************************************/

/* global variable */
struct parameter_str ps;

void parameters_getps(int *argcp, char ***argvp) {
  char emptystring[5];
  emptystring[0] = '\0';
  /* set defaults */
  ps.loops    = 10;
  ps.nchains  = 1;
  ps.temp     = 1;
  ps.mcmcflag = 0;
  ps.hypupdates = 0;  /* number of times to update hyperparameters.
		         by default, leave hyperparameters fixed.
		 	 ps.hypupdates = 5 seems to work OK */
  ps.alpha    = 1;
  ps.alphahyp = 1;
  ps.alphaupdate = 0; /* update alpha hyperparameter? */
  ps.betaprop = 1;    
  ps.betapropupdate= 0; /* update betaprop? */
  ps.betamag = 0.2; /* Beta(0.1, 0.1) prior by default */
  ps.betamagupdate= 0; /* update betamag? */

  ps.m = 0;
  ps.v = 5;
  ps.a = 2;
  ps.b = 0.3; 
  ps.datatype = BIN;

  strcpy(ps.outroot, emptystring);
  strcpy(ps.configfile, emptystring);
  strcpy(ps.graphname, emptystring);
  strcpy(ps.initfile, emptystring);

  optrega(&(ps.loops),	   OPT_INT,	'\0',  "loops", "number of loops");
  optrega(&(ps.nchains),   OPT_INT,	'\0',  "nchains", "number of chains");
  optrega(&(ps.temp),	   OPT_DOUBLE,  '\0',  "temp", "temperature for MC^3");
  optrega(&(ps.mcmcflag),  OPT_INT,	'\0',  "mcmcflag", "run MCMC?");
  optrega(&(ps.hypupdates),  OPT_INT,	'\0',  "hypupdates", "update hyperparameters?");
  optrega(&(ps.alpha),	   OPT_DOUBLE,  '\0',  "alpha", 
			   "concentration hyperparameter for the CRP");
  optrega(&(ps.alphahyp),  OPT_DOUBLE,  '\0',  "alphahyp", 
			   "prior on concentration hyperparameter for the CRP");
  optrega(&(ps.alphaupdate), OPT_INT,	'\0',  "alphaupdate", "change alpha?");
  optrega(&(ps.betaprop),  OPT_DOUBLE,  '\0',  "betaprop", 
  "ratio of the beta hyperparameters (from prior on entries in the eta matrix");
  optrega(&(ps.betapropupdate), OPT_INT,'\0',  "betapropupdate", "change betaprop?");
  optrega(&(ps.betamag),  OPT_DOUBLE,  '\0',  "betamag", 
			    "sum of the beta hyperparameters");
  optrega(&(ps.betamagupdate), OPT_INT,'\0',  "betamagupdate", "change betamag?");
  optrega(&(ps.m),  OPT_DOUBLE,  '\0',  "m", 
			    "mean for prior on mu (continuous data)");
  optrega(&(ps.v),  OPT_DOUBLE,  '\0',  "v", 
			    "mu has variance v*sigma^2(continuous data)");
  optrega(&(ps.a),  OPT_DOUBLE,  '\0',  "a", 
			    "sigma^2 distributed IG(a,b) (continuous data)");
  optrega(&(ps.b),  OPT_DOUBLE,  '\0',  "b", 
			    "sigma^2 distributed IG(a,b) (continuous data)");
  optrega(&(ps.datatype),  OPT_INT,  '\0',  "datatype", 
			    "binary, frequency or continuous data");
  optrega(&(ps.outroot),   OPT_CSTRING, '\0',  "outroot", 
			    "basename for output files");
  optrega(&(ps.configfile),OPT_CSTRING, '\0',  "configfile",
			    "file with domain and relation specifications ");
  optrega(&(ps.graphname), OPT_CSTRING, '\0',  "graphname",
			    "relation file");
  optrega(&(ps.initfile),  OPT_CSTRING, '\0',  "initfile", 
			    "file with initial class assignments");
  opt(argcp, argvp);
  opt_free();

  if (!ps.initfile[0]) {
    ps.outsideinit = 0;
  } else {
    ps.outsideinit = 1;
  }

  return;
}

