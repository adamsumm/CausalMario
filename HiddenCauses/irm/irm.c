#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "chaincolln.h"
#include "cokus.h"
#include "opt.h"
#include "parameters.h"

int main(int, char**);

int main(int argc, char **argv)
{
  chaincolln cc;

  /* seed for random number generator */
  seedMT(4357U);

  /* read parameters */
  parameters_getps(&argc, &argv);

  cc = chaincolln_readdata();

  if (ps.mcmcflag) {
    chaincolln_mcmc(cc, ps.loops); 
  } else {
    chaincolln_climb(cc, ps.loops); 
  }
  chaincolln_free(cc);
  return(0);
}
