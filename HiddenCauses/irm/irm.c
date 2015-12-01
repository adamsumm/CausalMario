#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "chaincolln.h"
#include "cokus.h"
#include "opt.h"
#include "parameters.h"
#include <time.h>

int main(int, char**);

int main(int argc, char **argv)
{
  chaincolln cc;

  srand(time(NULL));
  int r = rand();
  /* seed for random number generator */
  seedMT(r);

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
