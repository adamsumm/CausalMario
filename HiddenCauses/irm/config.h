#ifndef CONFIG_H 
#define CONFIG_H 

/* constants related to the search/sampling process */

#define RSCAN	      5		      /* number of restricted scans for
				         split-merge proposals*/
#define HYPUPDATES    5		      /* number of hyperparameter updates per 
					 gibbs or search pass */ 

/* constants used when updating hyperparameters */
#define MINALPHA       0	      /* minimum value for concentration
					 parameter */
#define MINBETAPROP    0	      /* minimum value for beta1/beta2 */
#define ALPHASD	      0.1	      /* st dev for proposals that update 
					 alpha */
#define BETASD	      0.1	      /* st dev for proposals that update 
					 betaprop */

/* miscellaneous constants */

#define MAXSTRING     300	      /* string length*/
#define MAXDOUBLE     100000000		
#define VERBOSE       0		      /* presently unused */
#define SPLIT	      1	
#define MERGE	      2		

/* data types */
#define BIN	      0		
#define FREQ 	      1		
#define CONT	      2		

#define DATAELSIZE    2		
#define DISTSIZE      4		


/* #define GSL   */	              /* compile with Gnu Scientific Library.
					 Currently used only for initializing 
					 cluster assignments in chain_colln.c */


#endif /* CONFIG_H */
