#ifndef PROC3OPTHEADER
#define PROC3OPTHEADER

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "tspmain.h"
#include "heap.h"

/******************************
CONSTANTS
******************************/

#define SCHEME_3_1   301  /* scheme r_1 in paper */
#define SCHEME_3_2   302  /* scheme r_2 in paper */
#define SCHEME_3_3   303  /* scheme r_3 in paper */
#define SCHEME_3_4   304  /* scheme r_4 in paper */
#define SCHEME_3_5   305  /* scheme r_5 in paper */
#define SCHEME_3_6   306  /* scheme r_6 in paper */
#define SCHEME_3_7   307  /* scheme r_7 in paper */
#define SCHEME_3_8   308  /* scheme r_8 in paper */
#define SCHEME_3_9   309  /* scheme r_9 in paper */


/*****************************************
   VARIABLES
 *****************************************/

extern int MAX_ALLOWED_COMB; // used to estimate BF moves time. 
                             // stop the enumeration after these many, and scale the time 
			     // in proportion to the # of combinations 


/*****************************************
   PROCEDURES
 *****************************************/

double estimateTime3OPTMoveBF( int* pi, int* schemelist, int totschemes );
void perform3OPTMove( int* pi, int* sel, int scheme );
long int numCompleteTrips( int N );
/* double find3OPTMoveBFOverall( int* pi, int* schlist, int totsch, int* bestsel, int* bestsch, 
                              long int* eval, int* nopt, long int* nEffective);
   double find3OPTMoveSFOverallSingleHeap( int* pi, int* schemes, int totschemes, int* bestsel, int* bscheme,
					int* heapSize, int* extract, double* timeHeapConstruction,
					double* timeHeapUsage, long int *eval, long int* evalGood );
*/
void find3optMoveSF( int *schemes, int nschemes, int* currentTour, OPTMOVE* move);

#endif
