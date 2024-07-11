#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tspmain.h"
#include "io4opt.h"
#include "proc3opt.h"
#include "heap.h"

#define SCHEME_3_1   301  /* scheme r_1 in paper */
#define SCHEME_3_2   302  /* scheme r_2 in paper */
#define SCHEME_3_3   303  /* scheme r_3 in paper */
#define SCHEME_3_4   304  /* scheme r_4 in paper */
#define SCHEME_3_5   305  /* scheme r_5 in paper */
#define SCHEME_3_6   306  /* scheme r_6 in paper */
#define SCHEME_3_7   307  /* scheme r_7 in paper */
#define SCHEME_3_8   308  /* scheme r_8 in paper */
#define SCHEME_3_9   309  /* scheme r_9 in paper */


/***************************************

CODE STARTS HERE

****************************************/

void perform3OPTMove( int* pi, int* sel, int scheme ) {

  int ii, jj, kk, hh;

  ii = sel[0];
  jj = sel[1];
  kk = sel[2];
  hh = -1;
  
  switch ( scheme ) {

  case SCHEME_3_1: // cut-and-paste <+1,+3,+2>
    reverse( pi, ii + 1, kk );                   // reverse (2,3) -> +1 -3 -2
    reverse( pi, ii + 1, ii + (kk - jj) );       // reverse -3    -> +1 +3 -2
    reverse( pi, ii + (kk - jj) + 1, kk );       // reverse -2    -> +1 +3 +2
    break;

  case SCHEME_3_2:
    reverse( pi, ii + 1, jj );                 // reverse (2)   -> +1 -2 +3
    reverse( pi, jj + 1, kk );                 // reverse (3)   -> +1 -2 -3
    break;

  case SCHEME_3_3:
    reverse( pi, ii + 1, kk );                 // reverse (2,3) -> +1 -3 -2
    reverse( pi, ii + 1, ii + (kk - jj) );     // reverse -3    -> +1 +3 -2
    break;

  case SCHEME_3_4: 
    reverse( pi, ii + 1, kk );                 // reverse (2,3) -> +1 -3 -2
    reverse( pi, kk - (jj - ii) + 1, kk );     // reverse (-2)  -> +1 -3 +2
    break;

  case SCHEME_3_5: // start +a 1 +b 2 +c 3
    hh = pi[ii];
    pi[ii] = pi[jj];
    pi[jj] = pi[kk];
    pi[kk] = hh;                          // end +a 2 +b 3 +c 1
    if ( ii == 0 )
      restore0( pi, kk );
    break;

  case SCHEME_3_6: // start +a 1 +b 2 +c 3
    hh = pi[ii];
    pi[ii] = pi[kk];
    pi[kk] = pi[jj];
    pi[jj] = hh;                          // end +a 3 +b 1 +c 2
    if ( ii == 0 )
      restore0( pi, jj );
    break;

  case SCHEME_3_7: // end +a 3 -b 2 -c 1
    hh = pi[ii];
    pi[ii] = pi[kk];
    pi[kk] = hh;
    reverse( pi, ii + 1, jj - 1 );
    reverse( pi, jj + 1, kk - 1 );
    if ( ii == 0 )
      restore0( pi, kk );
    break;

  case SCHEME_3_8: // end +a -c +b 3 1 2
    reverse( pi, ii, kk - 1 );
    reverse( pi, ii + kk - jj - 1, kk - 2 );

    hh = pi[kk];
    pi[kk] = pi[kk - 2];
    pi[kk - 2] = hh;
    if ( ii == 0 )
      restore0( pi, kk - 1);
    break;

  case SCHEME_3_9: // end +a 2 b+ 1 -c
    reverse( pi, jj + 1, kk );
    hh = pi[ii];
    pi[ii] = pi[jj];
    pi[jj] = hh;

    if ( ii == 0 )
      restore0( pi, jj);
    break;
    
  default:
    printf(" ERR perf move.  i %d  j %d  k %d  h %d scheme %d \n", ii, jj, kk, hh, scheme);
    exit (1);
  }  
}


static void setSchemes( int code, int* schemes, int* nsch )
{

  /* given code, assigns proper schemes in array schemes[]

    code = 0 -> 3opt "standard"  (schemes r_1, .., r_4  in paper)
    code = 1 -> 3opt "special 6opt" (schemes r_5...r_8 in paper)
    code = 2 -> 3opt "special 5opt" (scheme r_9 in paper)
    code = 3 -> all quadratic 3opt (schemes r_1...r_8 )
    code = 4 -> all 3opt (schemes r_1...r_9 )
  */
  
  int i, nschemes;

  switch ( code ) {
  case 0 :
    nschemes = 4;
    for ( i = 0; i < nschemes; i++ )
      schemes[i] = SCHEME_3_1 + i;
    break;
  case 1: 
    nschemes = 4;
    for ( i = 0; i < nschemes; i++ )
      schemes[i] = SCHEME_3_5 + i;
    break;
  case 2:
    nschemes = 1;
    schemes[0] = SCHEME_3_9;
    break;
  case 3:
    nschemes = 8;
    for ( i = 0; i < nschemes; i++ )
      schemes[i] = SCHEME_3_1 + i;
    break;
  case 4:
    nschemes = 9;
    for ( i = 0; i < nschemes; i++ )
      schemes[i] = SCHEME_3_1 + i;
    break;
  default :
    // ( code >= 0 && code <= 4) must be true
    nschemes = 0; // to shut up the compiler
  }
  *nsch = nschemes;
}
  
static double evaluate3OPTgeneral3( int* pi, int i, int j, int k, int nscheme ) {
  // evaluates cost of move given the 3 indexes

  double cost_out, cost_in;

  
  switch ( nscheme ) {
  case SCHEME_3_1:
    cost_out = cost[pi[i]][pi[SUCC(i)]] + cost[pi[j]][pi[SUCC(j)]] + cost[pi[k]][pi[SUCC(k)]];
    cost_in = cost[pi[i]][pi[SUCC(j)]] + cost[pi[k]][pi[SUCC(i)]] + cost[pi[j]][pi[SUCC(k)]];
    break;
    
  case SCHEME_3_2:
    cost_out = cost[pi[i]][pi[SUCC(i)]] + cost[pi[j]][pi[SUCC(j)]] + cost[pi[k]][pi[SUCC(k)]];
    cost_in = cost[pi[i]][pi[j]] + cost[pi[SUCC(j)]][pi[SUCC(k)]] + cost[pi[k]][pi[SUCC(i)]];
    break;

  case SCHEME_3_3:
    cost_out = cost[pi[i]][pi[SUCC(i)]] + cost[pi[j]][pi[SUCC(j)]] + cost[pi[k]][pi[SUCC(k)]];
    cost_in = cost[pi[i]][pi[SUCC(j)]] + cost[pi[j]][pi[k]] + cost[pi[SUCC(k)]][pi[SUCC(i)]];
    break;

  case SCHEME_3_4:
    cost_out = cost[pi[i]][pi[SUCC(i)]] + cost[pi[j]][pi[SUCC(j)]] + cost[pi[k]][pi[SUCC(k)]];
    cost_in = cost[pi[i]][pi[k]] + cost[pi[SUCC(i)]][pi[SUCC(j)]] + cost[pi[j]][pi[SUCC(k)]];
    break;
    
  case SCHEME_3_5:
    cost_out = cost[pi[i]][pi[SUCC(i)]] + cost[pi[j]][pi[SUCC(j)]] + cost[pi[k]][pi[SUCC(k)]] +
	           cost[pi[i]][pi[PRED(i)]] + cost[pi[j]][pi[PRED(j)]] + cost[pi[k]][pi[PRED(k)]];
    cost_in = cost[pi[i]][pi[SUCC(k)]] + cost[pi[i]][pi[PRED(k)]] +
              cost[pi[j]][pi[SUCC(i)]] + cost[pi[j]][pi[PRED(i)]] +
              cost[pi[k]][pi[SUCC(j)]] + cost[pi[k]][pi[PRED(j)]];

    break;
    
  case SCHEME_3_6:
    cost_out = cost[pi[i]][pi[SUCC(i)]] + cost[pi[j]][pi[SUCC(j)]] + cost[pi[k]][pi[SUCC(k)]] +
               cost[pi[i]][pi[PRED(i)]] + cost[pi[j]][pi[PRED(j)]] + cost[pi[k]][pi[PRED(k)]];
    cost_in = cost[pi[i]][pi[SUCC(j)]] + cost[pi[i]][pi[PRED(j)]] +
              cost[pi[j]][pi[SUCC(k)]] + cost[pi[j]][pi[PRED(k)]] +
              cost[pi[k]][pi[SUCC(i)]] + cost[pi[k]][pi[PRED(i)]];

    break;

  case SCHEME_3_7:
    cost_out = cost[pi[i]][pi[SUCC(i)]] + cost[pi[j]][pi[SUCC(j)]] + cost[pi[k]][pi[SUCC(k)]] +
               cost[pi[i]][pi[PRED(i)]] + cost[pi[j]][pi[PRED(j)]] + cost[pi[k]][pi[PRED(k)]];
    cost_in = cost[pi[i]][pi[SUCC(j)]] + cost[pi[i]][pi[SUCC(k)]] +
              cost[pi[j]][pi[SUCC(i)]] + cost[pi[j]][pi[PRED(k)]] +
              cost[pi[k]][pi[PRED(j)]] + cost[pi[k]][pi[PRED(i)]];

    break;
    
  case SCHEME_3_8:
    cost_out = cost[pi[i]][pi[SUCC(i)]] + cost[pi[j]][pi[SUCC(j)]] + cost[pi[k]][pi[SUCC(k)]] +
               cost[pi[i]][pi[PRED(i)]] + cost[pi[j]][pi[PRED(j)]] + cost[pi[k]][pi[PRED(k)]];
    cost_in = cost[pi[i]][pi[j]] + cost[pi[i]][pi[k]] +
              cost[pi[j]][pi[SUCC(k)]] + cost[pi[SUCC(j)]][pi[SUCC(i)]] +
              cost[pi[k]][pi[PRED(j)]] + cost[pi[PRED(k)]][pi[PRED(i)]];

    break;
    
  case SCHEME_3_9:
    cost_out = cost[pi[i]][pi[SUCC(i)]] + cost[pi[j]][pi[SUCC(j)]] + cost[pi[k]][pi[SUCC(k)]] +
               cost[pi[i]][pi[PRED(i)]] + cost[pi[j]][pi[PRED(j)]];
    cost_in = cost[pi[i]][pi[PRED(j)]] + cost[pi[i]][pi[k]] +
              cost[pi[j]][pi[SUCC(i)]] + cost[pi[j]][pi[PRED(i)]] +
              cost[pi[SUCC(k)]][pi[SUCC(j)]];

    break;
    
  default:
    printf("3OPT scheme %d unknown\n", nscheme);
    exit(1);
  }
  
  return cost_out - cost_in;
}

static double evaluate3OPTgeneral( int* pi, int* sel, int nscheme ) {
  // evaluates cost of move given the selection
  
  return evaluate3OPTgeneral3( pi, sel[0], sel[1], sel[2], nscheme );
}

/*****************************************************

  Functions used as keys to build the heaps for 
  various type of moves 

******************************************************/


static double tauplus( int* pi, int x, int y ) { 
  //  returns the value tau^+(x,y)
  return cost[pi[x]][pi[SUCC(x)]] - cost[pi[x]][pi[SUCC(y)]];
}

static double tauminus( int* pi, int x, int y ) {  
  //  returns the value tau^-(x,y)
  return cost[pi[x]][pi[PRED(x)]] - cost[pi[x]][pi[PRED(y)]];
}


double function3OPT(int* pi, int x, int y, int rscheme, int ndx0, int ndx1) {
  
  switch ( rscheme ) {
  case SCHEME_3_1 :
    if ( ndx0 == 1 ) // case j k 
      return tauplus(pi, x,y);
    else if ( ndx1 == 1 ) // case i j
      return tauplus(pi, x,y);
    else // case i k
      return tauplus(pi, y,x);
      
  case SCHEME_3_2 :
    if ( ndx0 == 1 ) // case j k 
      return tauminus(pi, SUCC(x),SUCC(SUCC(y)));
    else if ( ndx1 == 1 ) // case i j
      return tauplus(pi, x,PRED(y));
    else // case i k
      return tauplus(pi, y,x);
      
  case SCHEME_3_3 :
    if ( ndx0 == 1 ) // case j k 
      return tauplus(pi, x,PRED(y));
    else if ( ndx1 == 1 ) // case i j
      return tauplus(pi, x,y);
    else // case i k
      return tauminus(pi, SUCC(y),SUCC(SUCC(x)));
      
  case SCHEME_3_4 :
    if ( ndx0 == 1 ) // case j k 
      return tauplus(pi, x,y);
    else if ( ndx1 == 1 ) // case i j
      return tauminus(pi, SUCC(x),SUCC(SUCC(y)));
    else // case i k
      return tauplus(pi, y,PRED(x));
  
  case SCHEME_3_5 :
    if ( ndx0 == 1 ) // case j k 
      return tauplus(pi, y, x ) + tauminus(pi, y, x );
    else if ( ndx1 == 1 ) // case i j
      return tauplus(pi, y, x ) + tauminus(pi, y, x );
    else // case i k
      return tauplus(pi, x, y ) + tauminus(pi, x, y );
  
  case SCHEME_3_6 :
    if ( ndx0 == 1 ) // case j k 
      return tauplus(pi, x, y ) + tauminus(pi, x, y );
    else if ( ndx1 == 1 ) // case i j
      return tauplus(pi, x, y ) + tauminus(pi, x, y );
    else // case i k
      return tauplus(pi, y, x ) + tauminus(pi, y, x );
  
  case SCHEME_3_7 :
    if ( ndx0 == 1 ) // case j k 
      return tauminus(pi, x, y ) + tauminus(pi, y, x );
    else if ( ndx1 == 1 ) // case i j
      return tauplus(pi, x, y ) + tauplus(pi, y, x );
    else // case i k
      return tauplus(pi, PRED(x), PRED(y) ) + tauminus(pi, SUCC(y), SUCC(x) );

  case SCHEME_3_8 :
    if ( ndx0 == 1 ) // case j k
      {
	    return tauplus(pi, x, y ) + tauminus(pi, y, x );
      }
    else if ( ndx1 == 1 ) // case i j
      {
		return tauminus(pi, SUCC(x), SUCC(SUCC(y)) ) + tauminus(pi, y, SUCC(x) );
      }
    else // case i k
      {
        return tauplus(pi, PRED(x), PRED(PRED(y)) ) + tauplus(pi, y, PRED(x) );
      }

  case SCHEME_3_9 :
    if ( ndx0 == 1 ) // case j k 
      return tauminus(pi, SUCC(x), SUCC(SUCC(y)) );
    else if ( ndx1 == 1 ) // case i j
      return tauminus(pi, SUCC(x), SUCC(y) ) + tauplus(pi, PRED(x), PRED(y) ) + 
	         tauplus(pi, PRED(y), PRED(x) );
    else // case i k
      return tauplus(pi, y, PRED(x) );

  default:
    printf("Reins scheme %d for 3OPT unknown\n", rscheme );
    exit(1);
  }
}



/********************************************************
 ********************************************************

  Brute-force procedures

  ******************************************************
  ******************************************************/

static double find3OPTMoveBFOverall( int* pi, int* schemelist, int totschemes, int* bestsel,
                                     int* bestscheme, long int* eval, int* noptima, long int* nEffective ) {

  /* finds the best 3OPT move over all schemes of the array schemelist[],
     which has <totschemes> in it 

     Output: (bestsel[], bestscheme) -> the best move
             <noptima> -> the no. of global optima
  */

  double bestval, vv;
  int s, sch, numbest = 0;
  int _i0, _i1, _i2;  // the selection indices  
  long  int ncomb;
  long int nEff;
  
  bestval = -INFINITE;
  *bestscheme = -1;
  ncomb = 0;
  nEff=0;
  for ( s = 0; s < totschemes; s++ ) {
    sch = schemelist[s];
    for ( _i0 = 0; _i0 <= n - 5; _i0++ ) {
      for ( _i1 = _i0 + 2; _i1 <= n - 3; _i1++ ) {
	for ( _i2 = _i1 + 2; _i2 <= n - 1 - (_i0==0); _i2++ ) 
	  {
	    ncomb++;
	    dprintf("scheme %d %d %d %d\n", sch, _i0, _i1, _i2);
	    
	    vv = evaluate3OPTgeneral3( pi, _i0, _i1, _i2,  sch );
        
	    if ( (FINDBEST && (vv >= bestval - EPS) && (vv <= bestval + EPS)) )
	      numbest++;
	    
        
	    if ( (FINDBEST && vv > bestval + EPS) || (!FINDBEST && vv > EPS) ) 
        {
             nEff++;
             *bestscheme = sch;
	         bestval = vv;
	         bestsel[0] = _i0;
	         bestsel[1] = _i1;
	         bestsel[2] = _i2;
	         numbest = 1;

	         if ( !FINDBEST ) goto exit3BFloop;
	    }
	  }
      }
    }
  }

 exit3BFloop:
  
  *noptima = numbest;
  *eval = ncomb;
  *nEffective = nEff;
  return bestval;
}

static void find3optMoveBF( int* schemes, int nschemes, int* currentTour, OPTMOVE* move)
{

  static int bestsel[4];
  int bestscheme, i, noptima;
  double valmove, tbefore;
  long int evals;
  long int nEffective;
  
  tbefore = timer_elapsed_seconds();
  valmove = find3OPTMoveBFOverall( currentTour, schemes, nschemes, bestsel, &bestscheme,
                                   &evals, &noptima, &nEffective ); 
  move->seconds2Evaluate = timer_elapsed_seconds() - tbefore;

  move->found = (valmove > EPS);
  move->delta = valmove;
  move->moveEvaluations = evals; // shuold be = 4 * numCompleteTrips( n ) if best=true
  move->numberOfGlobalOptima = noptima; // defined only for BF w/best-improving
  move->effectiveExtractions = nEffective; 
  move->maxHeapSize = 0;          // undef for BF
  move->heapExtractions = 0;      // undef for BF 
  move->secondsToBuildHeap = 0;   // undef for BF 
  if ( move->found ) {  
    for ( i = 0; i < 3; i++ ) 
      move->selection[i] = bestsel[i];

    move->scheme = bestscheme;
  }
}

void find3optMove( int searchType, int code, int* currentTour, OPTMOVE* move)
{
  /*
  code = 0 -> 3opt "standard"  (schemes r_1, .., r_4  in paper)
  code = 1 -> 3opt "special 6opt" (schemes r_5...r_8 in paper)
  code = 2 -> 3opt "special 5opt" (scheme r_9 in paper)
  code = 3 -> all quadratic 3opt (schemes r_1...r_8 )
  code = 4 -> all 3opt (schemes r_1...r_9 )
  */

  int schemes[10], nschemes;

  if ( searchType == -1 ) { // undefined
    printf("ERR : 3OPT algorithm %d\n", searchType );
    assert( false );  
  }
  
  move->found = false;
  move->delta = -INFINITE;
  setSchemes( code, schemes, &nschemes );
  if ( searchType == SEARCH_3OPT_SF ) {
    find3optMoveSF( schemes, nschemes, currentTour, move);
  }
  else if ( searchType == SEARCH_3OPT_BF ) 
    find3optMoveBF( schemes, nschemes, currentTour, move);
  else
    assert( false );
}

/**************************************************************

    SMART-FORCE PROCEDURES

****************************************************************/

// +++++++++++++++++++

static int DIST( int x, int y ) {
  if ( y >= x )
    if ( x + n - y > y - x )
      return y-x;
    else
      return x + n - y;
  else
    if ( y + n - x > x - y )
      return x-y;
    else
      return y + n - x;
}


static void  randomSelectionFast(int* sel) {
  
  // fast generation of a random triple such that
  // dist(i,j), dist(j,k), dist(i,k) are all > 2
  
  int ii, kk, jj;

  ii = rand() % n;

  jj = rand() % n;
  while ( DIST(ii,jj) <= 1 )
    jj = rand() % n;

  kk = rand() % n;
  while ( DIST(ii, kk) <=1 || DIST(jj, kk) <= 1 ) 
    kk = rand() % n;

  if ( ii < jj  && jj < kk ) {
      sel[0] = ii;
      sel[1] = jj;
      sel[2] = kk;
  }
  else if ( ii < kk  && kk < jj ) {
    sel[0] = ii;
    sel[1] = kk;
    sel[2] = jj;
  }
  else if ( jj < ii  && ii < kk ) {
    sel[0] = jj;
    sel[1] = ii;
    sel[2] = kk;
  }
  else if ( jj < kk  && kk < ii ) {
    sel[0] = jj;
    sel[1]  = kk;
    sel[2]  = ii;
  }
  else if ( ii < jj  && jj < kk ) {
    sel[0] = ii;
    sel[1] = jj;
    sel[2] = kk;
  }
  else if ( kk < ii  && ii < jj ) {
    sel[0] = kk;
    sel[1] = ii;
    sel[2] = jj;
  }
  else {
    sel[0] = kk;
    sel[1] = jj;
    sel[2] = ii;
  }
}


static double getFirstChampion( int* pi, int* bsel, int* rscheme, long int* nevals ) {
  /* finds a starting solution of value > 0 or returns undefined, of value 0 
     returns the value of the starting solution.
     Counts the number of evaluations made and returns it in <nevals>
 */

  double bval, val;
  int r, try, i, nev;
  int sel[3];

  int NTRIES = 2 * n;   
  nev = 0;
  bval = 0;
  bsel[0] = bsel[1] = bsel[2] = -1;   // undef
  for ( r = SCHEME_3_1; r <= SCHEME_3_4; r++ )
    for ( try = 0; try < NTRIES; try++ ) {
      randomSelectionFast( sel ) ;
      val = evaluate3OPTgeneral ( pi, sel, r );
      nev++;
      if ( val > bval + EPS ) {
	bval = val;
	*rscheme = r;
	for ( i = 0; i < 3; i++ )
	  bsel[i] = sel[i];
	if ( !FINDBEST && bval > EPS )
	  goto exitChamp;
      }
    }

 exitChamp:
  
  *nevals = nev;
  return bval; 
}


double find3OPTMoveSFOverallSingleHeap( int* pi, int* schemes, int totschemes,
					int* bestsel, int* bscheme,
					int* heapSize, int* extract, double* timeHeapConstruction,
					double* timeHeapUsage, long int *eval, long int *evalGood ) {

  /* finds the best/first 3OPT move over all orbits.
     Uses a single heap as in the last version of the paper. Each selection/reins scheme
     is put on the heap for a partial move (see paper)

     Returns the selection (bestsel[]) and the move, both as a single code number identifying
	 the reinsertion scheme ( bscheme ) and a pair bestorbit,bestrot w/respect to the orbits
	 as in the paper of 4OPT

     in <heapSize> it returns the size of the starting heap
     in <extract[0]> it returns the number of pops that were done
     in <extract[1]> it returns the number of pops that led to an improvement of champion

  */
  
  double bestval, vv;
  int bestscheme, z, minz, maxz, missing;
  int sel[3], abortPreflop;
  long int neval, nevalGood;
  double tbefore;

  THEAP HEAP;
  TPAIR top;  // the max-element of the heap

  abortPreflop = false;
  bestval = 0;
  bestscheme = -1; // undef
  tbefore = timer_elapsed_seconds();
  nevalGood = 0;
  bestval = getFirstChampion( pi, bestsel, &bestscheme, &neval );
  if (!FINDBEST && bestval > EPS) {
    abortPreflop = true;
    *timeHeapConstruction = 0;
    *heapSize = 0;
    goto end3SFloop;
  }

  initializeSingleHeap( n, &HEAP, pi, totschemes, schemes, bestval / 3);
  *timeHeapConstruction = timer_elapsed_seconds() - tbefore;
  *heapSize = heapsize( HEAP.heap );

  tbefore = timer_elapsed_seconds();
  while ( heapsize( HEAP.heap ) > 0  ) {

    getMaxHeap( &top, HEAP.heap );
    if ( top.value < bestval / 3 + EPS ) 
      break; 

    dprintf("Pop el (x %d, y %d) val %g  r %d label %d\n", top.x, top.y, top.value, top.rscheme, top.paircode );

    switch ( top.paircode) {
    case 01:   // x=i, y= j, missing k
      sel[0] = top.x;
      sel[1] = top.y;
      missing = 2;
      minz = top.y + 2;
      maxz = n - 1 - (top.x==0); 
      break;

    case 12: // x=j, y= k, missing i
      sel[1] = top.x;
      sel[2] = top.y;
      missing = 0;
      minz = (top.y == n-1);
      maxz = top.x - 2; 
      break;

    case 02:  // x=i, y=k, missing j
      sel[0] = top.x;
      sel[2] = top.y;
      missing = 1;
      minz = top.x + 2;
      maxz = top.y - 2; 
      break;
    default:// to shut up compiler warnings
      missing=0; //any value
      minz=0;//any value
      maxz=0;//any value
    }

    for ( z = minz; z <= maxz; z++ ) {
      neval++;
      sel[missing]=z;
      vv = evaluate3OPTgeneral(pi, sel, top.rscheme );
      dprintf("EVAL  i %d j %d k %d  missing %d z %d  rsch %d val %g\n", 
	          sel[0], sel[1], sel[2], missing, z, top.rscheme, vv );
     if ( (FINDBEST && vv > bestval + EPS) || (!FINDBEST && vv > EPS) ) {
	nevalGood++;
	bestval = vv;
	bestscheme = top.rscheme;
	bestsel[0] = sel[0];
	bestsel[1] = sel[1];
	bestsel[2] = sel[2];
	if ( !FINDBEST )
	    goto end3SFloop;
      }
    }
  }

 end3SFloop:

  *timeHeapUsage = timer_elapsed_seconds() - tbefore;
  if ( abortPreflop ) { // heap was not created 
    extract[0] = 0;
    nevalGood = 0;
  }
  else {
    extract[0] = *heapSize - heapsize( HEAP.heap );  
	// startingSize - finalSize = elements popped
    freeHeapMem( &HEAP );
  }
  
  *bscheme = bestscheme;
  *eval = neval;
  *evalGood = nevalGood;

  return bestval;
}


void find3optMoveSF( int *schemes, int nschemes, int* currentTour, OPTMOVE* move)
{

/*
   A TSP for a complete graph is described by (cost[][], n), -global vars-
   where n is the number of graph nodes
   and cost[][] is a double[size][size] symmetric matrix of edge costs.

   A K-OPT move is described by K removed edges and K inserted edges
   (an edge can be both removed and inserted), each edge is a couple
   of indexes in [0,size-1].
   No order is required among edges nor between indexes in an edge.
*/

  int bestsel[4], scheme, heapsize, extractions[2], i;
  double valmove, time2build, time2eval;
  long int evals, evalsGood; 

  valmove = find3OPTMoveSFOverallSingleHeap( currentTour, schemes, nschemes, bestsel,
                                             &scheme, &heapsize, extractions, &time2build,
				             &time2eval, &evals, &evalsGood ); 

  move->found = (valmove > EPS);

  move->delta = valmove;
  move->maxHeapSize = heapsize;                // L value in our model 
  move->heapExtractions = extractions[0];      // M value in our model
  move->effectiveExtractions = extractions[1]; // number of heap extractions producing an update 
  move->secondsToBuildHeap = time2build;       // finding elements to insert into heap and building it 
  move->seconds2Evaluate = time2eval;          // extracting elements from heap and evaluating moves
  move->moveEvaluations = evals;               // approx nxM value in our model
  move->effectiveExtractions = evalsGood;      // a fraction of evaluations, possibly small

  if ( move->found ) {

    for ( i = 0; i < 3; i++ ) 
      move->selection[i] = bestsel[i];

    move->scheme = scheme;

  }
}



#ifdef TUTTODARIVEDERE

/****************************************************

  Global variables used by 3OPT procedures

*****************************************************/

static int MAX_ALLOWED_COMB = 32000000;    // used to estimate BF moves time. Aborts after
                                           // so many evaluations


static long int numCompleteTrips( int N ) {

  long a;
  a = (((long int)(N-2)*(N-3) / 2) * (N-4)) / 3;
  // must be a - n + 4 >= 0
  return a - (n - 4);
}




static double estimateTime3OPTMoveBF( int* pi, int* schemelist, int totschemes ) {

  /* estimates the time to find the best 3OPT move over all schemes of the array schemelist[], 
     which has <totschemes> in it 
  */

  double tbefore, lap;
  int s, sch;
  int _i0, _i1, _i2;  // the selection indices  
  double realTotMoves;
  long int ncomb;

  tbefore = timer_elapsed_seconds();
  ncomb = 0;
  for ( s = 0; s < totschemes; s++ ) {
    sch = schemelist[s];
    for ( _i0 = 0; _i0 <= n - 5; _i0++ ) {
      for ( _i1 = _i0 + 2; _i1 <= n - 3; _i1++ ) {
	for ( _i2 = _i1 + 2; _i2 <= n - 1 - (_i0==0); _i2++ ) 
	  {
	    ncomb++;
	    if ( ncomb == MAX_ALLOWED_COMB )
	      goto endofloop;
	    
	    evaluate3OPTgeneral3( pi, _i0, _i1, _i2,  sch );
	  }
      }
    }
  }

 endofloop:

  lap = timer_elapsed_seconds() - tbefore;
  realTotMoves = ((double)totschemes) * numCompleteTrips(n);  // risk of overflow?***
  return lap * (realTotMoves/(double)MAX_ALLOWED_COMB);
}

/********************************************************
 ********************************************************

  Smart-force procedures

  ******************************************************
  ******************************************************/



double estimateTime3optMoveBF( double** weights, int size, int* perm, int code ) {

  /*
    Returns estimated running time of a BF move over the specified perm
 */

  int schemes[10], nschemes;

  n = size;
  cost = weights;
  
  setSchemes( code, schemes, &nschemes );
  return estimateTime3OPTMoveBF( perm, schemes, nschemes );
}



void beforeRunning( int nnodes ) {

  /* allocates all data structures needed by 3OPT procedures, given number of nodes in graph */

  int c;

  bufferPerm = (int*) calloc( nnodes + 1, sizeof(int));
  SUCCV = (int*) calloc( nnodes + 1, sizeof(int));
  PREDV = (int*) calloc( nnodes + 1, sizeof(int));
  // (bufferPerm && SUCCV && PREDV) must be true
  for ( c = 0; c < nnodes; c++ ) {
    SUCCV[c] = (c+1) % nnodes;
    PREDV[c] = (c + nnodes - 1) % nnodes;
  }
}

void afterRunning() {

  free( bufferPerm );
  free( SUCCV );
  free( PREDV );
}

#endif
