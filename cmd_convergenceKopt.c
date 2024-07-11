/***************************************************************

command-line procedure to run a convergence to a local optimum

****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "tspmain.h"
#include "proc4opt.h"
#include "proc3opt.h"
#include "io4opt.h"


int n = -1;  // undefined
double** cost = NULL;
int OUTPUTON;
long int tottriples;
int bufferPerm[MAXN + 1];
TORBIT ORBIT3[2];
TORBIT ORBIT4[7];
int SUCCV[MAXN], PREDV[MAXN];
TPARAM param;
int FINDBEST;
int COSTI_EUCLIDEI; // se = 0, usa costi uniformi
int costReadFromFile;
int SHOWMOVE;    // mostra le mosse che trova lungo la LS

int SEEDTOUR, SEEDCOST;
int NUM_LS_STEPS;
double TIMELIMIT;
double globalTimeStart;
char startTourName[100];
char costFileName[100];
char outTourName[100];
char outLastTourName[100];
char outOptName[100];
char outMoveName[100];
char outTimeName[100];
char outValName[100];
char outTSPName[100];
double SWITCH2COEF;    /* switch to 2OPT CE at step = SWITCH2COEF * n */
double SWITCHIMPRCOEF;    /* change BI<->FI at step = SWITCHIMPRCOEF * n */

double BESTTVAL, AVGTVAL; // per le stampe, rimuovere poi

int SEARCHTYPE[5];  // searchtype[k] = algorithm used for k-opt

int NUMTESTS, NUMINSTANCES, PROJECTED, COMPUTE_FROM_STEP, DOPERTURB;
double MAX_NODE_U;   // each node i might get assigned 0 <= u[i] < MAX_NODE_U
                     // to then perturbe arcs as cost[i][j] + u[i] + u[j]
int orbitList[7], norbits, kloop[5];

int RANDOM_GRAPH_TYPE = -1;   // set initially to unknown
int PROJECTED = -1;   // set initially to undefined

// tpair* SORTED; /* array of sorted tour edges, used in a version of the 2opt search */
extern int sortedExists; /* semaphore saying that the array SORTED has already been built and can be used */

void allocateSorted(); // allocate mem for array of sorted tour edges, used in a version of the 2opt search
void disposeSorted();  // free above mem

typedef struct{

  int numsteps;
  double tottime;
  double totevals;
  double bestfound;
} TEXPSTAT;  // statistics for a single experiment

typedef struct{

  int numsteps;
  double tottime;
  double totevals;
  double bestfound;
} TCONVSTAT;  // statistics for a single convergence

typedef struct{

  double time;
  double nevals;
  double valfound;  // value after the move
} TSTEPSTAT;  // statistics for a single convergence

/********************************************
 some global vars used for the convergences
*********************************************/

int doReadTour, doWriteTour, doWriteLastTour;
int doWriteMove, doWriteTime, doWriteVal, doWriteOpt, doWriteTSP;
char fullname[200];
FILE *movelog = NULL, *timelog = NULL, *vallog = NULL;


void showMove( int neigh, OPTMOVE* mov ) {

  assert( !mov->found || ( mov->selection[0] < mov->selection[1] ));
  if ( neigh > 2)
    assert( !mov->found || ( mov->selection[1] < mov->selection[2] ));
  if ( neigh > 3)
    assert( !mov->found || ( mov->selection[2] < mov->selection[3] ));

  /*
  printf("Delta               : %lg\n", mov->delta );
  printf("Selection           : %d %d %d %d\n", mov->selection[0], mov->selection[1], mov->selection[2],  mov->selection[3]);
  printf("maxHeapSize         : %d\n", mov->maxHeapSize );
  printf("Scheme              : %d\n", mov->scheme );
  printf("OrbitNo             : %d\n", mov->orbitNo );
  printf("Phi                 : %d\n", mov->phi );
  printf("HeapExtractions     : %d\n", mov->heapExtractions );
  printf("effectExtractions   : %d\n", mov->effectiveExtractions);
  printf("MoveEvaluations     : %ld\n", mov->moveEvaluations );
  printf("timeToEvaluate      : %g\n", mov->secondsToEvaluate );
  */

  if ( SHOWMOVE ) {
    printf("\n");
    printf("#DELTA              : %9.5f\n", mov->delta );
    /**/
    printf("#I1                 : %d\n", mov->selection[0]);
    printf("#I2                 : %d\n", mov->selection[1]);
    if ( neigh > 2 )
      printf("#I3                 : %d\n", mov->selection[2]);
    if ( neigh > 3 ) {
      printf("#I4                 : %d\n", mov->selection[3]);
      printf("#SCHEME             : %d\n", mov->scheme );
    }

  /**/
  printf("#EVALS              : %ld\n", mov->moveEvaluations );
  printf("#TIME               : %g\n", mov->seconds2Evaluate );
  }
  
  //  printf("OrbitNo             : %d\n", mov->orbitNo );
  //  printf("Phi                 : %d\n", mov->phi );
  //  printf("maxHeapSize         : %d\n", mov->maxHeapSize );
  //  printf("HeapExtractions     : %d\n", mov->heapExtractions );
  //  printf("effectExtractions   : %d\n", mov->effectiveExtractions);
}

void showPerm(int* pi, int k)
{
  int i;

  //  if ( k > 20 ) return;

  printf("(");
  for (i=0; i<k; i++)
    {
       printf("%2d", pi[i]);
       if (i<k-1) printf(", "); 
    }
  printf(") %d\n", pi[k]);
}

void identity( int* pi, int n_ ) {

  int i;

  for ( i = 0; i < n_; i ++ )
    pi[i] = i;
  pi[n_] = 0;
}

void prepareSuccPred( int numnodes ) {

  // stores succ() and pred() in an array so as to make it faster to compute
  // them then with a function. numnodes shound be set to n
  
  int c;
  for ( c = 0; c < numnodes; c++ ) {
    SUCCV[c] = (c+1) % numnodes;
    PREDV[c] = (c + numnodes - 1) % numnodes;
  }
}



void allocateGlobalHeaps() {

  // allocate all the heaps that are used by our procedures

  int i;

  printf("  **** ALLOC GLOBAL HEAPS\n");

  // the heap for SF_SIMILAR3OPT
  heapSingle.heap = (TPAIR*) mycalloc( 1 + n * (n-1) / 2, sizeof(TPAIR));
  assert( heapSingle.heap );
  
  // the heaps for SF0 and SF1
  heapF12.heap = (TPAIR*) mycalloc( 1 + n * (n-1) / 2, sizeof(TPAIR));
  assert( heapF12.heap );

  heapF34.heap = (TPAIR*) mycalloc( 1 + n * (n-1) / 2, sizeof(TPAIR));
  assert( heapF34.heap );


  // the heaps sortedheap[0,1] are used as sorted arrays in the algorithms SF0, SF1
  for ( i = 0; i < 2; i++ ) {
    sortedheap[i] = (TPAIR*) mycalloc( 1 + n * (n-1) / 2, sizeof(TPAIR));
    assert( sortedheap[i] );
  }

  tmpHeap = (TPAIR*) mycalloc( 1 + n * (n-1) / 2, sizeof(TPAIR));
  assert( tmpHeap );
}

void disposeGlobalHeaps() {
  int i;
  
  for ( i = 0; i < 2; i++ ) 
    free( sortedheap[i] );
}

int nextMoveType( int currK ) {

  // determines which K-opt to search after the current K did not
  // find an improving move

  if ( currK == 4 )
    return 3;
  else if ( currK == 3 )
    return 4;
  else if ( currK == 2 )
    return 4;
  else {
    assert(false);
    return 0;
  }
}


int useSortedEdges2OPT( int alg ) {
  // checks if 2opt algorithm used is based on having the
  // array of sorted tour edges

  return ( alg == SEARCH_2OPT_SF_LEXI ||
	   alg == SEARCH_2OPT_SF_SORTPAIR ||
	   alg == SEARCH_2OPT_H_SORTPAIR );
}

void nearestNeighbor( int* nntour, int n ) {

  int unreached[MAXN];  // tour[] has size cnt, unriched[] has n-cnt
  int curr, cnt, i, next, first;
  int tour[MAXN];
  
  cnt = 0;
  first = rand() % n; // random first node
  curr = first;  // the tour is a path from first..curr of cnt nodes. we extend from curr
  tour[cnt++] = curr;
  for ( i = 1; i <= n - 1; i++ )
    unreached[i-1] = (curr + i) % n;
  while ( cnt < n ) {
    next = 0;
    for ( i = 1; i < n - cnt; i++ )
      if ( cost[curr][unreached[i]] < cost[curr][unreached[next]] - EPS )
	next = i;

    curr = unreached[next];
    tour[cnt] = curr;
    unreached[next] = unreached[n - cnt - 1];
    cnt++;
  }
  /* now makes sure tour starts at 0 ends at 0 */
  i = cnt = 0;
  while (tour[i] != 0) i++;
  printf("0 is at %d\n", i);
  for (cnt = 0; cnt < n; cnt++)
    nntour[cnt] = tour[(i+cnt)%n];
  nntour[n] = nntour[0];
  printf("Nearest neighbor tour found :\n");
  showPerm( nntour, n );
  printf("Value : %g\n", tourValue( nntour ));
  //  scanf("%d", &i);
}

int timesUp() {
  return timer_elapsed_seconds() - globalTimeStart >= TIMELIMIT;
}

void runSingleConvergence( int idexp, int idconv, TCONVSTAT* cstat, TCONVSTAT* cstatlap,
			   TSTEPSTAT* stepstat, int keepsteps, int* perm ) {

  // runs a single convergence, made of many steps, on the current instance
  // outputs global statistics on the convergence into <cstat>
  // and for the part relative on the steps from COMPUTE_FROM_STEP... into <cstatlap>
  // if <keepsteps> returns also statistics on each single step of the convergence
  // in the array <stepstat>, which is allocated by caller
  // the final tour permutation is returned in <perm[]>
  
  // <idexp> and <idconv> are the codes for the number of experiment and convergence within the experiment

  int step, i;
  long int totmoves = 0, totmovesLap = 0;
  int currK, nowLooking, nextLooking, doNotRepeat;
  double tourval;
  double netTime, netTimeLap = 0;
  OPTMOVE move;
  int ORG_SEARCHTYPE2;  // original algorithm for 2OPT. Might be changed to CE within the search
                        // and is restored at the end
  
  if ( doReadTour ) {
    if ( strcmp ( startTourName, "NN" ) == 0 )
      nearestNeighbor( perm, n );
    else {
      readTour( startTourName, perm, n );
      assert( NUMTESTS == 1 );
    }
  }
  else
    rndPerm0( perm, n );
    
  if ( doWriteMove ) {
    sprintf( fullname, "%s.%02d.%02d", outMoveName, idexp, idconv );  
    movelog = fopen( fullname, "w" );
  }

  if ( doWriteTime ) {
    sprintf( fullname, "%s.%02d.%02d", outTimeName, idexp, idconv );  
    timelog = fopen( fullname, "w" );
  }
    
  if ( doWriteVal ) {
    sprintf( fullname, "%s.%02d.%02d", outValName, idexp, idconv );  
    vallog = fopen( fullname, "w" );
  }
  
  if ( doWriteTour ) 
    writeTour( outTourName, perm, n, 0, idconv, false );

  if (  useSortedEdges2OPT( SEARCHTYPE[2] ) )
    sortedExists = false;
  
  ORG_SEARCHTYPE2 = SEARCHTYPE[2];

  tourval = tourValue( perm );

  /*
  int j;
  for ( i = 0; i < n && i < 10; i++ ) {
    for ( j = 0; j < n && j < 10; j++ ) 
      printf("%8.2g ", cost[i][j]);
    printf("\n");
  }
  printf("val %g reval %g\n", tourval, sumc);
  scanf("%d", &j);
  */
  
  totmoves = 0;
  step = 0;
  netTime = 0;
  
  nowLooking = 0;  // starting K is kloop[0]
  doNotRepeat = 0;  // if gets back to this finding move, stop
  
  while ( step < NUM_LS_STEPS ) {

    if ( step == (int) (SWITCHIMPRCOEF * n) + 1 )
      FINDBEST = !FINDBEST;
	
    if ( SEARCHTYPE[2] != 0 && step > SWITCH2COEF * n ) {
      printf(">>> Switching 2OPT from %d to CE at step %d..\n", SEARCHTYPE[2], step );
      SEARCHTYPE[2] = 0;
    }
    if ( step == COMPUTE_FROM_STEP ) {
      netTimeLap = netTime;
      totmovesLap = totmoves;
    }
    
    printf("\n[step %5d val %9.3f]", step++, tourval );
    if ( SHOWMOVE ) {
      printf(" TOUR : ");
      for ( i = 0; i <= n; i++ ) {
	printf("%d ", perm[i] );
	if ( i > 10 ) {
	  printf("...");
	  break;
	}
      }
    }
    printf("\n");
    
    do {
      currK = kloop[nowLooking];
      assert( currK >= 2 && currK <= 4 );
      printf("curr opt-%d\n", currK);
      if ( currK == 2 && SEARCHTYPE[2] != -1 ) 
	find2optMove( true, SEARCHTYPE[2], perm, &move);  
      else if ( currK == 3 && SEARCHTYPE[3] != -1 )
	find3optMove( SEARCHTYPE[3], 0, perm, &move);  // search all 4 reinsertion schemes
      else if ( currK == 4 && SEARCHTYPE [4] != -1 )
	find4optMove( SEARCHTYPE[4], norbits, orbitList, perm, &move);
      
      if ( !move.found ) {
	nextLooking =  (nowLooking + 1) % 3;
	if ( nextLooking == doNotRepeat ) {
	  printf("cycle completed. we are at loc opt\n");
	  break; // we are at local optimum
	}
	else {
	  printf("Switch from %d-OPT to %d-OPT\n", currK, kloop[nextLooking]);
	  nowLooking = nextLooking;
	}
      }
      else {
	//	  printf("found at k%d. do not return to %d", currK, currK );
	//	  scanf("%d", &i);
	doNotRepeat = nowLooking;
	break; // we found improving move
      }
    } while ( true );
      
    showMove( currK, &move );
    assert(  !move.found || move.delta > EPS );
    totmoves += move.moveEvaluations;
    netTime += move.seconds2Evaluate;
    
    if ( !move.found ) 
      break; // we are at local optimum
    
    // else perform the move	
    if ( currK == 4 )
      perform4OPTMove( perm, move.selection, move.orbitNo, move.phi );
    else if ( currK == 3 )
      perform3OPTMove( perm, move.selection, move.scheme );
    else if ( currK == 2 ) 
      perform2OPTMove( perm, move.selection, useSortedEdges2OPT( SEARCHTYPE[2] ) );
    
    else
      assert( false );
    
    tourval -= move.delta;
    
    if ( doWriteTour )
      writeTour( outTourName, perm, n, step, idconv, false );
    if ( doWriteMove )
      fprintf( movelog, "%d %ld\n", step, move.moveEvaluations );
    if ( doWriteTime )
      fprintf( timelog, "%d %8.5f\n", step, move.seconds2Evaluate );
    if ( doWriteVal )
      fprintf( vallog, "%d %10.4f\n", step, tourval );
    
    if ( keepsteps ) {
      assert( step < MAXCONVERGENCELEN );
      stepstat[step].time = move.seconds2Evaluate;
      stepstat[step].valfound = tourval;
      stepstat[step].nevals = move.moveEvaluations;
    }
  }
  
  cstat->numsteps = step;    
  cstat->tottime = netTime;    
  cstat->totevals = totmoves;    
  cstat->bestfound = tourval;

  cstatlap->numsteps = step - COMPUTE_FROM_STEP;    
  cstatlap->tottime = netTime - netTimeLap;    
  cstatlap->totevals = totmoves - totmovesLap ;    
  cstatlap->bestfound = tourval;

  if ( doWriteLastTour ) 
    writeTour( outLastTourName, perm, n, step, idconv, true );
    
  if ( movelog != NULL )
    fclose( movelog ) ;
  if ( timelog != NULL )
    fclose( timelog ) ;
  if ( vallog != NULL )
    fclose( vallog ) ;

  SEARCHTYPE[2] = ORG_SEARCHTYPE2; // restore 2-OPT algo, if it was changed to CE
}


void runSingleExperiment(int idexp, TEXPSTAT* expstat) {

  /*
    An experiment is a series of (-NT) convergences (to local opt or for -NS steps),
    each from a random tour or a given tour,
    on the same instance
  */

  int i, t;
  double bestTourval, totTourval;
  int perm[MAXN + 1], bestTour[MAXN + 1];
  //  double cumulativeSteps, cumulativeMoves, cumulativeNetTime;
  double cumulativeProjectedTime, cumulativeMovesLap = 0, cumulativeStepsLap = 0, cumulativeNetTimeLap = 0;
  // taken @ instant of LAP start
  TCONVSTAT cstat, cstatlap;
  TSTEPSTAT stepstat[MAXCONVERGENCELEN];


  doReadTour = strlen(startTourName) > 0;
  doWriteTour = strlen(outTourName) > 0;
  doWriteLastTour = strlen(outLastTourName) > 0;
  doWriteMove = strlen(outMoveName) > 0;
  doWriteTime = strlen(outTimeName) > 0;
  doWriteVal = strlen(outValName) > 0;
  doWriteOpt = strlen(outOptName) > 0;
  doWriteTSP = strlen(outTSPName) > 0;

  assert( !doReadTour || NUMINSTANCES == 1 );
  assert( !doWriteTour || NUMINSTANCES == 1 );
  assert( !doWriteLastTour || NUMINSTANCES == 1 );
  
  if ( doWriteTSP )
    writeTSP( outTSPName );

  if ( SEARCHTYPE[2] == -1 && SEARCHTYPE[3] == -1 &&  SEARCHTYPE[4] == -1 ) {
    printf("Error: No search algorithm specified\n");
    exit(1);
  }


  if (  useSortedEdges2OPT( SEARCHTYPE[2] ) )
    allocateSorted();
  
  if ( SEARCHTYPE[4] == SEARCH_SF_SIMILAR3OPT ||
       SEARCHTYPE[4] == SEARCH_SF0 ||
       SEARCHTYPE[4] == SEARCH_SF1 )
    allocateGlobalHeaps();
  
  if ( SEARCHTYPE[4] == SEARCH_GLOVER )
    allocateGloverMem();

  expstat->numsteps = 0;
  expstat->tottime = 0;
  expstat->totevals = 0;

  cumulativeProjectedTime = 0;

  bestTourval = INFINITE;
  totTourval = 0;

  for ( t = 0; t < NUMTESTS && !timesUp(); t++ ) {

    printf("*** CONVERGENCE N. %d ***\n\n", t); fflush(stdout);
    runSingleConvergence( idexp, t, &cstat, &cstatlap, stepstat, false, perm );
    
    printf("\n#CONVERGENCE_STATS for CONVERGENCE n %d:\n", t);
    printf("#NNODES**            %d\n", n);
    printf("#SEEDTOUR            %d\n", SEEDTOUR );
    printf("#SEEDCOST            %d\n", SEEDCOST );
    printf("#OPTVAL              %9.3f\n", cstat.bestfound );
    printf("#STEPS               %d\n", cstat.numsteps );
    printf("#STEPS (FROM LAP)    %d\n", cstatlap.numsteps );
    printf("#MOVES               %ld\n", (long int)cstat.totevals );
    printf("#MOVES (FROM LAP)    %ld\n", (long int)cstatlap.totevals );
    printf("#NETTIME             %9.3f\n", cstat.tottime );
    printf("#NETTIME (FROM LAP)  %9.3f\n", cstatlap.tottime );
    
    if ( cstat.bestfound < bestTourval ) {
      bestTourval = cstat.bestfound;
      for ( i = 0; i <= n; i++ )
	bestTour[i] = perm[i];
      
      /*
	FINDBEST = !FINDBEST;  // per alternare tra best e first -- togliere poi
        oppure si può fare +1 mod K per avere 1 volta ogni K una best (o una first...)
      */
    }

    totTourval += cstat.bestfound;
    
    
    if ( PROJECTED > 0 ) 
      printf("#PROJECTED on %4d steps    %9.3f\n", PROJECTED, (cstat.tottime * PROJECTED) / cstat.numsteps );
    
    expstat->tottime += cstat.tottime;
    cumulativeNetTimeLap += cstatlap.tottime;
    
    expstat->totevals += cstat.totevals;
    cumulativeMovesLap += cstatlap.totevals;
    cumulativeProjectedTime += (cstat.tottime * PROJECTED) / cstat.numsteps;
    expstat->numsteps += cstat.numsteps; 
    cumulativeStepsLap += cstatlap.numsteps;
  }

  // we are here either because all convergences on this instance have been completed
  // or the time was exceeded

  double DNUMTESTS;
  if ( timesUp() ) {
    DNUMTESTS = t;
    printf("\n\n\n");
    printf("* * * * TIME EXCEEDED * * *\n");
    printf("* * * * TIME EXCEEDED * * *\n");
    printf("* * * * TIME EXCEEDED * * *\n");
    printf("\n**** partial statistics under time limitations ****\n\n");
    if ( NUMINSTANCES > 1 ) {
      printf("Statistics valid only when -NINS = 1. Here it is %d\n", NUMINSTANCES);
      exit(1);
    }
    printf("Experiments completed      %d\n", (int) DNUMTESTS);
    printf("Best tour value found      %*.5f\n", 15, bestTourval);
    printf("Avg  tour value            %*.5f\n", 15, totTourval/DNUMTESTS );
    printf("Tot steps executed         %d\n", expstat->numsteps);
    printf("Avg steps per run          %*.2f\n", 10, (expstat->numsteps)/DNUMTESTS );
    printf("Tot moves evaluated        %*.1f\n", 10, expstat->totevals );
    printf("Avg move evals per step    %*.2f\n", 10, (expstat->totevals)/expstat->numsteps );
    printf("Avg time per run           %*.5g\n", 15, expstat->tottime / DNUMTESTS);
    printf("\n**** end of partial statistics ****\n");

    exit(1);
  }

  DNUMTESTS = NUMTESTS;
  printf("\n* *  * Overall avg results on all %d tests\n\n", NUMTESTS);
  printf("Average net time/LS %*.5g\n", 15, expstat->tottime / DNUMTESTS);
  printf("Average net time/mov %*.5g\n", 15, expstat->tottime / expstat->numsteps);
  printf("Average evaluations %*.5f\n", 15, expstat->totevals /  expstat->numsteps);
  printf("Average LS length   %*.5f\n", 15, expstat->numsteps / DNUMTESTS);
  printf("Best tour val found %*.5f\n", 15, bestTourval);
  printf("Average tour val    %*.5f\n", 15, totTourval/NUMTESTS );
  BESTTVAL = bestTourval;
  AVGTVAL = totTourval/NUMTESTS;
  printf("Convergences run/started %d out of %d\n", t, NUMTESTS);
  if ( PROJECTED > 0 ) {
    printf("Average projection/LS%*.5f\n", 15, cumulativeProjectedTime / DNUMTESTS );
    printf("Average projection/move%*.5f\n", 15, cumulativeProjectedTime / (DNUMTESTS * PROJECTED) );
  }
  if ( COMPUTE_FROM_STEP >= 1 ) {
    printf("**** SCORES FROM LAP %d ****\n", COMPUTE_FROM_STEP );
    printf("Average net time/LS %*.5g\n", 15, cumulativeNetTimeLap / DNUMTESTS);
    printf("Average net time/mov %*.5g\n", 15, cumulativeNetTimeLap / cumulativeStepsLap);
    printf("Average evaluations %*.5f\n", 15, cumulativeMovesLap / cumulativeStepsLap );
    printf("Average length      %*.5f\n", 15, cumulativeStepsLap / DNUMTESTS);
    if ( SEARCHTYPE[4] == SEARCH_GLOVER ) {
      printf("& %d & ", n);
      printf("%.2f & - & %.4f & - & ",
	      expstat->tottime / DNUMTESTS, 
	      expstat->tottime /  expstat->numsteps );
      printf("%.1f & ",  expstat->numsteps / DNUMTESTS );
      printf("%.2f & - & %.4f & - & ",   cumulativeNetTimeLap / DNUMTESTS, cumulativeNetTimeLap / cumulativeStepsLap );
      printf("%.1f & ",  cumulativeStepsLap / DNUMTESTS );
      printf("%.1f \\\\ \n ", expstat->tottime - cumulativeNetTimeLap );
    }
    else if  ( SEARCHTYPE[4] == SEARCH_SF1 ) {
      printf("& %d & ", n);
      printf("- & %.2f & - & %.4f & ",
	     expstat->tottime / DNUMTESTS, 
             expstat->tottime / expstat->numsteps );
      printf("%.1f & ",  expstat->numsteps / DNUMTESTS );
      printf("- & %.2f & - & %.4f & ",   cumulativeNetTimeLap / DNUMTESTS, cumulativeNetTimeLap / cumulativeStepsLap );
      printf("%.1f & ",  cumulativeStepsLap / DNUMTESTS );
      printf("%.1f \\\\ \n ",  cumulativeNetTimeLap );
    }
  }

  freeCost();

  if (  useSortedEdges2OPT( SEARCHTYPE[2] ) )
    disposeSorted();
  
  if ( SEARCHTYPE[4] == SEARCH_GLOVER )
    disposeGloverMem();

  if ( SEARCHTYPE[4] == SEARCH_SF_SIMILAR3OPT ||
       SEARCHTYPE[4] == SEARCH_SF0 ||
       SEARCHTYPE[4] == SEARCH_SF1 )
    disposeGlobalHeaps();


  printf("\nTotal time since start of all operations %lg\n", timer_elapsed_seconds() - globalTimeStart );

  if ( doWriteOpt ) 
    writeBestTour( outOptName, n, bestTour, bestTourval );

}

void runMultipleExperiments() {

  /* runs all the (-NINS) experiments, i.e., as many intstances and on each instance it
     will run (-NS) steps of local search, (-NT) times
  */

  int inst;
  TEXPSTAT estat;
  int totalsteps;
  double totalevals, totaltime;

  totalsteps = 0;
  totalevals = 0;
  totaltime = 0;
  
  timer_start();
  prepareOrbits4OPT();
  prepareSuccPred(n);
  
  globalTimeStart = timer_elapsed_seconds();

  OUTPUTON = 0;

  for ( inst = 0; inst < NUMINSTANCES; inst++ ) {

    printf("* * * EXP %d / %d\n", inst + 1, NUMINSTANCES ); fflush(stdout);

    if ( !costReadFromFile ) 
      setCosts();
    else
      assert( NUMINSTANCES == 1 );
    
    perturbCostsOnNodes();  // applies (if needed) perturb c(i,j)+u(i)+u(j)

    runSingleExperiment(inst, &estat);

    totalsteps += estat.numsteps;
    totalevals += estat.totevals;
    totaltime += estat.tottime;

    // #ifdef COSTS_ONFLY
    SEEDCOST++;  // cambia il seed, o tutte le matrici dei costi, generate al run con un seed fisso
                 // sarebbero uguali
    // #endif
  
  }
  printf("\n\nEND MULTIPLE (%d) EXPERIMENTS\n", NUMINSTANCES);
  //  printf("n %d A2 %d A3 %d A4 %d SWCOEF %g SIMPCOEF %g\n", n, SEARCHTYPE[2], SEARCHTYPE[3],	 SEARCHTYPE[4], SWITCH2COEF, SWITCHIMPRCOEF );
  printf("n %d A2 %d findbest %d SWCOEF %g SIMPCOEF %g NN %d FILE <%s>\n", n, SEARCHTYPE[2], FINDBEST,
	 SWITCH2COEF, SWITCHIMPRCOEF, strcmp(startTourName,"NN")==0,
	 costFileName );
  // togliere poi
  
  printf("Tour val %6.2f (best) %6.2f (avg)\n\n", BESTTVAL, AVGTVAL );

  printf("Total steps executed       %d\n", totalsteps);
  printf("AVG step/exper             %6.2f\n", ((double)totalsteps) / (NUMINSTANCES*NUMTESTS) );
  printf("Total move evals           %10lf\n", totalevals);
  printf("AVG  move evals/step       %5.2f\n", totalevals / totalsteps);
  printf("Total time                 %5.2f\n", totaltime);
  printf("AVG  time/exper            %5.2f\n", totaltime / (NUMINSTANCES*NUMTESTS) );
  printf("AVG  time/step             %9.4f\n", totaltime / totalsteps);

}




void Usage() {

  printf("Usage:\n\n");
  printf("./4optLS [parameters]\n");
  printf("\n");
  printf("Not all parameters are mandatory. Can be in any order.\n");
  printf("-A4 <algo>  (determines the algorithm for 4opt moves)\n");
  printf("   <algo>   : 0 basic SF, 1 SF w/o split, 2 SF w/split, 3 woeginger, 4 glover, 5 brute force\n");
  printf("   default <algo> = -1 (no 4-opt)\n\n");
  printf("-A3 <algo>  (determines the algorithm for 3opt moves)\n");
  printf("   <algo>   : %d brute force, %d basic SF\n", SEARCH_3OPT_BF, SEARCH_3OPT_SF );
  printf("   default <algo> = -1 (no 3-opt)\n\n");
  printf("-A2 <algo>  (determines the algorithm for 2opt moves)\n");
  printf("   <algo>   : %d brute f, %d basic SF, %d SF w/o dupl, %d SF w/LEXI, %d SF w/pps, %d H(delta) pps\n",
	                SEARCH_2OPT_BF, SEARCH_2OPT_SF, SEARCH_2OPT_SF_NODUP,
	                SEARCH_2OPT_SF_LEXI, SEARCH_2OPT_SF_SORTPAIR, SEARCH_2OPT_H_SORTPAIR);
  printf("   default <algo> = -1 (no 2-opt)\n\n");
  printf("-LOOP <k-order>  (determines in which order to run the k opts)\n");
  printf("   <k-order>   : permutation of {2,3,4}\n");
  printf("   default <k-order> = 234 \n\n");
  printf("-O4 <orbits>  (determines the 4opt orbits in the form of a string over 1..7)\n");
  printf("   <orbits> : e.g. 2 (only orbit 2) 135 (orbits 1,3 and 5) 1234567 (all orbits)\n");
  printf("   default <orbits> = 1234567\n\n");
  printf("-RU <n>  (creates a random uniform input on n nodes)\n");
  printf("   default <n> = 100\n\n");
  printf("-RE <n>  (creates a random euclidean input on n nodes)\n");
  printf("   default <n> = 100\n\n");
  printf("-SEEDT <seed>  (fixes seed for random number generator used for tours)\n");
  printf("   default <seed> = 1234\n\n");
  printf("-SEEDC <seed>  (fixes seed for random number generator used for costs)\n");
  printf("   default <seed> = 1234\n\n");
  printf("-NS <numsteps>  (how many steps of local search)\n");
  printf("   default <numsteps> = INFINITE (up to local opt)\n\n");
  printf("-NT <numtests> (runs more local searches of the same type on each instance)\n");
  printf("   default <numtests> = 1\n\n");
  printf("-NINS <numinst> (how many instances of this type)\n");
  printf("   default <numinst> = 1\n\n");
  printf("-ST <starting tour file | NN>  (file tour to start from)\n");
  printf("   if NN then starts from nearest neighbor tour\n");
  printf("   default = random tour\n\n");
  printf("-P <doperturb>  (perturb each cost by an epsilon  < 1/1000)\n\n");
  printf("   default = do not perturb. N.B: must put this before -F \n\n");
  printf("-F <filename>  (reads input graph from file)\n\n");
  printf("-DOFINDBEST <find best move>  (1=yes 0=no)\n");
  printf("   default = find the best move\n\n");
  printf("-DOSHOWMOVE <show each move>  (1=yes 0=no)\n");
  printf("   default = show each move\n\n");
  printf("-SWITCH2COEF <val>  (move to 2OPTCE after <val>*n steps of LS)\n");
  printf("   default = INFINITE\n\n");
  printf("-SWITCHIMPRCOEF <val>  (swap B.I. <-> F.I. at step <val>*n of LS)\n");
  printf("   default = INFINITE\n\n");
  printf("-WT <root_filename>  (writes intermediate tours on file)\n");
  printf("   default = do not write\n\n");
  printf("-WLT <root_filename>  (writes last tour of each convergence on file)\n");
  printf("   default = do not write\n\n");
  printf("-Wmove <filename>  (writes log of number of moves on file)\n");
  printf("   default = do not write\n\n");
  printf("-Wtime <filename>  (writes log of time for move on file)\n");
  printf("   default = do not write\n\n");
  printf("-Wval <filename>  (writes log of tours value on file)\n");
  printf("   default = do not write\n\n");
  printf("-WOPT <filename>  (writes best tour found overall on file)\n");
  printf("   default = do not write\n\n");
  printf("-PRO <projected steps>  (estimates time for projected steps given avg time for steps run)\n");
  printf("   default = do not write\n\n");
  printf("-LTF <lap time from step>  (takes time/stats also for the steps ltf,ltf+1,....)\n");
  printf("   default = step 0\n\n");
  printf("-PREC <p>  (sets the decimal precision for costs to p digits\n");
  printf("   default = 6 digits\n\n");
  printf("-CMIN <lbval>  random costs in [CMIN,CMAX)\n");
  printf("   default = 0.0\n\n");
  printf("-CMAX <ubval>  random costs in [CMIN,CMAX)\n");
  printf("   default = 1.0\n\n");
  printf("-TLIM <sec>  (time limit for overall procedure, over multiple runs, in secs\n");
  printf("   default = no limit\n\n");
  printf("-ALTER <maxU>  (changes costs into c[i,j]+u[i]+u[j] for random u[] <= maxU\n");
  printf("   default = maxU 0\n\n");
  printf("-TSPOUT <filename> (writes out the read or generated graph in TSP format)\n");
  printf("   default = do not write\n\n");

  
  exit(1);

}

int main( int argc, char** argv ){

  int i, c;
  
  if ( argc == 1 || (argc % 2 == 0 ))
    Usage();

  // defaults : 
  SEEDTOUR       = 1234;
  SEEDCOST       = 1234;
  n              = 100; 
  norbits        = 7;
  NUMTESTS       = 1;
  NUMINSTANCES   = 1;
  NUM_LS_STEPS   = INFINITE;
  TIMELIMIT      = INFINITE;
  DOPERTURB      = false;
  PROJECTED      = -1; // undef
  COMPUTE_FROM_STEP = 0;
  MAX_NODE_U     = 0;
  PRECISION      = 6;
  CMIN           = 0.0;
  CMAX           = 1.0;
  FINDBEST       = true;
  SHOWMOVE       = false;
  SWITCH2COEF    = INFINITE;
  SWITCHIMPRCOEF = INFINITE;
  costReadFromFile = false;
  COSTI_EUCLIDEI = 0;  // per ora suppone costi uniformi
  for ( c = 1; c <= 7; c++ )
    orbitList[c-1] = c;
  for ( c = 0; c < 3; c++ )
    kloop[c] = 2 + c;  // default loop = "234"
  for ( i = 2; i <= 4; i++ )
    SEARCHTYPE[i] = -1;   // undef
  strcpy( startTourName, "" );
  strcpy( costFileName, "" );
  strcpy( outTourName, "" );
  strcpy( outLastTourName, "" );
  strcpy( outMoveName, "" );
  strcpy( outTimeName, "" );
  strcpy( outValName, "" );
  strcpy( outOptName, "" );
  strcpy( outTSPName, "" );

  // parse the command line
  c = 1;
  while ( c < argc ) {
    // printf("curr arg a[%d] %s\n", c, argv[c] );
    if ( strcmp( argv[c], "-A4") == 0 ) {
      SEARCHTYPE[4] = atoi( argv[c+1] );
      if ( SEARCHTYPE[4] < 0 || SEARCHTYPE[4] > SEARCH_4OPT_TOTOPTIONS - 1 ) {
	printf("Error in 4OPT algorithm specification\n");
	Usage();
      }
      c += 2;
    }
    else if ( strcmp( argv[c], "-A3") == 0 ) {
      SEARCHTYPE[3] = atoi( argv[c+1] );
      if ( SEARCHTYPE[3] < 0 || SEARCHTYPE[3] > SEARCH_3OPT_TOTOPTIONS - 1  ) {
	printf("Error in 3OPT algorithm specification\n");
	Usage();
      }
      c += 2;
    }
    else if ( strcmp( argv[c], "-A2") == 0 ) {
      SEARCHTYPE[2] = atoi( argv[c+1] );
      if ( SEARCHTYPE[2] < 0 || SEARCHTYPE[2] > SEARCH_2OPT_TOTOPTIONS - 1 ) {
	printf("Error in 2OPT algorithm specification\n");
	Usage();
      }
      c += 2;
    }
    else if ( strcmp( argv[c], "-LOOP") == 0 ) {
      if ( strlen( argv[c+1] ) != 3 ) {
	printf("Error in permutation of K-opts\n");
	Usage();
      }
      for ( i = 0; i < 3; i++ )
	kloop[i] = argv[c+1][i] - '0';
      c+=2;
    }
    else if ( strcmp( argv[c], "-O4") == 0 ) {
      norbits = strlen( argv[c+1] );
      if ( norbits > 7 ) {
	printf("Error in number of 4OPT orbits\n");
	Usage();
      }
      for ( i = 0; i < norbits; i++ )
	orbitList[i] = argv[c+1][i] - '0';
      c+=2;
    }
    else if ( strcmp( argv[c], "-O3") == 0 ) {
      printf("Orbits for 3OPT moves still to implement\n");
      exit(1);
    }
    else if ( strcmp( argv[c], "-RU") == 0 ) {
      n = atoi( argv[c+1] );
      if ( n >= MAXN ) {
	printf("Error n too large\n");
	Usage();
      }
      RANDOM_GRAPH_TYPE = UNIFORM;
      c+=2;
    }
    else if ( strcmp( argv[c], "-RE") == 0 ) {
      n = atoi( argv[c+1] );
      if ( n >= MAXN ) {
	printf("Error n too large\n");
	Usage();
      }
      RANDOM_GRAPH_TYPE = EUCLIDEAN;
      c+=2;
    }
    else if ( strcmp( argv[c], "-P") == 0 ) {
      DOPERTURB = atoi( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-PREC") == 0 ) {
      PRECISION = atoi( argv[c+1] );
      if ( PRECISION > 8 || PRECISION < 0 ) {
	printf("Precision must be in range 0..8\n");
	Usage();
      }
      c+=2;
    }
    else if ( strcmp( argv[c], "-CMIN") == 0 ) {
      CMIN = atof( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-CMAX") == 0 ) {
      CMAX = atof( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-F") == 0 ) {
      strcpy( costFileName, argv[c+1] );
      readCostsFile( costFileName, DOPERTURB );
      costReadFromFile = true;
      c+=2;
    }
    else if ( strcmp( argv[c], "-DOFINDBEST") == 0 ) {
      FINDBEST = atoi( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-DOSHOWMOVE") == 0 ) {
      SHOWMOVE = atoi( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-SWITCH2COEF") == 0 ) {
      SWITCH2COEF = atof( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-SWITCHIMPRCOEF") == 0 ) {
      SWITCHIMPRCOEF = atof( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-SEEDT") == 0 ) {
      SEEDTOUR = atoi( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-SEEDC") == 0 ) {
      SEEDCOST = atoi( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-NS") == 0 ) {
      NUM_LS_STEPS = atoi( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-NT") == 0 ) {
      NUMTESTS = atoi( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-NINS") == 0 ) {
      NUMINSTANCES = atoi( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-TLIM") == 0 ) {
      TIMELIMIT = atof( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-ST") == 0 ) {
      strcpy( startTourName, argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-WT") == 0 ) {
      strcpy( outTourName, argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-WLT") == 0 ) {
      strcpy( outLastTourName, argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-Wmove") == 0 ) {
      strcpy( outMoveName, argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-Wtime") == 0 ) {
      strcpy( outTimeName, argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-Wval") == 0 ) {
      strcpy( outValName, argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-WOPT") == 0 ) {
      strcpy( outOptName, argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-TSPOUT") == 0 ) {
      strcpy( outTSPName, argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-PRO") == 0 ) {
      PROJECTED = atoi( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-LTF") == 0 ) {
      COMPUTE_FROM_STEP = atoi( argv[c+1] );
      c+=2;
    }
    else if ( strcmp( argv[c], "-ALTER") == 0 ) {
      MAX_NODE_U = atof( argv[c+1] );
      c+=2;
    }
    else {
      printf("Unknown parameter %s on usage string\n", argv[c] );
      Usage();
    }
  }
  
  COSTI_EUCLIDEI = ( RANDOM_GRAPH_TYPE == EUCLIDEAN );
#ifdef COSTS_ONFLY
  assert( !costReadFromFile );  // se letti da file, sono senz'altro in matrice
#endif
  

  if ( c > argc ) {
    printf("Too many arguments on usage string\n");  
    Usage();
  }
  
  if ( COMPUTE_FROM_STEP > 0 && PROJECTED >= 0 ) {
    printf("Cannot use both options -LTF and -PRO together\n");
    exit(1);
  }

  srand( SEEDTOUR );
  srand48( SEEDCOST );
  
  runMultipleExperiments();

  return 0;
}


