#include <assert.h>
#include <stdio.h>
#include "tspmain.h"
#include "io4opt.h"
#include "proc4opt.h"

/**************************************************************************

  Procedures to speed-up the search of most improving move in OPT4 for TSP

****************************************************************************/

// movesInOrbit[i][j] = j-th scheme of orbit i, listed in the order of
// the group action and numbered 1..25 as in paper.
// i = 1..7   j = 0..size(orbit_i)-1

int GOODQ = 0;

const  int movesInOrbit[8][10] = {
  {},
  { 1, 24, 23, 22 }, 
  { 2, 21 },
  { 3, 7, 13, 17 }, 
  { 4, 19, 11, 6 }, 
  { 5, 20, 14, 15, 12, 18, 9, 8 },
  { 10, 16 }, 
  { 25 }};


long int EVAL4COUNTER = 0; /* global variable, used as semaphor, to count how many heap evaluations
                              of 4opt moves there have been in a certain stage of the algorithm. It gets reset and 
                              tested from outside, and it is increased at each evaluation */

/* global heaps used by our smart force procedures */

THEAP heapSingle;  // heapSingle è usato da SF0 che funziona come 3opt

THEAP heapF12;  // heapF12 è usato da SF1 e SF2 per mettere le coppie in ordine di
// f1(p1(1),pi(2))+f2(pi(3),pi(4))

THEAP heapF34;  // heapF12 è usato da SF1 e SF2 per mettere le coppie in ordine di
// f1(p1(1),pi(2))+f2(pi(3),pi(4))

TPAIR* sortedheap[2]; // usati per trasformare gli heapF12 e heapF34 in vettori ordinati
// nella procedura master-servant di SF1, SF2

TPAIR* tmpHeap; // heap temporaneo usato per duplicazione in SF_ANALOGO3OPT

/* global variables for Glover's algorithm */

#define TYPEC  0
#define TYPED  1

#ifdef IMPLEMENTAZIONE_GLOVER_COMEPAPER

#define MAXNGLOVER 1000

/*
static double reach_c[MAXNGLOVER][MAXNGLOVER];
static double reach_d[MAXNGLOVER][MAXNGLOVER];
static double cross_c[MAXNGLOVER][MAXNGLOVER];
static double cross_d[MAXNGLOVER][MAXNGLOVER];
*/
static double bestreach[2][MAXNGLOVER][MAXNGLOVER];
static double bestcross[2][MAXNGLOVER][MAXNGLOVER];

static int predreach[2][MAXNGLOVER][MAXNGLOVER];
static int predcross[2][MAXNGLOVER][MAXNGLOVER][2];

#endif

double best_t;
int type[2];
int pred_t1[2];
int pred_t2[2];

double * memoizeA[2][MAXN];
double * memoizeB[2][MAXN];
int * memoizeBestA[2][MAXN]; // back pointer to best
int * memoizeBestB[2][2][MAXN]; // back pointers

int SWITCH_GLOVER_COST = 0; // if =1 then replace cost_c <-> cost_d in the recursions

void allocateGloverMem() {

  int i, j, k;

  /*
  printf("Allocate memory Glover for n %d\n", n);
  fflush( stdout );
  */
  
  assert ( n >= 0 && n < MAXN ) ;
  for ( i = 0; i < n; i++ )
    for ( j = 0; j < 2; j ++ ) {
      //      printf("***ALLOCA mamoizeA[%d][%d] : ", j, i);
      memoizeA[j][i] = (double*) mycalloc( n + 2, sizeof(double) );
      assert( memoizeA[j][i] );
      //      printf("***ALLOCA memoizeB[%d][%d] : ", j, i);
      memoizeB[j][i] = (double*) mycalloc( n + 2, sizeof(double) );
      assert( memoizeB[j][i] );
      //      printf("***ALLOCA memoizeBestA[%d][%d] : ", j, i);
      memoizeBestA[j][i] = (int*) mycalloc( n + 2, sizeof(int) );
      assert( memoizeA[j][i] );
      for ( k = 0; k < 2; k++ ) {
	//	printf("***ALLOCA mamoizeBestB[%d][%d][%d] : ", k, j, i);
	memoizeBestB[k][j][i] = (int*) mycalloc( n + 2, sizeof(int) );
	assert( memoizeBestB[k][j][i] );
      }
  }
}

void disposeGloverMem() {

  int i, j, k;

  /*
  printf("Release memory Glover\n");
  fflush( stdout );
  */
  
  for ( i = 0; i < n; i++ )
    for ( j = 0; j < 2; j ++ ) {
      //      printf("*** FREE memoizeA[%d][%d] : ", j, i);
      myfree ( memoizeA[j][i] );
      //      printf("*** FREE memoizeB[%d][%d] : ", j, i);
      myfree ( memoizeB[j][i] );
      //      printf("*** FREE memoizeBestA[%d][%d] : ", j, i);
      myfree ( memoizeBestA[j][i] );
      for ( k = 0; k < 2; k++ ) {
	//	printf("*** FREE memoizeBestB[%d][%d][%d] : ", k, j, i);
	myfree ( memoizeBestB[k][j][i] );
      }
    }
}

/*********************************************************
 *********************************************************

  Functions used as keys to build the heaps for various
  type of moves 

  *****************************************************
******************************************************/


double tauplus( int* pi, int x, int y ) { 
  //  returns the value tau^+(x,y)

  return cost[pi[x]][pi[SUCC(x)]] - cost[pi[x]][pi[SUCC(y)]];
}

double tauminus( int* pi, int x, int y ) {  
  //  returns the value tau^-(x,y)

  return cost[pi[x]][pi[PRED(x)]] - cost[pi[x]][pi[PRED(y)]];
}

double tauplus_00( int* pi, int x, int y ) { 
  //  returns the value tau^+(x,y)

  return cost[pi[x]][pi[SUCC(x)]] - cost[pi[x]][pi[SUCC(y)]];
}

double tauplus_01( int* pi, int x, int y ) { 
  //  returns the value tau^+(x,y-1)

  return cost[pi[x]][pi[SUCC(x)]] - cost[pi[x]][pi[y]];
}

double tauminus_11( int* pi, int x, int y ) {  
  //  returns the value tau^-(x+1,y+1)

  return cost[pi[SUCC(x)]][pi[x]] - cost[pi[SUCC(x)]][pi[y]];
}

double tauminus_12( int* pi, int x, int y ) {  
  //  returns the value tau^-(x+1,y+2)

  //  printf("TM_12 of %d %d is c(%d %d) - c(%d %d) is %g\n", x, y, pi[SUCC(x)], pi[x], pi[SUCC(x)], pi[SUCC(y)], cost[pi[SUCC(x)]][pi[x]] - cost[pi[SUCC(x)]][pi[SUCC(y)]]);
  return cost[pi[SUCC(x)]][pi[x]] - cost[pi[SUCC(x)]][pi[SUCC(y)]];
}

double tauplusInv( int* pi, int x, int y ) { 
  //  returns the value tau^+(y,x)

  return cost[pi[y]][pi[SUCC(y)]] - cost[pi[y]][pi[SUCC(x)]];
}

double ggA( int* pi, int x, int y ) { 

  // the function g^1() (and also g^2()) for the scheme r_25
  // the function g^2 for the scheme r_10

  return tauplus( pi, x, y ) + tauplus( pi, y, x );
}


double ggB( int* pi, int x, int y ) { 

  // the function g^1() (and also g^2()) for the scheme r_2
  // the function g^1 for the scheme r_10

  return tauplus_01( pi, x, y ) + tauminus_12( pi, y, x );
}


/********************************************************
 ********************************************************

  Defines the orbits for 4opt reinsertion schemes

  ******************************************************
  ******************************************************/


double evaluate4OPT( int* pi, int i, int j, int k, int h, int numorbit, int gelem );
double generalCaseOfAny4OPTMoveAnalogo3byOrbit( int* pi, int* bestsel, double* currBestVal,
						THEAP* theap, int numOrbit, int nphi );
double find4OPTWoegingerByOrbit( int* pi, int * sel, int* rot, int numOrbit, OPTMOVE* move );
double find4OPTGloverByOrbit( int* pi, int * sel, int* rot, int numOrbit, OPTMOVE* move );


void prepareOrbits4OPT() {


  int rho[4][4]; 
  int psirho[4][4]; 

  int rot, i;

  for ( rot = 0; rot < 4; rot++ )
    for ( i = 0; i < 4; i++ ) {
      rho[rot][i] = (i + rot) % 4;
      psirho[rot][i] = (6 - rot - i) % 4;
    }


  // 1st orbit. Scheme <+1,-2,-3,-4>

  ORBIT4[0].size  = 4;
  for ( rot = 0; rot < ORBIT4[0].size; rot++ ) {
    memcpy (ORBIT4[0].phi[rot], rho[rot], 4 * sizeof(int));
    ORBIT4[0].flip[rot] = false;
  }

  ORBIT4[0].fct[0] = &tauplus_01;
  ORBIT4[0].who[0][0] = 0;
  ORBIT4[0].who[0][1] = 1;
  ORBIT4[0].phiHeapCode[0] = 1;  // heap phi1

  ORBIT4[0].fct[1] = &tauminus_12;
  ORBIT4[0].who[1][0] = 3;
  ORBIT4[0].who[1][1] = 2;
  ORBIT4[0].phiHeapCode[1] = 17;  

  ORBIT4[0].fct[2] = &tauplus_00;
  ORBIT4[0].who[2][0] = 2;
  ORBIT4[0].who[2][1] = 0;
  ORBIT4[0].phiHeapCode[2] = 9;

  ORBIT4[0].fct[3] = &tauminus_11;
  ORBIT4[0].who[3][0] = 1;
  ORBIT4[0].who[3][1] = 3;
  ORBIT4[0].phiHeapCode[3] = 6;

  ORBIT4[0].g1 = NULL;
  ORBIT4[0].g2 = NULL;

  // 2nd orbit. Scheme <+1,-2,+3,-4>

  ORBIT4[1].size  = 2;
  for ( rot = 0; rot < ORBIT4[1].size; rot++ ) {
    memcpy (ORBIT4[1].phi[rot], rho[rot], 4 * sizeof(int));
    ORBIT4[1].flip[rot] = false;
  }

  ORBIT4[1].fct[0] = &tauplus_01;
  ORBIT4[1].who[0][0] = 0;
  ORBIT4[1].who[0][1] = 1;
  ORBIT4[1].phiHeapCode[0] = 1;

  ORBIT4[1].fct[1] = &tauplus_01;
  ORBIT4[1].who[1][0] = 2;
  ORBIT4[1].who[1][1] = 3;
  ORBIT4[1].phiHeapCode[1] = 10;

  ORBIT4[1].fct[2] = &tauminus_12;
  ORBIT4[1].who[2][0] = 1;
  ORBIT4[1].who[2][1] = 0;
  ORBIT4[1].phiHeapCode[2] = 5;

  ORBIT4[1].fct[3] = &tauminus_12;
  ORBIT4[1].who[3][0] = 3;
  ORBIT4[1].who[3][1] = 2;
  ORBIT4[1].phiHeapCode[3] = 17;

  ORBIT4[1].g1 = &ggB;
  ORBIT4[1].g2 = &ggB;


  // 3rd orbit. Scheme <+1,-2,+3,-4>

  ORBIT4[2].size  = 4;
  for ( rot = 0; rot < ORBIT4[2].size; rot++ ) {
    memcpy (ORBIT4[2].phi[rot], rho[rot], 4 * sizeof(int));
    ORBIT4[2].flip[rot] = false;
  }

  ORBIT4[2].fct[0] = &tauplus_01;
  ORBIT4[2].who[0][0] = 0;
  ORBIT4[2].who[0][1] = 1;
  ORBIT4[2].phiHeapCode[0] = 1;

  ORBIT4[2].fct[1] = &tauplus_00;
  ORBIT4[2].who[1][0] = 2;
  ORBIT4[2].who[1][1] = 3;
  ORBIT4[2].phiHeapCode[1] = 11;

  ORBIT4[2].fct[2] = &tauplus_00;
  ORBIT4[2].who[2][0] = 3;
  ORBIT4[2].who[2][1] = 0;
  ORBIT4[2].phiHeapCode[2] = 13;

  ORBIT4[2].fct[3] = &tauminus_12;
  ORBIT4[2].who[3][0] = 1;
  ORBIT4[2].who[3][1] = 2;
  ORBIT4[2].phiHeapCode[3] = 8;

  ORBIT4[2].g1 = NULL;
  ORBIT4[2].g2 = NULL;


  // 4th orbit. Scheme <+1,-2,+3,-4>

  ORBIT4[3].size  = 4;
  for ( rot = 0; rot < 2; rot++ ) {
    memcpy (ORBIT4[3].phi[rot], rho[rot], 4 * sizeof(int));
    ORBIT4[3].flip[rot] = false;
  }
  for ( rot = 0; rot < 2; rot++ ) {
    memcpy (ORBIT4[3].phi[2 + rot], psirho[rot], 4 * sizeof(int));
    ORBIT4[3].flip[2 + rot] = true;
  }

  ORBIT4[3].fct[0] = &tauplus_01;
  ORBIT4[3].who[0][0] = 0;
  ORBIT4[3].who[0][1] = 1;
  ORBIT4[3].phiHeapCode[0] = 1;

  ORBIT4[3].fct[1] = &tauplus_01;
  ORBIT4[3].who[1][0] = 3;
  ORBIT4[3].who[1][1] = 2;
  ORBIT4[3].phiHeapCode[1] = 15;

  ORBIT4[3].fct[2] = &tauminus_12;
  ORBIT4[3].who[2][0] = 2;
  ORBIT4[3].who[2][1] = 0;
  ORBIT4[3].phiHeapCode[2] = 12;

  ORBIT4[3].fct[3] = &tauminus_12;
  ORBIT4[3].who[3][0] = 1;
  ORBIT4[3].who[3][1] = 3;
  ORBIT4[3].phiHeapCode[3] = 7;

  ORBIT4[3].g1 = NULL;
  ORBIT4[3].g2 = NULL;

  // 5th orbit. Scheme <+1,-2,+4,+3>

  ORBIT4[4].size = 8; 

  for ( rot = 0; rot < 4; rot++ ) {
    memcpy (ORBIT4[4].phi[rot], rho[rot], 4 * sizeof(int));
    ORBIT4[4].flip[rot] = false;
  }
  for ( rot = 0; rot < 4; rot++ ) {
    memcpy (ORBIT4[4].phi[4 + rot], psirho[rot], 4 * sizeof(int));
    ORBIT4[4].flip[4 + rot] = true;
  }

  ORBIT4[4].fct[0] = &tauplus_01;
  ORBIT4[4].who[0][0] = 0;
  ORBIT4[4].who[0][1] = 1;
  ORBIT4[4].phiHeapCode[0] = 1;

  ORBIT4[4].fct[1] = &tauminus_11;
  ORBIT4[4].who[1][0] = 3;
  ORBIT4[4].who[1][1] = 2;
  ORBIT4[4].phiHeapCode[1] = 16;

  ORBIT4[4].fct[2] = &tauminus_12;
  ORBIT4[4].who[2][0] = 2;
  ORBIT4[4].who[2][1] = 0;
  ORBIT4[4].phiHeapCode[2] = 12;

  ORBIT4[4].fct[3] = &tauminus_11;
  ORBIT4[4].who[3][0] = 1;
  ORBIT4[4].who[3][1] = 3;
  ORBIT4[4].phiHeapCode[3] = 6;

  ORBIT4[4].g1 = NULL;
  ORBIT4[4].g2 = NULL;

  // 6th orbit. Scheme <+1,-3,-4,+2>

  ORBIT4[5].size  = 2;
  for ( rot = 0; rot < ORBIT4[5].size; rot++ ) {
    memcpy (ORBIT4[5].phi[rot], rho[rot], 4 * sizeof(int));
    ORBIT4[5].flip[rot] = false;
  }

  ORBIT4[5].fct[0] = &tauplus_01;
  ORBIT4[5].who[0][0] = 0;
  ORBIT4[5].who[0][1] = 2;
  ORBIT4[5].phiHeapCode[0] = 2;

  ORBIT4[5].fct[1] = &tauplus_00;
  ORBIT4[5].who[1][0] = 1;
  ORBIT4[5].who[1][1] = 3;
  ORBIT4[5].phiHeapCode[1] = 4;

  ORBIT4[5].fct[2] = &tauminus_12;
  ORBIT4[5].who[2][0] = 2;
  ORBIT4[5].who[2][1] = 0;
  ORBIT4[5].phiHeapCode[2] = 12;

  ORBIT4[5].fct[3] = &tauplus_00;
  ORBIT4[5].who[3][0] = 3;
  ORBIT4[5].who[3][1] = 1;
  ORBIT4[5].phiHeapCode[3] = 14;

  ORBIT4[5].g1 = &ggB;
  ORBIT4[5].g2 = &ggA;


  // 7th orbit. Scheme <+1,-3,-4,+2>

  ORBIT4[6].size  = 1;
  for ( rot = 0; rot < ORBIT4[6].size; rot++ ) {
    memcpy (ORBIT4[6].phi[rot], rho[rot], 4 * sizeof(int));
    ORBIT4[6].flip[rot] = false;
  }

  ORBIT4[6].fct[0] = &tauplus_00;
  ORBIT4[6].who[0][0] = 0;
  ORBIT4[6].who[0][1] = 2;
  ORBIT4[6].phiHeapCode[0] = 3;

  ORBIT4[6].fct[1] = &tauplus_00;
  ORBIT4[6].who[1][0] = 1;
  ORBIT4[6].who[1][1] = 3;
  ORBIT4[6].phiHeapCode[0] = 4;

  ORBIT4[6].fct[2] = &tauplus_00;
  ORBIT4[6].who[2][0] = 2;
  ORBIT4[6].who[2][1] = 0;
  ORBIT4[6].phiHeapCode[0] = 9;

  ORBIT4[6].fct[3] = &tauplus_00;
  ORBIT4[6].who[3][0] = 3;
  ORBIT4[6].who[3][1] = 1;
  ORBIT4[6].phiHeapCode[0] = 14;

  ORBIT4[6].g1 = &ggA;
  ORBIT4[6].g2 = &ggA;


  /*
  int j, h;
  for ( j = 0; j < 7; j++ ) {
    printf("Orbit %d size %d\n", j, ORBIT4[j].size);
    for ( rot = 0; rot < ORBIT4[j].size; rot++ ) {
      printf("Elem %d (flip %d) is ", rot, ORBIT4[j].flip[rot]);
      for ( h = 0; h < 4 ; h++ )
	printf("%d ", ORBIT4[j].phi[rot][h]);
      printf("\n");
    }
    printf("\n");
  }
  */
}


int removeReflection( int norb, int s )
{

  // norb is the number of an orbit, and s is the number of an element in  the orbit
  // if the orbit has also reflections, half elements are w/o reflections (i.e., rotations of
  // the representative) and half are w/reflections. With the trick of reversing pi[], the reflections
  // correspond to the moves of rotations, but for a different rotation code than s - orb(size)/2.
  // This procedure computes such rotation code

  if ( !ORBIT4[norb].flip[s] )
    return s;
  
  else switch (norb) {
    case 3:  // orbit no.4
      switch( s ) {
      case 0 ... 1:
	return s;
      case 2 ... 3:
	return 2 + (3 - s);
      }

    case 4:  // orbit no.5
      switch( s ) {
      case 0 ... 3:
	return s;
      case 4:
	return 7;
      case 5:
	return 4;
      case 6:
	return 5;
      case 7:
	return 6;
	/*
      case 4 ... 7:
	return 4 + ( 7 - s );
	*/
      }
    }

  assert(false);
}

void reflectSelection( int* sel ) {
  // given a selection, it reorgainizes it in a specular way, i.e.
  // as if it is read counterclockwise.
  // So, the first element becomes the last, ..., the last becomes
  // the first, as in    (i,j,k,h)<->(h,k,j,i)

  int i;

  for ( i = 0; i < 2; i++ ) {
    int tmp = sel[i];
    sel[i] = n - 1 - sel[3-i];
    sel[3-i] = n - 1 - tmp;
  }
}

int smalllex(int* sel, int i, int j, int k, int h ) {

  //*** needed for unique best move... remove after tests
  // cheks if (i,j,k,h) prec (s[0],s[1], s[2], s[3] )

  return false;  /* set this if finding any best move is ok */

  // continue here if we need the smallest lex order best move, i.e., unique

  // ***** c'è da controllare il codice dove usa il lex order perchè A2 A0 non funz...
  
  return
    (i < sel[0])  ||
    ( i == sel[0] && j < sel[1] ) ||
    ( i == sel[0] && j == sel[1] && k < sel[2] ) ||
    ( i == sel[0] && j == sel[1] && k == sel[2] && h < sel[3] ) ;
}
   

/*******************************************************
 *******************************************************

  Brute-force procedures

  ******************************************************
  ******************************************************/

double find4OPTMoveBFbyOrbit( double VAL2BEAT, int* pi, int* sel, int* rot, int numOrbit, OPTMOVE* move ) { 
  /* given permutation pi[], finds an improving (true) 4OPT move, within the orbit <numOrbit>
     i.e., 4 indices sel[0] < sel[1] < sel[2] < sel[3] defining a move, 
     plus the id of the group action element <*rot> corresponding to the best move in the orbit

     Here <orbitNum> goes from 1 to 7, with: 

(orbit 1) (-2, -3, -4) 

    edges (i,j) (i+1,k) (j+1,h) (k+1,h+1)

(orbit 2) (-2, +3, -4) 

    edges (i,j) (i+1,j+1) (k,h) (k+1,h+1)

- - - - -

  only moves that are better than VAL2BEAT are sought after

  if ( findbest ) then it returns quadruple such as delta  is as large as possible,
  otherwise it returns the first quadruple for which delta > VALUE2BEAT

  If delta <= VAL2BEAT, then we are at a local optimum. 
  It returns delta.

  This version uses Bruteforce, i.e. tries all quadruples.
  The true 4-OPT moves are such that no segment degenerates in a single point, i.e., 
  it must be j - i >=2 and k - j >= 2 and h - k >= 2

  */

  double delta, vv;
  int s, actuals, norb;
  int i, j, k, h;
  int revpi[MAXN];   // pi reversed
  int* actualpi[8];

  assert( numOrbit >= 1 && numOrbit <= 7 );
  norb = numOrbit - 1;   // seven orbits, 0..6 in c starts at 0
  *rot = -1; // undef

  if ( numOrbit == 4 || numOrbit == 5 ) // orbits w/reflection
    for ( i = 0; i <= n; i++ )
      revpi[i] = pi[n-i];

  for ( s = 0; s < ORBIT4[norb].size; s++ )
    if ( !ORBIT4[norb].flip[s] )
      actualpi[s] = pi;
    else
      actualpi[s] = revpi + 0;
  
  delta = VAL2BEAT;
  EVAL4COUNTER = 0;

  for ( s = 0; s < ORBIT4[norb].size; s++ ) {
    actuals = removeReflection( norb, s );
    dprintf("orb %d s %d actual %d\n", numOrbit, s, actuals );
    
    for ( i = 0; i <= n - 7; i++ ) {
      for ( j = i + 2; j <= n - 5; j++ ) {
	for ( k = j + 2; k <= n - 3; k++ ) {
	  for ( h = k + 2; h <= n - 1 - (i==0); h++ ) 
	    {
	      vv = evaluate4OPT(actualpi[s], i, j, k, h, numOrbit, s );
	      if ( (FINDBEST && vv > delta + EPS) || (!FINDBEST && vv > EPS) ||
		   ( FINDBEST && vv > delta - EPS && smalllex( sel, i, j, k, h )) ) {
		  delta = vv;
		  sel[0] = i;
		  sel[1] = j;
		  sel[2] = k;
		  sel[3] = h;
		  *rot = actuals;
		  		  
		  if ( !FINDBEST ) 
		    goto exitloop;
	      }
	    }
	}
      }
    }
  }

 exitloop:
  
  if ( *rot != -1 && ORBIT4[norb].flip[*rot] )
    reflectSelection( sel );

  dprintf("BF (OR%d G%d) i %4d  j %4d  k %4d  h %4d  val %g\n", numOrbit, *rot, sel[0], sel[1], sel[2], sel[3], delta );

  move->found = (delta > VAL2BEAT + EPS);
  move->delta = delta;
  for ( i = 0; i < 4; i++ )
    move->selection[i] = sel[i];

  move->maxHeapSize = 0;
  move->orbitNo = numOrbit;
  move->phi = *rot;
  move->heapExtractions = 0;
  move->effectiveExtractions = 0;
  move->moveEvaluations = EVAL4COUNTER;
  if ( *rot != -1 )
    move->scheme = movesInOrbit[numOrbit][*rot];
  else
    move->scheme = -1; // undef
  return delta;
}


double evaluate4OPT( int* pi, int i, int j, int k, int h, int numorbit, int gelem ) {

  // computes the value of selection sel[] for the move corresponding
  // to the <gelem> scheme in numorbit

  double cost_out, cost_in = 0;

  EVAL4COUNTER++;
  assert( numorbit >= 1 && numorbit <= 7 );

  cost_out = cost[pi[i]][pi[SUCC(i)]] + cost[pi[j]][pi[SUCC(j)]] +
             cost[pi[k]][pi[SUCC(k)]] + cost[pi[h]][pi[SUCC(h)]];

  switch ( numorbit ) {
  case 1:  
    switch ( gelem ) {
    case 0:      //      +1 -2 -3 -4
      cost_in = cost[pi[i]][pi[j]] + cost[pi[SUCC(k)]][pi[SUCC(h)]] +
	cost[pi[SUCC(j)]][pi[h]] + cost[pi[SUCC(i)]][pi[k]];
      break;

    case 1:    //       +1 +4 +3 -2
      cost_in = cost[pi[i]][pi[SUCC(k)]] + cost[pi[j]][pi[k]] +
	cost[pi[SUCC(i)]][pi[SUCC(h)]] + cost[pi[SUCC(j)]][pi[h]] ;
      break;

    case 2:    //       +1 +4 -3 +2
      cost_in = cost[pi[i]][pi[SUCC(k)]] + cost[pi[h]][pi[k]] +
	cost[pi[SUCC(i)]][pi[SUCC(j)]] + cost[pi[SUCC(h)]][pi[j]] ;
      break;
      
    case 3:    //       +1 -4 +3 +2
      cost_in = cost[pi[i]][pi[h]] + cost[pi[SUCC(i)]][pi[k]] +
	cost[pi[SUCC(k)]][pi[SUCC(j)]] + cost[pi[SUCC(h)]][pi[j]] ;
      break;

    default:
      assert(0);
    }
    break;


  case 2:  
    switch ( gelem ) {
    case 0:      //      +1 -2 +3 -4
      cost_in = cost[pi[i]][pi[j]] + cost[pi[k]][pi[h]] +
	cost[pi[SUCC(j)]][pi[SUCC(i)]] + cost[pi[SUCC(k)]][pi[SUCC(h)]];
      break;

    case 1:    //      +1 -4 +3 -2
      cost_in = cost[pi[i]][pi[h]] + cost[pi[j]][pi[k]] +
	cost[pi[SUCC(k)]][pi[SUCC(j)]] + cost[pi[SUCC(h)]][pi[SUCC(i)]] ;
      break;

    default:
      assert(0);
    }
    break;

  case 3:  
    switch ( gelem ) {
    case 0:      //      +1 -2 -4 +3
      cost_in = cost[pi[i]][pi[j]] + cost[pi[k]][pi[SUCC(h)]] +
	cost[pi[h]][pi[SUCC(i)]] + cost[pi[SUCC(j)]][pi[SUCC(k)]];
      break;

    case 1:    //      +1 +3 -2 -4
      cost_in = cost[pi[i]][pi[SUCC(j)]] + cost[pi[SUCC(k)]][pi[SUCC(h)]]
	+ cost[pi[h]][pi[SUCC(i)]] + cost[pi[j]][pi[k]] ;
      break;

    case 2:    //      +1 +3 -4 -2
      cost_in = cost[pi[k]][pi[h]] + cost[pi[j]][pi[SUCC(k)]] +
	cost[pi[i]][pi[SUCC(j)]] + cost[pi[SUCC(h)]][pi[SUCC(i)]] ;
      break;

    case 3:    //      +1 -4 -2 +3
      cost_in = cost[pi[i]][pi[h]] + cost[pi[j]][pi[SUCC(k)]] +
	cost[pi[k]][pi[SUCC(h)]] + cost[pi[SUCC(i)]][pi[SUCC(j)]] ;
      break;

    default:
      assert(0);
    }
    break;

  case 4:  
    switch ( gelem ) {
    case 0:      //      +1 -2 +4 -3
    case 2:
      cost_in = cost[pi[i]][pi[j]] + cost[pi[SUCC(i)]][pi[SUCC(k)]] +
	cost[pi[SUCC(j)]][pi[SUCC(h)]] + cost[pi[k]][pi[h]];
      dprintf("case2, edges in (%d %d) (%d %d) (%d %d) (%d %d) val %g\n", 
	     pi[i],pi[j],pi[SUCC(i)],pi[SUCC(k)],pi[SUCC(j)],pi[SUCC(h)],pi[k],pi[h],	cost_in);
      break;

    case 1:      //      +1 -4 +2 -3
    case 3:
      cost_in = cost[pi[i]][pi[h]] + cost[pi[j]][pi[k]] +
	cost[pi[SUCC(k)]][pi[SUCC(i)]] + cost[pi[SUCC(h)]][pi[SUCC(j)]] ;
      dprintf("case3, edges in (%d %d) (%d %d) (%d %d) (%d %d) val %g\n", 
	     pi[i],pi[h],pi[SUCC(k)],pi[SUCC(i)],pi[SUCC(h)],pi[SUCC(j)],pi[j],pi[k],	cost_in);
      break;

    default:
      assert(0);
    }
    break;

  case 5:
    switch ( gelem ) {
    case 0:    //      +1 -2 +4 +3
    case 4:    //      +1 -2 +4 +3
      cost_in = cost[pi[i]][pi[j]] + cost[pi[k]][pi[SUCC(h)]] +
	cost[pi[SUCC(i)]][pi[SUCC(k)]] + cost[pi[SUCC(j)]][pi[h]] ;
      break;
    case 1:    //      +1 +4 +2 -3
    case 5:    //      +1 +4 +2 -3
      cost_in = cost[pi[k]][pi[j]] + cost[pi[h]][pi[SUCC(i)]] +
	cost[pi[SUCC(k)]][pi[i]] + cost[pi[SUCC(j)]][pi[SUCC(h)]] ;
      break;
    case 2:    //      +1 +3 -4 +2
    case 6:    //      +1 +3 -4 +2
      cost_in = cost[pi[i]][pi[SUCC(j)]] + cost[pi[h]][pi[k]] +
	cost[pi[j]][pi[SUCC(h)]] + cost[pi[SUCC(k)]][pi[SUCC(i)]] ;
      break;
    case 3:    //      +1 -4 -2 -3
    case 7:    //      +1 -4 -2 -3
      cost_in = cost[pi[i]][pi[h]] + cost[pi[j]][pi[SUCC(k)]] +
	cost[pi[k]][pi[SUCC(i)]] + cost[pi[SUCC(h)]][pi[SUCC(j)]] ;
      break;
    }
    break;
    
  case 6:  
    switch ( gelem ) {
    case 0:      //      +1 -3 -4 +2
      cost_in = cost[pi[k]][pi[i]] + cost[pi[j]][pi[SUCC(h)]] +
	cost[pi[SUCC(k)]][pi[SUCC(i)]] + cost[pi[h]][pi[SUCC(j)]] ;
      break;

    case 1:     //      +1 +4 -2 -3
      cost_in = cost[pi[j]][pi[h]] + cost[pi[i]][pi[SUCC(k)]] +
	cost[pi[k]][pi[SUCC(i)]] + cost[pi[SUCC(h)]][pi[SUCC(j)]] ;
      break;

    default:
      assert(0);
    }
    break;

  case 7:    //     +1 +4 +3 +2
    assert( gelem == 0 );
    cost_in = cost[pi[i]][pi[SUCC(k)]] + cost[pi[j]][pi[SUCC(h)]] +
      cost[pi[k]][pi[SUCC(i)]] + cost[pi[SUCC(j)]][pi[h]] ;
    break;


  default:
    assert(false);
    break;
  }
  /*
  printf("...");
  scanf ("%d", &i);
  */

  dprintf("enter bis i %d  j %d  k %d  h %d  orb %d gelm %d  moveval %g\n", i, j, k, h,
	  numorbit, gelem, cost_out-cost_in );
  //  showPerm(pi, n+1);
  dprintf("edges out (%d %d) (%d %d) (%d %d) (%d %d) val %g\n", 
	 pi[i], pi[SUCC(i)], pi[j],pi[SUCC(j)], pi[k],pi[SUCC(k)],pi[h],pi[SUCC(h)], cost_out);
  

  return cost_out - cost_in;
}


/*************************************************

  SMART FORCE : 
  version analogous to 3OPT, with 4 phases

***********************************************/

double find4OPTMoveSFAnalogo3byOrbit( double VAL2BEAT, int* pi, int* sel, int* nphi, int nOrbit, OPTMOVE* move ) 
{ 
  /* given permutation pi[], finds an improving 4OPT move of the orbit nOrbit, i.e.
  sel[0] < sel[1] < sel[2] < sel[3] such that delta:= c(out) - c(in)
  is as large as possible.
  If it is non-positive, then we are at a local optimum. 
  It returns delta, if there exists delta > currBest, and it sets <*nphi> to the best group elem in 
  the orbit
  Else, it returns -1

  If FINDBEST = true, finds the best 4-opt move, else the first improving move

  Follows an algorithm similar to the 3OPT case in which we have 4 phases, one per heap

  <nOrbit> goes from 1 to 7

  returns: <move> with proper fields set, such as
  <heapsize> = size of the heaps built (total)
  <secondsToBuildHeap> = time spent in heap construction etc

  */

  
  //  THEAP HEAP;
  double delta, predelta;
  double tbefore,  time2build, time2use;
  int fs, s, i, phase, ndx1, ndx2;
  TORBIT* orbit;
  int revpi[MAXN], *actualpi;
  int hsize = 0;

  //  printf("find4optAnalogo3byOrbit VAL2BEAT %g nOrbit %d\n", VAL2BEAT, nOrbit);

  assert( nOrbit >= 1 && nOrbit <= 7 );

  hsize = 0;
  time2build = time2use = 0;
  orbit = ORBIT4 + (nOrbit - 1); // in c numbered from 0..6 rather than 1..7

  if ( nOrbit == 4 || nOrbit == 5 ) // orbits w/reflection
    for ( i = 0; i <= n; i++ )
      revpi[i] = pi[n-i];

  delta = VAL2BEAT;
  POPCOUNTER = 0;
  EVAL4COUNTER = 0;
  double tbeforeUse = timer_elapsed_seconds();
  for ( fs = 0; fs < orbit->size; fs++ ) {
    s =  ( orbit->flip[fs] ? fs - orbit->size / 2 : fs );
    actualpi =  ( orbit->flip[fs] ? revpi + 0 : pi );
    predelta = delta;
    for ( phase = 0; phase < 4; phase++ ) {

      ndx1 = orbit->who[phase][0];
      ndx2 = orbit->who[phase][1];

      dprintf("Fase %d crea heap per %d %d\n", phase, orbit->phi[s][ndx1], orbit->phi[s][ndx2]);
      tbefore = timer_elapsed_seconds();
      initializeHeap( false, &heapSingle, actualpi, orbit->phi[s][ndx1], orbit->phi[s][ndx2], false, orbit->fct[phase], false, BIGGER, delta / 4);
      time2build += timer_elapsed_seconds() - tbefore;
      hsize += heapsize( heapSingle.heap );
      dprintf("cerca la mossa  su orbita %d nphi %d\n", nOrbit, fs);
      tbefore = timer_elapsed_seconds();
      generalCaseOfAny4OPTMoveAnalogo3byOrbit( actualpi, sel, &delta, &heapSingle, nOrbit, s );
      if ( delta > VAL2BEAT + EPS && !FINDBEST ) {
	*nphi = removeReflection(nOrbit - 1, fs);
	printf("bbb\n");

	goto endprocedure;
      }
    }
    if ( delta > predelta + EPS ) {
      *nphi = removeReflection(nOrbit - 1, fs); // best sol was updated     
    //      *nphi = fs;
    }
  }
  
 endprocedure:

  //  printf("ddddd %d %d %d %d\n", sel[0], sel[1], sel[2], sel[3]);
  if (  delta > VAL2BEAT + EPS && orbit->flip[*nphi] ) 
    reflectSelection( sel );

  time2use = timer_elapsed_seconds() - tbeforeUse - time2build;

  dprintf("GEN i %4d  j %4d  k %4d  h %4d  bestrot %d val %g\n", sel[0], sel[1], sel[2], sel[3], *nphi, delta );

  move->found = (delta > VAL2BEAT + EPS);
  move->delta = delta;
  for ( i = 0; i < 4; i++ )
    move->selection[i] = sel[i];
  move->maxHeapSize = hsize;
  move->heapExtractions = POPCOUNTER;
  move->effectiveExtractions = UNDEF; //***//
  move->moveEvaluations = EVAL4COUNTER;
  if ( move->found ) {
    move->orbitNo = nOrbit;
    move->phi = *nphi;
    move->scheme = movesInOrbit[nOrbit][*nphi];
  }
  else {
    move->orbitNo = UNDEF;
    move->phi     = UNDEF;
    move->scheme  = UNDEF;
  }
  return delta;
}

double generalCaseOfAny4OPTMoveAnalogo3byOrbit( int* pi, int* bestsel, double* currBestVal,
						THEAP* theap, int numOrbit, int nphi ) {

  /* this is a generic case (case a,b,c,d) of a generic true 4opt move
     it will pick the pairs from the heap <heap>, which could be a real heap 
     or sorted array, to scan in order (the heap knows, object oriented-style)

     Each element of the heap identifies two of the 4 indexes. Which ones is specified by
     the booleans  <i_inheap>, <j_inheap>, <k_inheap>, <h_inheap> computed 
     from the flags stored in the theap

     These two values will then be offset by two offsets stored in the theap

     Then the search will be made on the third and fourth index (i.e., the ones for 
     which <x_inheap> is false

     The procedure updates the heap so that at the endf it becomes a sorted array and it only contains
     those elements which passed the test of potential candidates for the move (i.e., "good enough")

     <currBestVal> is the current value to beat. It also get updated to the new best value

     complete enum would be
     for ( i = 0; i <= n - 7; i++ )
       for ( j = i + 2; j <= n - 5; j++ ) 
         for ( k = j + 2; k <= n - 3; k++ ) 
           for ( h = k + 2; h <= n - 1 - (i==0); h++ ) ...
   */

  boole i_inheap, j_inheap, k_inheap, h_inheap;
  double delta, vv;
  int sel[4];
  TPAIR top;  // the max-element of the heap
  int tmpSize;
  TPAIR* heap;
  int i, j, k, h;
  int mini, maxi, minj, maxj, mink, maxk, minh, maxh;
  int BJ, BK, BH;
  

  heap = theap->heap;
  delta = *currBestVal;
#ifdef PRINTON
  printf("enter generic 4OPT case for delta %g with heap:\n", delta);
  showTHeap( theap );
#endif

  i_inheap = (theap->who1st == 0 || theap->who2nd == 0);
  j_inheap = (theap->who1st == 1 || theap->who2nd == 1);
  k_inheap = (theap->who1st == 2 || theap->who2nd == 2);
  h_inheap = (theap->who1st == 3 || theap->who2nd == 3);

  tmpSize = 0;
  while ( heapsize( heap ) > 0  ) {

    getMaxHeap( &top, heap );
    if ( top.value < delta / 4 + EPS )
      break; 

    tmpHeap[++tmpSize] = top;      // save it in the smaller "heap" (sorted array hereon)
    dprintf("Pop el (x %d, y %d) val %g\n", top.x, top.y, top.value );

    if ( i_inheap && j_inheap  ) { // x=i, y=j, z=k, w=h
      mini = maxi = top.x;
      minj = maxj = top.y;
      mink = top.y + 2;
      maxk = n - 3 - (top.x==0); 
      minh = top.y + 4;
      maxh = n - 1 - (top.x==0); 
    }
    else if ( i_inheap && k_inheap  ) {
      mini = maxi = top.x;
      mink = maxk = top.y;
      minj = top.x + 2;
      maxj = top.y - 2;
      minh = top.y + 2;
      maxh = n - 1 - (top.x==0); 
    }
    else if ( i_inheap && h_inheap  ) {
      mini = maxi = top.x;
      minh = maxh = top.y;
      minj = top.x + 2;
      maxj = top.y - 4;
      mink = top.x + 4;
      maxk = top.y - 2;
    }
    else if ( j_inheap && k_inheap  ) {
      minj = maxj = top.x;
      mink = maxk = top.y;
      mini = (top.y == n - 3);
      maxi = top.x - 2;
      minh = top.y + 2;
      maxh = n - 1;
    }
    else if ( j_inheap && h_inheap  ) {
      minj = maxj = top.x;
      minh = maxh = top.y;
      mini = (top.y == n - 1);
      maxi = top.x - 2;
      mink = top.x + 2;
      maxk = top.y - 2;
    }
    else if ( k_inheap && h_inheap  ) {
      mink = maxk = top.x;
      minh = maxh = top.y;
      mini = (top.y == n - 1);
      maxi = top.x - 4;
      minj = 2 + (top.y == n - 1);
      maxj = top.x - 2;
    }
    else
      assert( false );


    for ( i = mini; i <= maxi; i++ ) {
      sel[0] = i;
      BJ = MIN(maxj, n - 5 - (i==0));
      for ( j = MAX(minj, i+2); j <= BJ; j++ ) {
	sel[1] = j;
	BK = MIN(maxk, n - 3- (i==0));
	for ( k = MAX(mink, j+2); k <= BK; k++ ) {
	  sel[2] = k;
	  BH = MIN(maxh, n - 1 - (i==0));
	  for ( h = MAX(minh, k+2); h <= BH; h++ ) {
	    sel[3] = h;
	    dprintf("EVAL  i %d (%d)   j %d (%d)   k %d (%d)   h %d (%d)\n", i, i_inheap, j, j_inheap,
		    k, k_inheap, h, h_inheap);
	    vv = evaluate4OPT(pi, i, j, k, h, numOrbit, nphi );
	    if ( (FINDBEST && vv > delta + EPS) || (!FINDBEST && vv > EPS) ||
		 ( FINDBEST && vv > delta - EPS && smalllex( bestsel, sel[0], sel[1], sel[2], sel[3] )) ) {
	      dprintf("improvement at vv %g\n", vv);
	      delta = vv;
	      bestsel[0] = sel[0];
	      bestsel[1] = sel[1];
	      bestsel[2] = sel[2];
	      bestsel[3] = sel[3];
	      if ( !FINDBEST )
		goto endloop;
	    }
	  }
	}
      }
    }
  }
  
 endloop:
  
  setHeapSize( tmpHeap, tmpSize );
  duplicateHeap( heap, tmpHeap ); // updates the heap to its best elements
  setRealHeap( heap, false );
  
  *currBestVal = delta;
  return delta;
}

int goodQuadOld( int i, int j, int k, int h ) {

  if (( i < 0 ) || ( i > n - 7 ))
    return false;

  if (( j < i + 2 || j > n - 5 ))
    return false;

  if (( k < j + 2 || k > n - 3 ))
    return false;
  
  if (( h < k + 2 || h > n - 1 ))
    return false;

  if ( n - h + i < 2 )
    return false;

  return true;
}

int goodQuad( int i, int j, int k, int h ) {

  //  printf("test quad %d\n", GOODQ++);

  if ( j < i + 2 )
    return false;

  if ( k < j + 2 )
    return false;
  
  if ( h < k + 2 )
    return false;

  if (( i < 0 ) || ( h > n - 1 ))
    return false;

  if ( n - h + i < 2 )
    return false;


  return true;
}

int biggerlex(int i, int j, int k, int h, int ii, int jj, int kk, int hh ) {

  if ( i > ii )
    return true;
  if ( j > jj )
    return true;
  if ( k > kk )
    return true;

  return ( h > hh );
}


double generalCaseOfAny4OPTMove( int* pi, int* bestsel, double* currBestVal, THEAP* HEAP[2], int nOrbit, int nphi, boole perfectsplit ) {

  /* this is a generic case used to find either f1() + f2() > V/2 OR f3() + f4() > V/2
     (or, if <perfectsplit>, if (f1()+f3()) + (f2()+f4()) > V, i.e., g1() + g2() > V)
     
     The particular 4OPT orbits which "split perfectly" (orbits 2,6,7)
     have been renamed by adding 600 when we want to use this procedure.

     The input are two heaps, and the indices (i,j,k,h) are partitioned so that
     one of the two heaps refers to two of them and the other heap to the other two indices.

     The procedure has to identify one element of the first heap and one of the second heap
     (henceforth, a quadruple i,j,k,h) whose sum (of heap values) is at least Delta/DENOM
     (where Delta is the value to beat, and DENOM = 2 in !perfectsplit, while DENOM = 1 if perfectsplit)
     and compute the true value of the move identified by the 4 indices, possibly updating the
     bestval found.
 
     the loop stop when top(heap1) + top(heap2) do not realize at least Delta/DENOM
     (or when an improving move is found and we are not looking for the best move)

     When looking for pairs of heap elements, there are two (master) pointers in the heaps, one per heap.
     One of them  stays fixed (at the top) while the other (slave) runs through all the elements of
     the opposite heap which could still yield a good 4ple to try (i.e. the sum of heaps values is
     large enough). When this inner loop is over, we move on to decreasing one of the two
     master indices, and repeat the same process.

     Each element of the heap identifies two of the 4 indexes. Which ones is specified by
     the booleans  <i_inheap>, <j_inheap>, <k_inheap>, <h_inheap> computed 
     from the flags stored in the theap

     These two values will then be offset by two offsets stored in the theap

     The procedure updates the heap so that at the endf it becomes a sorted array and it only contains
     those elements which passed the test of potential candidates for the move (i.e., "good enough")

     <currBestVal> is the current value to beat. It also get updated to the new best value

     complete enum would be
     for ( i = 0; i <= n - 7; i++ )
       for ( j = i + 2; j <= n - 5; j++ ) 
         for ( k = j + 2; k <= n - 3; k++ ) 
           for ( h = k + 2; h <= n - 1 - (i==0); h++ ) ...


   */

  TPAIR* inheap[2];
  int size[2];
  int ptr[2]; // index within the heap
  int b, s;   // b = best heap (0 or 1)    s = other heap (=1-b)
  int cnt, i, j, k, h, overall = 0;
  double delta, vv;
  TPAIR top;
  int DENOM;

  delta = *currBestVal;

  if ( heapsize(HEAP[0]->heap) == 0 || heapsize(HEAP[1]->heap) == 0 )
    return delta;

  if ( perfectsplit )
    DENOM = 1;
  else
    DENOM = 2;
  
  for ( i = 0; i < 2; i++ ) {
    inheap[i] = HEAP[i]->heap;
    size[i] = 0;
    ptr[i] = 1;  // top of the heap
    sortedheap[i][0] = inheap[i][0];
  }

  //  #define CORREGGI
  
#ifdef CORREGGI
  printf("\ndelta %g\n", delta);
  printf("Enter. Inheap[0] is\n");
  showHeap( inheap[0]);
  printf("Inheap[1] is\n");
  showHeap( inheap[1]);
#endif
  
  getMaxHeap( &top, inheap[0] );
  size[0] += 1;
  sortedheap[0][size[0]] = top;
  do {
#ifdef CORREGGI
    printf("A) fix %g ) s0 %d  s1 %d  true %d\n", sortedheap[0][1].value, size[0], size[1], heapsize(inheap[1]) );
#endif
    getMaxHeap( &top, inheap[1] );
#ifdef CORREGGI
    printf("pop from H1  x %d  y %d v %g...", top.x, top.y, top.value);
#endif
    if ( top.value + sortedheap[0][1].value < delta / DENOM + EPS ) 
      break;
#ifdef CORREGGI
    printf("ok\n");
#endif
    size[1] += 1;
    sortedheap[1][size[1]] = top;
  } while( heapsize(inheap[1]) > 0 );
  
  do {
#ifdef CORREGGI
   printf("\nB) s0 %d  true %d  s1 %d  true %d\n", size[0], heapsize(inheap[0]), size[1], heapsize(inheap[1]) );
#endif

    getMaxHeap( &top, inheap[0] );
#ifdef CORREGGI
    printf("pop from H0  x %d  y %d v %g...", top.x, top.y, top.value);
#endif
    if ( top.value + sortedheap[1][1].value < delta / DENOM + EPS ) 
      break;

#ifdef CORREGGI
    printf("ok\n");
#endif
    size[0] += 1;
    sortedheap[0][size[0]] = top;
  } while( heapsize(inheap[0]) > 0 );

  for ( i = 0; i < 2; i++ ) {
    setHeapSize( sortedheap[i], size[i] );
    setRealHeap( sortedheap[i], false );
  }

#ifdef CORREGGI
  printf("After. Myheap[0] is\n");
  showHeap( sortedheap[0]);
  printf("Myheap[1] is\n");
  showHeap( sortedheap[1]);
 #endif
  while (  ptr[0] <= heapsize(sortedheap[0]) && ptr[1] <= heapsize(sortedheap[1]) && 
           valueAt(sortedheap[0],ptr[0]) + valueAt(sortedheap[1],ptr[1]) > delta/DENOM + EPS ) {  // rimettere > al posto di >=
    if ( valueAt(sortedheap[0],ptr[0]) > valueAt(sortedheap[1],ptr[1]))
      b=0;
    else
      b=1;
    s = 1 - b;
    //   printf("outer loop ptr %d\n", ptr[b]);
    cnt = ptr[s];
    while ( cnt <= heapsize(sortedheap[s]) &&
	    valueAt(sortedheap[s],cnt) + valueAt(sortedheap[b],ptr[b]) > delta/DENOM + EPS ) { // rimettere > al posto di >=
      //      printf("   inner loop cnt %d\n", cnt);
      i = j = k = h = -1;
      setProperIndices( HEAP[b], sortedheap[b]+ptr[b], &i, &j, &k, &h );
      setProperIndices( HEAP[s], sortedheap[s]+cnt, &i, &j, &k, &h );
#ifdef CORREGGI
      printf("   b %d  ptrb %d  valb %g\n", b, ptr[b], valueAt(sortedheap[b],ptr[b]));
      printf("   s %d  cnt  %d  vals %g\n", s, ptr[s], valueAt(sortedheap[s],cnt));
      printf("\n****** POTENTIAL i %d   j %d   k %d   h %d", i, j, k, h );
#endif
      overall++;
      if ( goodQuad(i,j,k,h) ) { // evaluate the quadruple

	vv = evaluate4OPT(pi, i, j, k, h, nOrbit, nphi );

#ifdef CORREGGI
	printf("****** VALUTAZ is %g\n", vv );
#endif
	if ( (FINDBEST && vv > delta + EPS) || (!FINDBEST && vv > EPS) ||
	     ( FINDBEST && vv > delta - EPS && smalllex( bestsel, i, j, k, h )) ) {
	      
	  delta = vv;
	  bestsel[0] = i;
	  bestsel[1] = j;
	  bestsel[2] = k;
	  bestsel[3] = h;
	  // 	 printf("update best i %d  j %d  k %d  h %d   %g\n", i, j, h, k, delta);
	  if ( !FINDBEST )
	    goto exitloop;
	}
      }
      cnt++;
    }
    //    printf("end inner\n");
    // now the inner loop is over. We determine the new best heap (the heap with the
    // largest value in the current top element) and re-iterate
    ptr[b] = ptr[b] + 1;
  }

exitloop:
  
  *currBestVal = delta;
  return delta;
}

double find4OPTMoveSFbyOrbit( double VAL2BEAT, int* pi, int* sel, int* rot, int numOrbit, OPTMOVE* move ) {
  
  /* 
     
  Given permutation pi[], finds an improving 4OPT move of value > VAL2BEAT
  If it is non-postive, then we are at a local optimum. 
  It returns delta, if there exists delta > currBest.
  Else, it returns -1
  Uses method with f()+f'() (see paper)
  If FINDBEST = true, finds the best 4-opt move, else the first improving move

  */

  THEAP* th[2];
  double delta, predelta, maxA;
  int fs, s, i, phase, norb, hsize = 0;
  TORBIT* orbit;
  int vv[2], uu[2];
  int* permphi;
  int revpi[MAXN], *actualpi;

  assert( numOrbit >= 1 && numOrbit <= 7 );

  norb = numOrbit - 1;  // starts from 0 in c
  orbit = ORBIT4 + norb;
  *rot = -1; // undefined
  
  if ( numOrbit == 4 || numOrbit == 5 ) // orbits w/reflection
    for ( i = 0; i <= n; i++ )
      revpi[i] = pi[n-i];

  delta = VAL2BEAT;
  POPCOUNTER = EVAL4COUNTER = 0;
  hsize = 0;
  //  printf("orbit %d size %d\n", numOrbit, orbit->size );
  //  exit(1);
  for ( fs = 0; fs < orbit->size; fs++ ) {
    s =  ( orbit->flip[fs] ? fs - orbit->size / 2 : fs );
    actualpi =  ( orbit->flip[fs] ? revpi + 0 : pi );
 
    predelta = delta;
    //    actuals = removeReflection( norb, s );
    permphi = orbit->phi[s];
    dprintf("PERMPHI: %d %d %d %d\n", permphi[0], permphi[1], permphi[2], permphi[3]);
    dprintf("\n\n************* NPHI %d\n", s);
    for ( phase = 0; phase < 2; phase++ ) {

      vv[0] = orbit->who[2* phase][0];
      vv[1] = orbit->who[2* phase][1];
      uu[0] = orbit->who[2* phase + 1][0];
      uu[1] = orbit->who[2* phase + 1][1];

      dprintf("Fase %d(a) crea heap per %d %d  vv0 %d  vv1 %d\n", phase, permphi[vv[0]], permphi[vv[1]], vv[0], vv[1]);
      initializeHeap( false, &heapF12, actualpi, permphi[vv[0]], permphi[vv[1]], false, orbit->fct[2 * phase], false, ALLPAIRS, -INFINITE);
      hsize += heapsize( heapF12.heap );
      maxA = getMaxHeapValue( &heapF12 );
      dprintf("delta %g maxA %g\n", delta, maxA);
      dprintf("Fase %d(b) crea heap per %d %d  uu0 %d  uu1 %d\n", phase, permphi[uu[0]], permphi[uu[1]], uu[0], uu[1]);
      initializeHeap( false, &heapF34, actualpi, permphi[uu[0]], permphi[uu[1]], false, orbit->fct[2 * phase + 1], false, BIGGER, delta / 2 - maxA);
      hsize += heapsize( heapF34.heap );

      th[0] = &heapF12;
      th[1] = &heapF34;
      dprintf("cerca la mossa  su orbita %d nphi %d\n", numOrbit, fs);
      generalCaseOfAny4OPTMove( actualpi, sel, &delta, th, numOrbit, s, false );
      // freeHeapMem( &heapF12 );
      //      freeHeapMem( &heapF34 );
      dprintf("After fase s %d parte %d delta is %g val2beat %g\n", s, phase, delta, VAL2BEAT);
      if ( delta > VAL2BEAT + EPS && !FINDBEST )  {	  
	*rot = removeReflection( norb, fs );
	goto endprocedure;
      }
    }
    if ( delta > predelta + EPS ) 
      *rot = removeReflection( norb, fs );  // best sol was updated
  }
  
 endprocedure:

  if ( *rot != -1 && orbit->flip[*rot] )
    reflectSelection( sel );

  dprintf("GENNEW (O%d rot%d) i %4d  j %4d  k %4d  h %4d  val %g\n", numOrbit, *rot, sel[0], sel[1], sel[2], sel[3], delta );


  move->found = (delta > VAL2BEAT + EPS);
  move->delta = delta;
  for ( i = 0; i < 4; i++ )
    move->selection[i] = sel[i];
  move->heapExtractions = POPCOUNTER;
  move->moveEvaluations = EVAL4COUNTER;
  move->effectiveExtractions = UNDEF; //***//
  move->orbitNo = numOrbit;
  move->phi = *rot;
  if ( *rot != -1 ) 
    move->scheme = movesInOrbit[numOrbit][*rot];
  else
    move->scheme = -1; // undef
  move->maxHeapSize = hsize;

  return delta;
}


double find4OPTMovePerfectSplit( double VAL2BEAT, int* pi, int* sel, int* rot, int numOrbit, OPTMOVE* move ) {
  /* given permutation pi[], finds an improving 4OPT move.
  If it is non-postive, then we are at a local optimum. 
  It returns delta, if there exists delta > currBest.  Else, it returns -1
  Uses method with g^1()+g^2() for moves with perfect split (see paper)


  If FINDBEST = true, finds the best 4-opt move, else the first improving move


  This version applies to those moves in which ijkh are splitted in the same way in f_1() and f_2()
  so that we can group them in a different manner than before.
  These are orbits O2, O6 and O7 
  */

  THEAP* th[2];
  double delta, predelta, maxA;
  int s, norb, hsize = 0, i;
  TORBIT* orbit;
  int vv[2], uu[2];
  int* permphi;

  assert( numOrbit == 2 || numOrbit == 6 || numOrbit == 7 );

  norb = numOrbit - 1;  // because c numbers from 0
  orbit = ORBIT4 + norb;
  *rot = -1; // undef

  delta = VAL2BEAT;

  POPCOUNTER = EVAL4COUNTER = 0;
  for ( s = 0; s < orbit->size; s++ ) {
    predelta = delta;
    permphi = orbit->phi[s];
    dprintf("\n\n**PS*********** NPHI %d\n", s);

    vv[0] = orbit->who[0][0];
    vv[1] = orbit->who[0][1];
    uu[0] = orbit->who[1][0];
    uu[1] = orbit->who[1][1];

    dprintf("Fase g1 crea heap per %d %d   (where vv0=%d  vv1=%d)\n", permphi[vv[0]], permphi[vv[1]], vv[0], vv[1]);
    initializeHeap( false, &heapF12, pi, permphi[vv[0]], permphi[vv[1]], false, orbit->g1, false, ALLPAIRS, -INFINITE);
    hsize += heapsize( heapF12.heap );
    maxA = getMaxHeapValue( &heapF12 );
    dprintf("delta %g maxA %g\n", delta, maxA);
    dprintf("Fase g2 crea heap per %d %d\n", permphi[uu[0]], permphi[uu[1]]);
    initializeHeap( false, &heapF34, pi, permphi[uu[0]], permphi[uu[1]], false, orbit->g2, false, BIGGER, delta - maxA);
    hsize += heapsize( heapF34.heap );

    th[0] = &heapF12;
    th[1] = &heapF34;
    dprintf("PS: cerca la mossa  su orbita %d nphi %d\n", numOrbit, s);
    generalCaseOfAny4OPTMove( pi, sel, &delta, th, numOrbit, s, true );
    if ( delta > VAL2BEAT + EPS && !FINDBEST ) {
      *rot = s;
      goto endprocedure;
    }
    if ( delta > predelta + EPS ) *rot = s;  // best sol was updated
  }
  
 endprocedure:

  dprintf("4PARTIC (O%d rot%d) i %4d  j %4d  k %4d  h %4d  val %g\n", numOrbit, *rot, sel[0], sel[1], sel[2], sel[3], delta );

  move->found = (delta > VAL2BEAT + EPS);
  move->delta = delta;
  for ( i = 0; i < 4; i++ )
    move->selection[i] = sel[i];
  move->orbitNo = numOrbit;
  move->phi = *rot;
  if ( *rot != -1 )
    move->scheme = movesInOrbit[numOrbit][*rot];
  else
    move->scheme = -1;
  move->heapExtractions = POPCOUNTER;
  move->moveEvaluations = EVAL4COUNTER;
  move->maxHeapSize = hsize;
  move->effectiveExtractions = UNDEF; //***//

  return delta;
}


/********************************************************

PERFORMING THE MOVES

*********************************************************/

void perform4OPTMoveByScheme( int* pi, int* sel, int nscheme ) {

  int ii, jj, kk, hh;

  ii = sel[0];
  jj = sel[1];
  kk = sel[2];
  hh = sel[3]; 

  assert( nscheme >= 1 && nscheme <= 25 );
  switch ( nscheme ) {

    //ORBIT 1
  case 1: // +1 -2 -3 -4
      reverse( pi, ii + 1, jj );                            // reverse (2)     -> +1 -2 +3 +4
      reverse( pi, jj + 1, kk );                            // reverse (3)     -> +1 -2 -3 +4
      reverse( pi, kk + 1, hh );                            // reverse (4)     -> +1 -2 -3 -4
      break;

  case 24:    // +1 +4 +3 -2
      reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2
      reverse( pi, ii + (hh - kk) + 1, hh - (jj - ii));     // reverse (-3)    -> +1 -4 +3 -2
      reverse( pi, ii + 1, ii + (hh - kk));                 // reverse (-4)    -> +1 +4 +3 -2
      break;

  case 23:    // +1 +4 -3 +2
      reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2
      reverse( pi, ii + 1, ii + (hh - kk));                 // reverse (-4)    -> +1 +4 -3 -2
      reverse( pi, hh - (jj - ii) + 1, hh );                // reverse (-2)    -> +1 +4 -3 +2
      break;

  case 22:    // +1 -4 +3 +2 
      reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2
      reverse( pi, ii + (hh - kk) + 1, hh - (jj - ii));     // reverse (-3)    -> +1 -4 +3 -2
      reverse( pi, hh - (jj - ii) + 1, hh );                // reverse (-2)    -> +1 -4 +3 +2
      break;

      // ORBIT 2
  case 2:   // +1 -2 +3 -4 
      reverse( pi, ii + 1, jj );                            // reverse (2)     -> +1 -2 +3 +4
      reverse( pi, kk + 1, hh );                            // reverse (4)     -> +1 -2 +3 -4
      break;

  case 21:   // +1 -4 +3 -2 
      reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2
      reverse( pi, ii + (hh - kk) + 1, hh - (jj - ii));     // reverse (-3)    -> +1 -4 +3 -2
      break;

      // ORBIT 3
  case 3:   // +1 -2 -4 +3 
      reverse( pi, ii + 1, jj );                            // reverse (2)     -> +1 -2 +3 +4
      reverse( pi, jj + 1, kk );                            // reverse (3)     -> +1 -2 -3 +4
      reverse( pi, jj + 1, hh );                            // reverse (-3,4)  -> +1 -2 -4 +3
      break;
      
  case 7:  // +1 +3 -2 -4
      reverse( pi, ii + 1, kk );                            // reverse (2,3)   -> +1 -3 -2 +4
      reverse( pi, ii + 1, ii + (kk - jj) );                // reverse (-3)    -> +1 +3 -2 +4
      reverse( pi, kk + 1, hh );                            // reverse (4)     -> +1 +3 -2 -4
      break;

  case 13: // +1 +3 -4 -2
      reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2 
      reverse( pi, ii + 1, ii + (hh - jj) );                // reverse (-4,-3) -> +1 +3 +4 -2
      reverse( pi, ii + (kk - jj) + 1, hh - (jj - ii) );    // reverse (+4)    -> +1 +3 -4 -2
      break;

  case 17: // +1 -4 -2 +3
      reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2 
      reverse( pi, ii + (hh - kk) + 1, hh );                // reverse (-3-2)  -> +1 -4 +2 +3 
      reverse( pi, ii + (hh - kk) + 1, hh -(kk - jj) );     // reverse (+2)    -> +1 -4 -2 +3 
      break;

      // ORBIT 4
  case 4: // +1 -2 +4 -3
      reverse( pi, ii + 1, jj );                            // reverse (2)     -> +1 -2 +3 +4
      reverse( pi, jj + 1, hh );                            // reverse (3,4)   -> +1 -2 -4 -3
      reverse( pi, jj + 1, jj + (hh - kk));                 // reverse (-4)    -> +1 -2 +4 -3
      break;

  case 19: // +1 -4 +2 -3
      reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2 
      reverse( pi, ii + (hh - kk) + 1, hh );                // reverse (-3,-2) -> +1 -4 +2 +3
      reverse( pi, hh - (kk - jj) + 1, hh);                 // reverse (+3)    -> +1 -4 +2 -3
      break;

  case 11: // +1 -3 +4 -2
      reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2 
      reverse( pi, ii + 1, ii + (hh - jj) );                // reverse (-4,-3) -> +1 +3 +4 -2 
      reverse( pi, ii + 1, ii + (kk - jj) );                // reverse (+3)    -> +1 -3 +4 -2 
      break; 

  case 6: // +1 -3 +2 -4
      reverse( pi, ii + 1, jj );                            // reverse (2)     -> +1 -2 +3 +4
      reverse( pi, ii + 1, kk );                            // reverse (-2,+3) -> +1 -3 +2 +4
      reverse( pi, kk + 1, hh );                            // reverse (4)     -> +1 -3 +2 -4
      break;
      
      // ORBIT 5
  case 5:   // +1 -2 +4 +3   (r_5)
      reverse( pi, ii + 1, jj );                            // reverse (2)     -> +1 -2 +3 +4
      reverse( pi, jj + 1, hh );                            // reverse (3,4)   -> +1 -2 -4 -3
      reverse( pi, jj + 1, jj + (hh - kk));                 // reverse (-4)    -> +1 -2 +4 -3
      reverse( pi, jj + (hh - kk) + 1, hh);                 // reverse (-3)    -> +1 -2 +4 +3
      break;

  case 20:  // +1 +4 +2 -3   (r_20)
      reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2 
      reverse( pi, ii + (hh - kk) + 1, hh );                // reverse (-3,-2) -> +1 -4 +2 +3
      reverse( pi, hh - (kk - jj) + 1, hh);                 // reverse (+3)    -> +1 -4 +2 -3
      reverse( pi, ii + 1, ii + (hh - kk));                 // reverse (-4)    -> +1 +4 +2 -3
    break;

  case 14:   // +1 +3 -4 +2  (r_14)
      reverse( pi, ii + 1, kk );                            // reverse (2,3)   -> +1 -3 -2 +4
      reverse( pi, ii + (kk - jj) + 1, hh );                // reverse (-2,+4) -> +1 -3 -4 +2
      reverse( pi, ii + 1, ii + (kk - jj) );                // reverse (-3)    -> +1 +3 -4 +2
      break;
      
  case 15:   // +1 -4 -2 -3  (r_15)
      reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2 
      reverse( pi, ii + (hh - kk) + 1, hh );                // reverse (-3,-2) -> +1 -4 +2 +3
      reverse( pi, ii + (hh - kk) + 1, hh -(kk - jj) );     // reverse (+2)    -> +1 -4 -2 +3
      reverse( pi, hh - (kk - jj) + 1, hh);                 // reverse (+3)    -> +1 -4 -2 -3
      break;
      
  case 12:   // +1 -3 +4 +2  (r_12)
      reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2 
      reverse( pi, ii + 1, ii + (hh - jj) );                // reverse (-4,-3) -> +1 +3 +4 -2
      reverse( pi, ii + 1, ii + (kk - jj) );                // reverse (+3)    -> +1 -3 +4 -2
      reverse( pi, hh - (jj - ii) + 1, hh );                // reverse (-2)    -> +1 -3 +4 +2
      break;
      
  case 18:   // +1 +4 -2 +3  (r_18)
      reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2 
      reverse( pi, ii + (hh - kk) + 1, hh );                // reverse (-3,-2) -> +1 -4 +2 +3
      reverse( pi, ii + (hh - kk) + 1, hh -(kk - jj) );     // reverse (+2)    -> +1 -4 -2 +3
      reverse( pi, ii + 1, ii + (hh - kk));                 // reverse (-4)    -> +1 +4 -2 +3
      break;
      
  case 9:   // +1 -3 -4 -2  (r_9)
      reverse( pi, ii + 1, kk );                            // reverse (2,3)   -> +1 -3 -2 +4
      reverse( pi, ii + (kk - jj) + 1, kk );                // reverse (-2)    -> +1 -3 +2 +4
      reverse( pi, ii + (kk - jj) + 1, hh );                // reverse (+2,+4) -> +1 -3 -4 -2
      break;

  case 8:   // +1 +3 +2 -4  (r_8)
      reverse( pi, ii + 1, kk );                            // reverse (2,3)   -> +1 -3 -2 +4
      reverse( pi, ii + 1, ii + (kk - jj) );                // reverse (-3)    -> +1 +3 -2 +4
      reverse( pi, ii + (kk - jj) + 1, kk );                // reverse (-2)    -> +1 +3 +2 +4
      reverse( pi, kk + 1, hh );                            // reverse (4)     -> +1 +3 +2 -4
      break;
      
      // ORBIT 6
  case 10:   // +1 -3 -4 +2 
      reverse( pi, ii + 1, kk );                            // reverse (2,3)   -> +1 -3 -2 +4
      reverse( pi, ii + (kk - jj) + 1, hh );                // reverse (-2,4)  -> +1 -3 -4 +2
      break;

  case 16:   // +1 +4 -2 -3 
      reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2 
      reverse( pi, ii + (hh - kk) + 1, hh );                // reverse (-3-2)  -> +1 -4 +2 +3 
      reverse( pi, ii + (hh - kk) + 1, hh -(kk - jj) );     // reverse (+2)    -> +1 -4 -2 +3
      reverse( pi, hh - (kk - jj) + 1, hh);                 // reverse (+3)    -> +1 -4 -2 -3
      reverse( pi, ii + 1, ii + (hh - kk));                 // reverse (-4)    -> +1 +4 -2 -3
      break;

      // ORBIT 7
  case 25:   // +1 +4 +3 +2
    reverse( pi, ii + 1, hh );                            // reverse (2,3,4) -> +1 -4 -3 -2 
    reverse( pi, ii + (hh - kk) + 1, hh - (jj - ii));     // reverse (-3)    -> +1 -4 +3 -2 
    reverse( pi, ii + 1, ii + (hh - kk));                 // reverse (-4)    -> +1 +4 +3 -2
    reverse( pi, hh - (jj - ii) + 1, hh );                // reverse (-2)    -> +1 +4 +3 +2
    break;


  default:
    printf(" ERR perf move.  i %d  j %d  k %d  h %d scheme %d\n", ii, jj, kk, hh, nscheme);
    assert ( false );
  }  
}

void perform4OPTMove( int* pi, int* sel, int norbit, int nphi ) {

  int ii, jj, kk, hh;

  ii = sel[0];
  jj = sel[1];
  kk = sel[2];
  hh = sel[3]; 

  assert( norbit >= 1 && norbit <= 7 );
  switch ( norbit ) {

  case 1:
    switch ( nphi ) {
    case 0: // +1 -2 -3 -4
      perform4OPTMoveByScheme( pi, sel, 1 );
      break;

    case 1:    // +1 +4 +3 -2
      perform4OPTMoveByScheme( pi, sel, 24 );
      break;

    case 2:    // +1 +4 -3 +2
      perform4OPTMoveByScheme( pi, sel, 23 );
      break;

    case 3:    // +1 -4 +3 +2 
      perform4OPTMoveByScheme( pi, sel, 22 );
      break;
      
    default:
      assert(false);
    }
    break;

  case 2:
    switch ( nphi ) {
    case 0:   // +1 -2 +3 -4
      perform4OPTMoveByScheme( pi, sel, 2 );
      break;

    case 1:   // +1 -4 +3 -2 
      perform4OPTMoveByScheme( pi, sel, 21 );
      break;

    default:
      assert(false);
    }
    break;

  case 3:
    switch ( nphi ) {
    case 0:   // +1 -2 -4 +3 
      perform4OPTMoveByScheme( pi, sel, 3 );
      break;
      
    case 1:  // +1 +3 -2 -4
      perform4OPTMoveByScheme( pi, sel, 7 );
      break;

    case 2: // +1 +3 -4 -2
      perform4OPTMoveByScheme( pi, sel, 13 );
      break;

    case 3: // +1 -4 -2 +3
      perform4OPTMoveByScheme( pi, sel, 17 );
      break;

    default:
      assert(false);
    }
    break;

  case 4:
    switch ( nphi ) {
    case 0: // +1 -2 +4 -3
      perform4OPTMoveByScheme( pi, sel, 4 );
      break;

    case 1: // +1 -4 +2 -3
      perform4OPTMoveByScheme( pi, sel, 19 );
      break;

    case 2: // +1 -3 +4 -2
      perform4OPTMoveByScheme( pi, sel, 11 );
      break; 

    case 3: // +1 -3 +2 -4
      perform4OPTMoveByScheme( pi, sel, 6 );
      break;
      
    default:
      assert(false);
    }
    break;
    
  case 5:
    switch ( nphi ) {
    case 0:   // +1 -2 +4 +3   (r_5)
      perform4OPTMoveByScheme( pi, sel, 5 );
      break;

    case 1:  // +1 +4 +2 -3   (r_20)
      perform4OPTMoveByScheme( pi, sel, 20 );
    break;

    case 2:   // +1 +3 -4 +2  (r_14)
      perform4OPTMoveByScheme( pi, sel, 14 );
      break;
      
    case 3:   // +1 -4 -2 -3  (r_15)
      perform4OPTMoveByScheme( pi, sel, 15 );
      break;
      
    case 4:   // +1 -3 +4 +2  (r_12)
      perform4OPTMoveByScheme( pi, sel, 12 );
      break;
      
    case 5:   // +1 +4 -2 +3  (r_18)
      perform4OPTMoveByScheme( pi, sel, 18 );
      break;
      
    case 6:   // +1 -3 -4 -2  (r_9)
      perform4OPTMoveByScheme( pi, sel, 9 );
      break;

    case 7:   // +1 +3 +2 -4  (r_8)
      perform4OPTMoveByScheme( pi, sel, 8 );
      break;
      
    default:
      assert(false);
    }
    break;
    


  case 6:
    switch ( nphi ) {
    case 0:   // +1 -3 -4 +2 
      perform4OPTMoveByScheme( pi, sel, 10 );
      break;

    case 1:   // +1 +4 -2 -3 
      perform4OPTMoveByScheme( pi, sel, 16 );
      break;

    default:
      assert(false);
    }
    break;

  case 7:   // +1 +4 +3 +2
      perform4OPTMoveByScheme( pi, sel, 25 );
    break;


  default:
    printf(" ERR perf move.  i %d  j %d  k %d  h %d orb %d  phi %d\n", ii, jj, kk, hh, norbit, nphi);
    assert ( false );
  }  
}



/*****************************************************

  Combinatorial functions to count complete selections

******************************************************/

long int numCompleteQuads( int N ) {

  long a, b;
  a = ((long int)(N-3)*(N-4) / 2) * ((long int)(N-5)*(N-6) / 2) / 6;
  b = (long int)(N-5)*(N-6) / 2;
  return a - b;
}

long int numTotQuads( int N ) {

  return ((long int)(N)*(N-1) / 2)  * ((long int)(N-2)*(N-3) / 12);
}

/****************************************************************************************/

/****************************************************************************************/

/****************************************************************************************/


void resetSolution( OPTMOVE* move ) {

  //  zeroes all fields of a solution

  int i;

  /* first, the part identifying the move */

  move->delta = -INFINITE;
  move->found = false;
  for ( i = 0; i < 4; i++ )
    move->selection[i] = UNDEF;
  move->orbitNo = UNDEF;
  move->phi = UNDEF;
  move->scheme = UNDEF;


  /* then the fields that are counters which can be incremented by a sequence
     of such moves */
  
  move->maxHeapSize = 0;
  move->heapExtractions = 0;
  move->effectiveExtractions = 0; 
  move->moveEvaluations = 0;
  move->seconds2Evaluate = 0; 
}

void copyOntoSolution( OPTMOVE* m1, OPTMOVE* m2 ) {

  // given a move m2, copies the part of data regarding the optimal solution
  // (i.e., selection, scheme, value, etc) onto m1

  int i;

  /* first, the part identifying the move */

  m1->delta = m2->delta;
  m1->found = m2->found;
  for ( i = 0; i < 4; i++ )
    m1->selection[i] = m2->selection[i];
  m1->orbitNo = m2->orbitNo;
  m1->phi = m2->phi;
  m1->scheme = m2->scheme;
}

void addToSolution( OPTMOVE* move1, OPTMOVE* move2 ) {

  // given a move move2, adds the incremental data, such as running times, number of
  // evaluations etc, to the same fields of move1

  move1->maxHeapSize +=  move2->maxHeapSize;
  move1->heapExtractions +=  move2->heapExtractions;
  move1->effectiveExtractions += move2->effectiveExtractions; 
  move1->moveEvaluations += move2->moveEvaluations;
}

void find4optMove(int searchType, int norbits, int* orbitcode, int* currentTour, OPTMOVE* move)
{

/*
   A TSP for a complete graph is described by (cost, n), -global vars-
   where n is the number of graph nodes
   and cost is a double[n][n] symmetric matrix of edge costs.

   A K-OPT move is described by K removed edges and K inserted edges
   (an edge can be both removed and inserted), each edge is a couple
   of indexes in [0,n-1].
   No order is required among edges nor between indexes in an edge.

   Input: 
   - <searchType> is a code, i.e.: 
       = SEARCH_4OPT_BF -> search is brute force 
       = SEARCH_SF0/_SF1 -> search is smart force with master/servant heap 
       = SEARCH_SF_SIMILAR3OPT -> search is smart force with 4 heaps same as 3OPT, slower than prec
       = SEARCH_WOEGINGER  -> search is woeginger dyn prog
       = SEARCH_GLOVER  -> search is glover dyn prog // added jun 22
       = -1  -> undefined, skip search
   - <orbitcode[]> are the codes of the orbits to try (1..7), in total <norbits>
     (the move will be searched over all elements of the orbit)
   - <currentTour> the tour before the best move is applied

   Output
   - <newTour> the tour after the move is applied  (TOLTO)
   - <move> info e statistics on the best move found
*/

  int bestsel[4],  o, orbit;
  //  int heapsize, extractions[2];
  // double valmove, time2build, time2eval;
  //  long int evals, evalsGood;
  int bestrot;
  OPTMOVE tmpmove;


  bestsel[0] = bestsel[1] = bestsel[2] = bestsel[3] = -1; // undef
  
  if ( searchType == -1 ) { // undefined
    printf("ERR : 4OPT algorithm undefined %d\n", searchType );
    assert( false );  
  }
  resetSolution( move );
  double tbefore = timer_elapsed_seconds();
  move->delta = 0; // value to beat is at least 0
  for ( o = 0; o < norbits; o++ ) {
    if ( !FINDBEST && move->delta > EPS )
      goto endfor;
    
    orbit = orbitcode[o];
    printf("IN %14.2f orbit %d OUT ", move->delta, orbitcode[o] );
    fflush(stdout);
    if ( searchType == SEARCH_4OPT_BF ) {
      find4OPTMoveBFbyOrbit( move->delta, currentTour, bestsel, &bestrot, orbit, &tmpmove );
    }
    else if ( searchType == SEARCH_SF_SIMILAR3OPT ) {
      find4OPTMoveSFAnalogo3byOrbit( move->delta, currentTour, bestsel, &bestrot, orbit, &tmpmove );
    }
    else if ( searchType == SEARCH_SF1 && (orbit == 2 || orbit == 6 || orbit == 7) ) {
      find4OPTMovePerfectSplit( move->delta, currentTour, bestsel, &bestrot, orbit, &tmpmove );
    }
    else if ( searchType == SEARCH_SF0 || searchType == SEARCH_SF1 ) {
      find4OPTMoveSFbyOrbit( move->delta, currentTour, bestsel, &bestrot, orbit, &tmpmove );
      
    }
    else if ( searchType == SEARCH_WOEGINGER ) {
      find4OPTWoegingerByOrbit( currentTour, bestsel, &bestrot, orbit, &tmpmove );
    }
    else if ( searchType == SEARCH_GLOVER ) {
      find4OPTGloverByOrbit( currentTour, bestsel, &bestrot, orbit, &tmpmove );
    }
    else {
      printf("Search type %d not valid\n", searchType );
      assert( false );
    }
    if ( tmpmove.delta > move->delta + EPS ) 
      copyOntoSolution( move, &tmpmove) ;

    addToSolution( move, &tmpmove );
    printf("%14.2f\n", move->delta);
    fflush(stdout);
  }

 endfor:

  move->seconds2Evaluate =  timer_elapsed_seconds() - tbefore;
}


/*-----------------------------------------------
parte relativa alle mosse di Glover, aggiunta per testing in 06/22
alle procedure "dopo ODS" del 2019
--------------------------------------------------------------*/

static double cost_c( int* pi, int i, int j ) {

  // computes the cost of a 2-OPT move which re-creates the tour (i.e., a "regular"
  // move. Use the same notation as in Glover's paper (but we index from 0, while
  // he indexes from 1)
  // the arcs removed are (pi(i),pi(i+1)) and ( pi(j), pi(j+1) )
  
  double cost_out, cost_in;

  cost_out = cost[pi[i]][pi[SUCC(i)]] + cost[pi[j]][pi[SUCC(j)]];
  cost_in = cost[pi[i]][pi[j]] + cost[pi[SUCC(i)]][pi[SUCC(j)]];
  
  return cost_out - cost_in;
}

static double cost_d( int* pi, int i, int j ) {

  // computes the cost of a 2-OPT move which disconnectes the tour (i.e., not a valid 2OPT
  // move. Use the same notation as in Glover's paper (but we index from 0, while
  // he indexes from 1)
  
  double cost_out, cost_in;

  cost_out = cost[pi[i]][pi[SUCC(i)]] + cost[pi[j]][pi[SUCC(j)]];
  cost_in = cost[pi[i]][pi[SUCC(j)]] + cost[pi[SUCC(i)]][pi[j]];
  /*
  printf("costd %d %d out %d in %d tot %d\n", i, j, (int) cost_out, (int) cost_in, (int) (cost_out-cost_in));
  scanf("%d",&i);
  */
  return cost_out - cost_in;
}


#ifdef IMPLEMENTAZIONE_GLOVER_COMEPAPER

static double cost_x( int*pi, int typex, int i, int j ) {


  if ( typex == TYPEC )
    return cost_c( pi, i, j );
  else if ( typex == TYPED )
    return cost_d( pi, i, j );
  else
    assert( false );
}


double find4OPTGlover(int* pi, int * ii, int* jj, int* kk, int* hh ) { 
  /* given permutation pi[], finds the best 4OPT move, of type 1 or 2 (as in Glover's paper)
     i.e., 4 indices i < j < k < h defining a move of type 1 or 2, which in our paper are
     either orbit R25 or R10 

- - - - -

  Uses O(n^2) method by Glover

  We implement the graph S with four matrices of size N x N of which only about half the
  cells are used. 

*/

  double sinkval;
  int i, j;
  int back[4];      // back[.] are the indices of the pivots
  int backtype[4];  // backtype refers to which submatrix (0,1,2,3) the backpointer owes its value

  
  assert ( n < MAXNGLOVER );
  
  /* initialization */

  for ( i = 0; i < n; i++ )
    for  ( j = i + 2; j < n; j++ )
      reach_c[i][j] = reach_d[i][j] = cross_c[i][j] = cross_d[i][j] = -INFINITE;
  
  for ( i = 0; i <= n - 4; i++ )
    for ( j = i + 2; j <= n - 2; j++ )
      reach_c[i][j] = cost_c(pi, i, j);

  for ( i = 0; i <= n - 4; i++ )
    for ( j = i + 2; j <= n - 2; j++ )
      reach_d[i][j] = cost_c(pi, i, j);

  sinkval = -INFINITE;

  // computes the longest path

  // 1. reach_c
  for ( i = 1; i <= n - 4; i++ )
    for ( j = i + 2; j <= n - 2; j++ )
      if ( reach_c[i][j] < reach_c[i-1][j]  ) 
	reach_c[i][j] = reach_c[i-1][j]; 

  // 2. reach_d
  for ( i = 1; i <= n - 4; i++ )
    for ( j = i + 2; j <= n - 2; j++ )
      if ( reach_d[i][j] < reach_d[i-1][j] ) 
	reach_d[i][j] = reach_d[i-1][j]; 

  // 3. cross_c
  for ( i = 0; i <= n - 4; i++ )
    for ( j = i + 2; j <= n - 2; j++ ) {
      if ( cross_c[i][j] < reach_c[i][j] ) 
	cross_c[i][j] = reach_c[i][j];
      if ( cross_c[i][j] < cross_c[i][j-1] ) 
	cross_c[i][j] = cross_c[i][j-1];

      // try to update the sink value

      if ( sinkval < cross_c[i][j] + cost_d(pi, i+1, j+1) )
	sinkval = cross_c[i][j] + cost_d(pi, i+1, j+1);
      if ( sinkval < cross_c[i][j] + cost_c(pi, i+1, j+1) )
	sinkval = cross_c[i][j] + cost_c(pi, i+1, j+1);
    }

  // 4. cross_d
  for ( i = 0; i <= n - 4; i++ )
    for ( j = i + 2; j <= n - 2; j++ ) {
      if ( cross_d[i][j] < reach_d[i][j] ) 
	cross_d[i][j] = reach_d[i][j];
      if ( cross_d[i][j] < cross_d[i][j-1] ) 
	cross_d[i][j] = cross_d[i][j-1];

      // try to update the sink value

      if ( sinkval < cross_d[i][j] + cost_d(pi, i+1, j+1) )
	sinkval = cross_d[i][j] + cost_d(pi, i+1, j+1);
    }

  *ii = -1;
  *jj = -1;
  *kk = -1;
  *hh = -1;
  return sinkval;
}

double find4OPTGloverAsPaper(int* pi, int * ii, int* jj, int* kk, int* hh, int* myorbit, int* myrho ) { 
  /* given permutation pi[], finds the best 4OPT move, of type 1 or 2 (as in Glover's paper)
     i.e., 4 indices i < j < k < h defining a move of type 1 or 2, which in our paper are
     either orbit R25 or R10 

- - - - -

  Uses O(n^2) method by Glover

*/

  int i, j, x, frst, scnd;
  int* pfirst, *psecond;
  
  assert ( n < MAXNGLOVER );
  
  /* initialization */

  best_t = -INFINITE;

  for ( x = TYPEC; x <= TYPED; x++ ) 
    for ( j = 2; j <= n-2; j++ ) {
      bestreach[x][0][j] = cost_x(pi, x, 0, j);
      predreach[x][0][j] = 0;
      simpleTerminalNodeUpdate( pi, 0, j );
    }

  for ( x = TYPEC; x <= TYPED; x++ ) 
    for ( i = 1; i <= n - 3; i++ ) {
      bestcross[x][i-1][i+1] = bestreach[x][i-1][i+1];
      predcross[x][i-1][i+1][0] = predreach[x][i-1][i+1];
      predcross[x][i-1][i+1][1] = i+1;
      for ( j = i + 2; j <= n - 1; j++ ) {
	bestCrossUpdate( pi, x, i, j );
	bestReachUpdate( pi, x, i, j );
	bestTerminalNodeUpdate( pi, x, i, j );	  
	simpleTerminalNodeUpdate( pi, i, j );
      }
    }

  // recovers solution


  if ( type[0] == TYPEC && type[1] == TYPED ) {

    // best move degenerates in simple 2OPT

    *ii = pred_t2[0];
    *jj = pred_t2[1];
    *kk = -1;
    *hh = -1;
    *myorbit = -1;
  }
  else {

    // faccio in modo che pred_t1[first][0] sia l'indice minore dei pivot

    if (pred_t1[0] < pred_t2[0]) {
      frst = 0;
      scnd = 1;
      pfirst = &(pred_t1[0]);
      psecond = &(pred_t2[0]);
    }
    else {
      frst = 1;
      scnd = 0;
      pfirst = &(pred_t2[0]);
      psecond = &(pred_t1[0]);
    }
    
    *ii = pfirst[0];
    *kk = pfirst[1];
    *jj = psecond[0];
    *hh = psecond[1];

    if ( type[frst] == TYPED && type[scnd] == TYPED ) {
      *myorbit = 7;
      *myrho = 0;
    }
    else {
      *myorbit = 6;
      if ( type[frst] == TYPEC )
	*myrho = 0;
      else
	*myrho = 1;
    }
  }
  
  return best_t;
}


void bestCrossUpdate( int*pi, int x, int i, int startj ) {

  int j;

  for ( j = startj; j <= n - 2; j++ )
    if ( bestreach[x][i-1][j] > bestcross[x][i-1][j-1] ) {
      bestcross[x][i-1][j] = bestreach[x][i-1][j]; 
      predcross[x][i-1][j][0] = predreach[x][i-1][j]; 
      predcross[x][i-1][j][1] = j;
    }
    else {
      bestcross[x][i-1][j] = bestcross[x][i-1][j-1]; 
      predcross[x][i-1][j][0] = predcross[x][i-1][j-1][0]; 
      predcross[x][i-1][j][1] = predcross[x][i-1][j-1][1]; 
    }
}

void bestReachUpdate( int*pi, int x, int i, int startj ) {

  int j;

  for ( j = startj; j <= n - 2; j++ )
    if ( cost_x( pi, x, i , j ) > bestreach[x][i-1][j] ) {
      bestreach[x][i][j] = cost_x( pi, x, i, j ); 
      predreach[x][i][j] = i;
    }
    else {
      bestreach[x][i][j] = bestreach[x][i-1][j]; 
      predreach[x][i][j] = predreach[x][i-1][j];
    }
}


void bestTerminalNodeUpdate( int* pi, int x, int i, int j ) {

  int y, miny, maxy;

  maxy = TYPED;
  if ( x == TYPED )
    miny = TYPEC;
  else
    miny = TYPED;

  for ( y = miny; y <= maxy; y++ ) {
    if ( cost_x( pi, x, i, j ) + bestcross[y][i-1][j-1] > best_t ) {
      best_t = cost_x( pi, x, i, j ) + bestcross[y][i-1][j-1];
      pred_t1[0] = predcross[y][i-1][j-1][0];
      pred_t1[1] = predcross[y][i-1][j-1][1];
      pred_t2[0] = i;
      pred_t2[1] = j;
      type[0] = x;
      type[1] = y;
    }
  }
}

void simpleTerminalNodeUpdate( int*pi, int i, int j ) {

  if ( cost_c( pi, i, j) > best_t ) {
    best_t = cost_c(pi, i, j);
    pred_t2[0] = i;
    pred_t2[1] = j;
    type[0] = TYPEC;
    type[1] = TYPEC;
  }
}

#endif

// *************** mia implementazione Glover *************

double A(int* pi, int x, int y, int minx, int* p ) {

  // returns the value of a best possible 2-OPT disconnecting move
  // which removes an edge (p,p+1) with minx <= p <= x
  // and the edge (y, y+1)
  // minx =  0 or = 1

  double removeX, keepX;

  assert( y >= x );
  if ( x < minx || y > n - 1 ) 
    return -INFINITE;

  if ( memoizeA[minx][x][y] == INFINITE ) { // computes the cell once and forall

    if ( x == PRED(y) || x == SUCC(y) )
      removeX = -INFINITE;
    else
      if ( !SWITCH_GLOVER_COST )
	removeX = cost_d( pi, x, y );
      else
	removeX = cost_c( pi, x, y );

    if ( x > minx )
      keepX = A(pi, x-1, y, minx, p);
    else
      keepX = -INFINITE;

    //  printf("At A(%d %d) removeX %d keepX %d\n",x, y, (int)removeX, (int)keepX);
    if ( keepX > removeX ) {
      memoizeBestA[minx][x][y] = *p;
      memoizeA[minx][x][y] = keepX;
    }
    else {
      memoizeBestA[minx][x][y]= x;
      memoizeA[minx][x][y] = removeX;
    }
  }

  // cell has already been comupted

  *p = memoizeBestA[minx][x][y];
  return memoizeA[minx][x][y];
}

double B(int* pi, int x, int y, int minx, int*p, int*q ) {

  // returns the value of a best possible 2-OPT disconnecting move
  // which removes an edge (p,p+1) with minx <= p <= x - 2
  // and an edge (q, q+1) with x+2 <= q <= y - 2
  // minx can be either 0 or 1

  double removeY, keepY;
  int pr;

  assert( y >= x );
  if ( x - 2 < minx || y < x + 4 )
    return -INFINITE;

  if ( memoizeB[minx][x][y] == INFINITE ) { // computes the cell once and forall

    removeY = A(pi, x - 2, y - 2, minx, &pr);
    if ( y - 2 > x + 2 )
      keepY = B(pi, x, y - 1, minx, p, q);
    else
      keepY = -INFINITE;

    if ( keepY > removeY ) {
      memoizeBestB[0][minx][x][y] = *p;
      memoizeBestB[1][minx][x][y] = *q;
      memoizeB[minx][x][y] = keepY;
    }
    else {
      memoizeBestB[0][minx][x][y] = pr;
      memoizeBestB[1][minx][x][y] = y - 2;
      memoizeB[minx][x][y] = removeY;
    }
  }
  
  // cell has already been comupted

  *p = memoizeBestB[0][minx][x][y];
  *q = memoizeBestB[1][minx][x][y];
  return memoizeB[minx][x][y];
}


double find4OPTGloverDynProgOrbit7(int* pi, int * ii, int* jj, int* kk, int* hh, int* rot ) {

  // finds the best move of orbit7 using dynamic programming

  int i, j, k, h;
  int besti = -1, bestj = -1, besth = -1, bestk = -1;
  double opt, v;


  SWITCH_GLOVER_COST = false;
  
  for ( i = 0; i < n - 1; i++ )
    for ( j = i; j < n; j++ ) {
      memoizeA[0][i][j] = INFINITE; // dummy value that cannot be achieved, represents unset
      memoizeA[1][i][j] = INFINITE; // dummy value that cannot be achieved, represents unset
      memoizeB[0][i][j] = INFINITE; // dummy value that cannot be achieved, represents unset
      memoizeB[1][i][j] = INFINITE; // dummy value that cannot be achieved, represents unset
    }
  
  opt = - INFINITE;

  // #define RIMUOVI
  #ifdef RIMUOVI

  int minx;
  int p, q;

  for (minx=0; minx<=1; minx++ ) {
    
    for ( i = minx; i < n -2; i++ ) {
      for ( j = i + 2; j < n; j++ ) {
	printf("A[%d](%d %d) = %d  ", minx, i, j, (int)A(pi, i, j, minx, &k));
	printf("[%d] ", k);
	printf("check\n");
	/*
	while ( true ){
	  scanf("%d", &p);
	  if ( p < 0 )
	    break;
	  scanf("%d", &q );
	  if ( q < 0 )
	    break;
	  printf("cost(%d %d) = %d\n", p, q, (int)cost_d(pi, p,q));
	}
	*/
      }
      printf("\n");
    }
    printf("\n");
  }


  for ( minx = 0; minx <= 1; minx++ )
    for ( i = 0; i < n -1; i++ ) {
      for ( j = i + 1; j < n; j++ ) {
	printf("B[%d](%d %d) = %d  ", minx, i, j, (int)B(pi, i, j, minx, &k, &h));
	printf("[%d %d] ", k, h);
      }
      printf("\n");
    }

  #endif
  
  for ( j = 2; j <= n-5; j++ )
    for ( h = j + 4; h <= n - 1; h++ ) {
      v = cost_d( pi, j, h ) + B( pi, j, h, (h==n-1), &i, &k );
      if ( opt < v ) {
	//	printf("update opt i %d j %d k %d h %d at %d\n", i, j, k, h, (int) v);
	bestj = j;
	bestk = k;
	besti = i;
	besth = h;
	opt = v;
      }
    }

  *ii = besti;
  *jj = bestj;
  *kk = bestk;
  *hh = besth;
  *rot = 0;
  return opt;
}

double find4OPTGloverDynProgOrbit6(int* pi, int * ii, int* jj, int* kk, int* hh, int* rot ) {

  // finds the best move of orbit6 using dynamic programming

  int i, j, k, h;
  int besti = -1, bestj = -1, besth = -1, bestk = -1; 
  double opt, v;

  // finds best move for r_16
  
  SWITCH_GLOVER_COST = false;

  for ( i = 0; i < n - 1; i++ )
    for ( j = i; j < n; j++ ) {
      memoizeA[0][i][j] = INFINITE; // dummy value that cannot be achieved, represents unset
      memoizeA[1][i][j] = INFINITE; // dummy value that cannot be achieved, represents unset
      memoizeB[0][i][j] = INFINITE; // dummy value that cannot be achieved, represents unset
      memoizeB[1][i][j] = INFINITE; // dummy value that cannot be achieved, represents unset
    }
  
  opt = - INFINITE;
  
  for ( j = 2; j <= n-5; j++ )
    for ( h = j + 4; h <= n - 1; h++ ) {
      v = cost_c( pi, j, h ) + B( pi, j, h, (h==n-1), &i, &k );
      if ( v > opt + EPS ) {
	//	printf("update opt i %d j %d k %d h %d at %d\n", i, j, k, h, (int) v);
	bestj = j;
	bestk = k;
	besti = i;
	besth = h;
	opt = v;
	*rot = 1;
      }
    }

  //  printf("At the end of orbit6 R_16 best is %d %d %d %d val %d\n", besti, bestk, bestk, besth, (int)opt);

  // finds best move for r_10
  
  SWITCH_GLOVER_COST = true;
  
  for ( i = 0; i < n - 1; i++ )
    for ( j = i; j < n; j++ ) {
      memoizeA[0][i][j] = INFINITE; // dummy value that cannot be achieved, represents unset
      memoizeA[1][i][j] = INFINITE; // dummy value that cannot be achieved, represents unset
      memoizeB[0][i][j] = INFINITE; // dummy value that cannot be achieved, represents unset
      memoizeB[1][i][j] = INFINITE; // dummy value that cannot be achieved, represents unset
    }
    
  for ( j = 2; j <= n-5; j++ )
    for ( h = j + 4; h <= n - 1; h++ ) {
      v = cost_d( pi, j, h ) + B( pi, j, h, (h==n-1), &i, &k );
      if ( v > opt + EPS ) {
	//	printf("update opt i %d j %d k %d h %d at %d\n", i, j, k, h, (int) v);
	bestj = j;
	bestk = k;
	besti = i;
	besth = h;
	opt = v;
	*rot = 0;
      }
    }

  //  printf("At the end of orbit6 R_10 best is %d %d %d %d val %d\n", besti, bestk, bestk, besth, (int)opt);
  
  *ii = besti;
  *jj = bestj;
  *kk = bestk;
  *hh = besth;
  return opt;
}

double find4OPTGloverByOrbit( int* pi, int * sel, int* rot, int orbit, OPTMOVE* move ) {

  double delta;
  int i, j, k, h;

  assert( orbit == 6 || orbit == 7 );

  /*
  for ( i = 0; i < n ; i++ ) {
    for ( j = 0; j < n; j++ ) {
      printf("%4d ", (int)cost[i][j] );
      if ( (j+1) % 5 == 0 )
	printf(" | ");
    }
    if ( (i+1) % 5 == 0 )
      printf("\n----");
    printf("\n");
  }
  */

  if ( orbit == 6 )
    delta = find4OPTGloverDynProgOrbit6( pi, &i, &j, &k, &h, rot );
  else
    delta = find4OPTGloverDynProgOrbit7( pi, &i, &j, &k, &h, rot );

  sel[0] = i;
  sel[1] = j;
  sel[2] = k;
  sel[3] = h;

  move->found = (delta > EPS);
  move->delta = delta;
  for ( i = 0; i < 4; i++ )
    move->selection[i] = sel[i];
  move->heapExtractions = 0;
  move->moveEvaluations = 0;
  move->effectiveExtractions = 0;
  move->orbitNo = orbit;
  move->phi = *rot;
  move->scheme = movesInOrbit[orbit][*rot];
  move->maxHeapSize = 0;

  return delta;
}

