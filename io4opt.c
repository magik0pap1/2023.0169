#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "tspmain.h"
#include "io4opt.h"
#include "proc4opt.h"

int PRECISION;  // set decimal precision in costs
int POWERPREC = 1; // 10^PRECISION
double CMIN = 0.0;  // costs are in [CMIN,CMAX)
double CMAX = 1.0;  // costs are in [CMIN,CMAX)

double xcoor[MAXN], ycoor[MAXN];

double u[MAXN];  /* used to perturb costs as c[i][j] + u[i] + u[j] -2sum_i u[i]/n */

double putPrecision( double v ) {

  // resets the value in v so as its precision becomes PRECISION digits
  // e.g., if PRECISION = 2,  1.02356 -> 1.02, 0.17533 -> 0.17, etc
  // value is truncated, not rounded

  return ((int)(v * POWERPREC)) / (double) POWERPREC;
}

void setPowerPrec() {

  // sets POWERPREC to be 10^PRECISION

  int i;

  POWERPREC = 1;
  for ( i = 1; i <= PRECISION; i++ )
    POWERPREC = POWERPREC * 10;

  /*
  printf("PREC %d POW %d\n", PRECISION, POWERPREC );  double v;
  while (1) {

    printf("v : ");
    scanf("%lg", &v);
    v = putPrecision(v);
    printf("becomes %lg\n", v);
  }
  */
  
}


void* mycalloc( int numelem, int sizeelem ) {
  void* p;
  
   p = (void*) calloc( numelem, sizeelem );
   assert(p);
   //    printf("\nAllocate %d elem of size %d at addr %ld\n", numelem, sizeelem, (long int) p );
   //   int i;
   //   scanf("%d", &i);
  return p;
}

void myfree( void* p ) {
  //  printf("\nfree addr %ld\n", (long int) p );
  //  int i;
  //  scanf("%d", &i);
  free( p );
}

/*********************************

input/output procedures 

**********************************/

double tourValue( int* p ) { 
  // computes the value of a tour

  double v = 0;
  int i;

  assert( p[0] == 0 );
  assert( p[n] == 0 );
  for ( i = 0; i < n; i++ )
    v += edgecost(p[i], p[i+1]);
  return v;
}

void randomTple( int* sel, int n_, int k ) {

  //  generates a random k-ple of numbers in 0..n_-1

  int i, j, left, tmp;
  int* rem;

  rem = (int*) mycalloc(n_ + 1, sizeof(int));

  for ( i = 0; i < n_; i++ )
    rem[i] = i;
  
  left = n_;
  
  for ( i = 0; i < k; i++ ) {
    j = rand() % left;
    sel[i] = rem[j];
    rem[j] = rem[left-1];
    left--;
  }
  
  myfree(rem);  

  // ordina la selezione
  
  for ( i = 0; i < k - 1; i++ )
    for ( j = i + 1; j < k; j ++ )
      if ( sel[i] > sel[j] ) {
	tmp = sel[i];
	sel[i] = sel[j];
	sel[j] = tmp;
      }
}

void rndPerm0( int *pi, int n_ ) 
/* generates pi (of elements in 0..n-1) uniformly at random,
   but with pi[0] = 0.
   memory for pi is allocated by caller. Also pi[n] = 0 to close the cycle*/
{
  int i, j, left;
  int* rem;

  rem = (int*) mycalloc(n_ + 1, sizeof(int));

  pi[0] = pi[n] = 0;
  for ( i = 1; i <= n_ - 1; i++ )
    rem[i-1] = i;
  left = n_ - 1;
  
  for ( i = 1; i < n_; i++ ) {
    j = rand() % left;
    pi[i] = rem[j];
    rem[j] = rem[left-1];
    left--;
  }

  myfree(rem);
}

double randomCost() {

  double v;
  
  // random in [CMIN,CMAX) :
  v = CMIN + drand48() * (CMAX - CMIN);
  return v;
}

double edgecost( int i, int j )
{
  if ( i == j )
    return INFINITE;

#ifndef COSTS_ONFLY
  return cost[i][j];
#else
  if ( COSTI_EUCLIDEI ) {
    // printf("x(%d) %g y(%d) %g x(%d) %g y(%d) %g\n", i, xcoor[i], i, ycoor[i], j, xcoor[j], j, ycoor[j]);
    return euclDist(xcoor[i], ycoor[i], xcoor[j], ycoor[j]);
  }
  else /* costi uniformi in (0,1) */ {
    if ( i < j ) 
      srand48( ((i*n) + j) * 100 + SEEDCOST );
    else 
      srand48( ((j*n) + i) * 100 + SEEDCOST );
    return drand48();
  }
#endif
}

void showCosts() {
  int  i, j;
  
  for ( i = 0; i < n && i < 10; i++ ) {
    for ( j = 0; j < n && j < 10; j++ ) 
      printf("%8.2g ", edgecost(i, j));
    printf("\n");
  }
}



void freeCost() { // release memory
  int i;

  if ( cost == NULL )
    return ;
  
  for ( i = 0; i < n; i++ ) {
    free( cost[i] );
    cost[i] = NULL;
  }
  free (cost);
  cost = NULL;
}

void setUniformCosts() { 
  int i, j;

  freeCost();
  setPowerPrec();
  cost = (double**) calloc( n, sizeof(double*));
  for ( i = 0; i < n; i++ )
    cost[i] = (double*) calloc( n, sizeof(double));
  
  for ( i = 0; i < n; i++ ) {
    cost[i][i] = 0;
    for ( j = i + 1; j < n; j++ ) { 
      cost[i][j] = cost[j][i] = putPrecision( randomCost() );
    }
  }

  if ( OUTPUTON )
    showCosts();
}

double euclDist( double xi, double yi, double xj, double yj ) {

  //  float xx = (xcoor[i] - xcoor[j]) * (xcoor[i] - xcoor[j]); 
  //  float yy = (ycoor[i] - ycoor[j]) * (ycoor[i] - ycoor[j]); 
  double xx = (xi - xj) * (xi - xj); 
  double yy = (yi - yj) * (yi - yj); 
  return sqrt( xx + yy );
}

void setEuclideanCosts() { 
  int i, j;

  freeCost();
  setPowerPrec();
  cost = (double**) mycalloc( n, sizeof(double*));
  for ( i = 0; i < n; i++ )
    cost[i] = (double*) mycalloc( n, sizeof(double));

  for ( i = 0; i < n; i++ ) {
    xcoor[i] = drand48();
    ycoor[i] = drand48();
  }
    
  for ( i = 0; i < n; i++ ) {
    cost[i][i] = 0;
    for ( j = i + 1; j < n; j++ ) { 
      cost[i][j] = cost[j][i] = putPrecision( euclDist( xcoor[i], ycoor[i], xcoor[j], ycoor[j] ) );
    }
  }

  if ( OUTPUTON )
    showCosts();
}

void perturbCostsOnNodes() {

  if ( MAX_NODE_U == 0 )
    return;

  int i, j;
  double USUM, minCost;

  USUM = 0;
  for ( i = 0; i < n; i++ ) {
    u[i] = drand48() * MAX_NODE_U;
    USUM += u[i];
  }

  USUM = (2*USUM) / n;  // per non dovere fare ogni volta questo calcolo
  minCost = INFINITE;
  for ( i = 0; i < n - 1; i++ )
    for ( j = i + 1; j < n; j++ ) {
      cost[i][j] = cost[j][i] = cost[i][j] + u[i] + u[j] - USUM;
      if ( cost[i][j] < minCost - EPS )
	minCost = cost[i][j];
    }

  if ( minCost < -EPS ) {
    printf("There are negative arcs %g\n", minCost );
    assert(false);
  }
}

void setCosts( int uniform ) {

  if ( RANDOM_GRAPH_TYPE == UNIFORM )
    setUniformCosts();
  else  if ( RANDOM_GRAPH_TYPE == EUCLIDEAN )
    setEuclideanCosts();
  else {
    printf("Unknown graph type %d\n", RANDOM_GRAPH_TYPE );
    assert( false );
  }
}


void readCostsFile( char* fname, int doPerturb ) { 

  // reads the cost matrix from a file

  FILE* f;
  int i, j;
  double v;

  f = fopen( fname, "r" );
  if ( f == NULL ) {
    printf("File %s not found\n", fname );
    exit(1);
  }

  printf("Reading file %s...", fname);
  fflush(stdout);
  
  fscanf( f, "%d", &n );

  freeCost();
  cost = (double**) mycalloc( n, sizeof(double*));
  for ( i = 0; i < n; i++ )
    cost[i] = (double*) mycalloc( n, sizeof(double));

  for ( i = 0; i < n; i++ ) 
    cost[i][i] = 0;

  for ( i = 0; i < n - 1; i++ )
    for ( j = i + 1; j < n; j++ )  
      cost[i][j] = cost[j][i] = INFINITE;

  while ( fscanf(f, "%d %d", &i, &j) != EOF ) {
    assert( i >=0 && i < j && j < n );
    fscanf(f, "%lg", &v);
    if ( !doPerturb )
      cost[i][j] = cost[j][i] = v;
    else 
      cost[i][j] = cost[j][i] = v + drand48() * 0.0001;
  }

  if ( OUTPUTON )
    showCosts();

  fclose(f);
  printf("Done\n");
  fflush(stdout);
}

void readTour( char* fname, int* pi, int n ) {

  // reads a tour of n nodes from a file
  // the number of nodes is already known

  FILE* f;
  int nnod, i, el;

  printf("Read starting tour <%s> ... ", fname );
  fflush( stdout );
  f = fopen( fname, "r" );
  if ( f == NULL ) {
    printf("File not found\n");
    exit(1);    
  }
  (void) fscanf( f, "%d", &nnod );
  if ( nnod != n ) {
    printf("Nodes %d are not %d as expected\n", nnod, n);
    exit(1);    
  }

  for ( i = 0; i < n; i++ ) {
    fscanf( f, "%d", &el );
    if ( i == 0 && el != 0 ) {
      printf("Tour starts from %d and not 0 as expected\n", el);
      exit(1);    
    }
    pi[i] = el;
  }
  pi[n] = 0;
  fclose(f);
  printf("done\n");
}

void writeTour( char* fname, int* pi, int n, int step, int numRun, int islast ) {

  // writes the tour resulting after <step> moves of the <numRun>-th
  // local search have been made
  
  FILE* f;
  int i;
  char fullname[100];

  if ( islast )
    sprintf( fullname, "%s.%02d.LAST", fname, numRun );  
  else
    sprintf( fullname, "%s.%02d.%04d", fname, numRun, step );  
  printf("write curr tour on <%s> ... ", fullname );
  fflush( stdout );
  f = fopen( fullname, "w" );
  fprintf( f, "%d\n", n );
  for ( i = 0; i <= n; i++ ) 
    fprintf( f, "%d ", pi[i] );
  
  fprintf(f, "\n");
  printf("done\n");
  fclose( f );
  fflush(stdout);
}

void writeBestTour( char* fname, int n, int* pi, double tval ) {
  // writes the tour resulting after <step> moves of the <numRun>-th
  // local search have been made
  
  FILE* f;
  int i;

  printf("write opt tour on <%s> ... ", fname );
  fflush( stdout );
  f = fopen( fname, "w" );
  fprintf( f, "%d\n", n );
  for ( i = 0; i <= n; i++ ) 
    fprintf( f, "%d ", pi[i] );
  
  fprintf(f, "\n");
  fprintf(f, "TOURVAL\n%lg\n", tval );
  printf("done\n");
  fclose( f );
  fflush(stdout);
}

void writeTSP( char* fname ) {

  // writes the read or generated graph in a special TSP format (compatible wih LKH software)
  FILE* f;
  int i,j;
  
  printf("write graph on <%s> ... ", fname );
  fflush( stdout );
  f = fopen( fname, "w" );
  fprintf(f, "TYPE: TSP\n");
  fprintf(f, "DIMENSION: %d\n", n);
  fprintf(f, "EDGE_WEIGHT_TYPE: EXPLICIT\n");
  fprintf(f, "EDGE_WEIGHT_FORMAT: UPPER_ROW\n");
  fprintf(f, "EDGE_WEIGHT_SECTION\n");
  for ( i = 0; i < n; i++ )
    for (j = i+1; j < n; j++)
      {
	fprintf(f, "%d", (int)round(edgecost(i, j)) );
	if (j != n-1)
	  fprintf(f, " ");
	else
	  fprintf(f, "\n");
      }
  printf("done\n");
  fclose( f );
  fflush(stdout);
}


/***********************************************************
  PERMUTATIONS UPDATES AND PROCS
**********************************************************/

void reverse( int* pi, int a, int b ) {

  // reversal of segment between <a> and <b>, (with a < b), both included

  int tmp;

  while ( a < b ) {
    tmp = pi[a];
    pi[a++] = pi[b];
    pi[b--] = tmp;
  }

#ifdef OLDVER
  int i;
  int len = b - a + 1;
  for ( i = 0; i < len / 2; i++ ) {
    tmp = pi[a + i];
    pi[a + i] = pi[b - i];
    pi[b - i] = tmp;
  }
#endif
}

void reverseExternal( int* pi, int a, int b ) {

  // reversal of segment between <a> and <b>, both included
  // here a > b and the part between [b+1,...,a-1] remains the same
  // while the rest is reversed
  // In total, n - a + b + 1 elements get reversed for a total of <nswaps> swaps.
  
  int tmp, i, nswaps;

  nswaps = ( n - (a - b - 1) ) / 2;
  //  printf("a %d  b %d nswaps %d\n", a, b, nswaps );
  //  scanf("%d", &tmp);
  for ( i = 0; i < nswaps; i++ ) {
    //    printf("swap %d and %d\n", a, b );
    //    scanf("%d", &tmp);
    tmp = pi[a];
    pi[a] = pi[b];
    pi[b] = tmp;
    a = (a + 1) % n;
    b = (b - 1 + n) % n;
  }
}


void restore0( int* pi, int pos0 ) {

  // pi[] is a permutation such that pi[pos0] = 0
  // it shifts it so that pi[0] = 0
  
  int i, j;
  
  // (pi[pos0] == 0) must be true
  i = pos0;
  for ( j = 0; j < n; j++ ) {
    bufferPerm[j] = pi[i];
    i = (i + 1) % n;
  }
  for ( j = 0; j < n; j++ )
    pi[j] = bufferPerm[j];
  pi[n] = 0;
}


