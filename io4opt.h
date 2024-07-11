#ifndef IO4OPT
#define IO4OPT

extern int PRECISION;  // set decimal precision in costs
extern double CMIN, CMAX;    // costs in [CMIN,CMAX)
extern double xcoor[MAXN], ycoor[MAXN];
extern int COSTI_EUCLIDEI; /* se = 0 usa costi uniformi */



/*********************************

input/output procedures 

**********************************/

double tourValue( int* p );
// computes the value of a tour

double randomCost();

void showCosts();

void freeCost();
// release memory

void setCosts();

double edgecost( int i, int j );
// returns cost[i][j] either from matrix or computed

void readCostsFile( char* fname, int doPerturb );

void perturbCostsOnNodes();  // applies (if needed) perturb c(i,j)+u(i)+u(j)


void readTour( char* fname, int* pi, int n );
// reads a tour of n nodes from a file
// the number of nodes is already known

void writeTour( char* fname, int* pi, int n, int step, int numLS, int islast );
// writes the tour resulting after <step> moves of the <numLS>-th
// local search have been made
// <islast>=true if it is the last of a convergence

void writeBestTour( char* fname, int n, int* pi, double tval );
// writes the best tour found

void writeTSP( char* fname );
// writes the instance in TSPLIB format

// **** PERM MANIPULATION

void reverse( int* pi, int a, int b );
 // reversal of segment between <a> and <b>, w/ a < b, both included

void reverseExternal( int* pi, int a, int b );
 // reversal of segment between <a> and <b>, w/ a > b, both included

void restore0( int* pi, int pos0 );
// makes pi[0] = pi[n] = 0

void rndPerm0( int *pi, int n_ ) ;
/* generates pi (of elements in 0..n-1) uniformly at random,
   but with pi[0] = 0.
   memory for pi is allocated by caller. Also pi[n] = 0 to close the cycle*/

void randomTple( int* sel, int n_, int k );
// generates a random ordered selection of k numbers in 0..n_

double euclDist( double xi, double yi, double xj, double yj );

#endif
