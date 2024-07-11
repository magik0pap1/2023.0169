#ifndef TSPMAIN
#define TSPMAIN

/*
   A TSP for a complete graph is described by (weights, size),
   where size is the number of graph nodes
   and weights is a double[size][size] symmetric matrix of edge costs.
   
   A tour is described by a (size+1)array of indexes belonging to [0,size-1],
   with tour[0]=tour[size]=0.
   
   A K-OPT move is described by K removed edges and K inserted edges
   (an edge can be both removed and inserted), each edge is a couple
   of indexes in [0,size-1].
   No order is required among edges nor between indexes in an edge.
*/

#include "stopwatch.h"

/******************************************

        CONSTANTS AND MACROS

*******************************************/

#define MAXN 34000                /* max n. of nodes */
#define MAXPOSSIBLEMOVES 1000

#define EPS         0.00000001
#define INFINITE    1000000000
#define true        1
#define false       0
#define UNIFORM     0
#define EUCLIDEAN   1

#define UNDEF -9999  /* code for values which should be undefined */

// code for various k-opt algos
#define SEARCH_4OPT_BF           0    /* search is brute force */
#define SEARCH_SF0               1    /* search is smart force with master/servant heap, no perfect split */
#define SEARCH_SF1               2    /* search is smart force with master/servant heap, yes perfect split */
#define SEARCH_WOEGINGER         3    /* search is woeginger dyn prog */
#define SEARCH_GLOVER            4    /* search is glover dyn prog */
#define SEARCH_SF_SIMILAR3OPT    5    /* search is smart force with 4 heaps same as 3OPT */
#define SEARCH_4OPT_TOTOPTIONS   6    /* how many variants for 4OPT */

#define SEARCH_3OPT_BF           0
#define SEARCH_3OPT_SF           1
#define SEARCH_3OPT_TOTOPTIONS   2    /* how many variants for 3OPT */

#define SEARCH_2OPT_BF           0
#define SEARCH_2OPT_SF           1   /* org AG (pivot+expand) */
#define SEARCH_2OPT_SF_NODUP     2   /* org AG (pivot+expand) w/o duplicate moves */
#define SEARCH_2OPT_SF_LEXI      3   /* A_L in resubm with lexi on sorted edges */
#define SEARCH_2OPT_SF_SORTPAIR  4   /* A_g in resubm with heap for pps */
#define SEARCH_2OPT_H_SORTPAIR   5   /* alg H(delta) in resubmission  */
#define SEARCH_2OPT_TOTOPTIONS   6    /* how many variants for 2OPT */



// #define PRINTON

#ifdef PRINTON
  #define dprintf printf
#else
  #define dprintf(...) ((void)0) //strip out PRINT instructions from code
#endif

#define dequal( x, y ) (( x  <= y + EPS ) && ( x => y - EPS ))

//#define PRED(x) ( x > 0 ? x - 1 : n - 1)
//#define SUCC(x) ( x+1 < n ? x + 1 : 0 )

#define PRED(x) ( PREDV[x] )
#define SUCC(x) ( SUCCV[x] )




/*****************************************

           TYPES

******************************************/

typedef struct {
  int found;                // 0 means not found
  double delta;             // the value of the move
  int selection[4];         // meaningless if !found
  int scheme;               // meaningless if !found
                            // in 4OPT it is in 1-1 correspondence with (orbit,phi), namely
                            // each element of orbit has also a name as a scheme (see paper, names 1..25)
  int orbitNo;              // used only in 4OPT
  int phi;                  // action on the representative to get current elem of orbit
  //  int* tour;                // meaningless if !found  NON DIREI CHE SERVE QUI. Non capisco cos'è, è
  // il tour prima della mossa?
  int maxHeapSize;          // L value in our model
  int heapExtractions;      // M value in our model
  long int moveEvaluations; // If Best-Impr: approx nxM value in SF model, all triples in BF.
  // If first-impr: it could be anything
  int effectiveExtractions; // number of heap extractions producing an update
  double seconds2Evaluate; // extracting elements from heap and evaluating moves
  double secondsToBuildHeap; // used in 3opt sometimes
  int numberOfGlobalOptima; // used in 3opt sometimes
} OPTMOVE;

typedef struct{

  char TSPFILENAME[50];
  int PERM[MAXPOSSIBLEMOVES];
  int nperm;
  int DETERMINISTIC;
  int NUMSTARTS;
  int BRUTEFORCE;
  int PERTURB;
  int REFRESHRATE;
  double RESHUFFLEPERC;
  int ANALOGOA3OPT;
} TPARAM;



/*****************************************

           EXTERNAL GLOBAL VARS

******************************************/

extern int n;
extern double** cost;
extern int OUTPUTON;
extern long int ncombinations;
extern int bufferPerm[MAXN + 1];
extern float TIMERATIO;  // used for estimating running time of "long" BF moves
extern int RANDOM_GRAPH_TYPE;  /* used when graph is generated at random */

extern double MAX_NODE_U;   // each node i might get assigned 0 <= u[i] < MAX_NODE_U
                            // to then perturbe arcs as cost[i][j] + u[i] + u[j]

extern double TIMELIMIT;            // overall time limit in seconds 
extern double globalTimeStart;   // clock at start of overall procedures

extern TPARAM param;
extern int FINDBEST;

extern int SUCCV[MAXN];
extern int PREDV[MAXN];

extern int SEEDTOUR;  // seed to generate tours and permutations
extern int SEEDCOST;  // seed for random edge costs

/*****************************************

           PROCEDURES

******************************************/

/* 
   Fills tour uniformly at random using seed,
   memory is allocated by caller.
*/
void generateRandomTour(int *tour, int size, int seed);

/*
   Selects either 3OPT or 4OPT using code.
*/
void findMove(double** weights, int size, int best, int smartforce,
		          int code, int* currentTour, int* newTour, OPTMOVE* move);
				  
/*
   Looks for the FIRST/BEST 3opt move in TSP described by (weights, size)
   starting from currentTour and producing a better newTour, filling move struct
   (everything allocated by caller).
   If best == 0 looks for the FIRST move, otherwise the BEST one.
   If the function sets move->found == 0, newTour is not used.

   code = 0 -> 3opt "standard"  (schemes r_1, .., r_4  in paper)
   code = 1 -> 3opt "special 6opt" (schemes r_5...r_8 in paper)
   code = 2 -> 3opt "special 5opt" (scheme r_9 in paper)
   code = 3 -> all quadratic 3opt (schemes r_1...r_8 )
   code = 4 -> all 3opt (schemes r_1...r_9 )

   If <smartforce>=true it uses Smart Force (Heaps) 
   else it uses Brute Force (complete enumeration)
*/

void perform2OPTMove( int* pi, int* sel, int updateSorted );

void find2optMove(int best, int searchType, int* currentTour, OPTMOVE* move);

/*
   Looks for the FIRST/BEST 2opt move in TSP described by globals (cost[][], n)
   starting from currentTour and producing a better newTour, filling move struct
   (everything allocated by caller).
   If best == 0 looks for the FIRST move, otherwise the BEST one.

   If <searchType>=SEARCH_2OPT_SF it uses Smart Force (Heaps) 
   else if = SEARCH_2OPT_BF it uses Brute Force (complete enumeration)
*/

void perform3OPTMove( int* pi, int* sel, int scheme );

void find3optMove( int searchType, int code, int* currentTour, OPTMOVE* move);

/*
   Looks for the FIRST/BEST 3opt move in TSP described by globals (cost[][], n)
   starting from currentTour and producing a better newTour, filling move struct
   (everything allocated by caller).

   code = 0 -> 3opt "standard"  (schemes r_1, .., r_4  in paper)
   code = 1 -> 3opt "special 6opt" (schemes r_5...r_8 in paper)
   code = 2 -> 3opt "special 5opt" (scheme r_9 in paper)
   code = 3 -> all quadratic 3opt (schemes r_1...r_8 )
   code = 4 -> all 3opt (schemes r_1...r_9 )

   If <searchType>=SEARCH_3OPT_SF it uses Smart Force (Heaps) 
   else if = SEARCH_3OPT_BF it uses Brute Force (complete enumeration)
*/

void find4optMove(int searchType, int norbits, int* orbitcode,int* currentTour, OPTMOVE* move);

/*
   Looks for the FIRST/BEST 4opt move in TSP described by (cost, n) -global vars-
   --where n is the number of graph nodes  and cost is a 
     double[n][n] symmetric matrix of edge costs--

   starting from currentTour and filling move struct
   (everything allocated by caller).
   If the function sets move->found == 0, newTour is not used.

   Input: 
   - <searchType> is a code, i.e.: 
       = SEARCH_BF -> search is brute force 
       = SEARCH_SF0/_SF1 -> search is smart force with master/servant heap 
       = SEARCH_SF_SIMILAR3OPT -> search is smart force with 4 heaps same as 3OPT, slower than prec
       = SEARCH_WOEGINGER  -> search is woeginger dyn prog
       = SEARCH_GLOVER  -> search is glover dyn prog // added jun 22
   - <orbitcode[]> are the codes of the orbits to try (1..7), in total <norbits>
     (the move will be searched over all elements of the orbit)
   - <currentTour> the tour before the best move is applied

   Output
   - <newTour> the tour after the move is applied  (TOLTO)
   - <move> info e statistics on the best move found

*/

/* Preliminary operations to set-up for data structures and other bootstrap
   operations needed by 3OPT procedures, given number of nodes in graph */
void beforeRunning( int nnodes );

/* Un-do whatever done by beforeRunning(). Free memory */
void afterRunning();

/* with code described above, returns the estimated time of a move specified
   by code over the permutation perm
*/
double estimateTime3optMoveBF( double** weights, int size, int* perm, int code );

/* with code described above, returns the estimated time of a move specified
   by code over the permutation perm
*/
double estimateTime4optMoveBF( double** weights, int size, int* perm, int code );

void describeTime( double secs, char* stringa );
  /* Transforms a number <secs> of seconds into a human readable format, e.g.
     112 -> 1m 52s
     3602 ->1h 2s
     3665 -> 1h 1m 5s
     etc.
  */

void find4optMoveNewStrategy(int* currentTour, int orbitcode, OPTMOVE* move);

#endif
