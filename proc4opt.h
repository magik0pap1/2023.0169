#ifndef PROC4OPTHEADER
#define PROC4OPTHEADER

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "stopwatch.h"
#include "heap.h"


// movesInOrbit[i][j] = j-th scheme of orbit i, listed in the order of
// the group action and numbered 1..25 as in paper

extern  const  int movesInOrbit[8][10];

/*****************************************
   TYPES
 *****************************************/


typedef struct{
  double time;
  double moveval;
  int orbit;
  int phi;
  boole bestimpr;
} RECORDMOVE;

typedef struct{
  // una mossa aggregata, ossia una sequenza di mosse della lista, tutte dalla
  // stessa soluzione iniziale.
  // se findbest=false, può essersi fermata a una qualsiasi dei tipi di mossa
  // altrimenti li prova tutti e sceglie il migliore
  
  double time; // tempo complessivo
  RECORDMOVE* list; // lista dei record delle mosse provate
  int listlen; // quanti elementi in list
  double moveval; // valore della mossa migliore nella lista
  int orbit;  // chi è la mossa migliore (orbit,phi)
  int phi;
  boole bestimpr; // era un bestimp o firstimp
} RECORDAGGREGATE;

#define MAXCONVERGENCELEN 30000  /* in 30000 steps we must have found a local optimum, or we abort */

typedef struct{
  double startval;
  int startpi[MAXN];
  double endval;
  int endpi[MAXN];
  double time;
  int numsteps;
  RECORDAGGREGATE step[MAXCONVERGENCELEN];
} RECORDCONVERGENCE;

typedef struct{

  char datafilename[100];
  int numconvergences;
  RECORDCONVERGENCE* pconv;
  int* moveperm;
  int nperm;
  boole findbest;
  int rndseed;
  boole bruteforce;
  boole deterministic;
  int bestpi[MAXN + 1];
  double bestval;
  double time;
} RECORDRUN;

typedef struct{

  int size;
  int phi[8][4];       // the permutations for each element of the orbit
  boole flip[8];       // the permutation includes a reflection or not
  fgenericTau fct[5];  // functions f^1,f^2,f^3 in paper
  int phiHeapCode[5];  // aggiunto ott 2022, per ridurre gli heap usati a 17 casi complessivi
  int who[5][2];       // who[f][0] is the id of the 1st argument of the f-th function (0=i, 1=j, 2=k, 3=h)
                       // who[f][1]   "              2nd   "
  fgenericTau g1, g2;  // used for orbits with perfect split
} TORBIT;


/*****************************************
   VARIABLES
 *****************************************/

extern TORBIT ORBIT4[7];


/* global heaps used by our smart force procedures */

extern THEAP heapF12, heapF34, heapSingle;
extern TPAIR* sortedheap[2];
extern TPAIR* tmpHeap;



/*****************************************
   PROCEDURES
 *****************************************/

void* mycalloc( int numelem, int sizeelem );
void myfree( void* p );

void rndPerm0( int *pi, int n_ );
double val( int* p );
double tauplus( int* pi, int x, int y );
double tauplus_00( int* pi, int x, int y );
double tauplus_01( int* pi, int x, int y );
double tauminus( int* pi, int x, int y );
double tauminus_11( int* pi, int x, int y );
double tauminus_12( int* pi, int x, int y );
double function3OPT(int* pi, int x, int y, int rscheme, int ndx0, int ndx1);
void prepareOrbits4OPT();
void perform4OPTMove( int* pi, int* sel, int norbit, int nphi );
double evaluate4OPT( int* pi, int i, int j, int k, int h, int norbit, int nphi );
int removeReflection( int norb, int s );
void reflectSelection( int* sel );
void perform4OPTMoveByScheme( int* pi, int* sel, int nscheme );

void freeCost();
double edgelist(int* pi, int* sel, int norbit, int rot );

double find4OPTGloverDynProgOrbit7(int* pi, int * ii, int* jj, int* kk, int* hh, int* rot );
double find4OPTGloverDynProgOrbit6(int* pi, int * ii, int* jj, int* kk, int* hh, int* rot );
void allocateGloverMem();  // allocate memory for the D.P. tables of glover algorithm
void disposeGloverMem();



#endif
