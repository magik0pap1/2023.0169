#ifndef HEAPSHEADER
#define HEAPSHEADER

// general macros

#define ALLPAIRS 0
#define BIGGER 1
#define PAIRPOSITIVE 2


#define MOVECLOCK(x,off) ( (x+off) % n )
#define MOVECOUNTER(x,off) ( (x+n-off) % n )
#define MOVE(x,off) ( off >= 0? (x+off) % n : (x+n+off) % n )

#define MIN(x,y) ( x > y ? y : x )
#define MAX(x,y) ( x > y ? x : y )

void print(const char *fmt, ...);

// data types

typedef unsigned char boole;

typedef struct{
  double value;
  int x; // derived from one of i,j,k,h with possibly an offset
  int y; // derived from one of i,j,k,h with possibly an offset
  int paircode;  // 01, 02, 03, 12, 13, 23 a number which specifies the two index types, in order
  int rscheme;   // identifies the reinsertion scheme
} TPAIR;  


typedef double (*fgenericTau)(int*, int, int);

typedef struct{
  TPAIR* heap;

  /**** i seguenti campi sono usati solo nella versione vecchia dell'algoritmo, in cui
        c'erano 3 heaps. Adesso, se l'heap è singolo come nel nuovo paper,  sono ignorati.

        update 27/5/19 . adesso che abbiamo messo insieme 3 opt e 4 opt, la versione 
        di 4 opt che usa i 4 heap come quella di 3 opt ne usava 3 non è stata aggiornata
        ad usare uno heap solo, come invece abbiamo fatto per 3 opt. Quindi reesta ancora
        questa parte ai fini di creare statistiche nel paper. Siccome il metodo non impostato
        come il vecchio 3 opt è meglio (ossia il metodo master/servant) alla fine questa parte
        si potrà togliere, così come le procedure che la usano

  *****/
  
  int who1st;    // 0 = i,  1 = j,  2 = k,  3 = h
  int who2nd;    // 0 = i,  1 = j,  2 = k,  3 = h
  int paircode;  // 01, 02, 03, 12, 13, 23 a number which specifies the two index types, in order
  int threeopt;  // heap used for 3OPT (1) or for 4OPT (0)
  fgenericTau f;
  int rscheme;   // identifies the reinsertion scheme
} THEAP;

extern int POPCOUNTER; /* global variable, used as semaphor, to count how many heap extractions there
                    have been in a certain stage of the algorithm. It gets reset and tested from outside,
                    it is increased at each extraction */

void initializeHeap( int doAlloc, THEAP* H, int* pi, int w1, int w2, int opt3, fgenericTau fct, boole sorted, int filter, double threshold );
void initializeSingleHeap( int doAlloc, THEAP* H, int* pi, int totschemes, int* schemes, double threshold );
double valueAt( TPAIR* heap, int ndx );
int isRealHeap( TPAIR* heap );
int heapsize( TPAIR* heap );
int getFirstElem( TPAIR* heap );
void nextFirstElem( TPAIR* heap );
void decrHeapSize( TPAIR* heap );
void setHeapSize( TPAIR* heap, int size );
void heapify( TPAIR* heap, int i );
void setRealHeap( TPAIR* heap, boole real );
void getMaxHeap( TPAIR* maxel, TPAIR* heap );
double getMaxHeapValue( THEAP* heap );
void showHeap( TPAIR* heap );
void showTHeap( THEAP* H );
void freeHeapMem( THEAP* H );
void duplicateHeap( TPAIR* heap1, TPAIR* heap2 );
void setProperIndices( THEAP* theap, TPAIR* pnt, int* pi, int* pj, int* pk, int* ph );
void initialize2OPTHeap( THEAP* H, int* pi, float (*keyval)(int*, int));
void getMax2OPTHeap( double* value, int* index, TPAIR* heap );


#endif
