/*  Heap.c

    procedures for managing heaps / sorted arrays

    New version, for use with orbits

    last check : 27 may 2019 good for 3-opt and 4-opt
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdarg.h>
#include "tspmain.h"
// #include "proc3opt.h"
#include "proc4opt.h"

#define STRONGDEBUG 1

/***** 
       can be both a true heap or a sorted array.
       if a heap is indeed a sorted array, then it first unused element (the current "top")
       is stored in heap[0].y. The size of the heap in is heap[0].x
       The fact that the heap is a real heap or not is in heap[0].value (=1 real heap, =0 sorted array)

 *****/

int POPCOUNTER = 0;  /* global variables, used as semaphor, to count how many heap extractions there
                    have been in a certain stage of the algorithm. It gets reset and tested from outside,
                    it is increased at each extraction */

int NHEAP = 0; // how many heaps have been created

double valueAt( TPAIR* heap, int ndx ) {

  // returns the value field of the element heap[ndx]

  return heap[ndx].value;
}

void setTopElement( TPAIR* heap, int val ) {

  // when heap is a sorted array, it tells where the current first element is
  heap[0].y = val;
}

void setRealHeap( TPAIR* heap, boole real ){
  heap[0].value = real;
  if ( !real )
    setTopElement( heap, 1 );
}

int isRealHeap( TPAIR* heap ){
  return (int)heap[0].value;
}

int heapsize( TPAIR* heap ) {

  return heap[0].x;  // field x of H[0] used to store the size of the heap
}

int posFirstElem( TPAIR* heap ) {

  return heap[0].y;  // field y of H[0] used to store the position of the 1st elem of sorted array
}

void nextFirstElem( TPAIR* heap ) {

  heap[0].y +=1;  // moves to the next first elem. 
}

void  decrHeapSize( TPAIR* heap ) {
  // decreases the size of the heap

  heap[0].x -= 1;
}

void  setHeapSize( TPAIR* heap, int size ) {
  // sets the size of the heap

  heap[0].x = size;
}

void heapify( TPAIR* heap, int i ) {  

  // enforces the heap property in the subtree rooted at i

  int l, r, largest;
  TPAIR tmp;

  l = 2 * i; // left(i)
  r = 2 * i + 1; // right(i)
  if ( (l <= heapsize(heap)) && (heap[l].value > heap[i].value ) )
    largest = l;
  else 
    largest = i;

  if ( (r <= heapsize(heap)) && (heap[r].value > heap[largest].value ) )
    largest = r;

  if ( largest != i) {
    tmp = heap[i];
    heap[i] = heap[largest];
    heap[largest] = tmp;
    heapify( heap, largest );
  }
}

void getMaxHeap( TPAIR* maxel, TPAIR* heap ) {


  // if <trueHeap> we have an actual heap. Then:
  // extracts the max-element (the root) of the heap
  // and rebuilds the heap

  // if <trueHeap>=false, then the "heap" is in fact a sorted array, and we just pick
  // its first unused element, which is addressed through the field currFirst

  
  if ( isRealHeap(heap) ) {
    //  printf("get at pos %d on %d elem of value %g\n", 1, heapsize(heap), heap[1].value );

    POPCOUNTER++;
    maxel[0] = heap[1];
    heap[1] = heap[heapsize(heap)];
    decrHeapSize( heap );
    heapify( heap, 1 );
  }
  else {
    maxel[0] = heap[posFirstElem(heap)];
    nextFirstElem( heap );
    decrHeapSize( heap );
  }
}

double getMaxHeapValue( THEAP* H ) {

  // returns the value of the largest element of the heap

  // if <trueHeap> we have an actual heap. 
  // if <trueHeap>=false, then the "heap" is in fact a sorted array, and we just look at
  // its first unused element, which is addressed through the field currFirst

  TPAIR* h;

  h = H->heap;

  if ( isRealHeap( h ) ) 
    return h[1].value;
  else 
    return h[posFirstElem(h)].value;
}

void duplicateHeap( TPAIR* heap1, TPAIR* heap2 ) {
  // makes a copy of heap2 onto heap1

  int i; 

  for ( i = 0; i <= heapsize(heap2); i++ )
    heap1[i] = heap2[i];
}

void showHeap( TPAIR* heap ) {

  int i;

  if ( isRealHeap(heap) )
    printf("REAL HEAP - size %d\n", heapsize(heap));
  else 
    printf("SORTED ARRAY - size %d top at %d\n", heapsize(heap), heap[0].y );
  for ( i = 1; i <= heapsize( heap ); i++ )
    printf("HEAP el %d is  x %d  y %d  val %g  r %d  label %d\n", i, heap[i].x, heap[i].y, heap[i].value,
	                                                          heap[i].rscheme, heap[i].paircode );

  printf("...<ENDHEAP>..."); 
  //  scanf("%d", &i);
}

void setelement( TPAIR* heap, int a, int b, double v, int cnt ) {

  //printf("insert (%d %d %g) at %d\n", a, b, v, cnt);
  heap[cnt].x = a;
  heap[cnt].y = b;
  heap[cnt].value = v;
}

void setelementSingleHeap( TPAIR* heap, int a, int b, double v, int pcode, int rsch, int cnt ) {

  //  printf("insert (%d %d %g) at %d\n", a, b, v, cnt);
  heap[cnt].x = a;
  heap[cnt].y = b;
  heap[cnt].value = v;
  heap[cnt].paircode = pcode;
  heap[cnt].rscheme = rsch;
}

void filterElements( TPAIR* heap ) {

  // (this should be called when an heap was created with the PAIRPOSITIVE option)
  // given an array of pairs (elem 0 contains info only) removes
  // all elements E such that it is impossible that val(E)+val(X) >0
  // for any X in the array. Namely if M is the max val(X), removes
  // E if val(E) + M <= 0

  double M;
  int i, cnt;

  M = -INFINITE;
  for ( i = 1; i <= heapsize(heap); i++ )
    if ( heap[i].value > M )
      M = heap[i].value;

  cnt = heapsize( heap );
  i = 1;
  while ( i <= cnt ) {
    if ( heap[i].value + M <= EPS ) {
      heap[i] = heap[cnt];
      cnt--;
    }
    else
      i++;
  }
  setHeapSize( heap, cnt );
}

void initializeSingleHeap( int doAlloc, THEAP* H, int* pi, int totschemes, int* schemes, double threshold ) {
  /* Creates the partial moves heap

    INPUT : totschemes = how many schemes are considered
                        (for each scheme we will run over all pairs (x,y) of type 01, 12, 02
            schemes[] : the specific schemes


  // if doAlloc = TRUE
  // the memory for the heap is allocated here and can be released by calling freeHeapMem()
  // if doAlloc = FALSE
  // the memory for heap is allocated by caller

  // the actual heap is stored in the array H->heap[] starting at index 1.
  // Left(i) = 2*i, right(i) = 2*i +1
  // father(i) = i/2
  // the size of the heap is in the element 0 and accessed through heapsize()


  */


  int i, j, k, nr, r, cnt;
  double fval;

  if ( doAlloc ) {
    H->heap = (TPAIR*) mycalloc( 1 + totschemes * 3 * n * (n-1) / 2, sizeof(TPAIR));
    assert(H->heap);
  }

  cnt = 0;
  for ( nr = 0; nr < totschemes; nr++ ) {
    r = schemes[nr];

    // enumera le coppie (i,j)

    for ( i = 0; i <= n - 5; i++ )
      for ( j = i + 2; j <= n - 3 - (i==0); j++ ) {
	fval = function3OPT(pi, i, j, r, 0, 1);
	if ( fval > threshold + EPS ) 
	  setelementSingleHeap( H->heap, i, j, fval, 01, r, ++cnt);
      }

    // enumera le coppie (j,k)
    
    for ( j = 2; j <= n - 3; j++ )
      for ( k = j + 2; k <= n - 1 - (j==2); k++ ) {
	fval = function3OPT(pi, j, k, r, 1, 2);
	if ( fval > threshold + EPS ) 
	  setelementSingleHeap( H->heap, j, k, fval, 12, r, ++cnt);
      }

    // enumera le coppie (i,k)
    
    for ( i = 0; i <= n - 5; i++ )
      for ( k = i + 4; k <= n - 1 - (i==0); k++ ) {
	fval = function3OPT(pi, i, k, r, 0, 2);
	if ( fval > threshold + EPS ) 
	  setelementSingleHeap( H->heap, i, k, fval, 02, r, ++cnt);
      }
  }

  setHeapSize( H->heap, cnt ); // uses field H[0].x to store size

  // makes it into a heap
  for ( i = cnt / 2; i >= 1; i-- ) 
    heapify( H->heap,  i );
  
  setRealHeap( H->heap, true );

  //  showHeap( H->heap );
}

int compare (const void * a, const void * b) //++++
{
  if ( ((TPAIR*)a)->value <  ((TPAIR*)b)->value ) return 1;
  if ( ((TPAIR*)a)->value >  ((TPAIR*)b)->value ) return -1;
  
  return 0;
}


void initializeHeap( int doAlloc, THEAP* H, int* pi, int w1, int w2, int opt3, fgenericTau fct, boole sorted, int filter, double threshold ){

  // creates an heap to store pairs in which <w1> identifies the
  // type of 1st index (0=i, 1=j, 2=k, 3=h) argument of the function <fct>
  // and <w2> of the second.

  // if doAlloc = TRUE
  // the memory for the heap is allocated here and can be released by calling freeHeapMem()
  // if doAlloc = FALSE
  // the memory for heap is allocated by caller

  // the actual heap is stored in the array heap[] starting at index 1.
  // Left(i) = 2*i, right(i) = 2*i +1
  // father(i) = i/2
  // the size of the heap is in the element 0 and accessed through heapsize()

  // the value of the pairs in the heap is computed by using the
  // function fct which is passed to the procedure.
  // (for instance fct can be tau^+, tau^- etc)

  // the parameter <sorted> if true, makes the input into a sorted array
  // if <sorted>=false, the heap is initially unsorted and has simply the structure of a max-heap

  // the parameter <filter> can be 
  // ALLELEM (0)     : put all elements in heap
  // BIGGER  (1)    : put in heap only elements with value > threshold
  // PAIRPOSITIVE (2): put in heap only elements which can yield a sum > 0 in pairs. I.e., only elements s.t.
  //    value + maxvalue > 0
  
  // if <opt3> the heap will be used for 3opt moves, i.e. with
  // pairs taken over 3 indices (i,j,k), otherwise for 4OPT, and pairs are
  // taken over 4 indices, (i,j,k,h)


  int a, b, i, cnt, wdiff;
  double value;
  int min1, max1, max2, minW12, maxW12;
  int ordered;

  //  printf("crea heap n. %d\n", ++NHEAP );


  dprintf("creating heap w1 %d  w2 %d  filter %d  thres %g\n", w1, w2, filter, threshold);
  minW12 = MIN(w1,w2);
  maxW12 = MAX(w1,w2);
  ordered = ( w1 < w2 );

  H->paircode = minW12 * 10 + maxW12;
  H->who1st = w1;  // 0 = i,  1 = j,  2 = k,  3 = h
  H->who2nd = w2;  // 0 = i,  1 = j,  2 = k,  3 = h

  H->threeopt = opt3;
  H->f = fct;
  wdiff = maxW12 - minW12;

  if ( doAlloc ) {
    //  printf("*** In initializeHeap ALLOCATE H->heap : ");
    H->heap = (TPAIR*) mycalloc( 1 + n * (n-1) / 2, sizeof(TPAIR));
    assert(H->heap);
  }
  else {
    assert(H->heap);
  }
  cnt = 0;

  switch ( minW12 ) {
  case 0:
    min1 = 0;
    max1 = n - 5 - 2*(!opt3);
    break;
  case 1:
    min1 = 2;
    max1 = n - 3 -2*(!opt3);
    break;
  case 2:
    min1 = 4;
    max1 = n - 1 -2*(!opt3);
    break;
  default:
    assert(false);
  }

  switch ( maxW12 ) {
  case 1:
    max2 = n - 3 -2*(!opt3);
    break;
  case 2:
    max2 = n - 1 -2*(!opt3);
    break;
  case 3:
    assert(!opt3 && n>=8);
    max2 = n - 1;
    break;
  default:
    printf("Error at MIN %d MAX %d\n", minW12, maxW12);
    assert(false);
  }

  dprintf("List all values of elem %d from min %d to max %d\n", minW12, min1, max1);
  dprintf("List all values of elem %d to max %d\n", maxW12, max2);

  for ( a = min1; a <= max1; a++ )
    for ( b = a + 2 * wdiff; b <= max2 - (a==min1); b++ ) {
      if ( ordered )
        value = (*fct)( pi, a, b );
      else
	value = (*fct)( pi, b, a );
	
      dprintf("Loop: a %d  b %d  val %g\n", a, b, value);
      if ( filter != 1 || value > threshold + EPS )
	setelement( H->heap, a, b, value, ++cnt);
    }
  
  setHeapSize( H->heap, cnt ); // uses field H[0].x to store size

  if ( filter == 2 ) {
    //    printf("before filter size %d\n", heapsize(heap));
    filterElements( H->heap );
    //    printf("after filter  size %d\n", heapsize(H->heap));
    //    scanf("%d", &i);
  }

  if (!sorted) {
    // makes it into a heap
    for ( i = cnt / 2; i >= 1; i-- ) 
      heapify( H->heap,  i );
    
    setRealHeap( H->heap, true );
  }
  else {
    qsort ( H->heap + 1, cnt, sizeof(TPAIR), compare);
    setRealHeap( H->heap, false );
  }
}

void freeHeapMem( THEAP* H ) {
  if ( H->heap ) {
    myfree( H->heap );
  }
  H->heap = NULL;
}


void showTHeap( THEAP* H ) {
  char letter[4] = {'I', 'J', 'K', 'H'};
  char let1 =  letter[H->who1st];
  char let2 =  letter[H->who2nd];
    
  printf("\n***(3OPT %d) Heap structure for pair %02d ", H->threeopt, H->paircode);
  printf(" ( %c", let1 );
  /*
  if ( off1 > 0 )
    printf("+%d ", off1 );
  else if ( off1 < 0 )
    printf("-%d ", -off1 );
  */
  printf("  %c", let2 );
  /*
  if ( off2 > 0 )
    printf("+%d ", off2 );
  else if ( off2 < 0 )
    printf("-%d ", -off2 );
  */
  printf(")\n");
  showHeap( H->heap );
}

void setProperIndices( THEAP* theap, TPAIR* pnt, int* pi, int* pj, int* pk, int* ph ) {

  // looks at element HEAP[cnt] and sets two of the 4 indices to the corresponding
  // values. The heap knows which indices it concerns
  // The values of the other two indices is left untouched (it could have been set previously)
  // It is assumed that indices are initialized at -1 outside the call

  TPAIR top;

  /*
  printf("beg set indices at cnt %d on  i %d j%d k %d h %d\n", cnt, *pi, *pj, *pk, *ph);
  showTHeap( theap );
  */

  top = *pnt;
  switch ( theap->paircode ) {
  case 01:
    assert( *pi == -1 && *pj == -1 );
    *pi = MIN( top.x, top.y );
    *pj = MAX( top.x, top.y );
    break;

  case 02:
    assert( *pi == -1 && *pk == -1 );
    *pi = MIN( top.x, top.y );
    *pk = MAX( top.x, top.y );
    break;

  case 03:
    assert( *pi == -1 && *ph == -1 );
    *pi = MIN( top.x, top.y );
    *ph = MAX( top.x, top.y );
    break;

  case 12:
    assert( *pj == -1 && *pk == -1 );
    *pj = MIN( top.x, top.y );
    *pk = MAX( top.x, top.y );
    break;

  case 13:
    assert( *pj == -1 && *ph == -1 );
    *pj = MIN( top.x, top.y );
    *ph = MAX( top.x, top.y );
    break;

  case 23:
    assert( *pk == -1 && *ph == -1 );
    *pk = MIN( top.x, top.y );
    *ph = MAX( top.x, top.y );
    break;

  default:
    assert(false);
  }

  /*
  printf("end set indices at cnt %d on  i %d j%d k %d h %d\n", cnt, *pi, *pj, *pk, *ph);
  */
}

void initialize2OPTHeap( THEAP* H, int* pi, float (*keyval)(int*, int)) {
  /* crea l'heap per il nuovo 2 opt

    Usa i vecchi heap del 3 e 4 opt, ma solo riempiendo alcuni campi, non tutti.

  */


  int i, cnt;

  H->heap = (TPAIR*) calloc( 1 + n, sizeof(TPAIR));
  assert(H->heap);

  cnt = 0;
  for ( i = 0; i < n; i++ ) {
    //    setelementSingleHeap( H->heap, i, -1, edgecost(pi[i], pi[SUCC(i)]), -1, -1, ++cnt);
    setelementSingleHeap( H->heap, i, -1, (*keyval)(pi, i), -1, -1, ++cnt);
  }

  setHeapSize( H->heap, cnt ); // uses field H[0].x to store size

  // makes it into a heap
  for ( i = cnt / 2; i >= 1; i-- ) 
    heapify( H->heap,  i );
  
  setRealHeap( H->heap, true );

  //  showHeap( H->heap );
}

void getMax2OPTHeap( double* val, int* index, TPAIR* heap ) {

  TPAIR maxel;
  getMaxHeap( &maxel, heap );
  *val = maxel.value;
  *index = maxel.x;
}
