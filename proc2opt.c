#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "tspmain.h"
#include "io4opt.h"
#include "heap.h"

float minIncident[MAXN];
int ABORT_AT_SWITCH; /* Usato nella funzione searchhbeta, per test statistici
		       al fine di determinare il miglior beta. Andrebbe tolto una volta trovato
		       beta e fatte le tabelle */

int DETAIL_LS_STEPS = false;
int jolly; // used for output and printing options per statistiche

TPAIR* SORTED; /* array of sorted tour edges, used in a version of the 2opt search */
int sortedExists; /* semaphore saying that the array SORTED has already been built and can be used */

int LASTI;

int neverExp[MAXN];  // used to avoid double evaluations (see paper)


double ABS( double x ) {

  if ( x >= -EPS )
    return x;
  else
    return -x;
}


#ifdef VECCHIAVERSIONE

/* questa versione di perform2OPTMove faceva il reverse della regione più corta, per risparmiare
   tempo. In questo modo, però, non si mantiene l'invariante pi0[0]=0 che in alcuni casi sembra
   essere usata da qualche procedura (mi pare dalle 4opt). Quindi la versione viene sostituita con
   una in cui l'elemento 0 non si muove (ossia si inverte soltanto la zona tra i e j , con i < j
*/


void perform2OPTMove( int* pi, int* sel ) {
  
  // effettua il 2-opt specificato da sel[0], sel[1]
  // Dovendo rovesciare una parte tra x e y (x<y) si può altrettanto lasciare la parte tra x e y come è
  // e rovesciare il resto (da 0 a x e da y alla fine). Scegliamo di rovesciare quello che
  // tra i due è al max n/2
    
  int ii, jj;
  
  ii = sel[0];
  jj = sel[1];
  if ( ii > jj ) {
    ii = sel[1];
    jj = sel[0];
  }

  if (( jj - ii ) <= n/2) 
    reverse( pi, ii + 1, jj );
  else {
    //    printf("EXTERNAL %d %d  len %d\n", (jj+1)%n, ii, n - jj + ii);
    //    fflush(stdout);
    reverseExternal( pi, (jj+1) % n, ii);
  }
  
  pi[n] = pi[0];
  /*  
  assert( pi[0] == 0 );
  if ( ii == 0 )
    restore0( pi, ii );  // controllare**
  */
}

#else

void perform2OPTMove( int* pi, int* sel, int useSortedArray ) {

  
  // performs the 2-opt move corresponding to sel[0], sel[1]
  // if <useSortedArray> it also updates the sorted array of tour edges

  
  int ii, jj;
  int a, b, c;
  double val[MAXN];  // cost of tour edges
  int    tail[MAXN]; // tail of tour edges
  double cmp[3];     // cost of the two entering edges
  int    etail[2];
  

  ii = sel[0];
  jj = sel[1];
  if ( ii > jj ) {
    ii = sel[1];
    jj = sel[0];
  }


  if ( useSortedArray ) {

    if ( edgecost(pi[ii], pi[jj]) > edgecost(pi[SUCC(ii)], pi[SUCC(jj)]) ) {
      etail[0] = ii;
      cmp[0] = edgecost(pi[ii], pi[jj]);
      etail[1] = jj;  // after the reversal, this will be jj-th node
      cmp[1] = edgecost(pi[SUCC(ii)], pi[SUCC(jj)]);
    }
    else {
      etail[0] = jj;
      cmp[0] = edgecost(pi[SUCC(ii)], pi[SUCC(jj)]);
      etail[1] = ii;
      cmp[1] = edgecost(pi[ii], pi[jj]);
    }
    cmp[2] = -INFINITE;  // dummy value used as a trick to stop after inserting 2 edges

  
    /*
    
    printf("eseguo la 2 opt %d (%d) %d (%d)\n", ii, pi[ii], jj, pi[jj] );
    printf("entering (%d %d) cost %g  and (%d %d) cost %g\n", pi[ii], pi[jj], cost[pi[ii]][pi[jj]],
	   pi[SUCC(ii)], pi[SUCC(jj)], cost[pi[SUCC(ii)]][pi[SUCC(jj)]]);
    printf("TOUR BEF :\n");
    for ( a=0; a<n; a++ )
      printf("Edge n.%d (%d %d) cost %g\n", a, pi[a], pi[SUCC(a)], cost[pi[a]][pi[SUCC(a)]]); 
    printf("\n");

    int i, j;
    printf("sorted edges BEF :\n");
    for ( a=1; a<=n; a++ ) {
      i = SORTED[a].x;
      j = SUCC(i);
      printf("SORT Edge n.%d (%d %d) cost %g %g (NDX %d %d)\n", a, pi[i], pi[j], cost[pi[i]][pi[j]], SORTED[a].value, i, j);
    }
    printf("\n");
    */
  }
  

  reverse( pi, ii + 1, jj );
  
  pi[n] = pi[0];

  if ( useSortedArray ) {

    /*
    printf("\n****\n\n");
    printf("TOUR AFT :\n");
    for ( a=0; a<n; a++ )
      printf("Edge n.%d (%d %d) cost %g\n", a, pi[a], pi[SUCC(a)], cost[pi[a]][pi[SUCC(a)]]); 
    printf("\n");
    */

    a = 1;  // counter on the old sorted array (goes from 1(!) to n)
    b = 0;  // counter on the entering edges
    c = 0;  // counter on the new sorted array
    while ( c < n ) {
      if ( a <= n && (SORTED[a].x == ii || SORTED[a].x == jj) ) { // skip removed edges
	a++;
	continue;
      }
      /*
      if ( a > n ) 
	printf("no cmp S[%d]. simply insert %g\n", a, cmp[b]);
      else
	printf("cmp S[%d] val %g with insert %g\n", a, SORTED[a].value, cmp[b]);
      */

      if ( a <= n && SORTED[a].value > cmp[b] ) {
	tail[c] = SORTED[a].x;
	val[c] = SORTED[a].value;
	//	printf("keep tail[%d]=%d\n", c, tail[c]);
	if ( tail[c] > ii && tail[c] < jj ) {

	  // printf("Should reverse at tail %d (node %d) \n", tail[c], pi[tail[c]] );
	  tail[c] = ii + ABS( tail[c] - jj );
	}
	a++;
      }
      else {
	// printf("ins here\n");
	tail[c] = etail[b];
	val[c] = cmp[b];
	b++;
      }
      c++;
      //      scanf("%d", &i);
    }
    
    //    printf("sorted edges AFT :\n");
    for ( a=0; a<n; a++ ) {
      SORTED[a+1].x = tail[a];
      SORTED[a+1].value = val[a];
      //      ii = SORTED[a+1].x;
      //      jj = SUCC(ii);
      //      printf("SORT Edge n.%d (%d %d) cost %g (NDX %d %d)\n", a, pi[ii], pi[jj], cost[pi[ii]][pi[jj]], ii, jj);
    }
    //    printf("\n");
  }

}

#endif

void find2OPTMoveBF(int* pi, int* sel, double* Delta, int* nmoves ) {
  
  int i, j, ii, jj, offset, besti = -1, bestj = -1;
  double delta, best;
  int evaluations = 0; 
  /* offset is used not to start always looking from the same point
     because otherwise, when First Improvement was not fould e.g. for all
     moves up to (i,j) it still tries all moves up to (i,j) in
     the same order. So it starts looking from where it left the step earlier
  */
  
  //  best = 0;  // must find a move of value > 0
  best = -INFINITE;  // facciamogliela trovare comunque, anche se non è improving
  //  if ( !FINDBEST )  offset = lrand48() % n;
  offset = LASTI;
    //    printf("offset %d\n", offset);
    //  scanf("%d", &i);
    
  for ( ii = 0; ii < n - 2; ii++ ) 
    for ( jj = ii + 2; jj < n - (ii == 0) ; jj++ ) {  
      evaluations++;
      i = (offset + ii) % n;
      j = (offset + jj) % n;
      // printf("i %d j %d\n", i, j);
      delta = edgecost(pi[i], pi[SUCC(i)]) + edgecost(pi[j], pi[SUCC(j)]) -
	(edgecost(pi[i], pi[j]) + edgecost(pi[SUCC(i)], pi[SUCC(j)]));
      if ( delta > best + EPS ) {
	//	printf("new best %f\n", best);
	best = delta;
	besti = i;
	bestj = j;
	if ( !FINDBEST && delta > EPS ) {
	  LASTI = i;
	  goto exit2BFloop;
	}
      }
    }

 exit2BFloop:

  if ( besti < bestj ) {
    sel[0] = besti;
    sel[1] = bestj;
  }
  else {
    sel[0] = bestj;
    sel[1] = besti;
  }
    
  *Delta = best;
  *nmoves = evaluations;
}


float tauBasic( int* pi, int i ) {

  // versione "standard" del valore per verificare se un nodo è good e va espanso,
  // ossia c(i,i+1)
  
  /*
  printf("i   %d pi(%d) = %d\n", i, i, pi[i] );
  printf("i+1 %d pi(%d) = %d\n", SUCC(i), SUCC(i), pi[SUCC(i)] );
  printf("c(%d,%d) %g\n", pi[i], pi[SUCC(i)], edgecost(pi[i], pi[SUCC(i)]));
  */

  return edgecost(pi[i], pi[SUCC(i)]);
}

float tauStrong( int* pi, int i ) {

  // versione più forte in cui c'è anche il contributo del minimo arco incidente
 
  /*
  printf("i   %d pi(%d) = %d\n", i, i, pi[i] );
  printf("i+1 %d pi(%d) = %d\n", SUCC(i), SUCC(i), pi[SUCC(i)] );
  printf("mininc(%d) %g\n", pi[i], minIncident[pi[i]]);
  printf("c(%d,%d) %g\n", pi[i], pi[SUCC(i)], edgecost(pi[i], pi[SUCC(i)]));
  printf("taunew(%d) %g\n\n", i, edgecost(pi[i], pi[SUCC(i)]) - 0.5 * (minIncident[pi[i]] + minIncident[pi[SUCC(i)]]));
  */

  //  printf("tau new for %d (%d)  : edg (%d %d) %g  min(%d) %g  min(%d) %g\n", i, pi[i], pi[i], pi[SUCC(i)], edgecost(pi[i], pi[SUCC(i)]), pi[i], minIncident[pi[i]], pi[SUCC(i)],  minIncident[pi[SUCC(i)]]);
	 
  return edgecost(pi[i], pi[SUCC(i)]) - 0.5 * (minIncident[pi[i]] + minIncident[pi[SUCC(i)]] );
}

int min2(int a, int b) {
  if ( a < b )
    return a;
  return b;
}

int max2(int a, int b) {
  if ( a > b )
    return a;
  return b;
}

void find2OPTMoveSFHeap(int* pi, int* sel, double* Delta, int* nmoves, int strongForm ) {
  /* implementazione dell'algoritmo A_g del paper, ver JOC subm

     input: pi (il tour) strongForm (0=normale 1=versione con un tau più forte)
     output: sel (la mossa) Delta (il suo valore)  nmoves (quante mosse valutate)

     old: senza l'uso del vettore neverExp[]
   */
  
  THEAP HEAP;
  TPAIR* heap;
  int i, j, besti = -1, bestj = -1;
  double delta, best;
  double Tau;
  int evaluations = 0; 


  //  double tbefore = timer_elapsed_seconds();

  evaluations = 0;
  
  if ( strongForm )
    initialize2OPTHeap( &HEAP, pi, tauStrong );
  else
    initialize2OPTHeap( &HEAP, pi, tauBasic );
  best = -INFINITE;
  heap = HEAP.heap;
  while ( heapsize( heap ) > 0  ) {

    getMax2OPTHeap( &Tau, &i, heap );
    //    printf("popped (%d %d) cost %f thre %f\n", pi[i], pi[SUCC(i)], Tau, best/2 );
    if ( Tau <= best / 2 + EPS ) 
      break; 

    //    printf("expand popped (%d %d) cost %f thre %f\n", pi[i], pi[SUCC(i)], Tau, best/2 );
    for ( j = 0; j < n; j++ )
      if ( j < i - 1 || j > i + 1 ) { 
	evaluations++;
        delta = edgecost(pi[i], pi[SUCC(i)]) + edgecost(pi[j], pi[SUCC(j)]) -
               (edgecost(pi[i], pi[j]) + edgecost(pi[SUCC(i)], pi[SUCC(j)]));
	/*
	  printf("Evaluating i %d j %d delta %g best %g tau %g thres %g\n", i, j, delta, best, Tau, best/2);
	  int h;
	  scanf("%d", &h);
	  printf("i %d si %d j %d sj%d\n", i, SUCC(i), j, SUCC(j));
	*/
	if ( delta > best + EPS ) {
	  best = delta;
	  besti = i;
	  bestj = j;
	  //	  printf("found move %d (%d %d) %d (%d %d) val %g\n", i, pi[i], pi[SUCC(i)], j, pi[j], pi[SUCC(j)], delta);
	  if ( !FINDBEST && delta > EPS )
	    goto exit2SFloop;
	}
      }
  }

 exit2SFloop:
  
  sel[0] = min2(besti, bestj);
  sel[1] = max2(besti, bestj);
  freeHeapMem( &HEAP );
  *Delta = best;
  *nmoves = evaluations;
  /*
  printf("Evaluated %d pairs to best %g\n", *nmoves, best );
  */
}


void find2OPTMoveSFHeapNoDup(int* pi, int* sel, double* Delta, int* nmoves, int strongForm ) {
  /* implementazione dell'algoritmo A_g del paper, JOC subm

     input: pi (il tour) strongForm (0=normale 1=versione con un tau più forte)
     output: sel (la mossa) Delta (il suo valore)  nmoves (quante mosse valutate)

     questa versione utilizza in vettore neverExp[] per evitare che delle mosse vengano valutate 2 volte (vedi paper)
   */
  
  THEAP HEAP;
  TPAIR* heap;
  int i, j, k, besti = -1, bestj = -1, neSize;
  double delta, best;
  double Tau;
  int evaluations = 0; 


  //  double tbefore = timer_elapsed_seconds();

  evaluations = 0;
  
  if ( strongForm )
    initialize2OPTHeap( &HEAP, pi, tauStrong );
  else
    initialize2OPTHeap( &HEAP, pi, tauBasic );

  for ( i = 0; i < n; i++ )
    neverExp[i] = i;
  neSize = n;  // size of neverExp[]
  best = -INFINITE;
  heap = HEAP.heap;
  while ( heapsize( heap ) > 0  ) {

    getMax2OPTHeap( &Tau, &i, heap );
    //    printf("popped (%d %d) cost %f thre %f\n", pi[i], pi[SUCC(i)], Tau, best/2 );
    if ( Tau <= best / 2 + EPS ) 
      break; 

    //    printf("expand popped (%d %d) cost %f thre %f\n", pi[i], pi[SUCC(i)], Tau, best/2 );
    //    printf("size %d\n", neSize);
    k = 0;
    while ( k < neSize ) {
      j = neverExp[k];
      //      printf("in loop k %d nesize %d\n", k, neSize);
      //      printf("j is %d\n", j);
      if ( j == i ) {
	neverExp[k] = neverExp[neSize-1];
	neSize--;
      }
      else {
	k++;
	if ( j < i - 1 || j > i + 1 ) {
	  evaluations++;
	  delta = edgecost(pi[i], pi[SUCC(i)]) + edgecost(pi[j], pi[SUCC(j)]) -
	    (edgecost(pi[i], pi[j]) + edgecost(pi[SUCC(i)], pi[SUCC(j)]));
	  /*
	  printf("Evaluating i %d j %d delta %g best %g tau %g thres %g\n", i, j, delta, best, Tau, best/2);
	  int h;
	  scanf("%d", &h);
	  // printf("i %d si %d j %d sj%d\n", i, SUCC(i), j, SUCC(j));
	  */
	  if ( delta > best + EPS ) {
	    best = delta;
	    besti = i;
	    bestj = j;
	    //	  printf("found move %d (%d %d) %d (%d %d) val %g\n", i, pi[i], pi[SUCC(i)], j, pi[j], pi[SUCC(j)], delta);
	    if ( !FINDBEST && delta > EPS )
	      goto exit2SFloop;
	  }
	}
      }
    }
  }
  
 exit2SFloop:
  
  sel[0] = min2(besti, bestj);
  sel[1] = max2(besti, bestj);
  freeHeapMem( &HEAP );
  *Delta = best;
  *nmoves = evaluations;
  /*
  printf("Evaluated %d pairs to best %g\n", *nmoves, best );
  */
}
  
void allocateSorted() { // allocate mem for array of sorted tour edges, used in a version of the 2opt search

  printf("allocate SORTED[.]\n");

  SORTED = (TPAIR*) calloc( 1 + n, sizeof(TPAIR));
  assert(SORTED);
  setRealHeap(SORTED, false);
  setHeapSize(SORTED, n);
  sortedExists = false; /* array is not sorted yet */
}

void disposeSorted() { // free above mem

  free( SORTED );
}

void find2OPTMoveSFSorted(int* pi, int* sel, double* Delta, int* nmoves, int strongForm ) {
  /* implementazione dell'algoritmo A_g nella versione con gli archi del tour sempre ordinati

     Corrisponde all'algoritmo A_L (LEXI) nella versione JOC resubm
     
     input: pi (il tour) strongForm (0=normale 1=versione con un tau più forte)
     output: sel (la mossa) Delta (il suo valore)  nmoves (quante mosse valutate)

   */

  THEAP HEAP;
  int i, j, f, l, besti = -1, bestj = -1;
  double delta, best;
  int evaluations = 0; 

  //  double tbefore = timer_elapsed_seconds();

  evaluations = 0;
 
  if ( strongForm )
    assert( false ); // per ora non implementata la forma forte
  else {
    if ( !sortedExists ) { // per ora fa il sorting sempre a costo n logn. poi si può fare update di sorted a costo (n)
      printf("sorting missinng: we do sort...\n");
      initialize2OPTHeap( &HEAP, pi, tauBasic );
      //      showHeap( HEAP.heap );
      for ( i = 1; i <= n; i++ ) 
	getMaxHeap( SORTED + i, HEAP.heap );
      sortedExists = true; 
      freeHeapMem( &HEAP );
    }
    else {
      printf("sorting not needed...\n");
    }
  }

  /*
    printf("sorted tour edges:\n");
    showHeap( SORTED );
  */
  
  best = -INFINITE;
  f = 1;
  while ( f < n &&
	  edgecost(pi[SORTED[f].x], pi[SUCC(SORTED[f].x)]) +
	  edgecost(pi[SORTED[f+1].x], pi[SUCC(SORTED[f+1].x)])
	  > best + EPS ) {
    l = f+1;
    do {
      evaluations++;
      i = SORTED[f].x;
      j = SORTED[l].x;
      /*
	printf("in loop f %d (%g) l %d (%g) cost %g  best %g\n", f, edgecost(pi[i], pi[SUCC(i)]),
	                                                       l, edgecost(pi[j], pi[SUCC(j)]),
	     edgecost(pi[i], pi[SUCC(i)]) + edgecost(pi[j], pi[SUCC(j)]) -
	     (edgecost(pi[i], pi[j]) + edgecost(pi[SUCC(i)], pi[SUCC(j)])), best);
        int t; printf("..."); scanf("%d", &t);

      */
      delta = edgecost(pi[i], pi[SUCC(i)]) + edgecost(pi[j], pi[SUCC(j)]) -
	(edgecost(pi[i], pi[j]) + edgecost(pi[SUCC(i)], pi[SUCC(j)]));
      if ( delta > best + EPS ) {
	best = delta;
	besti = i;
	bestj = j;
	// printf("found move %d (%d %d) %d (%d %d) val %g\n", i, pi[i], pi[SUCC(i)], j, pi[j], pi[SUCC(j)], delta);
	if ( !FINDBEST && delta > EPS )
	  goto exit2SFSORTloop;
      }

      l++;
    } while ((l <= n) &&
	     (edgecost(pi[SORTED[f].x], pi[SUCC(SORTED[f].x)]) +
	      edgecost(pi[SORTED[l].x], pi[SUCC(SORTED[l].x)])
	      > best + EPS ));
    f++;
  }

 exit2SFSORTloop:

  /*
    printf("exit loop\n");
    printf("..."); scanf("%d", &i);
  */
  
  sel[0] = min2(besti, bestj);
  sel[1] = max2(besti, bestj);
  *Delta = best;
  *nmoves = evaluations;

  /*
  printf("Evaluated %d pairs to best %g\n", *nmoves, best );
  */
}

void creaPairHeap( THEAP* H, int* pi ) {

  int i, ii, jj;

  H->heap = (TPAIR*) calloc( 1 + n, sizeof(TPAIR));
  assert(H->heap);

  for ( i = 1; i <= n - 1; i++ ) {
    H->heap[i].x = i;
    H->heap[i].y = i + 1;
    ii = SORTED[i].x;
    jj = SORTED[i + 1].x;
    H->heap[i].value = edgecost(pi[ii], pi[SUCC(ii)]) + edgecost(pi[jj], pi[SUCC(jj)]);
    //    setelementSingleHeap( H->heap, i, i+1, c[i] + c[i+1], -1, -1, ++cnt);
  }

  setHeapSize( H->heap, n-1 ); // uses field H[0].x to store size

  // makes it into a heap
  for ( i = (n-1) / 2; i >= 1; i-- ) 
    heapify( H->heap,  i );
  
  setRealHeap( H->heap, true );

  //  showHeap( H->heap );
}

void find2OPTMoveSFPairs(int* pi, int* sel, double* Delta, int* nmoves, int strongForm, int heur, double delta_n ) {
  /* implementazione dell'algoritmo A_g nella versione con gli archi del tour ordinati
     e in cui vengono scelte le coppie (i,j) in ordine decrescente di
     c(i,i+1) + c(j,j+1)

     
     input: pi (il tour) strongForm (0=normale 1=versione con un tau più forte)
            heur (T/F) è T se si usa l'euristica, ossia un confronto con delta_n fissato
     output: sel (la mossa) Delta (il suo valore)  nmoves (quante mosse valutate)

   */

  THEAP EDGEHEAP, MOVEHEAP;
  TPAIR* heap;
  // int PR=0;
  int i, j, ii, jj, besti = -1, bestj = -1;
  double mval, best, thre;
  int evaluations = 0; 

  //  double tbefore = timer_elapsed_seconds();

  evaluations = 0;

  if ( strongForm )
    assert( false ); // per ora non implementata la forma forte
  else {
    if ( !sortedExists ) {
      initialize2OPTHeap( &EDGEHEAP, pi, tauBasic );
      //  showHeap( EDGEHEAP.heap );
      for ( i = 1; i <= n; i++ ) 
	getMaxHeap( SORTED + i, EDGEHEAP.heap );

      sortedExists = true; 
      freeHeapMem( &EDGEHEAP );
    }
  }

  /*
    printf("sorted tour edges:\n");
    showHeap( SORTED );
  */

  //  int PR=0;
  creaPairHeap( &MOVEHEAP, pi );
  heap = MOVEHEAP.heap;
  best = -INFINITE;
  if ( heur )
    thre = delta_n;
  else
    thre = best;
  /*  printf("thre : ");
  scanf("%lg", &thre);
  */
  while ( heap[1].value > thre + EPS ) {

    i = heap[1].x;
    j = heap[1].y;
    ii = SORTED[i].x;
    jj = SORTED[j].x;
    //    printf("** pair %d (%d, %d) in sorted (%d, %d) val %g  %g  thre %g\n", PR++, ii, jj, i, j, edgecost(pi[ii], pi[SUCC(ii)]) + edgecost(pi[jj], pi[SUCC(jj)]), heap[1].value, best );

    evaluations++;
    /*
      printf("in loop f %d (%g) l %d (%g) cost %g  best %g\n", f, edgecost(pi[i], pi[SUCC(i)]),
      l, edgecost(pi[j], pi[SUCC(j)]),
      edgecost(pi[i], pi[SUCC(i)]) + edgecost(pi[j], pi[SUCC(j)]) -
      (edgecost(pi[i], pi[j]) + edgecost(pi[SUCC(i)], pi[SUCC(j)])), best);
      int t; printf("..."); scanf("%d", &t);
      
    */

    mval = edgecost(pi[ii], pi[SUCC(ii)]) + edgecost(pi[jj], pi[SUCC(jj)]) - (edgecost(pi[ii], pi[jj]) + edgecost(pi[SUCC(ii)], pi[SUCC(jj)]));
    if ( mval > best + EPS ) {
      best = mval;
      if ( !heur )
	thre = best;
      besti = ii;
      bestj = jj;
      //  printf("found move %d (%d %d) %d (%d %d) val %g\n", ii, pi[ii], pi[SUCC(ii)], jj, pi[jj], pi[SUCC(jj)], mval);
      if ( !FINDBEST && mval > EPS )
	goto exit2SFSORTPAIRloop;
    }

    heap[1].y = heap[1].y + 1;
    if ( heap[1].y == n + 1 )
      heap[1].value = -INFINITE;
    else {
      jj = SORTED[j+1].x;
      heap[1].value = edgecost(pi[ii], pi[SUCC(ii)]) + edgecost(pi[jj], pi[SUCC(jj)]);
    }
    
    heapify( heap, 1 );
    // showHeap( MOVEHEAP.heap );
    
  }

exit2SFSORTPAIRloop:

  /*
    printf("exit loop\n");
    printf("..."); scanf("%d", &i);
  */
  
  sel[0] = min2(besti, bestj);
  sel[1] = max2(besti, bestj);
  *Delta = best;
  *nmoves = evaluations;
  freeHeapMem( &MOVEHEAP );

  /*
  printf("Evaluated %d pairs to best %g\n", *nmoves, best );
  */
}

double getDelta_n() {

  // returns a value for delta_n to be used on the current instance, in the
  // heuristic version of A_g, as of JOC resubmission

  if ( RANDOM_GRAPH_TYPE == UNIFORM ) {
    return 2 - 3.8 / sqrt(n);
  }

  else if ( RANDOM_GRAPH_TYPE == EUCLIDEAN ) {
    return 2 * sqrt(2) - 5 * pow( n, -0.25 );
  }
  else {
    printf("unknown random graph type %d\n", RANDOM_GRAPH_TYPE);
    exit(1);
  }
  
  return 0;
}

void find2optMove(int best, int searchType, int* currentTour, OPTMOVE* move) {

  int sel[4];
  double Delta, tbefore;
  int nmoves;

  tbefore = timer_elapsed_seconds();   
  if ( searchType == SEARCH_2OPT_BF )
    find2OPTMoveBF( currentTour, sel, &Delta, &nmoves );
  else if ( searchType == SEARCH_2OPT_SF )  // greedy A_g basic
    find2OPTMoveSFHeap( currentTour, sel, &Delta, &nmoves, false );
  else if ( searchType == SEARCH_2OPT_SF_NODUP )  // greedy A_g w/o duplicate evaluations
    find2OPTMoveSFHeapNoDup( currentTour, sel, &Delta, &nmoves, false );
  else if ( searchType == SEARCH_2OPT_SF_LEXI )  // greedy based on sorted edges 
    find2OPTMoveSFSorted(currentTour, sel, &Delta, &nmoves, false );
  else if ( searchType == SEARCH_2OPT_SF_SORTPAIR )  // greedy based on sorted pairs of edges 
    find2OPTMoveSFPairs(currentTour, sel, &Delta, &nmoves, false, false, 0 );
  else if ( searchType == SEARCH_2OPT_H_SORTPAIR )  // heur based on sorted pairs of edges w/delta_n
    find2OPTMoveSFPairs(currentTour, sel, &Delta, &nmoves, false, true, getDelta_n() );

  else {
    printf("ERR : 2OPT algorithm %d\n", searchType );
    assert( false );   // other searches to implement, of algorithm for 2OPT unspecified, i.e., = -1
  }
  
  move->seconds2Evaluate = timer_elapsed_seconds() - tbefore;
  
  move->delta = Delta;
  move->found = (Delta > EPS);
  move->selection[0] = sel[0];
  move->selection[1] = sel[1];
  move->selection[2] = -1;
  move->selection[3] = -1;
  move->scheme = -1;
  move->orbitNo = -1;
  move->phi = -1;
  move->maxHeapSize = 0;
  move->heapExtractions = 0;   
  move->numberOfGlobalOptima = -1;
  move->effectiveExtractions = 0;
  move->moveEvaluations = nmoves; 
  move->secondsToBuildHeap = 0;
  //  printf("move evals %ld\n", move->moveEvaluations);
  //  int i;
  // scanf("%d", &i);
}

