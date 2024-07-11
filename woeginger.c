#include <assert.h>
#include "tspmain.h"
#include "proc4opt.h"

/*
  Procedures for the dynamic program by Woeginger et al. for 4OPT
*/


//#define PRINTON


double costWoeg( int* pi, int nmove, int type, int tail, int pt1, int pt2 ) {

  switch ( nmove ) {

  case 1: // +1 -2 -3 -4
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1]] - cost[pi[tail + 1]][pi[pt2]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1+1]] - cost[pi[tail + 1]][pi[pt2+1]]; 

  case 2: // +1 -2 +3 -4
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1]] - cost[pi[tail + 1]][pi[pt1+1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2]] - cost[pi[tail + 1]][pi[pt2+1]]; 

  case 3: // +1 -2 -4 +3
  case 12: // +1 -3 +4 +2
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1]] - cost[pi[tail + 1]][pi[pt2+1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1+1]] - cost[pi[tail + 1]][pi[pt2]]; 

  case 4: // +1 -2 +4 -3
  case 11: // +1 -3 +4 -2
      if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1]] - cost[pi[tail + 1]][pi[pt2+1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2]] - cost[pi[tail + 1]][pi[pt1+1]]; 

  case 5: // +1 -2 +4 +3
  case 9: // +1 -3 -4 -2
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1]] - cost[pi[tail + 1]][pi[pt2]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2+1]] - cost[pi[tail + 1]][pi[pt1+1]]; 

  case 6: // +1 -3 +2 -4
  case 19: // +1 -4 +2 -3
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2]] - cost[pi[tail + 1]][pi[pt1+1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1]] - cost[pi[tail + 1]][pi[pt2+1]]; 

  case 7: // +1 +3 -2 -4
  case 8: // +1 +3 +2 -4
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2]] - cost[pi[tail + 1]][pi[pt1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1+1]] - cost[pi[tail + 1]][pi[pt2+1]]; 

  case 10: // +1 -3 -4 +2
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1]] - cost[pi[tail + 1]][pi[pt1+1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2+1]] - cost[pi[tail + 1]][pi[pt2]]; 
    
  case 13: // +1 +3 -4 -2
  case 14: // +1 +3 -4 +2
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2+1]] - cost[pi[tail + 1]][pi[pt1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2]] - cost[pi[tail + 1]][pi[pt1+1]]; 
    
  case 15: // +1 -3 +2 +4
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1+1]] - cost[pi[tail + 1]][pi[pt2]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1]] - cost[pi[tail + 1]][pi[pt2+1]]; 

  case 16: // +1 +4 -2 -3
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1+1]] - cost[pi[tail + 1]][pi[pt1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2]] - cost[pi[tail + 1]][pi[pt2+1]]; 

  case 17: // +1 -4 -2 +3
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2+1]] - cost[pi[tail + 1]][pi[pt1+1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1]] - cost[pi[tail + 1]][pi[pt2]]; 

  case 18: // +1 +4 -2 +3
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2]] - cost[pi[tail + 1]][pi[pt1+1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2+1]] - cost[pi[tail + 1]][pi[pt1]]; 

  case 20: // +1 +4 +2 -3
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2]] - cost[pi[tail + 1]][pi[pt1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1+1]] - cost[pi[tail + 1]][pi[pt2+1]]; 

  case 21: // +1 -4 +3 -2
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2]] - cost[pi[tail + 1]][pi[pt2+1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1]] - cost[pi[tail + 1]][pi[pt1+1]]; 

  case 22: // +1 -4 +3 +2
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1+1]] - cost[pi[tail + 1]][pi[pt2+1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1]] - cost[pi[tail + 1]][pi[pt2]]; 

  case 23: // +1 +4 -3 +2
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2+1]] - cost[pi[tail + 1]][pi[pt1+1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2]] - cost[pi[tail + 1]][pi[pt1]]; 
    
  case 24: // +1 +4 +3 -2
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2]] - cost[pi[tail + 1]][pi[pt1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2+1]] - cost[pi[tail + 1]][pi[pt1+1]]; 

  case 25: // +1 +4 +3 +2
    if ( type == 0 ) 
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt1+1]] - cost[pi[tail + 1]][pi[pt1]]; 
    else
      return cost[pi[tail]][pi[tail+1]] - cost[pi[tail]][pi[pt2+1]] - cost[pi[tail + 1]][pi[pt2]]; 

  default:
    assert(0);
  }
}


double dynProg(int *pi, int movetype, int fix0, int fix1, int* pa, int* pb) 
{

  // Dynamic programming to place two missing edges and form a 4-opt move
  // it is assumed that two of the removed edges (a set R1, i.e. (fix0,fix0+1) and (fix1,fix1+1)) 
  // of the move have already been
  // chosen, and two (a set R2) are missing. The missing ones will carry the weight also of the
  // inserted edges in the cost function, and hence are such that these edges have
  // endpoints only to the edges in R1

  // We will determine the starting, a and b, of the edges in R2 (i.e. the edges will
  // be (a,a+1) and (b,b+1). These two endpoints have to be placed in the interval
  // [from,...,to], with a < b

  int a, b;
  double dpval[2][MAXN];
  int from0, to0, from1, to1;

  switch (movetype) {
  case 9 ... 12:  
  case 15 ... 16:  
  case 19 ... 22:
  case 24 ... 25: // fix0=I fix1=J
    from0 = fix1 + 2;   to0 = n-1-(fix0==0);
    from1 = fix1 + 2;   to1 = n-1-(fix0==0);
    break;

  case 2:
  case 3:
  case 7:
  case 13: 
  case 17:  // fix0=I fix1=K
    from0 = fix0 + 2;   to0 = fix1 - 2;
    from1 = fix1 + 2;   to1 = n-1-(fix0==0);
    break;

  case 1:
  case 4 ... 6:
  case 8:
  case 14:
  case 18:   
  case 23: // fix0=I fix1=H
    from0 = fix0 + 2;   to0 = fix1 - 4;
    from1 = fix0 + 4;   to1 = fix1 - 2;
    break;

  default:
    assert(false);
  }

  assert( from0 >= 2 && from1 >= from0 );
  dpval[0][from0 - 1] = dpval[0][from0 - 2] = -INFINITE;
  dpval[1][from1 - 1] = dpval[1][from1 - 2] = -INFINITE;

  dpval[0][from0] = costWoeg( pi, movetype, 0, from0, fix0, fix1 );
  for ( a = from0 + 1; a <= to0; a++ ) {
    dpval[0][a] = MAX( dpval[0][a-1], costWoeg( pi, movetype, 0, a, fix0, fix1 ));
    //    printf("set d[%d %d] to %g\n", 0, a, dpval[0][a] );
  }

  for ( b = from1; b <= to1; b++ ) {
    dpval[1][b] = MAX( dpval[1][b-1], costWoeg( pi, movetype, 1, b, fix0, fix1 ) + dpval[0][MIN(b-2,to0)] );
    //    printf("set d[%d %d] to %g\n", 1, b, dpval[1][b] );
  }

  //#define PPP
#ifdef PPP
  printf("TABELLA WOEG DP [FROM %d to %d] per ii %d  jj %d\n", from0, to1, fix0, fix1);
  for ( a = from0 -2; a <= to1; a++ )
    printf("%11d", a);
  printf("\n\n");
  for ( b = from0 -2; b <= to0; b++ )
    printf("%11d", (int)dpval[0][b]);
  printf("\n");
  for ( b = from0 -2; b <= to1; b++ )
    if ( b < from1 ) 
      printf("%11d", -9999);
    else
      printf("%11d", (int)dpval[1][b]);
  printf("\n...");
  scanf("%d", &a);
#endif

  b = to1;
  while ( dpval[1][b-1] == dpval[1][b] ) b--;
  a = MIN(b-2, to0);
  while ( dpval[0][a-1] == dpval[0][a] ) a--;

  dprintf("OPT a %d b %d \n\n", a, b);
  *pa = a;
  *pb = b;
  return dpval[1][to1];
}


void retrieve4tuple( int x, int y, int u, int z, int* ii, int* jj, int* kk, int* hh ) {
  // given x < y and u < z, sorts them so that ii is the min, jj the second,
  // kk the third and hh the max

  if ( y < u ) {
    *ii = x;    *jj = y;    *kk = u;   *hh = z;
  }
  else if ( y < z && x < u ) {
    *ii = x;    *jj = u;    *kk = y;   *hh = z;
  }
  else if ( y < z && u < x ) {
    *ii = u;    *jj = x;    *kk = y;   *hh = z;
  }
  else if ( x < u ) {
    *ii = x;    *jj = u;    *kk = z;   *hh = y;
  }
  else if ( x < z ) {
    *ii = u;    *jj = x;    *kk = z;   *hh = y;
  }
  else {
    *ii = u;    *jj = z;    *kk = x;   *hh = y;
  }
}


double find4OPTWoeginger( int* pi, int * ii, int* jj, int* kk, int* hh, int movetype ) { 
  /* given permutation pi[], finds an improving (true) 4OPT move, of type <movetype>
     i.e., 4 indices i < j < k < h defining a move, 

- - - - -

  Uses O(n^3) method by woeginger

  */

  int x, y, u, z;
  int minx, maxx, maxy, incr;
  double delta, vv;
  double partout, restofcost;

  dprintf("WOEG %d *****\n", movetype);

  assert( movetype >= 1 && movetype <= 25 );
  switch ( movetype ) {
  case 9 ... 12:  
  case 15 ... 16:  
  case 19 ... 22:
  case 24 ... 25: // fix0=I fix1=J
    minx = 0;    maxx = n - 7;
    incr = 2;    maxy = n - 5;
    break;

  case 2:
  case 3:
  case 7:
  case 13: 
  case 17:  // fix0=I fix1=K
    minx = 0;    maxx = n - 7;
    incr = 4;    maxy = n - 3;
    break;

  case 1:
  case 4 ... 6:
  case 8:
  case 14:
  case 18:   
  case 23: // fix0=I fix1=H
    minx = 0;    maxx = n - 7;
    incr = 6;    maxy = n - 1;
    break;

  default:
    assert(false);
    break;
  }

  delta = -INFINITE;
  for ( x = minx; x <= maxx; x++ )
    for ( y = x + incr;  y <= maxy - (x==minx); y++ ) {
      partout = cost[pi[x]][pi[x+1]] + cost[pi[y]][pi[y+1]];
      restofcost = dynProg(pi, movetype, x, y, &u, &z); 
      vv = partout + restofcost;
      dprintf("Val per x %d y %d part %g rest %g totvv %g\n", x, y, partout, restofcost, vv);
      if ( (FINDBEST && vv > delta + EPS ) || (!FINDBEST && vv > EPS ) ) { 
        retrieve4tuple( x, y, u, z, ii, jj, kk, hh);
	delta = vv;
	dprintf("**** UPDATE i %d j %d k %d h %d for val %g\n", *ii, *jj, *kk, *hh, vv);
	if ( !FINDBEST ) 
	  break;
      }
    }

  dprintf("WOEG (M%d) i %4d  j %4d  k %4d  h %4d  val %g\n", movetype, *ii, *jj, *kk, *hh, delta );
  return delta;
}


double find4OPTWoegingerByOrbit( int* pi, int * sel, int* rot, int numOrbit, OPTMOVE* move )
{

  double delta, vv;

  int s, actuals, norb;
  int i, j, k, h, scheme;
  int revpi[MAXN];   // pi reversed
  int* actualpi[8];

  assert( numOrbit >= 1 && numOrbit <= 7 );
  norb = numOrbit - 1;   // seven orbits, 0..6

  if ( numOrbit == 4 || numOrbit == 5 ) // orbits w/reflection
    for ( i = 0; i <= n; i++ )
      revpi[i] = pi[n-i];

  delta = -INFINITE;

  for ( s = 0; s < ORBIT4[norb].size; s++ ) {
    if ( !ORBIT4[norb].flip[s] ) {
      actualpi[s] = pi;
      scheme = movesInOrbit[numOrbit][s];
    }
    else {
      actualpi[s] = revpi + 0;
      scheme = movesInOrbit[numOrbit][s - ORBIT4[norb].size/2 ];
    }

    actuals = removeReflection( norb, s );
    dprintf("orb %d s %d actual %d\n", numOrbit, s, actuals );

    vv = find4OPTWoeginger( actualpi[s], &i, &j, &k, &h, scheme );

    if ( (FINDBEST && vv > delta + EPS ) || (!FINDBEST && vv > EPS ) ) { 
      delta = vv;
      sel[0] = i;
      sel[1] = j;
      sel[2] = k;
      sel[3] = h;
      *rot = actuals;
      
      if ( !FINDBEST ) 
	break;
    }
  }
  
  if ( ORBIT4[norb].flip[*rot] )
    reflectSelection( sel );

  dprintf("WOEGORB (OR%d G%d) i %4d  j %4d  k %4d  h %4d  val %g\n", numOrbit, *rot, sel[0], sel[1], sel[2], sel[3], delta );

  move->found = (delta > EPS);
  move->delta = delta;
  for ( i = 0; i < 4; i++ )
    move->selection[i] = sel[i];
  move->heapExtractions = 0;
  move->moveEvaluations = UNDEF;
  move->effectiveExtractions = UNDEF;
  move->orbitNo = numOrbit;
  move->phi = *rot;
  move->scheme = movesInOrbit[numOrbit][*rot];
  move->maxHeapSize = 0;

  return delta;
}
