#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

/*

  legge un file di tipo TSPLIB e lo trasforma in un formato con la
  lista degli archi e relativi costi, i.e.,

  <n>
  <i1> <j1> <cost1>
  ...
  <i_m> <j_m> <cost_m>


  il tipo TSPLIB può essere EUC_2D, EXPLICIT (LOWER O UPPER TRIANGULAR)
*/

#define true          1
#define false         0
#define MAXN      34000

float xcoo[MAXN], ycoo[MAXN];

#define MAXTOKS           1000     /* max tokens per line */
#define MAXLINELEN       10000

typedef char* tokenarray[MAXTOKS];
typedef char  line[MAXLINELEN];

FILE* inf, * outf;
tokenarray token;
int ntokens;
int n;


int _getline(FILE* f, char* s)
/* reads a line from a file. Returns TRUE if eof */
{
  char c;
  int i, eof;
  
  i = 0;
  while (1) {
    c = getc(f);
    if ( c == '\r' ) continue;
    eof = ( c == EOF );
    if ( eof || c =='\n')
      c = '\0';
    if ( i == MAXLINELEN ) {
      printf("Error: line too long in getline\n");
      exit(1);
    }
    s[i++] = c;
    if ( c == '\0')
      break;
  }

  // deletes trailing blanks
  i--;
  assert( s[i] == '\0');
  while (i>0 && s[i-1] == ' ')
    s[--i] = '\0';
  
  printf(">>%s<<\n", s);
  return eof;
}

int  line_tokens(FILE* inf, tokenarray token, int* ntokens)
/* reads a line from the file inf and decomposes it into *ntokens tokens. */
{

   char line[MAXLINELEN];
   char* strptr;
   int eof;

   eof = _getline(inf, line);
   strptr = line;
   *ntokens = 0;
   while ( (token[*ntokens] = strtok(strptr, " \t\n")) != NULL ) {
       printf("Tok %d is %s\n", *ntokens, token[*ntokens]); 
      strptr = NULL;
      *ntokens += 1;
      if ( *ntokens == MAXTOKS ) {
        printf("Error: too many tokens in line_tokens\n");
        exit(1);
      }
   }
   return eof;
}

int nint( double x ) {
  
  return (int) (x+0.5);
}

void usage() {

  printf("Usage: convTSPLIB <infile> <outfile>\n");
  exit(1);
}

void readLOWERDIAG() {
  printf("read LOWER_DIAG to do\n");
  assert(false);
}

void readUPPERDIAG() {

  int  i, j, cij, esito;

  while ( true ) {
    
    esito = line_tokens(inf, token, &ntokens); 
    printf("dd %d\n", ntokens); fflush(stdout);
    
    if ( esito )
      break;

    if ( strcmp( token[0], "EDGE_WEIGHT_SECTION" ) == 0 )
      break;
  }

  printf("saving upper row format\n");
  fprintf( outf, "%d\n", n );
  for ( i = 0; i < n; i++ )
    for ( j = i; j < n; j++ ) {
      fscanf( inf, "%d", &cij );
      if ( j > i )
	fprintf(outf, "%d %d %d\n", i, j, cij );
    }
}

void readFULLMATRIX() {
  printf("read FULL_MATRIX to do\n");
  assert(false);
}

void readUPPERROW() {

  int  i, j, cij, esito;

  while ( true ) {
    
    esito = line_tokens(inf, token, &ntokens); 
    
    if ( esito )
      break;

    if ( strcmp( token[0], "EDGE_WEIGHT_SECTION" ) == 0 )
      break;
  }

  printf("saving upper row format\n");
  fprintf( outf, "%d\n", n );
  for ( i = 0; i < n - 1; i++ )
    for ( j = i + 1; j < n; j++ ) {
      fscanf( inf, "%d", &cij );
      fprintf(outf, "%d %d %d\n", i, j, cij );
    }
}	




void readEXPLICIT() {
  
  int i, esito;
  char sss[100];

  while ( true ) {

    esito = line_tokens(inf, token, &ntokens); 
   
    if ( esito )
      break;

    strcpy( sss, token[0] );
    if ( strcmp( sss, "EDGE_WEIGHT_FORMAT:" ) == 0 ) {
      if ( strcmp( token[1], "LOWER_DIAG_ROW" ) == 0 ) {
	readLOWERDIAG();
	break;
      }
      else if ( strcmp( token[1], "UPPER_DIAG_ROW" ) == 0 ) {
	readUPPERDIAG();
	break;
      }
      else if ( strcmp( token[1], "UPPER_ROW" ) == 0 ) {
	readUPPERROW();
	break;
      }
      else if ( strcmp( token[1], "FULL_MATRIX" ) == 0 ) {
	readFULLMATRIX();
	break;
      }
      else
	assert(false);
    }
  
    else if ( strcmp( sss, "EDGE_WEIGHT_FORMAT" ) == 0 ) {
      if ( strcmp( token[2], "LOWER_DIAG_ROW" ) == 0 ) {
	readLOWERDIAG();
	break;
      }
      else if ( strcmp( token[2], "LOWER_DIAG_ROW" ) == 0 ) {
	readUPPERDIAG();
	break;
      }
      else if ( strcmp( token[2], "UPPER_ROW" ) == 0 ) {
	readUPPERROW();
	break;
      }
      else if ( strcmp( token[2], "FULL_MATRIX" ) == 0 ) {
	readFULLMATRIX();
	break;
      }
      else
	assert(false);
    }
  }
}

void readEUC_2D() {

  // prosegue con la lettura di un file di tipo EUC_2D

  double xd, yd;
  int i, j, esito;
  char sss[100];
  
  while ( true ) {

    esito = line_tokens(inf, token, &ntokens); 
   
    if ( esito )
      break;

    /*
    printf("dd %d\n", ntokens); fflush(stdout);
    printf(">>");
    for ( i = 0; i < ntokens; i++ )
      printf("(%d [%s])", i, token[i]);
    printf("\n");
    */
    
    strcpy( sss, token[0] );
    if ( strcmp( sss, "DIMENSION" ) == 0 ) {
      n = atoi( token[2] );
    }
    else if ( strcmp( sss, "DIMENSION:" ) == 0 ) {
      printf("comp [%s] with [%s] result %d\n", sss, "DIMENSION",  strcmp( sss, "DIMENSION:" ));
      n = atoi( token[1] );
    }
    else if ( strcmp( sss, "NODE_COORD_SECTION" ) == 0 ) {
      for ( i = 0; i < n; i++ ) {
	fscanf(inf, "%d", &j);
	assert( j == i + 1 );
	printf("j is %d\n", j);
	fflush(stdout);
	fscanf(inf, "%g %g\n", &(xcoo[i]), &(ycoo[i]));
      }
    }
  }

  printf("ecco i nodi\n");
  for ( i = 0; i < n; i++ )
    printf("node %d  x %g y %g\n", i, xcoo[i], ycoo[i] );

  
  fprintf( outf, "%d\n", n );
  for ( i = 0; i < n - 1; i++ )
    for ( j = i + 1; j < n; j++ ) {
      xd = xcoo[i] - xcoo[j];
      yd = ycoo[i] - ycoo[j];
      fprintf( outf, "%d %d %d\n", i, j, nint( sqrt( xd*xd + yd*yd) ) );
    }
}


int main( int argc, char ** argv ) {

  int i, esito;
  char sss[100];
  
  if ( argc != 3 )
    usage();

  inf = fopen( argv[1], "r" );
  assert( inf );

  outf = fopen( argv[2], "w" );


  while ( true ) {

    esito = line_tokens(inf, token, &ntokens); 
   
    if ( esito )
      break;
    /*
    printf("dd %d\n", ntokens); fflush(stdout);
    printf(">>");
    for ( i = 0; i < ntokens; i++ )
      printf("(%d [%s])", i, token[i]);
    printf("\n");
    */
    strcpy( sss, token[0] );
    if ( strcmp( sss, "EDGE_WEIGHT_TYPE:" ) == 0 ) {
      if ( strcmp( token[1], "EXPLICIT" ) == 0 ) {
	readEXPLICIT();
	break;
      }
      else if ( strcmp( token[1], "EUC_2D" ) == 0 ) {
	readEUC_2D();
	break;
      }
      else
	assert(false);
    }
    else if ( strcmp( sss, "EDGE_WEIGHT_TYPE" ) == 0 ) {
      if ( strcmp( token[2], "EXPLICIT" ) == 0 ) {
	readEXPLICIT();
	break;
      }
      else if ( strcmp( token[2], "EUC_2D" ) == 0 ) {
	readEUC_2D();
	break;
      }
      else
	assert(false);
    }
    else if ( strcmp( sss, "DIMENSION" ) == 0 ) {
      n = atoi( token[2] );
    }
    else if ( strcmp( sss, "DIMENSION:" ) == 0 ) {
      printf("comp [%s] with [%s] result %d\n", sss, "DIMENSION",  strcmp( sss, "DIMENSION:" ));
      n = atoi( token[1] );
    }
    else if ( strcmp( sss, "NODE_COORD_SECTION" ) == 0 ) {
      assert(false);
    }
  }

  fclose( inf);
  fclose(outf);
}
