#include <sys/resource.h>
/*
   Questa funzione va invocata nel punto in cui si vuole "far partire il cronometro":
   il tempo (di utilizzo della CPU) verrà misurato a partire da qui.
   (va invocata, ad esempio, all'inizio del programma o all'inizio della sezione di
   programma di cui si vuole misurare il tempo d'esecuzione).
   Resettando il cronometro (cioè invocando di nuovo questa funzione) si possono
   misurare diversi intervalli di tempo, per diverse sezioni del programma.
   Esempio:
   
   timer_start();
   ...
   double tempo_dal_primo_reset = timer_elapsed_seconds();
   ...
   double tempo_dal_primo_reset_2 = timer_elapsed_seconds();
   ...
   timer_start();
   ...
   double tempo_dal_secondo_reset = timer_elapsed_seconds();
   ...
   
   Le varie invocazioni possono essere nel main e/o in funzioni diverse, perché l'istante
   in cui è avvenuto il reset più recente viene memorizzato in una variabile static esterna,
   cioè condivisa da tutte le funzioni.
*/
void timer_start();
/*
   Questa funzione va invocata nel punto in cui si vuole "leggere il valore del cronometro"
   e restituisce il numero di secondi di utilizzo della CPU (in formato double).
   Questa funzione NON resetta il cronometro, quindi la si può invocare più volte e
   restituisce di volta in volta il tempo trascorso dall'ultimo reset del cronometro.
*/
double timer_elapsed_seconds();
/*
   Le due funzioni seguenti consentono un utilizzo più generale,
   con cronometri annidati. Il chiamante deve creare una struttura di tipo rusage
   e fornirla come argomento a entrambe le funzioni, in questo modo non c'è
   interferenza con altri cronometri né con il cronometro principale gestito
   dalle due funzioni precedenti.
*/
void timer_start_nestable(struct rusage* rs);
double timer_elapsed_seconds_nestable(struct rusage* rs);


