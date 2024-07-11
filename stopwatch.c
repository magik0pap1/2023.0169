#include "stopwatch.h"

// comments in stopwatch.h

static struct rusage timer_started; // used by timer_start() and timer_elapsed_seconds()

static const double SECONDS_PER_USEC = 0.000001;

void timer_start()
{	getrusage(RUSAGE_SELF, &timer_started);
}

static double elapsed_seconds(struct timeval* to, struct timeval* from)
{
    return (to->tv_sec + SECONDS_PER_USEC * to->tv_usec
            - (from->tv_sec + SECONDS_PER_USEC * from->tv_usec));
}

double timer_elapsed_seconds()
{	struct rusage now;
	getrusage(RUSAGE_SELF, &now);
	return (elapsed_seconds(&now.ru_utime, &timer_started.ru_utime));     // user-mode time
            //+ elapsed_seconds(&now.ru_stime, &timer_started.ru_stime)); // system-mode time
}

void timer_start_nestable(struct rusage* ts)
{	getrusage(RUSAGE_SELF, ts);
}

double timer_elapsed_seconds_nestable(struct rusage* ts)
{	struct rusage now;
	getrusage(RUSAGE_SELF, &now);
	return (elapsed_seconds(&now.ru_utime, &ts->ru_utime));    // user-mode time
            //+ elapsed_seconds(&now.ru_stime, &ts->ru_stime)); // system-mode time
}
