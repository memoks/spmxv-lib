
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <string.h>

#include "include/config.h"
#include "include/timer/custom_timer.h"

#ifdef ARCH_GENERIC

struct timeval start;
void startTimer(void)
{
	gettimeofday(&start, NULL);
}

void stopTimer(double* timeElapsed)
{
	struct timeval end;
	gettimeofday(&end, NULL);

	// convert to seconds
	*timeElapsed = ((double) (end.tv_sec - start.tv_sec)) +
			(end.tv_usec - start.tv_usec) * 0.00000001;
}

#else

#include <mkl.h>

/*
double timeStart;
void startTimer(void)
{
	// do a dummy dsecnd() call to improve accuracy of timing
	dsecnd();
	dSecondTimeStart = dsecnd();
}

void stopTimer(void)
{
	double timeEnd = dsecnd();
	PRINTF("done in %3.3f seconds.\n", timeEnd - timeStart);
}
*/

double dSecondTimeStart;
void startTimer(void)
{
	// do a dummy dsecnd() call to improve accuracy of timing
	dsecnd();
	dSecondTimeStart = dsecnd();
}

void stopTimer(double* timeElapsed)
{
	double dSecondTimeEnd = dsecnd();
	*timeElapsed = dSecondTimeEnd - dSecondTimeStart;
}

#endif

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
