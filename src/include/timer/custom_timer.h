
#ifndef CUSTOMTIMER_H_
#define CUSTOMTIMER_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>

/**
 * Time is measured in seconds.
 */
extern void startTimer(void);
extern void stopTimer(double* timeElapsed_out);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif


#endif /* CUSTOMTIMER_H_ */
