
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>

#include "include/config.h"
#include "include/scheduler/ring_sched.h"

void ring_sched_init(int blockId, int numBlocks, void* ringSched)
{
	((ring_sched_t*) ringSched)->stealRange = 0;
	((ring_sched_t*) ringSched)->victimId = blockId;
}

void ring_sched_terminate(void* ringSched)
{
}

ring_sched_t* ring_sched_new(int blockId, int numBlocks)
{
	ring_sched_t* ringSched = (ring_sched_t*) malloc(sizeof(ring_sched_t));
	ring_sched_init(blockId, numBlocks, ringSched);
	return ringSched;
}

void ring_sched_delete(ring_sched_t* ringSched)
{
	ring_sched_terminate(ringSched);
	free(ringSched);
}

void ring_sched_reset(int blockId, void* ringSched)
{
	((ring_sched_t*) ringSched)->stealRange = 0;
	((ring_sched_t*) ringSched)->victimId = blockId;
}

void ring_sched_print(ring_sched_t* ringSched)
{
	PRINTF("Victim-id: %d, steal-range: %d\n",
			ringSched->victimId, ringSched->stealRange);
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
