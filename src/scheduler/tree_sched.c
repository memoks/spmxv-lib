
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>

#include "include/scheduler/tree_sched.h"

void tree_sched_init(int blockId, int numBlocks, void* treeSched)
{
	int victimCount = numBlocks - 1;

	((tree_sched_t*) treeSched)->stealCount = 0;
	((tree_sched_t*) treeSched)->victimIndex = 0;
	((tree_sched_t*) treeSched)->lastVictim = blockId;
	((tree_sched_t*) treeSched)->victimCount = victimCount;

	((tree_sched_t*) treeSched)->victims = (int*) malloc(sizeof(int) * victimCount);
	int i;
	for(i = 0; i < victimCount; ++i)
		((tree_sched_t*) treeSched)->victims[i] = -1;
}

void tree_sched_terminate(void* treeSched)
{
	if(((tree_sched_t*) treeSched)->victims != NULL)
		free(((tree_sched_t*)treeSched)->victims);
}

tree_sched_t* tree_sched_new(void)
{
	tree_sched_t* treeSched = (tree_sched_t*) malloc(sizeof(tree_sched_t));
	treeSched->stealCount = 0;
	treeSched->victimIndex = 0;
	treeSched->lastVictim = 0;
	treeSched->victims = NULL;
	treeSched->victimCount = 0;
	return treeSched;
}

void tree_sched_delete(tree_sched_t* treeSched)
{
	tree_sched_terminate(treeSched);
	free(treeSched);
}

void tree_sched_reset(int blockId, void* treeSched)
{
	((tree_sched_t*) treeSched)->stealCount = 0;
	((tree_sched_t*) treeSched)->victimIndex = 0;
	((tree_sched_t*) treeSched)->lastVictim = blockId;
}

void tree_sched_print(tree_sched_t* treeSched)
{
	if(treeSched->victims == NULL)
		return;

	PRINTF("Victim-Id: %d, victimIndex: %d, victimCount: %d, steal-count: %d\n",
			treeSched->lastVictim, treeSched->victimIndex,
			treeSched->victimCount, treeSched->stealCount);
	PRINTF("Victim Arr-%d:", treeSched->victimCount);

	int i;
	for(i = 0; i < treeSched->victimCount; ++i)
		PRINTF(" %d", treeSched->victims[i]);
	PRINTF("\n");
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
