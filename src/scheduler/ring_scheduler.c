
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <omp.h>

#include "include/config.h"
#include "include/timer/custom_timer.h"
#include "include/scheduler/ring_sched.h"
#include "include/scheduler/block_info.h"

#include "include/scheduler/ring_scheduler.h"

// Some helper functions
// ------------------------------------------------------------------------------------------------------

static void __ring_scheduler_calculatePrevNext(block_info_t* threadInfo, int numBlocks,
		int* prevId_out, int* nextId_out);

// ------------------------------------------------------------------------------------------------------

void ring_scheduler_init(block_info_t** blockInfos, int numBlocks, int jobBatchCount)
{
	#pragma omp parallel num_threads(numBlocks)
	{
		int blockId = omp_get_thread_num();
		ring_sched_init(blockId, numBlocks, &blockInfos[blockId]->sched.ringSched);
	}
}

void ring_scheduler_terminate(block_info_t** blockInfos, int numBlocks)
{
	#pragma omp parallel num_threads(numBlocks)
	{
		int blockId = omp_get_thread_num();
		ring_sched_terminate(&blockInfos[blockId]->sched.ringSched);
	}
}

int ring_scheduler_findVictim(block_info_t* myBlock, block_info_t** blockInfos, int numBlocks)
{
	int prevId;
	int nextId;

	__ring_scheduler_calculatePrevNext(myBlock, numBlocks, &prevId, &nextId);

	int carryOn = TRUE;
	int victimId = myBlock->sched.ringSched.victimId;
	while(carryOn)
	{
		if(victimId == nextId)
		{
#if PROGRAM_MODE >= TRACE_MODE
			PRINTF("findVictim: P%d increasing steal_range: %d to %d.\n",
					myBlock->id, myBlock->sched.ringSched.stealRange, myBlock->sched.ringSched.stealRange + 1);
#endif
			++myBlock->sched.ringSched.stealRange;
			if((2 * myBlock->sched.ringSched.stealRange) > numBlocks)
				carryOn = FALSE;

			__ring_scheduler_calculatePrevNext(myBlock, numBlocks, &prevId, &nextId);
			continue;
		}
		else if(victimId == prevId)
		{
			victimId = nextId;
#if PROGRAM_MODE >= TRACE_MODE
			PRINTF("findVictim: P%d (steal_range: %d) trying for victim P%d(%d)...\n",
					myBlock->id, myBlock->sched.ringSched.stealRange, victimId, blockInfos[victimId]->queue.batchCount);
#endif
		}
		else
		{
			victimId = prevId;
#if PROGRAM_MODE >= TRACE_MODE
			PRINTF("findVictim: P%d (steal_range: %d) trying for victim P%d(%d)...\n",
					myBlock->id, myBlock->sched.ringSched.stealRange, victimId, blockInfos[victimId]->queue.batchCount);
#endif
		}

		if(blockInfos[victimId]->queue.executableBatchCount > 0)
		{
			myBlock->sched.ringSched.victimId = victimId;
#if PROGRAM_MODE >= TRACE_MODE
			PRINTF("findVictim: P%d (steal_range: %d) chosen P%d(%d) as a victim.\n",
					myBlock->id, myBlock->sched.ringSched.stealRange, victimId, blockInfos[victimId]->queue.batchCount);
#endif

			return TRUE;
		}
	}

	return FALSE;
}

// Some helper functions
// --------------------------------------------------------------------------------------------------------------

static void __ring_scheduler_calculatePrevNext(
		block_info_t* threadInfo, int numBlocks,
		int* prevId_out, int* nextId_out)
{
	*prevId_out = (threadInfo->id + numBlocks - threadInfo->sched.ringSched.stealRange) % numBlocks;
	*nextId_out = (threadInfo->id + threadInfo->sched.ringSched.stealRange) % numBlocks;
}


#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
