
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>

#include "include/control_unit/cu_options.h"
#include "include/scheduler/block_info.h"


void block_init(block_info_t* blockInfo, int blockId, int numBlocks, int stealTreshold, int stealingScheme)
{
	blockInfo->id = blockId;
	blockInfo->schedType = stealingScheme;
	job_queue_init(&blockInfo->queue, stealTreshold);

	if(stealingScheme == STEALING_SCHEME_RING)
	{
		blockInfo->scheduler_findVictim = ring_scheduler_findVictim;
		blockInfo->sched_init = ring_sched_init;
		blockInfo->sched_reset = ring_sched_reset;
		blockInfo->sched_terminate = ring_sched_terminate;

	}
	else // if(stealingScheme == STEALING_SCHEME_TREE)
	{
		blockInfo->scheduler_findVictim = tree_scheduler_findVictim;
		blockInfo->sched_init = tree_sched_init;
		blockInfo->sched_reset = tree_sched_reset;
		blockInfo->sched_terminate = tree_sched_terminate;
	}

	blockInfo->sched_init(blockId, numBlocks, &blockInfo->sched);
}

void block_deleteNonPtr(block_info_t* blockInfo)
{
	job_queue_deleteNonPtr(&blockInfo->queue);
}

block_info_t* block_new(int blockId, int stealTreshold, int stealingScheme, int numBlocks)
{
	block_info_t* blockInfo = (block_info_t*) malloc(sizeof(block_info_t));
	block_init(blockInfo, blockId, numBlocks, stealTreshold, stealingScheme);

	return blockInfo;
}

block_info_t** block_newMultiple(int numBlocks, int stealTreshold, int stealingScheme)
{
	block_info_t** blockInfos = (block_info_t**) malloc(sizeof(block_info_t*) * numBlocks);

	int i;
	for(i = 0; i < numBlocks; ++i)
		blockInfos[i] = block_new(i, stealTreshold, stealingScheme, numBlocks);

	return blockInfos;
}

void block_delete(block_info_t* blockInfo)
{
	blockInfo->sched_terminate((void*) &blockInfo->sched);
	block_deleteNonPtr(blockInfo);
	free(blockInfo);
}

void block_deleteMultiple(block_info_t** blockInfos, int numBlocks)
{
	int i;
	for(i = 0; i < numBlocks; ++i)
		block_delete(blockInfos[i]);

	free(blockInfos);
}

inline block_info_t* block_getBlockFromRingSched(ring_sched_t* ringSched)
{
	return block_entry((sched_t*) ringSched, block_info_t, sched);
	// return block_entry(ringSched, block_info_t, ringSched);
}

inline block_info_t* block_getBlockFromTreeSched(tree_sched_t* treeSched)
{
	return block_entry((sched_t*) treeSched, block_info_t, sched);
	// return block_entry(treeSched, block_info_t, treeSched);
}

void block_reset(block_info_t* blockInfo)
{
	job_queue_reset(&blockInfo->queue);
	blockInfo->sched_reset(blockInfo->id, (void*) &blockInfo->sched);
}

void block_resetMultiple(block_info_t** blockInfos, int numBlocks)
{
	int i;
	for(i = 0; i < numBlocks; ++i)
		block_reset(blockInfos[i]);
}

void block_reform(block_info_t* blockInfo, lg_t* executedJobBatchList)
{
	block_reset(blockInfo);
	job_queue_emptyQueue(&blockInfo->queue);
	job_queue_addJobBatchListWithoutLocking(&blockInfo->queue, executedJobBatchList);
}

void block_reformMultiple(block_info_t** blockInfos, lg_t** executedJobBatchListPerBlock, int numBlocks)
{
	int i;
	for(i = 0; i < numBlocks; ++i)
	{
		block_reset(blockInfos[i]);
		job_queue_emptyQueue(&blockInfos[i]->queue);
	}

	for(i = 0; i < numBlocks; ++i)
	{
		job_queue_addJobBatchListWithoutLocking(&blockInfos[i]->queue, executedJobBatchListPerBlock[i]);
	}
}

void block_print(block_info_t* blockInfo)
{
	PRINTF("--------------------------- BLOCK-%d ---------------------------\n", blockInfo->id);

	if(blockInfo->schedType == STEALING_SCHEME_RING)
		ring_sched_print(&blockInfo->sched.ringSched);
	else
		tree_sched_print(&blockInfo->sched.treeSched);

	job_queue_printAttributes(&blockInfo->queue);
	// PRINTF("-----------------------------------------------------------------\n");
}

void block_printDetailed(block_info_t* blockInfo, void* spm)
{
	PRINTF("--------------------------- BLOCK-%d ---------------------------\n", blockInfo->id);

	if(blockInfo->schedType == STEALING_SCHEME_RING)
		ring_sched_print(&blockInfo->sched.ringSched);
	else
		tree_sched_print(&blockInfo->sched.treeSched);

	job_queue_printDetailed(&blockInfo->queue, spm);
	// PRINTF("-----------------------------------------------------------------\n");
}

void block_printMultiple(block_info_t** blockInfos, int numBlocks)
{
	int i;
	for(i = 0; i < numBlocks; ++i)
		block_print(blockInfos[i]);
}

void block_printDetailedMultiple(block_info_t** blockInfos, int numBlocks, void* spm)
{
	int i;
	for(i = 0; i < numBlocks; ++i)
		block_printDetailed(blockInfos[i], spm);
}

void block_convertToJDSPartialMultiple(
		spm_cmp_t* spmCsrCounterpart, DECIMAL* permutation,
		block_info_t** blockInfos_inout, int numBlocks)
{
	#pragma omp parallel num_threads(numBlocks)
	{
		int blockId = omp_get_thread_num();
		block_info_t* blockInfo = blockInfos_inout[blockId];
		job_queue_convertToJDSPartial(spmCsrCounterpart, permutation, &blockInfo->queue);
	}
}

void block_convertFromPartialJDSToHybridJDSCSRMultiple(
		block_info_t** blockInfos_inout, int numBlocks)
{
	#pragma omp parallel num_threads(numBlocks)
	{
		int blockId = omp_get_thread_num();
		block_info_t* blockInfo = blockInfos_inout[blockId];

		// convert job queue in block info from partial JDS to hybrid_JDS_CSR
		// and store it in a new queue
		job_queue_t hybridJdsCsrQueue;
		job_queue_init(&hybridJdsCsrQueue, blockInfo->queue.stealTreshold);
		job_queue_convertFromPartialJDSToHybridJDSCSR(&blockInfo->queue, &hybridJdsCsrQueue);

		// empty queue in block info
		job_queue_deleteNonPtr(&blockInfo->queue);

		// fill queue in block info from hybrid_JDS_CSR queue
		job_queue_copy(&hybridJdsCsrQueue, &blockInfo->queue);

		// no need to clean up hybrid_JDS_CSR queue since the pointers
		// it contains will be swept when block info is destroyed.
	}
}

int block_localityAwareStealHalf(block_info_t* myBlock, block_info_t** blockInfos, int numBlocks)
{
	// try to steal from current victim
	job_batch_t* stolenBatchStart = NULL;
	int stealCount = 0;

	int victimId = myBlock->sched.ringSched.victimId;
	block_info_t* victimThread = blockInfos[victimId];

	// steal from the last victim if victim is not the current thread
	if(myBlock->id != victimId)
		job_queue_getLastBatchHalf(&victimThread->queue, &stolenBatchStart, &stealCount);

	while(stealCount <= 0)
	{
	#if PROGRAM_MODE >= TRACE_MODE
			PRINTF("localityAwareStealHalf: P%d (steal_range: %d) tried"	\
					"to steal from P%d but nothing to steal.\n",
					myBlock->id, myBlock->sched.ringSched.stealRange, victimId);
	#endif

			if(!myBlock->scheduler_findVictim(myBlock, blockInfos, numBlocks))
				break;

			victimId = myBlock->sched.ringSched.victimId;
			victimThread = blockInfos[victimId];

			job_queue_getLastBatchHalf(&victimThread->queue, &stolenBatchStart, &stealCount);
		}

		// Nothing to steal
		if(stealCount <= 0)
		{
	#if PROGRAM_MODE >= TRACE_MODE
			PRINTF("Queue-%d => nothing to steal. Finished for now.\n", myBlock->id);
	#endif
			return FALSE;
		}

		int i;
		list_head_t* head = &stolenBatchStart->head;
		for(i = 0; i < stealCount; ++i)
		{
			job_queue_addStolenBatchWithoutLocking(&myBlock->queue, job_batch_getBatch(head));
			head = head->next;
		}

	#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Queue-%d stole %d from Queue-%d(%d)\n", myBlock->id, stealCount,
				victimId, victimThread->queue.executableBatchCount);
	#endif

		return TRUE;
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
