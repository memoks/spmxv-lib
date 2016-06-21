
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "include/config.h"
#include "include/util/utility.h"
#include "include/scheduler/job_batch.h"
#include "include/scheduler/job_queue.h"
#include "include/scheduler/ring_scheduler.h"
#include "include/scheduler/tree_scheduler.h"
#include "include/task_decomposition/partitioning.h"
#include "inc/spmxv_gws.h"

#include "inc/spmxv_dws.h"

// Core
// ---------------------------------------------------------------------------------------------------------------------

// Row Parallel
// ---------------------------------------------------------------------------------------------------------------------

void spmxv_dws_rowWise_CSR(dws_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	int numWarps = args->numWarps;
	block_info_t** warps = args->warps;
	spm_cmp_t* spmCsr = args->spmCsr;

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("\n\n");
	PRINTF("---------------- Distribute Work Stealing algorithm execution trace ----------------\n");
#endif

	lg_t** execHistoryPerWarp = (lg_t**) malloc(sizeof(lg_t*) * numWarps);

	#pragma omp parallel num_threads(numWarps)
	{
		int warpId = omp_get_thread_num();
		block_info_t* myBlock = warps[warpId];
		lg_t* warpExecHistHead = lg_new();
		execHistoryPerWarp[warpId] = warpExecHistHead;

		job_batch_t* currBatch = job_queue_getFirstBatch(&myBlock->queue);
		while(currBatch != NULL)
		{
			JOB_BATCH_EXECUTE_SUBMTX_CSR(currBatch, spmCsr, x, y_inout);

#if PROGRAM_MODE >= TRACE_MODE
			char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
			job_batch_toStringDetailed(currBatch, (void*) spmCsr, temp);
			PRINTF("Warp-%d executed %s\n", warpId, temp);
#endif

			lg_addTailData((void*) currBatch, warpExecHistHead);
			currBatch = job_queue_getFirstBatch(&myBlock->queue);
		}

#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Warp %d finished execution, looking for victims...\n", myBlock->id);
#endif

		int goon = block_localityAwareStealHalf(myBlock, warps, numWarps);
		while(goon)
		{
			while(job_queue_hasStolenBatch(&myBlock->queue))
			{
				currBatch = job_queue_removeStolenBatchWithoutLocking(&myBlock->queue);
				JOB_BATCH_EXECUTE_SUBMTX_CSR(currBatch, spmCsr, x, y_inout);

#if PROGRAM_MODE >= TRACE_MODE
				char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
				job_batch_toStringDetailed(currBatch, (void*) spmCsr, temp);
				PRINTF("Warp-%d executed STOLEN-BATCH: %s\n", warpId, temp);
#endif
				lg_addTailData((void*) currBatch, warpExecHistHead);
			}

			goon = block_localityAwareStealHalf(myBlock, warps, numWarps);
		}
	}

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("---------------------------------- End of execution trace -----------------------------------\n");
	PRINTF("\n");
#endif

	// return values
	args->execHistoryPerWarp = execHistoryPerWarp;
}

void spmxv_dws_sharedBlock_rowWise_CSR(dws_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	int numThreads = args->staticArgs.numThreads;
	int threadsPerBlock = args->staticArgs.threadsPerBlock;
	int numWarps = args->numWarps;
	block_info_t** warps = args->warps;
	spm_cmp_t* spmCsr = args->spmCsr;

	lg_t** execHistoryPerThread = (lg_t**) malloc(sizeof(lg_t*) * numThreads);
	lg_t** execHistoryPerWarp = (lg_t**) malloc(sizeof(lg_t*) * numWarps);

	int i;
	for(i = 0; i < numWarps; ++i)
		execHistoryPerWarp[i] = lg_new();

	#pragma omp parallel num_threads(numThreads)
	{
		int globalThreadId = omp_get_thread_num();
		int threadId = globalThreadId % threadsPerBlock;
		int warpId = globalThreadId / threadsPerBlock;

		lg_t* execHistHead = lg_new();
		lg_t* blockExecHistHead = execHistoryPerWarp[warpId];

		block_info_t* myWarp = warps[warpId];
		job_batch_t* currBatch = NULL;

		// Execute local jobs (loop works as long as there is a local job in block's job queue or a detached job left)
		// ------------------------------------------------------------------------------------------------------------
		while(myWarp->queue.state == JOB_QUEUE_STATE_INITIAL || currBatch != NULL)
		{
			omp_set_lock(&myWarp->queue.writeLock);
			// check queue state to see whether there are any jobs or not
			if(myWarp->queue.state == JOB_QUEUE_STATE_INITIAL)
			{
				currBatch = job_queue_getFirstBatchWithoutLocking(&myWarp->queue);

				// If queue is empty then set its state to EMPTY.
				if(currBatch == NULL)
				{
					myWarp->queue.state = JOB_QUEUE_STATE_EMPTY;
				}
				else
				{
					// record block scheduling info
					lg_addTailData(currBatch, blockExecHistHead);

#if PROGRAM_MODE >= TRACE_MODE
					char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
					job_batch_toStringDetailed(currBatch, NULL, temp);
					PRINTF("Warp-%d, Thread-%d taken local batch %s.\n", warpId, threadId, temp);
#endif
				}
			}
			omp_unset_lock(&myWarp->queue.writeLock);

			// Execute detached job batch
			if(currBatch != NULL)
			{
				// record thread scheduling info
				lg_addTailData(currBatch, execHistHead);

				// Execute job batch
				JOB_BATCH_EXECUTE_SUBMTX_CSR(currBatch, spmCsr, x, y_inout);

#if PROGRAM_MODE >= TRACE_MODE
				char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
				job_batch_toStringDetailed(currBatch, NULL, temp);
				PRINTF("Warp-%d, Thread-%d executed local batch %s.\n", warpId, threadId, temp);
#endif

				// reset current batch
				currBatch = NULL;
			}
		}

		int currThreadLookingForSteals = FALSE;

#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Warp-%d Thread-%d done local execution. " \
				"Moving on (Global_id: %d). Queue State: %d\n",
				warpId, threadId, globalThreadId, myWarp->queue.state);
#endif

		// Execute stolen jobs (loop works until there are no steal attempts or detached jobs left)
		// ------------------------------------------------------------------------------------------------------------
		while(myWarp->queue.state != JOB_QUEUE_STATE_DONE || currBatch != NULL)
		{
			// Execute job batch if current thread has any detached jobs
			// Then record scheduling information for thread itself and its block (latter requires locking)
			if(currBatch != NULL)
			{
				// record thread execution history
				lg_addTailData(currBatch, execHistHead);

				// record block execution histoyr
				omp_set_lock(&myWarp->queue.writeLock);
				lg_addTailData(currBatch, blockExecHistHead);
				omp_unset_lock(&myWarp->queue.writeLock);

				// Execute job batch
				JOB_BATCH_EXECUTE_SUBMTX_CSR(currBatch, spmCsr, x, y_inout);

#if PROGRAM_MODE >= TRACE_MODE
				char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
				job_batch_toStringDetailed(currBatch, NULL, temp);
				PRINTF("Warp-%d Thread-%d executed stolen batch %s.\n", warpId, threadId, temp);
#endif


				// reset current batch pointer
				currBatch = NULL;
			}

			omp_set_lock(&myWarp->queue.writeLock);
			{
				if(myWarp->queue.state == JOB_QUEUE_STATE_EMPTY)
				{
					currBatch = job_queue_removeStolenBatchWithoutLocking(&myWarp->queue);

					// Prepare block for stealing if there are no stolen batches left
					if(currBatch == NULL)
					{
						myWarp->queue.state = JOB_QUEUE_STATE_STEALING;
						currThreadLookingForSteals = TRUE;

#if PROGRAM_MODE >= TRACE_MODE
						PRINTF("Warp-%d Thread-%d CurrBatch is NULL. Queue State is " \
								"block for stealing (queue_state: %d Global_id: %d).\n",
								warpId, threadId, myWarp->queue.state, globalThreadId);
#endif
					}
#if PROGRAM_MODE >= TRACE_MODE
					else
					{
						PRINTF("Warp-%d Thread-%d detached a stolen job (global_id: %d).\n", warpId, threadId, globalThreadId);
					}
#endif
				}
			}
			omp_unset_lock(&myWarp->queue.writeLock);

			if(currThreadLookingForSteals == TRUE)
			{
#if PROGRAM_MODE >= TRACE_MODE
				PRINTF("Warp-%d Thread-%d looking to steal jobs (Global-thread-id %d)...\n",
						warpId, threadId, globalThreadId);
#endif

				int goon = block_localityAwareStealHalf(myWarp, warps, numWarps);

				// secure a batch for stealing thread first
				currBatch = job_queue_removeStolenBatchWithoutLocking(&myWarp->queue);

				// TODO Beware (syntactic sugar) !! If STATE values change for some reason, do not leave it like this.
				omp_set_lock(&myWarp->queue.writeLock);
				myWarp->queue.state = goon;
				omp_unset_lock(&myWarp->queue.writeLock);

				// Safer alternative for the block above
				// omp_set_lock(&myBlock->queue.writeLock);
				// {
				//	if(goon == FALSE)
				//		myBlock->queue.state = JOB_QUEUE_STATE_DONE;
				//	else
				//		myBlock->queue.state = JOB_QUEUE_STATE_EMPTY;
				// }
				// omp_unset_lock(&myBlock->queue.writeLock);

				currThreadLookingForSteals = FALSE;

#if PROGRAM_MODE >= TRACE_MODE
				if(myWarp->queue.state == JOB_QUEUE_STATE_EMPTY)
					PRINTF("Warp-%d Thread-%d (global-thread-id %d) steal successful. Detached a stolen job.\n",
							warpId, threadId, globalThreadId);
				else
					PRINTF("Warp-%d Thread-%d (global-thread-id %d) Nothing to steal. DONE.\n",
							warpId, threadId, globalThreadId);
#endif
			}
		}

#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Warp-%d Thread-%d literally DONE (global_id: %d). Queue state: %d.\n",
				warpId, threadId, globalThreadId, myWarp->queue.state);
#endif

		// add recored scheduling info
		execHistoryPerThread[globalThreadId] = execHistHead;
	}

	// return values
	args->execHistoryPerThread = execHistoryPerThread;
	args->execHistoryPerWarp = execHistoryPerWarp;
}

// ---------------------------------------------------------------------------------------------------------------------

void spmxv_dws_rowWise_JDS(dws_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("\n\n");
	PRINTF("---------------- Distribute Work Stealing algorithm execution trace ----------------\n");
#endif

	int numWarps = args->numWarps;
	block_info_t** warps = args->warps;

	lg_t** execHistoryPerWarp = (lg_t**) malloc(sizeof(lg_t*) * numWarps);

	#pragma omp parallel num_threads(numWarps)
	{
		int warpId = omp_get_thread_num();
		block_info_t* myWarp = warps[warpId];
		lg_t* execHistHead = lg_new();
		execHistoryPerWarp[warpId] = execHistHead;

		job_batch_t* currBatch = job_queue_getFirstBatch(&myWarp->queue);
		while(currBatch != NULL)
		{
			JOB_BATCH_EXECUTE_PARTIAL_JDS(currBatch, x, y_inout);

#if PROGRAM_MODE >= TRACE_MODE
			char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
			job_batch_toStringDetailed(currBatch, NULL, temp);
			PRINTF("Warp-%d executed %s\n", warpId, temp);
#endif

			lg_addTailData(currBatch, execHistHead);
			currBatch = job_queue_getFirstBatch(&myWarp->queue);
		}

#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Warp-%d finished execution, looking for victims...\n", myWarp->id);
#endif

		int goon = block_localityAwareStealHalf(myWarp, warps, numWarps);
		while(goon)
		{
			while(job_queue_hasStolenBatch(&myWarp->queue))
			{
				currBatch = job_queue_removeStolenBatchWithoutLocking(&myWarp->queue);
				JOB_BATCH_EXECUTE_PARTIAL_JDS(currBatch, x, y_inout);

#if PROGRAM_MODE >= TRACE_MODE
				char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
				job_batch_toStringDetailed(currBatch, NULL, temp);
				PRINTF("Warp-%d executed STOLEN-BATCH: %s\n", warpId, temp);
#endif

				lg_addTailData(currBatch, execHistHead);
			}

			goon = block_localityAwareStealHalf(myWarp, warps, numWarps);
		}
	}

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("---------------------------------- End of execution trace -----------------------------------\n");
	PRINTF("\n");
#endif

	// return values
	args->execHistoryPerWarp = execHistoryPerWarp;
}

void spmxv_dws_sharedBlock_rowWise_JDS(dws_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	int numThreads = args->staticArgs.numThreads;
	int numWarps = args->numWarps;
	int threadsPerWarp = args->staticArgs.threadsPerBlock;
	block_info_t** warps = args->warps;

	lg_t** execHistoryPerThread = (lg_t**) malloc(sizeof(lg_t*) * numThreads);
	lg_t** execHistoryPerWarp = (lg_t**) malloc(sizeof(lg_t*) * numWarps);

	int i;
	for(i = 0; i < numWarps; ++i)
		execHistoryPerWarp[i] = lg_new();

	#pragma omp parallel num_threads(numThreads)
	{
		int globalThreadId = omp_get_thread_num();
		int threadId = globalThreadId % threadsPerWarp;
		int warpId = globalThreadId / threadsPerWarp;

		lg_t* threadExecHistHead = lg_new();
		lg_t* blockExecHistHead = execHistoryPerWarp[warpId];

		block_info_t* myWarp = warps[warpId];
		job_batch_t* currBatch = NULL;

		// Execute local jobs (loop works as long as there is a local job in block's job queue or a detached job left)
		// ------------------------------------------------------------------------------------------------------------
		while(myWarp->queue.state == JOB_QUEUE_STATE_INITIAL || currBatch != NULL)
		{
			omp_set_lock(&myWarp->queue.writeLock);
			// check queue state to see whether there are any jobs or not
			if(myWarp->queue.state == JOB_QUEUE_STATE_INITIAL)
			{
				currBatch = job_queue_getFirstBatchWithoutLocking(&myWarp->queue);

				// If queue is empty then set its state to EMPTY.
				if(currBatch == NULL)
				{
					myWarp->queue.state = JOB_QUEUE_STATE_EMPTY;
				}
				else
				{
					// record block scheduling info
					lg_addTailData((void*) currBatch, blockExecHistHead);

#if PROGRAM_MODE >= TRACE_MODE
					char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
					job_batch_toStringDetailed(currBatch, NULL, temp);
					PRINTF("Warp-%d, Thread-%d taken local batch %s.\n", warpId, threadId, temp);
#endif

				}
			}
			omp_unset_lock(&myWarp->queue.writeLock);

			// Execute detached job batch
			if(currBatch != NULL)
			{
				// record thread scheduling info
				lg_addTailData(currBatch, threadExecHistHead);

				// Execute job batch
				JOB_BATCH_EXECUTE_PARTIAL_JDS(currBatch, x, y_inout);

#if PROGRAM_MODE >= TRACE_MODE
				char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
				job_batch_toStringDetailed(currBatch, NULL, temp);
				PRINTF("Warp-%d, Thread-%d executed local batch %s.\n", warpId, threadId, temp);
#endif

				// reset current batch
				currBatch = NULL;
			}
		}

		int currThreadLookingForSteals = FALSE;

#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Warp-%d Thread-%d done local execution. " \
				"Moving on (Global_id: %d). Queue State: %d\n",
				warpId, threadId, globalThreadId, myWarp->queue.state);
#endif

		// Execute stolen jobs (loop works until there are no steal attempts or detached jobs left)
		// ------------------------------------------------------------------------------------------------------------
		while(myWarp->queue.state != JOB_QUEUE_STATE_DONE || currBatch != NULL)
		{
			// Execute job batch if current thread has any detached jobs
			// Then record scheduling information for thread itself and its block (latter requires locking)
			if(currBatch != NULL)
			{
				// record thread execution history
				lg_addTailData((void*) currBatch, threadExecHistHead);

				// record block execution histoyr
				omp_set_lock(&myWarp->queue.writeLock);
				lg_addTailData(currBatch, blockExecHistHead);
				omp_unset_lock(&myWarp->queue.writeLock);

				// Execute job batch
				JOB_BATCH_EXECUTE_PARTIAL_JDS(currBatch, x, y_inout);

#if PROGRAM_MODE >= TRACE_MODE
				char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
				job_batch_toStringDetailed(currBatch, NULL, temp);
				PRINTF("Warp-%d Thread-%d executed stolen batch %s.\n", warpId, threadId, temp);
#endif

				// reset current batch pointer
				currBatch = NULL;
			}

			omp_set_lock(&myWarp->queue.writeLock);
			{
				if(myWarp->queue.state == JOB_QUEUE_STATE_EMPTY)
				{
					currBatch = job_queue_removeStolenBatchWithoutLocking(&myWarp->queue);

					// Prepare block for stealing if there are no stolen batches left
					if(currBatch == NULL)
					{
						myWarp->queue.state = JOB_QUEUE_STATE_STEALING;
						currThreadLookingForSteals = TRUE;

#if PROGRAM_MODE >= TRACE_MODE
						PRINTF("Warp-%d Thread-%d CurrBatch is NULL. Queue State is " \
								"block for stealing (queue_state: %d Global_id: %d).\n",
								warpId, threadId, myWarp->queue.state, globalThreadId);
#endif

					}
#if PROGRAM_MODE >= TRACE_MODE
					else
					{
						PRINTF("Warp-%d Thread-%d detached a stolen job (global_id: %d).\n", warpId, threadId, globalThreadId);
					}
#endif
				}
			}
			omp_unset_lock(&myWarp->queue.writeLock);

			if(currThreadLookingForSteals == TRUE)
			{
#if PROGRAM_MODE >= TRACE_MODE
				PRINTF("Warp-%d Thread-%d looking to steal jobs (Global-thread-id %d)...\n",
						warpId, threadId, globalThreadId);
#endif

				int goon = block_localityAwareStealHalf(myWarp, warps, numWarps);

				// secure a batch for stealing thread first
				currBatch = job_queue_removeStolenBatchWithoutLocking(&myWarp->queue);

				// TODO Beware (syntactic sugar) !! If STATE values change for some reason, do not leave it like this.
				omp_set_lock(&myWarp->queue.writeLock);
				myWarp->queue.state = goon;
				omp_unset_lock(&myWarp->queue.writeLock);

				// Safer alternative for the block above
				// omp_set_lock(&myBlock->queue.writeLock);
				// {
				//	if(goon == FALSE)
				//		myBlock->queue.state = JOB_QUEUE_STATE_DONE;
				//	else
				//		myBlock->queue.state = JOB_QUEUE_STATE_EMPTY;
				// }
				// omp_unset_lock(&myBlock->queue.writeLock);

				currThreadLookingForSteals = FALSE;

#if PROGRAM_MODE >= TRACE_MODE
				if(myWarp->queue.state == JOB_QUEUE_STATE_EMPTY)
					PRINTF("Warp-%d Thread-%d (global-thread-id %d) steal successful. Detached a stolen job.\n",
							warpId, threadId, globalThreadId);
				else
					PRINTF("Warp-%d Thread-%d (global-thread-id %d) Nothing to steal. DONE.\n",
							warpId, threadId, globalThreadId);
#endif

			}
		}

#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Warp-%d Thread-%d literally DONE (global_id: %d). Queue state: %d.\n",
				warpId, threadId, globalThreadId, myWarp->queue.state);
#endif

		// add recored scheduling info
		execHistoryPerThread[globalThreadId] = threadExecHistHead;
	}

	// return values
	args->execHistoryPerThread = execHistoryPerThread;
	args->execHistoryPerWarp = execHistoryPerWarp;
}

void spmxv_dws_rowWise_hybrid_JDS_CSR(dws_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("\n\n");
	PRINTF("---------------- Distribute Work Stealing algorithm execution trace ----------------\n");
#endif

	int numWarps = args->numWarps;
	block_info_t** warps = args->warps;

	lg_t** execHistoryPerWarp = (lg_t**) malloc(sizeof(lg_t*) * numWarps);

	#pragma omp parallel num_threads(numWarps)
	{
		int warpId = omp_get_thread_num();
		block_info_t* myWarp = warps[warpId];
		lg_t* warpExecHistHead = lg_new();
		execHistoryPerWarp[warpId] = warpExecHistHead;

		job_batch_t* currBatch = job_queue_getFirstBatch(&myWarp->queue);
		while(currBatch != NULL)
		{
			lg_addTailData(currBatch, warpExecHistHead);

#if PROGRAM_MODE >= TRACE_MODE
			char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
			job_batch_toStringDetailed(currBatch, NULL, temp);
			PRINTF("Warp-%d executing: %s\n", warpId, temp);
#endif

			JOB_BATCH_EXECUTE_HYBRID_JDS_CSR(currBatch, x, y_inout);
			currBatch = job_queue_getFirstBatch(&myWarp->queue);
		}

#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Warp-%d finished execution, looking for victims...\n", myWarp->id);
#endif

		int goon = block_localityAwareStealHalf(myWarp, warps, numWarps);
		while(goon)
		{
			while(job_queue_hasStolenBatch(&myWarp->queue))
			{
				currBatch = job_queue_removeStolenBatchWithoutLocking(&myWarp->queue);
				lg_addTailData(currBatch, warpExecHistHead);

#if PROGRAM_MODE >= TRACE_MODE
				char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
				job_batch_toStringDetailed(currBatch, NULL, temp);
				PRINTF("Warp-%d stole: %s\n", warpId, temp);
#endif

				JOB_BATCH_EXECUTE_HYBRID_JDS_CSR(currBatch, x, y_inout);
			}

			goon = block_localityAwareStealHalf(myWarp, warps, numWarps);
		}
	}

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("________________________________EXECUTION_HISTORY_PER_WARP___________________________________\n");
		lg_job_batch_printMultiple(execHistoryPerWarp, numWarps);
	PRINTF("---------------------------------- End of execution trace -----------------------------------\n");
	PRINTF("\n");
#endif

	// return values
	args->execHistoryPerWarp = execHistoryPerWarp;
}

void spmxv_dws_sharedBlock_rowWise_hybrid_JDS_CSR(
		dws_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	int numThreads = args->staticArgs.numThreads;
	int numWarps = args->numWarps;
	int threadsPerWarp = args->staticArgs.threadsPerBlock;
	block_info_t** warps = args->warps;

	lg_t** execHistoryPerThread = (lg_t**) malloc(sizeof(lg_t*) * numThreads);
	lg_t** execHistoryPerWarp = (lg_t**) malloc(sizeof(lg_t*) * numWarps);

	int i;
	for(i = 0; i < numWarps; ++i)
		execHistoryPerWarp[i] = lg_new();

	#pragma omp parallel num_threads(numThreads)
	{
		int globalThreadId = omp_get_thread_num();
		int threadId = globalThreadId % threadsPerWarp;
		int warpId = globalThreadId / threadsPerWarp;

		lg_t* threadExecHistHead = lg_new();
		lg_t* warpExecHistHead = execHistoryPerWarp[warpId];

		block_info_t* myWarp = warps[warpId];
		job_batch_t* currBatch = NULL;

		// Execute local jobs (loop works as long as there is a local job in block's job queue or a detached job left)
		// ------------------------------------------------------------------------------------------------------------
		while(myWarp->queue.state == JOB_QUEUE_STATE_INITIAL || currBatch != NULL)
		{
			omp_set_lock(&myWarp->queue.writeLock);
			// check queue state to see whether there are any jobs or not
			if(myWarp->queue.state == JOB_QUEUE_STATE_INITIAL)
			{
				currBatch = job_queue_getFirstBatchWithoutLocking(&myWarp->queue);

				// If queue is empty then set its state to EMPTY.
				if(currBatch == NULL)
				{
					myWarp->queue.state = JOB_QUEUE_STATE_EMPTY;
				}
				else
				{
					// record block scheduling info
					lg_addTailData(currBatch, warpExecHistHead);

#if PROGRAM_MODE >= TRACE_MODE
					char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
					job_batch_toStringDetailed(currBatch, NULL, temp);
					PRINTF("Warp-%d, Thread-%d taken local batch %s.\n", warpId, threadId, temp);
#endif

				}
			}
			omp_unset_lock(&myWarp->queue.writeLock);

			// Execute detached job batch
			if(currBatch != NULL)
			{
				// record thread scheduling info
				lg_addTailData(currBatch, threadExecHistHead);

				// Execute job batch
				JOB_BATCH_EXECUTE_HYBRID_JDS_CSR(currBatch, x, y_inout);

#if PROGRAM_MODE >= TRACE_MODE
				char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
				job_batch_toStringDetailed(currBatch, NULL, temp);
				PRINTF("Warp-%d, Thread-%d executed local batch %s.\n", warpId, threadId, temp);
#endif

				// reset current batch
				currBatch = NULL;
			}
		}

		int currThreadLookingForSteals = FALSE;

#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Warp-%d Thread-%d done local execution. " \
				"Moving on (Global_id: %d). Queue State: %d\n",
				warpId, threadId, globalThreadId, myWarp->queue.state);
#endif

		// Execute stolen jobs (loop works until there are no steal attempts or detached jobs left)
		// ------------------------------------------------------------------------------------------------------------
		while(myWarp->queue.state != JOB_QUEUE_STATE_DONE || currBatch != NULL)
		{
			// Execute job batch if current thread has any detached jobs
			// Then record scheduling information for thread itself and its block (latter requires locking)
			if(currBatch != NULL)
			{
				// record thread execution history
				lg_addTailData(currBatch, threadExecHistHead);

				// record block execution histoyr
				omp_set_lock(&myWarp->queue.writeLock);
				lg_addTailData(currBatch, warpExecHistHead);
				omp_unset_lock(&myWarp->queue.writeLock);

				// Execute job batch
				JOB_BATCH_EXECUTE_HYBRID_JDS_CSR(currBatch, x, y_inout);

#if PROGRAM_MODE >= TRACE_MODE
				char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
				job_batch_toStringDetailed(currBatch, NULL, temp);
				PRINTF("Warp-%d Thread-%d executed stolen batch %s.\n", warpId, threadId, temp);
#endif

				// reset current batch pointer
				currBatch = NULL;
			}

			omp_set_lock(&myWarp->queue.writeLock);
			{
				if(myWarp->queue.state == JOB_QUEUE_STATE_EMPTY)
				{
					currBatch = job_queue_removeStolenBatchWithoutLocking(&myWarp->queue);

					// Prepare block for stealing if there are no stolen batches left
					if(currBatch == NULL)
					{
						myWarp->queue.state = JOB_QUEUE_STATE_STEALING;
						currThreadLookingForSteals = TRUE;

#if PROGRAM_MODE >= TRACE_MODE
						PRINTF("Warp-%d Thread-%d CurrBatch is NULL. Queue State is " \
								"block for stealing (queue_state: %d Global_id: %d).\n",
								warpId, threadId, myWarp->queue.state, globalThreadId);
#endif

					}
#if PROGRAM_MODE >= TRACE_MODE
					else
					{
						PRINTF("Warp-%d Thread-%d detached a stolen job (global_id: %d).\n", warpId, threadId, globalThreadId);
					}
#endif
				}
			}
			omp_unset_lock(&myWarp->queue.writeLock);

			if(currThreadLookingForSteals == TRUE)
			{
#if PROGRAM_MODE >= TRACE_MODE
				PRINTF("Warp-%d Thread-%d looking to steal jobs (Global-thread-id %d)...\n",
						warpId, threadId, globalThreadId);
#endif

				int goon = block_localityAwareStealHalf(myWarp, warps, numWarps);

				// secure a batch for stealing thread first
				currBatch = job_queue_removeStolenBatchWithoutLocking(&myWarp->queue);

				// TODO Beware (syntactic sugar) !! If STATE values change for some reason, do not leave it like this.
				omp_set_lock(&myWarp->queue.writeLock);
				myWarp->queue.state = goon;
				omp_unset_lock(&myWarp->queue.writeLock);

				// Safer alternative for the block above
				// omp_set_lock(&myBlock->queue.writeLock);
				// {
				//	if(goon == FALSE)
				//		myBlock->queue.state = JOB_QUEUE_STATE_DONE;
				//	else
				//		myBlock->queue.state = JOB_QUEUE_STATE_EMPTY;
				// }
				// omp_unset_lock(&myBlock->queue.writeLock);

				currThreadLookingForSteals = FALSE;

#if PROGRAM_MODE >= TRACE_MODE
				if(myWarp->queue.state == JOB_QUEUE_STATE_EMPTY)
					PRINTF("Warp-%d Thread-%d (global-thread-id %d) steal successful. Detached a stolen job.\n",
							warpId, threadId, globalThreadId);
				else
					PRINTF("Warp-%d Thread-%d (global-thread-id %d) Nothing to steal. DONE.\n",
							warpId, threadId, globalThreadId);
#endif

			}
		}

#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Warp-%d Thread-%d literally DONE (global_id: %d). Queue state: %d.\n",
				warpId, threadId, globalThreadId, myWarp->queue.state);
#endif

		// add recored scheduling info
		execHistoryPerThread[globalThreadId] = threadExecHistHead;
	}

	// return values
	args->execHistoryPerThread = execHistoryPerThread;
	args->execHistoryPerWarp = execHistoryPerWarp;
}


// Pre-multiplication routines / data preparations
// ---------------------------------------------------------------------------------------------------------------------

void spmxv_dws_prepFromSubMtxTree_rowWise_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromSubMtxTree_rowWise_CSR\n");

	int numWarps = args->numWarps;
	int stealTreshold = args->stealTreshold;
	int stealingScheme = args->stealingScheme;
	sub_mtx_tree_t* subMtxHead = args->staticArgs.subMtxHead;
	spm_cmp_t* spmCsr = args->spmCsr;

	block_info_t** warps = block_newMultiple(numWarps, stealTreshold, stealingScheme);

	// Find uppermost nodes and put everything under them to different queues
	tree_node_t** nodesPerWarp = sub_mtx_tree_getJobDistribution(subMtxHead, numWarps);

#if PROGRAM_MODE >= TRACE_MODE
	int j;
	for(j = 0; j < numWarps; ++j)
	{
		if(nodesPerWarp[j] == NULL)
			PRINTF("nodesPerWarp[%d] is NULL\n", j);
		else
			PRINTF("nodesPerWarp[%d] child_count: %d\n",
					j, tree_node_getLeafCount(nodesPerWarp[j]));
	}
#endif

	#pragma omp parallel num_threads(numWarps)
	{
		int warpId = omp_get_thread_num();
		block_info_t* currThreadInfo = warps[warpId];

		// Put every child node under corresponding assignment to queue
		if(nodesPerWarp[warpId] != NULL)
		{
			job_queue_fillFromSubMtxTree_rowWiseColNet_CSR(
					spmCsr, &currThreadInfo->queue, nodesPerWarp[warpId]);
		}
	}

	// local clean-up
	free(nodesPerWarp);

	// return values
	args->warps = warps;
}

void spmxv_dws_prepFromSubMtxTreeScatter_rowWise_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromSubMtxTreeScatter_rowWise_CSR\n");

	int numWarps = args->numWarps;
	int numBlocks = args->staticArgs.numBlocks;
	int numThreads = args->staticArgs.numThreads;
	int threadsPerWarp = args->staticArgs.threadsPerBlock;
	int stealTreshold = args->stealTreshold;
	int stealingScheme = args->stealingScheme;
	sub_mtx_tree_t* subMtxHead = args->staticArgs.subMtxHead;
	spm_cmp_t* spmCsr = args->spmCsr;

	block_info_t** warps = block_newMultiple(numWarps, stealTreshold, stealingScheme);

	// Find uppermost nodes and put everything under them to different queues
	tree_node_t** nodesPerWarp = sub_mtx_tree_getJobDistribution(subMtxHead, numBlocks);

#if PROGRAM_MODE >= TRACE_MODE
	int j;
	for(j = 0; j < numWarps; ++j)
	{
		if(nodesPerWarp[j] == NULL)
			PRINTF("nodesPerWarp[%d] is NULL\n", j);
		else
			PRINTF("nodesPerWarp[%d] child_count: %d\n",
					j, tree_node_getLeafCount(nodesPerWarp[j]));
	}
#endif

	#pragma omp parallel num_threads(numBlocks)
	{
		int blockId = omp_get_thread_num();
		int firstThreadGlobalId = blockId * threadsPerWarp;
		job_queue_t* jobQueuesPerWarp[threadsPerWarp];

		// collect job queues in this block
		int i;
		for(i = 0; i < threadsPerWarp; ++i)
			jobQueuesPerWarp[i] = &warps[firstThreadGlobalId + i]->queue;

		// Put every child node under corresponding assignment
		// into given queue pointer array in a scatter fashion
		job_queue_fillFromSubMtxTreeScatter_rowWiseColNet_CSR(
				spmCsr, jobQueuesPerWarp, threadsPerWarp, nodesPerWarp[blockId]);
	}

	// local clean-up
	free(nodesPerWarp);

	// return values
	args->warps = warps;
}

void spmxv_dws_prepFromSubMtxList_rowWise_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromSubMtxList_rowWise_CSR\n");

	int numWarps = args->numWarps;
	int stealTreshold = args->stealTreshold;
	int stealingScheme = args->stealingScheme;
	lg_t* subMtxList = args->staticArgs.subMtxList;

	block_info_t** warps = block_newMultiple(numWarps, stealTreshold, stealingScheme);

	// convert sub-matrix list to job batch list
	lg_t* jobBatchList = NULL;
	job_batch_extractJobBatchListFromSubMtxList_rowWiseColNet_CSR(subMtxList, &jobBatchList);

	// divide single job batch list to "numBlocks" number of lists
	lg_t** jobBatchListPerWarp = NULL;
	lg_split(jobBatchList, numWarps, &jobBatchListPerWarp);

	#pragma omp parallel num_threads(numWarps)
	{
		int warpId = omp_get_thread_num();
		lg_t* warpJobBatchList = jobBatchListPerWarp[warpId];
		block_info_t* threadInfo = warps[warpId];

		// add everything in corresponding job batch list into threads's job queue
		job_queue_addJobBatchListWithoutLocking(&threadInfo->queue, warpJobBatchList);

		// delete thread job batch list (which is empty by now)
		lg_deleteShallow(warpJobBatchList);
	}

	// global clean-up
	free(jobBatchListPerWarp);
	lg_deleteShallow(jobBatchList);

	// return values
	args->warps = warps;
}

void spmxv_dws_prepFromPartitionVector_rowWise_colNet_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromPartitionVector_rowWise_colNet_CSR\n");

	int numWarps = args->numWarps;
	int stealTreshold = args->stealTreshold;
	int stealingScheme = args->stealingScheme;
	partitioning_vector_t* rowPartitioning = args->rowPartitioning;
	int partitionCount = rowPartitioning->length - 1;
	spm_cmp_t* spmCsr = args->spmCsr;

	block_info_t** warps = block_newMultiple(numWarps, stealTreshold, stealingScheme);

	int minJobBatchCountPerBlock = partitionCount / numWarps;
	int remainingJobBatchCount = partitionCount % numWarps;
	int* jobBatchCountAcc = (int*) malloc(sizeof(int) * (numWarps + 1));
	jobBatchCountAcc[0] = 0;

	#pragma omp parallel num_threads(numWarps)
	{
		int threadId = omp_get_thread_num();
		int blockJobBatchCount = minJobBatchCountPerBlock;
		if(remainingJobBatchCount > threadId)
			++blockJobBatchCount;

		jobBatchCountAcc[threadId + 1] = blockJobBatchCount;

		#pragma omp barrier
		#pragma omp single
		{
			accumulate(jobBatchCountAcc, numWarps + 1);
		}
		#pragma omp barrier

		int jobBatchStartIndex = jobBatchCountAcc[threadId];
		int jobBatchEndIndex = jobBatchCountAcc[threadId + 1];
		job_queue_fillFromParitioningVector_rowWiseColNet_CSR(
				spmCsr, &warps[threadId]->queue,
				rowPartitioning, jobBatchStartIndex, jobBatchEndIndex);
	}

	// clean up
	free(jobBatchCountAcc);

	// return values
	args->warps = warps;
}

void spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_CSR\n");

	int numWarps = args->numWarps;
	int numBlocks = args->staticArgs.numBlocks;
	int threadsPerBlock = args->staticArgs.threadsPerBlock;
	int stealTreshold = args->stealTreshold;
	int stealingScheme = args->stealingScheme;
	partitioning_vector_t* rowPartitioning = args->rowPartitioning;
	int partitionCount = rowPartitioning->length - 1;

	spm_cmp_t* spmCsr = args->spmCsr;

	block_info_t** warps = block_newMultiple(numWarps, stealTreshold, stealingScheme);

	// calculate minimum sub-matrix count per block
	int minJobBatchCountPerBlock = partitionCount / numBlocks;
	int remainingJobBatchCount = partitionCount % numBlocks;
	int* blockJobBatchCountAcc = (int*)	malloc(sizeof(int) * (numBlocks + 1));
	blockJobBatchCountAcc[0] = 0;

	#pragma omp parallel num_threads(numBlocks)
	{
		// calculate thread IDs
		int blockId = omp_get_thread_num();
		int blockFirstThreadId = blockId * threadsPerBlock;

		// create a job queue pointer array for each block
		job_queue_t* jobQueuePtrArr[threadsPerBlock];
		int j;
		for(j = 0; j < threadsPerBlock; ++j)
			jobQueuePtrArr[j] = &warps[blockFirstThreadId + j]->queue;

		// calculate sub-matrix count for current block
		int blockJobBatchCount = minJobBatchCountPerBlock;
		if(remainingJobBatchCount > blockId)
			++blockJobBatchCount;

		blockJobBatchCountAcc[blockId + 1] = blockJobBatchCount;

		#pragma omp barrier
		#pragma omp single
		{
			accumulate(blockJobBatchCountAcc, numBlocks + 1);
		}
		#pragma omp barrier

		// calculate partition vector start and end index for current block
		int blockJobBatchStartIndex = blockJobBatchCountAcc[blockId];
		int blockJobBatchEndIndex = blockJobBatchCountAcc[blockId + 1];

		// fill job queues in jobQueuePtrArr in scatter fashion
		job_queue_fillFromParitioningVectorScatter_rowWiseColNet_CSR(
				spmCsr, jobQueuePtrArr, threadsPerBlock, rowPartitioning,
				blockJobBatchStartIndex, blockJobBatchEndIndex);
	}

	// clean up
	free(blockJobBatchCountAcc);

	// return values
	args->warps = warps;
}

void spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_CSR\n");

	int numWarps = args->numWarps;
	int stealTreshold = args->stealTreshold;
	int stealingScheme = args->stealingScheme;
	partitioning_vector_t* rowPartitioning = args->rowPartitioning;
	int partitionCount = rowPartitioning->length - 1;
	spm_cmp_t* spmCsr = args->spmCsr;

	block_info_t** warps = block_newMultiple(numWarps, stealTreshold, stealingScheme);

	DECIMAL* subMtxCountPerWarp = partitioning_calculateSubMtxCountPerBlockPowerOf2(numWarps, partitionCount);

	accumulate(subMtxCountPerWarp, numWarps);

	#pragma omp parallel num_threads(numWarps)
	{
		int blockId = omp_get_thread_num();
		block_info_t* block = warps[blockId];

		int jobBatchStartIndex = (blockId == 0) ? 0 : subMtxCountPerWarp[blockId - 1];
		int jobBatchEndIndex = subMtxCountPerWarp[blockId];

		job_queue_fillFromParitioningVector_rowWiseColNet_CSR(
				spmCsr, &block->queue, rowPartitioning, jobBatchStartIndex, jobBatchEndIndex);
	}

	// clean up
	free(subMtxCountPerWarp);

	// return values
	args->warps = warps;
}

void spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_CSR\n");

	int numWarps = args->numWarps;
	int numBlocks = args->staticArgs.numBlocks;
	int threadsPerBlock = args->staticArgs.threadsPerBlock;
	int stealTreshold = args->stealTreshold;
	int stealingScheme = args->stealingScheme;
	partitioning_vector_t* rowPartitioning = args->rowPartitioning;
	int partitionCount = rowPartitioning->length - 1;
	spm_cmp_t* spmCsr = args->spmCsr;

	block_info_t** warps = block_newMultiple(numWarps, stealTreshold, stealingScheme);

	// calculate minimum sub-matrix count per block
	int minJobBatchCountPerBlock = partitionCount / numBlocks;
	int remainingJobBatchCount = partitionCount % numBlocks;
	int* blockJobBatchCountAcc = (int*)	malloc(sizeof(int) * (numBlocks + 1));
	blockJobBatchCountAcc[0] = 0;

	#pragma omp parallel num_threads(numBlocks)
	{
		// calculate thread IDs
		int blockId = omp_get_thread_num();
		int blockFirstThreadId = blockId * threadsPerBlock;

		// create a job queue pointer array for each block
		job_queue_t* jobQueuePtrArr[threadsPerBlock];
		int j;
		for(j = 0; j < threadsPerBlock; ++j)
			jobQueuePtrArr[j] = &warps[blockFirstThreadId + j]->queue;

		// calculate sub-matrix count for current block
		int blockJobBatchCount = minJobBatchCountPerBlock;
		if(remainingJobBatchCount > blockId)
			++blockJobBatchCount;

		blockJobBatchCountAcc[blockId + 1] = blockJobBatchCount;

		#pragma omp barrier
		#pragma omp single
		{
			accumulate(blockJobBatchCountAcc, numBlocks + 1);
		}
		#pragma omp barrier

		// calculate partition vector start and end index for current block
		int blockJobBatchStartIndex = blockJobBatchCountAcc[blockId];
		int blockJobBatchEndIndex = blockJobBatchCountAcc[blockId + 1];

		// fill job queues in jobQueuePtrArr in scatter fashion
		job_queue_fillFromParitioningVectorScatter_rowWiseColNet_CSR(
				spmCsr, jobQueuePtrArr, threadsPerBlock, rowPartitioning,
				blockJobBatchStartIndex, blockJobBatchEndIndex);
	}

	int jobBatchCount = rowPartitioning->length;

	// clean up
	free(blockJobBatchCountAcc);

	// return values
	args->warps = warps;
}

// ----------------------------------------------------------------------------------------------------------------------

void spmxv_dws_prepFromSubMtxTree_rowWise_JDS(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromSubMtxTree_rowWise_JDS\n");

	spmxv_dws_prepFromSubMtxTree_rowWise_CSR(args);

	spm_cmp_t* spmCsr = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertToJDSPartialMultiple(spmCsr, permutation, warps, numWarps);
}

void spmxv_dws_prepFromSubMtxTreeScatter_rowWise_JDS(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromSubMtxTreeScatter_rowWise_JDS\n");

	spmxv_dws_prepFromSubMtxTreeScatter_rowWise_CSR(args);

	spm_cmp_t* spmCsr = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertToJDSPartialMultiple(spmCsr, permutation, warps, numWarps);
}

void spmxv_dws_prepFromSubMtxList_rowWise_JDS(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromSubMtxList_rowWise_JDS\n");

	spmxv_dws_prepFromSubMtxList_rowWise_CSR(args);

	spm_cmp_t* spmCsr = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertToJDSPartialMultiple(spmCsr, permutation, warps, numWarps);
}

void spmxv_dws_prepFromPartitionVector_rowWise_colNet_JDS(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromPartitionVector_rowWise_colNet_JDS\n");

	spmxv_dws_prepFromPartitionVector_rowWise_colNet_CSR(args);

	spm_cmp_t* spmCsr = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertToJDSPartialMultiple(spmCsr, permutation, warps, numWarps);
}

void spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_JDS(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_JDS\n");

	spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_CSR(args);

	spm_cmp_t* spmCsr = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertToJDSPartialMultiple(spmCsr, permutation, warps, numWarps);
}

void spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_JDS(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_JDS\n");

	spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_CSR(args);

	spm_cmp_t* spmCsr = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertToJDSPartialMultiple(spmCsr, permutation, warps, numWarps);
}

void spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_JDS(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_JDS\n");

	spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_CSR(args);

	spm_cmp_t* spmCsr = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertToJDSPartialMultiple(spmCsr, permutation, warps, numWarps);
}

// ----------------------------------------------------------------------------------------------------------------------

void spmxv_dws_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR\n");

	spmxv_dws_prepFromSubMtxTree_rowWise_JDS(args);

	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertFromPartialJDSToHybridJDSCSRMultiple(warps, numWarps);
}

void spmxv_dws_prepFromSubMtxTreeScatter_rowWise_hybrid_JDS_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromSubMtxTreeScatter_rowWise_hybrid_JDS_CSR\n");

	spmxv_dws_prepFromSubMtxTreeScatter_rowWise_JDS(args);

	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertFromPartialJDSToHybridJDSCSRMultiple(warps, numWarps);
}

void spmxv_dws_prepFromSubMtxList_rowWise_hybrid_JDS_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromSubMtxList_rowWise_hybrid_JDS_CSR\n");

	spmxv_dws_prepFromSubMtxList_rowWise_JDS(args);

	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertFromPartialJDSToHybridJDSCSRMultiple(warps, numWarps);
}

void spmxv_dws_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR\n");

	spmxv_dws_prepFromPartitionVector_rowWise_colNet_JDS(args);

	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertFromPartialJDSToHybridJDSCSRMultiple(warps, numWarps);
}

void spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_hybrid_JDS_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_hybrid_JDS_CSR\n");

	spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_JDS(args);

	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertFromPartialJDSToHybridJDSCSRMultiple(warps, numWarps);
}

void spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_hybrid_JDS_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_hybrid_JDS_CSR\n");

	spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_JDS(args);

	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertFromPartialJDSToHybridJDSCSRMultiple(warps, numWarps);
}

void spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_hybrid_JDS_CSR(dws_args_t* args)
{
	DEBUG("spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_hybrid_JDS_CSR\n");

	spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_JDS(args);

	block_info_t** warps = args->warps;
	int numWarps = args->numWarps;
	block_convertFromPartialJDSToHybridJDSCSRMultiple(warps, numWarps);
}

// Data structure clean-up
// ----------------------------------------------------------------------------------------------------------------------

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
