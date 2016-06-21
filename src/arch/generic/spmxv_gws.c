
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "include/timer/custom_timer.h"

#include "inc/spmxv_gws.h"

// Core
// ---------------------------------------------------------------------------------------------------------------------

// Row Parallel
// ---------------------------------------------------------------------------------------------------------------------

void spmxv_gws_multiply_CSR(gws_args_t* args, vector_real_t* x, vector_real_t* y_inout)
//		spm_cmp_t* spmCsr, vector_real_t* x, vector_real_t* y_inout,
//		job_queue_t* globalQueue, int numThreads,
//		lg_t*** execHistoryPerBlock_out)
{
	spm_cmp_t* spmCsr = args->staticArgs.spmCsr;
	job_queue_t* globalQueue = args->globalQueue;
	int numThreads = args->staticArgs.numThreads;

	lg_t** execHistoryPerThread = (lg_t**) malloc(sizeof(lg_t*) * numThreads);

	#if PROGRAM_MODE >= TRACE_MODE
		int* jobCountPerThread = (int*) malloc(sizeof(int) * numThreads);
		int i;
		for(i = 0; i < numThreads; ++i)
			jobCountPerThread[i] = 0;

		PRINTF("\n\n");
	#endif

	#pragma omp parallel num_threads(numThreads)
	{
		int threadId = omp_get_thread_num();

		execHistoryPerThread[threadId] = lg_new();
		lg_t* threadExecHistory = execHistoryPerThread[threadId];

		job_batch_t* currBatch = NULL;
		do
		{
			currBatch = job_queue_getFirstBatch(globalQueue);
			if(currBatch == NULL)
				break;

			lg_addTailData((void*) currBatch, threadExecHistory);
			JOB_BATCH_EXECUTE_SUBMTX_CSR(currBatch, spmCsr, x, y_inout);

			#if PROGRAM_MODE >= TRACE_MODE
				++jobCountPerThread[threadId];
			#endif

		} while(TRUE);
	}

	#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Printing job count executed by each thread;\n");
		for(i = 0; i < numThreads; ++i)
			if(jobCountPerThread[i] > 0)
				PRINTF("Thread-%d: %d jobs.\n", i, jobCountPerThread[i]);
		free(jobCountPerThread);

		PRINTF("\n");
	#endif

	// return values
	args->execHistoryPerThread = execHistoryPerThread;
}

void spmxv_gws_multiply_JDS(
		gws_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	job_queue_t* globalQueue = args->globalQueue;
	int numThreads = args->staticArgs.numThreads;

	lg_t** execHistoryPerThread = (lg_t**) malloc(sizeof(lg_t*) * numThreads);

	#if PROGRAM_MODE >= TRACE_MODE
		int* jobCountPerThread = (int*) malloc(sizeof(int) * numThreads);
		int i;
		for(i = 0; i < numThreads; ++i)
			jobCountPerThread[i] = 0;

		PRINTF("\n\n");
	#endif

	#pragma omp parallel num_threads(numThreads)
	{
		int threadId = omp_get_thread_num();

		execHistoryPerThread[threadId] = lg_new();
		lg_t* threadExecHistory = execHistoryPerThread[threadId];

		job_batch_t* currBatch = NULL;
		do
		{
			currBatch = job_queue_getFirstBatch(globalQueue);
			if(currBatch == NULL)
				break;

			lg_addTailData((void*) currBatch, threadExecHistory);
			JOB_BATCH_EXECUTE_PARTIAL_JDS(currBatch, x, y_inout);

			#if PROGRAM_MODE >= TRACE_MODE
				++jobCountPerThread[threadId];
			#endif

		} while(TRUE);
	}

	#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Printing job count executed by each thread;\n");
		for(i = 0; i < numThreads; ++i)
			if(jobCountPerThread[i] > 0)
				PRINTF("Thread-%d: %d jobs.\n", i, jobCountPerThread[i]);
		free(jobCountPerThread);

		PRINTF("\n");
	#endif

	// return values
	args->execHistoryPerThread = execHistoryPerThread;
}

void spmxv_gws_multiply_hybrid_JDS_CSR(
		gws_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	job_queue_t* globalQueue = args->globalQueue;
	int numThreads = args->staticArgs.numThreads;

	#if PROGRAM_MODE >= TRACE_MODE
		int* jobCountPerThread = (int*) malloc(sizeof(int) * numThreads);
		int i;
		for(i = 0; i < numThreads; ++i)
			jobCountPerThread[i] = 0;

		PRINTF("\n\n");
	#endif

		lg_t** execHistoryPerThread = (lg_t**) malloc(sizeof(lg_t*) * numThreads);

	#pragma omp parallel num_threads(numThreads)
	{
		int threadId = omp_get_thread_num();
		execHistoryPerThread[threadId] = lg_new();
		lg_t* threadExecHistory = execHistoryPerThread[threadId];

		job_batch_t* currBatch = NULL;
		do
		{
			currBatch = job_queue_getFirstBatch(globalQueue);
			if(currBatch == NULL)
				break;

			lg_addTailData((void*) currBatch, threadExecHistory);
			JOB_BATCH_EXECUTE_HYBRID_JDS_CSR(currBatch, x, y_inout);

			#if PROGRAM_MODE >= TRACE_MODE
				++jobCountPerThread[threadId];
			#endif

		} while(TRUE);
	}

	#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Printing job count executed by each thread;\n");
		for(i = 0; i < numThreads; ++i)
			if(jobCountPerThread[i] > 0)
				PRINTF("Thread-%d: %d jobs.\n", i, jobCountPerThread[i]);
		free(jobCountPerThread);

		PRINTF("\n");
	#endif

	// return values
	args->execHistoryPerThread = execHistoryPerThread;
}

// Intermediate data preparation functions
// ---------------------------------------------------------------------------------------------------------------------

void spmxv_gws_prepFromSubMtxList_colNet_CSR(gws_args_t* args)
{
	DEBUG("spmxv_gws_prepFromSubMtxList_colNet_CSR\n");

	spm_cmp_t* spmCsr = args->staticArgs.spmCsr;
	lg_t* subMtxList = args->staticArgs.subMtxList;
	int stealTreshold = 0; // not needed for GWS routine

	lg_t* jobBatchList = NULL;
	job_batch_extractJobBatchListFromSubMtxList_rowWiseColNet_CSR(subMtxList, &jobBatchList);

	job_queue_t* globalQueue = job_queue_new(stealTreshold);
	job_queue_fillFromJobBatchList(globalQueue, jobBatchList);

	// clean up
	lg_deleteShallow(jobBatchList);

	// return values
	args->globalQueue = globalQueue;
}

void spmxv_gws_prepFromSubMtxTree_colNet_CSR(gws_args_t* args)
{
	DEBUG("spmxv_gws_prepFromSubMtxTree_colNet_CSR\n");

	spm_cmp_t* spmCsr = args->staticArgs.spmCsr;
	sub_mtx_tree_t* head = args->staticArgs.subMtxHead;
	int stealTreshold = 0; // not needed for GWS routine

	job_queue_t* globalQueue = job_queue_new(stealTreshold);
	job_queue_fillFromSubMtxTree_rowWiseColNet_CSR(spmCsr, globalQueue, &head->node);

	// return values
	args->globalQueue = globalQueue;
}

void spmxv_gws_prepFromPartitionVector_colNet_CSR(gws_args_t* args)
{
	DEBUG("spmxv_gws_prepFromPartitionVector_colNet_CSR\n");

	spm_cmp_t* spmCsr = args->staticArgs.spmCsr;
	sub_mtx_tree_t* head = args->staticArgs.subMtxHead;
	partitioning_vector_t* rowPartitioning = args->staticArgs.rowPartitioning;
	int partitionCount = rowPartitioning->length - 1;
	int stealTreshold = 0; // not needed for GWS routine

	job_queue_t* globalQueue = job_queue_new(0);
	job_queue_fillFromParitioningVector_rowWiseColNet_CSR(
			spmCsr, globalQueue, rowPartitioning, 0, partitionCount);

	// return values
	args->globalQueue = globalQueue;
}

// ---------------------------------------------------------------------------------------------------------------------

void spmxv_gws_prepFromSubMtxList_colNet_JDS(gws_args_t* args)
{
	DEBUG("spmxv_gws_prepFromSubMtxList_colNet_JDS\n");

	spmxv_gws_prepFromSubMtxList_colNet_CSR(args);

	spm_cmp_t* spmCsrCounterpart = args->staticArgs.spmCsr;
	DECIMAL* permutation = args->staticArgs.permutation;
	job_queue_t* globalQueue = args->globalQueue;

	job_queue_convertToJDSPartial(spmCsrCounterpart, permutation, globalQueue);

	// return values
	args->globalQueue = globalQueue;
}

void spmxv_gws_prepFromSubMtxTree_colNet_JDS(gws_args_t* args)
{
	DEBUG("spmxv_gws_prepFromSubMtxTree_colNet_JDS\n");

	spmxv_gws_prepFromSubMtxTree_colNet_CSR(args);

	spm_cmp_t* spmCsrCounterpart = args->staticArgs.spmCsr;
	DECIMAL* permutation = args->staticArgs.permutation;
	job_queue_t* globalQueue = args->globalQueue;

	job_queue_convertToJDSPartial(spmCsrCounterpart, permutation, globalQueue);

	// return values
	args->globalQueue = globalQueue;
}

void spmxv_gws_prepFromPartitionVector_colNet_JDS(gws_args_t* args)
{
	DEBUG("spmxv_gws_prepFromPartitionVector_colNet_JDS\n");

	spmxv_gws_prepFromPartitionVector_colNet_CSR(args);

	spm_cmp_t* spmCsrCounterpart = args->staticArgs.spmCsr;
	DECIMAL* permutation = args->staticArgs.permutation;
	job_queue_t* globalQueue = args->globalQueue;

	job_queue_convertToJDSPartial(spmCsrCounterpart, permutation, globalQueue);

	// return values
	args->globalQueue = globalQueue;
}

// ---------------------------------------------------------------------------------------------------------------------

void spmxv_gws_prepFromSubMtxList_colNet_hybrid_JDS_CSR(gws_args_t* args)
{
	DEBUG("spmxv_gws_prepFromSubMtxList_colNet_hybrid_JDS_CSR\n");

	spmxv_gws_prepFromSubMtxList_colNet_JDS(args);

	job_queue_t* partialJdsQueue = args->globalQueue;
	job_queue_t* globalQueue = job_queue_new(0);
	job_queue_convertFromPartialJDSToHybridJDSCSR(partialJdsQueue, globalQueue);

	// clean up
	job_queue_delete(partialJdsQueue);

	// return values
	args->globalQueue = globalQueue;
}


void spmxv_gws_prepFromSubMtxTree_colNet_hybrid_JDS_CSR(gws_args_t* args)
{
	DEBUG("spmxv_gws_prepFromSubMtxTree_colNet_hybrid_JDS_CSR\n");

	spmxv_gws_prepFromSubMtxTree_colNet_JDS(args);

	job_queue_t* partialJdsQueue = args->globalQueue;
	job_queue_t* globalQueue = job_queue_new(0);
	job_queue_convertFromPartialJDSToHybridJDSCSR(partialJdsQueue, globalQueue);

	// clean up
	job_queue_delete(partialJdsQueue);

	// return values
	args->globalQueue = globalQueue;
}

void spmxv_gws_prepFromPartitionVector_colNet_hybrid_JDS_CSR(gws_args_t* args)
{
	DEBUG("spmxv_gws_prepFromPartitionVector_colNet_hybrid_JDS_CSR\n");

	spmxv_gws_prepFromPartitionVector_colNet_JDS(args);

	job_queue_t* partialJdsQueue = args->globalQueue;
	job_queue_t* globalQueue = job_queue_new(partialJdsQueue->stealTreshold);
	job_queue_convertFromPartialJDSToHybridJDSCSR(partialJdsQueue, globalQueue);

	// clean up
	job_queue_delete(partialJdsQueue);

	// return values
	args->globalQueue = globalQueue;
}
