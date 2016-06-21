
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/io/converter.h"
#include "include/util/utility.h"
#include "include/data_structure/sub_mtx.h"
#include "include/task_decomposition/partitioning.h"

#include "arch/generic/inc/spmxv_static.h"

// Row Parallel
// ----------------------------------------------------------------------------------------------------------------

void spmxv_static_rowWise_CSR(
		static_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	int numThreads = args->numThreads;
	spm_cmp_t* spmCsr = args->spmCsr;
	job_batch_t** jobBatchArrPerThread = args->jobBatchArrPerThread;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;


#if PROGRAM >= TRACE_MODE
	PRINTF("\n\n");
#endif

	#pragma omp parallel num_threads(numThreads)
	{
		int threadId = omp_get_thread_num();
		job_batch_t* myJobBatchArr = jobBatchArrPerThread[threadId];
		int length = jobBatchCountPerThread[threadId];

#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Thread-%d batch-count: %d\n", threadId, length);
#endif

		int i;
		for(i = 0; i < length; ++i)
		{

#if PROGRAM_MODE >= TRACE_MODE
			char jobBatchStr[DEFAULT_STR_BUFF_SIZE] = "\0";
			job_batch_toStringDetailed(&myJobBatchArr[i], (void*) spmCsr, jobBatchStr);
			PRINTF("Thread-%d executing %s\n", threadId, jobBatchStr);
#endif

			job_batch_t* currBatch = &myJobBatchArr[i];
			JOB_BATCH_EXECUTE_SUBMTX_CSR(currBatch, spmCsr, x, y_inout);
		}
	}
}

void spmxv_static_rowWise_JDS(
		static_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	int numThreads = args->numThreads;
	job_batch_t** jobBatchArrPerThread = args->jobBatchArrPerThread;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("\n\n");
#endif

	#pragma omp parallel num_threads(numThreads)
	{
		int threadId = omp_get_thread_num();
		job_batch_t* myJobBatchArr = jobBatchArrPerThread[threadId];
		int length = jobBatchCountPerThread[threadId];

#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Thread-%d batch-count: %d\n", threadId, length);
#endif

		int i;
		for(i = 0; i < length; ++i)
		{

#if PROGRAM_MODE >= TRACE_MODE
			char jobBatchStr[DEFAULT_STR_BUFF_SIZE] = "\0";
			job_batch_toStringDetailed(&myJobBatchArr[i], NULL, jobBatchStr);
			PRINTF("Thread-%d executing %s\n", threadId, jobBatchStr);
#endif

			job_batch_t* currBatch = &myJobBatchArr[i];
			JOB_BATCH_EXECUTE_PARTIAL_JDS(currBatch, x, y_inout);
		}
	}
}

void spmxv_static_rowWise_hybrid_JDS_CSR(
		static_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	int numThreads = args->numThreads;
	job_batch_t** jobBatchArrPerThread = args->jobBatchArrPerThread;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("\n\n");
#endif

	#pragma omp parallel num_threads(numThreads)
	{
		int threadId = omp_get_thread_num();
		int length = jobBatchCountPerThread[threadId];

		job_batch_t* myJobBatchArr = jobBatchArrPerThread[threadId];


#if PROGRAM_MODE >= TRACE_MODE
		PRINTF("Thread-%d batch-count: %d\n", threadId, length);
#endif

		int i;
		for(i = 0; i < length; ++i)
		{

#if PROGRAM_MODE >= TRACE_MODE
			char jobBatchStr[DEFAULT_STR_BUFF_SIZE] = "\0";
			job_batch_toStringDetailed(&myJobBatchArr[i], NULL, jobBatchStr);
			PRINTF("Thread-%d executing %s\n", threadId, jobBatchStr);
#endif

			job_batch_t* currBatch = &myJobBatchArr[i];
			JOB_BATCH_EXECUTE_HYBRID_JDS_CSR(currBatch, x, y_inout);
		}
	}
}

// Pre-multiplication routines / data preparations
// ----------------------------------------------------------------------------------------------------------------

void spmxv_static_prepFromSubMtxTree_rowWise_CSR(static_args_t* args)
{
	DEBUG("Prep: spmxv_static_prepFromSubMtxTree_rowWise_CSR\n");

	int numThreads = args->numThreads;
	spm_cmp_t* spmCsr = args->spmCsr;
	sub_mtx_tree_t* subMtxHead = args->subMtxHead;

	job_batch_t** jobBatchArrPerThread = (job_batch_t**) malloc(sizeof(job_batch_t*) * numThreads);
	int* jobBatchCountPerThread = (int*) malloc(sizeof(int) * numThreads);

	tree_node_t** treeNodePerThread = sub_mtx_tree_getJobDistribution(subMtxHead, numThreads);

	#pragma omp parallel num_threads(numThreads)
	{
		int blockId = omp_get_thread_num();
		tree_node_t* threadTreeNode = treeNodePerThread[blockId];
		lg_t* threadJobBatchList = NULL;
		job_batch_t* threadBatchArr = NULL;
		int blockBatchCount = 0;

		job_batch_extractJobBatchListFromSubMtxTree_rowWiseColNet_CSR(spmCsr, threadTreeNode, &threadJobBatchList);

		// ASSUMPTION job-batches used in this function are also created in this function.
		// So there is no shallow copy. We can do whatever we want with batch-lists, no need to copy.
		// merge batches if possible
		// TODO implement merging and come back here
		// if(mergeJobBatches == TRUE)
		//	batch_list_mergeBatchList(&blockBatchList);

		lg_job_batch_toArray(threadJobBatchList, &threadBatchArr, &blockBatchCount);

		jobBatchArrPerThread[blockId] = threadBatchArr;
		jobBatchCountPerThread[blockId] = blockBatchCount;

		// thread clean up
		lg_job_batch_deleteDeep(threadJobBatchList);
	}

	// clean up
	free(treeNodePerThread);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
	args->jobBatchCountPerThread = jobBatchCountPerThread;
}

void spmxv_static_prepFromSubMtxTreeScatter_rowWise_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromSubMtxTreeScatter_rowWise_CSR\n");

	spm_cmp_t* spmCsr = args->spmCsr;
	sub_mtx_tree_t* subMtxHead = args->subMtxHead;
	int numBlocks = args->numBlocks;
	int threadsPerBlock = args->threadsPerBlock;
	int numThreads = args->numThreads;

	job_batch_t** jobBatchArrPerThread = (job_batch_t**) malloc(sizeof(job_batch_t*) * numThreads);
	int* jobBatchCountPerThread = (int*) calloc(numThreads, sizeof(int));

	tree_node_t** treeNodePerBlock = sub_mtx_tree_getJobDistribution(subMtxHead, numBlocks);

	#pragma omp parallel num_threads(numBlocks)
	{
		int blockId = omp_get_thread_num();
		tree_node_t* blockTreeNode = treeNodePerBlock[blockId];

		lg_t* blockJobBatchList = NULL;
		job_batch_extractJobBatchListFromSubMtxTree_rowWiseColNet_CSR(spmCsr, blockTreeNode, &blockJobBatchList);

		// scatter block job batch list for each thread and delete empty block list
		lg_t** scatteredJobBatchList = NULL;
		job_batch_scatterJobBatchList(blockJobBatchList, threadsPerBlock, &scatteredJobBatchList);
		lg_deleteShallow(blockJobBatchList);

		// Convert job batch list to dense array structure for each thread in this block
		for(int i = 0; i < threadsPerBlock; ++i)
		{
			lg_job_batch_toArray(scatteredJobBatchList[i],
					&jobBatchArrPerThread[threadsPerBlock * blockId + i],
					&jobBatchCountPerThread[threadsPerBlock * blockId + i]);
		}

		// delete job batch list per thread
		lg_job_batch_deleteDeepMultiple(scatteredJobBatchList, threadsPerBlock);
	}

	// clean up
	free(treeNodePerBlock);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
	args->jobBatchCountPerThread = jobBatchCountPerThread;
}

void spmxv_static_prepFromSubMtxList_rowWise_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromSubMtxList_rowWise_CSR\n");

	spm_cmp_t* spmCsr = args->spmCsr;
	lg_t* subMtxList = args->subMtxList;
	int numThreads = args->numThreads;

	job_batch_t** jobBatchArrPerThread = (job_batch_t**) malloc(sizeof(job_batch_t*) * numThreads);
	int* jobBatchCountPerThread = (int*) malloc(sizeof(int) * numThreads);

	lg_t** subMtxListPerBlock = NULL;
	lg_sub_mtx_decomposeDeep(subMtxList, numThreads, &subMtxListPerBlock);

	#pragma omp parallel num_threads(numThreads)
	{
		int blockId = omp_get_thread_num();

		jobBatchArrPerThread[blockId] = NULL;
		jobBatchCountPerThread[blockId] = 0;

		lg_t* jobBatchList = NULL;
		job_batch_extractJobBatchListFromSubMtxList_rowWiseColNet_CSR(
				subMtxListPerBlock[blockId], &jobBatchList);

		lg_job_batch_toArray(
				jobBatchList,
				&jobBatchArrPerThread[blockId],
				&jobBatchCountPerThread[blockId]);

		// local clean up
		lg_sub_mtx_deleteDeep(subMtxListPerBlock[blockId]);
		lg_job_batch_deleteDeep(jobBatchList);
	}

	// clean up
	free(subMtxListPerBlock);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
	args->jobBatchCountPerThread = jobBatchCountPerThread;
}

void spmxv_static_prepFromSubMtxListScatter_rowWise_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromSubMtxListScatter_rowWise_CSR\n");

	spm_cmp_t* spmCsr = args->spmCsr;
	lg_t* subMtxList = args->subMtxList;
	int numBlocks = args->numBlocks;
	int threadsPerBlock = args->threadsPerBlock;
	int numThreads = args->numThreads;

	job_batch_t** jobBatchArrPerThread = (job_batch_t**) malloc(sizeof(job_batch_t*) * numThreads);
	int* jobBatchCountPerThread = (int*) malloc(sizeof(int) * numThreads);

	lg_t** subMtxListPerBlock = NULL;
	lg_sub_mtx_scatterDeep(subMtxList, numThreads, &subMtxListPerBlock);

	#pragma omp parallel num_threads(numThreads)
	{
		int blockId = omp_get_thread_num();

		jobBatchArrPerThread[blockId] = NULL;
		jobBatchCountPerThread[blockId] = 0;

		lg_t* jobBatchList = NULL;
		job_batch_extractJobBatchListFromSubMtxList_rowWiseColNet_CSR(
				subMtxListPerBlock[blockId], &jobBatchList);

		lg_job_batch_toArray(
				jobBatchList,
				&jobBatchArrPerThread[blockId],
				&jobBatchCountPerThread[blockId]);

		// local clean up
		lg_sub_mtx_deleteDeep(subMtxListPerBlock[blockId]);
		lg_job_batch_deleteDeep(jobBatchList);
	}

	// clean up
	free(subMtxListPerBlock);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
	args->jobBatchCountPerThread = jobBatchCountPerThread;
}

void spmxv_static_prepFromPartitionVector_rowWise_colNet_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromPartitionVector_rowWise_colNet_CSR\n");

	spm_cmp_t* spmCsr = args->spmCsr;
	int numThreads = args->numThreads;
	partitioning_vector_t* rowPartitioning = args->rowPartitioning;
	int partitionCount = rowPartitioning->length - 1;

	job_batch_t** jobBatchArrPerThread = (job_batch_t**) malloc(sizeof(job_batch_t*) * numThreads);
	int* jobBatchCountPerThread = (int*) malloc(sizeof(int) * numThreads);
	int* jobBatchCountAcc = (int*) malloc(sizeof(int) * (numThreads + 1));
	jobBatchCountAcc[0] = 0;

	// calculate minimum batch count per block
	int minBatchCount = partitionCount / numThreads;
	int remainingBatchCount = partitionCount % numThreads;

	#pragma omp parallel num_threads(numThreads)
	{
		int threadId = omp_get_thread_num();

		// calculate job batch count for the current block
		jobBatchCountPerThread[threadId] = minBatchCount;
		if(remainingBatchCount > threadId)
			++jobBatchCountPerThread[threadId];

		jobBatchCountAcc[threadId + 1] = jobBatchCountPerThread[threadId];

		#pragma omp barrier
		#pragma omp single
		{
			accumulate(jobBatchCountAcc, numThreads + 1);
		}
		#pragma omp barrier

		// initialize for healthy memory clean up
		jobBatchArrPerThread[threadId] = NULL;

		// extract job batches from partitioning vector in the form of lists
		lg_t* threadJobBatchList = NULL;
		int threadBatchStartIndex = jobBatchCountAcc[threadId];
		int threadBatchEndIndex = jobBatchCountAcc[threadId + 1];
		job_batch_extractJobBatchListFromPartitionVector_rowWiseColNet_CSR(
				spmCsr, rowPartitioning, threadBatchStartIndex, threadBatchEndIndex, &threadJobBatchList);

		// ASSUMPTION job-batches used in this function are also created in this function.
		// So there is no shallow copy. We can do whatever we want with batch-lists, no need to copy.
		// merge batches if possible
		// TODO implement merge and come back here
		// if(mergeJobBatches == TRUE)
		//	batch_list_mergeBatchList(&blockBatchList);

		// convert to dense job batch array for better performance
		lg_job_batch_toArray(threadJobBatchList, &jobBatchArrPerThread[threadId], &jobBatchCountPerThread[threadId]);

		// partitioning vector start index for next block
		threadBatchStartIndex = threadBatchEndIndex;

		// per thread local clean up
		lg_job_batch_deleteDeep(threadJobBatchList);
	}

	// clean up
	free(jobBatchCountAcc);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
	args->jobBatchCountPerThread = jobBatchCountPerThread;
}

void spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_CSR\n");

	spm_cmp_t* spmCsr = args->spmCsr;
	int numBlocks = args->numBlocks;
	int threadsPerBlock = args->threadsPerBlock;
	int numThreads = args->numThreads;
	partitioning_vector_t* rowPartitioning = args->rowPartitioning;
	int partitionCount = args->rowPartitioning->length - 1;

	job_batch_t** jobBatchArrPerThread = (job_batch_t**) malloc(sizeof(job_batch_t*) * numThreads);
	int* jobBatchCountPerThread = (int*) malloc(sizeof(int) * numThreads);
	int* jobBatchCountAcc = (int*) malloc(sizeof(int) * (numBlocks + 1));
	jobBatchCountAcc[0] = 0;

	// initialize for healthy memory clean up
	int i;
	for(i = 0; i < numThreads; ++i)
	{
		jobBatchArrPerThread[i] = NULL;
		jobBatchCountPerThread[i] = 0;
	}

	// calculate minimum batch count per block
	int minBatchCountPerBlock = partitionCount / numBlocks;
	int remainingBatchCount = partitionCount % numBlocks;

	#pragma omp parallel num_threads(numBlocks)
	{
		int blockId = omp_get_thread_num();

		// calculate batch count per block
		int blockJobBatchCount = minBatchCountPerBlock;
		if(remainingBatchCount > blockId)
			++blockJobBatchCount;

		jobBatchCountAcc[blockId + 1] = blockJobBatchCount;

		#pragma omp barrier
		#pragma omp single
		{
			accumulate(jobBatchCountAcc, numBlocks + 1);
		}
		#pragma omp barrier


		// extract job batches from partitioning vector in the form of lists
		int blockBatchStartIndex = jobBatchCountAcc[blockId];
		int blockBatchEndIndex = jobBatchCountAcc[blockId + 1];
		lg_t* blockJobBatchList = NULL;
		job_batch_extractJobBatchListFromPartitionVector_rowWiseColNet_CSR(
				spmCsr, rowPartitioning, blockBatchStartIndex, blockBatchEndIndex, &blockJobBatchList);

		// job batches are distributed between each thread in a single block
		// first divide block list to threadsPerBlock number of lists and fill those in scatter fashion
		// convert job batch to dense format in a scatter fashion
		lg_t** jobBatchListPerThread = NULL;
		job_batch_scatterJobBatchList(blockJobBatchList, threadsPerBlock, &jobBatchListPerThread);

		// convert each thread's workload to dense format
		int j;
		for(j = 0; j < threadsPerBlock; ++j)
		{
			lg_job_batch_toArray(jobBatchListPerThread[j],
					&jobBatchArrPerThread[blockId * threadsPerBlock + j],
					&jobBatchCountPerThread[blockId * threadsPerBlock + j]);
		}

		// block local clean up
		lg_job_batch_deleteDeepMultiple(jobBatchListPerThread, threadsPerBlock);
		lg_deleteShallow(blockJobBatchList);
	}

	// clean up
	free(jobBatchCountAcc);


	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
	args->jobBatchCountPerThread = jobBatchCountPerThread;
}

// TODO test
void spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_CSR\n");

	spm_cmp_t* spmCsr = args->spmCsr;
	int numThreads = args->numThreads;
	partitioning_vector_t* rowPartitioning = args->rowPartitioning;
	int partitionCount = args->rowPartitioning->length - 1;

	job_batch_t** jobBatchArrPerThread = (job_batch_t**) malloc(sizeof(job_batch_t*) * numThreads);
	int* jobBatchCountPerThread = (int*) malloc(sizeof(int) * numThreads);

	// calculate job count per thread
	DECIMAL* subMtxCountPerThread = partitioning_calculateSubMtxCountPerBlockPowerOf2(numThreads, rowPartitioning->length);

	#pragma omp parallel num_threads(numThreads)
	{
		int threadId = omp_get_thread_num();
		jobBatchCountPerThread[threadId] = subMtxCountPerThread[threadId];

		#pragma omp barrier
		#pragma omp single
		{
			accumulate(subMtxCountPerThread, numThreads);
		}
		#pragma omp barrier

		// initialize for healthy memory clean up afterwards
		jobBatchArrPerThread[threadId] = NULL;

		// calculate starting & ending indices from accumulated array
		int jobBatchStartIndex = (threadId == 0) ? 0 : subMtxCountPerThread[threadId - 1];
		int jobBatchEndIndex = subMtxCountPerThread[threadId];

		// extract job batches from partitioning vector in the form of lists
		lg_t* threadJobBatchList = NULL;
		job_batch_extractJobBatchListFromPartitionVector_rowWiseColNet_CSR(
				spmCsr, rowPartitioning, jobBatchStartIndex, jobBatchEndIndex, &threadJobBatchList);

		// ASSUMPTION job-batches used in this function are also created in this function.
		// So there is no shallow copy. We can do whatever we want with batch-lists, no need to copy.
		// merge batches if possible
		// TODO implement merge and come back here
		// if(mergeJobBatches == TRUE)
		//	batch_list_mergeBatchList(&blockBatchList);

		// convert list to dense job batch array for better performance
		lg_job_batch_toArray(threadJobBatchList, &jobBatchArrPerThread[threadId], &jobBatchCountPerThread[threadId]);

		// per thread local clean up: delete list
		lg_job_batch_deleteDeep(threadJobBatchList);
	}

	// clean up
	free(subMtxCountPerThread);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
	args->jobBatchCountPerThread = jobBatchCountPerThread;
}

// TODO test
void spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_CSR\n");

	spm_cmp_t* spmCsr = args->spmCsr;
	int numBlocks = args->numBlocks;
	int threadsPerBlock = args->threadsPerBlock;
	int numThreads = args->numThreads;
	partitioning_vector_t* rowPartitioning = args->rowPartitioning;
	int partitionCount = rowPartitioning->length - 1;

	job_batch_t** jobBatchArrPerThread = (job_batch_t**) malloc(sizeof(job_batch_t*) * numThreads);
	int* jobBatchCountPerThread = (int*) malloc(sizeof(int) * numThreads);

	// initialize for healthy memory clean up
	int i;
	for(i = 0; i < numThreads; ++i)
	{
		jobBatchArrPerThread[i] = NULL;
		jobBatchCountPerThread[i] = 0;
	}

	// calculate job count per block
	int* subMtxCountPerBlock = partitioning_calculateSubMtxCountPerBlockPowerOf2(numBlocks, partitionCount);

	#pragma omp parallel num_threads(numBlocks)
	{
		int blockId = omp_get_thread_num();

		// calculate batch count per block
		int blockJobBatchCount = subMtxCountPerBlock[blockId];

		#pragma omp barrier
		#pragma omp single
		{
			accumulate(subMtxCountPerBlock, numBlocks);
		}
		#pragma omp barrier


		// extract job batches from partitioning vector in the form of lists
		int blockBatchStartIndex = (blockId == 0) ? 0: subMtxCountPerBlock[blockId - 1];
		int blockBatchEndIndex = subMtxCountPerBlock[blockId];
		lg_t* blockJobBatchList = NULL;
		job_batch_extractJobBatchListFromPartitionVector_rowWiseColNet_CSR(
				spmCsr, rowPartitioning, blockBatchStartIndex, blockBatchEndIndex, &blockJobBatchList);

		// job batches are distributed between each thread in a single block
		// first divide block list to threadsPerBlock number of lists and fill those in scatter fashion
		// convert job batch to dense format in a scatter fashion
		lg_t** jobBatchListPerThread = NULL;
		job_batch_scatterJobBatchList(blockJobBatchList, threadsPerBlock, &jobBatchListPerThread);

		// convert each thread's workload to dense format
		int j;
		for(j = 0; j < threadsPerBlock; ++j)
		{
			lg_job_batch_toArray(jobBatchListPerThread[j],
					&jobBatchArrPerThread[blockId * threadsPerBlock + j],
					&jobBatchCountPerThread[blockId * threadsPerBlock + j]);
		}

		// block local clean up
		lg_job_batch_deleteDeepMultiple(jobBatchListPerThread, threadsPerBlock);
		lg_deleteShallow(blockJobBatchList);
	}

	// clean up
	free(subMtxCountPerBlock);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
	args->jobBatchCountPerThread = jobBatchCountPerThread;
}


// ----------------------------------------------------------------------------------------------------------------

void spmxv_static_prepFromSubMtxTree_rowWise_JDS(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromSubMtxTree_rowWise_JDS\n");

	spmxv_static_prepFromSubMtxTree_rowWise_CSR(args);

	spm_cmp_t* spmCsrCounterpart = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	job_batch_t** subMtxArrPerThread = args->jobBatchArrPerThread;
	job_batch_t** jobBatchArrPerThread = NULL;
	int numThreads = args->numThreads;

	converter_extractPartialJDSArrMultiple(
			spmCsrCounterpart, numThreads, subMtxArrPerThread,
			jobBatchCountPerThread, permutation, &jobBatchArrPerThread);

	// local cleanup
	job_batch_deleteDenseArrMultiple(subMtxArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromSubMtxTreeScatter_rowWise_JDS(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromSubMtxTreeScatter_rowWise_JDS\n");

	spmxv_static_prepFromSubMtxTreeScatter_rowWise_CSR(args);

	spm_cmp_t* spmCsrCounterpart = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	job_batch_t** subMtxArrPerThread = args->jobBatchArrPerThread;
	job_batch_t** jobBatchArrPerThread = NULL;
	int numThreads = args->numThreads;

	converter_extractPartialJDSArrMultiple(
			spmCsrCounterpart, numThreads, subMtxArrPerThread,
			jobBatchCountPerThread, permutation, &jobBatchArrPerThread);

	// local cleanup
	job_batch_deleteDenseArrMultiple(subMtxArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromSubMtxList_rowWise_JDS(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromSubMtxList_rowWise_JDS\n");

	spmxv_static_prepFromSubMtxList_rowWise_CSR(args);

	spm_cmp_t* spmCsrCounterpart = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	job_batch_t** subMtxArrPerThread = args->jobBatchArrPerThread;
	job_batch_t** jobBatchArrPerThread = NULL;
	int numThreads = args->numThreads;


	converter_extractPartialJDSArrMultiple(
			spmCsrCounterpart, numThreads, subMtxArrPerThread,
			jobBatchCountPerThread, permutation, &jobBatchArrPerThread);

	// local cleanup
	job_batch_deleteDenseArrMultiple(subMtxArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromSubMtxListScatter_rowWise_JDS(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromSubMtxListScatter_rowWise_JDS\n");

	spmxv_static_prepFromSubMtxListScatter_rowWise_CSR(args);

	spm_cmp_t* spmCsrCounterpart = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	job_batch_t** subMtxArrPerThread = args->jobBatchArrPerThread;
	job_batch_t** jobBatchArrPerThread = NULL;
	int numThreads = args->numThreads;


	converter_extractPartialJDSArrMultiple(
			spmCsrCounterpart, numThreads, subMtxArrPerThread,
			jobBatchCountPerThread, permutation, &jobBatchArrPerThread);

	// local cleanup
	job_batch_deleteDenseArrMultiple(subMtxArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromPartitionVector_rowWise_colNet_JDS(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromPartitionVector_rowWise_colNet_JDS\n");

	spmxv_static_prepFromPartitionVector_rowWise_colNet_CSR(args);

	spm_cmp_t* spmCsrCounterpart = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	job_batch_t** subMtxArrPerThread = args->jobBatchArrPerThread;
	job_batch_t** jobBatchArrPerThread = NULL;
	int numThreads = args->numThreads;

	converter_extractPartialJDSArrMultiple(
			spmCsrCounterpart, numThreads, subMtxArrPerThread,
			jobBatchCountPerThread, permutation, &jobBatchArrPerThread);

	// local cleanup
	job_batch_deleteDenseArrMultiple(subMtxArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_JDS(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_JDS\n");

	spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_CSR(args);

	spm_cmp_t* spmCsrCounterpart = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	job_batch_t** subMtxArrPerThread = args->jobBatchArrPerThread;
	job_batch_t** jobBatchArrPerThread = NULL;
	int numThreads = args->numThreads;

	converter_extractPartialJDSArrMultiple(
			spmCsrCounterpart, numThreads, subMtxArrPerThread,
			jobBatchCountPerThread, permutation, &jobBatchArrPerThread);

	// local cleanup
	job_batch_deleteDenseArrMultiple(subMtxArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_JDS(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_JDS\n");

	spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_CSR(args);

	spm_cmp_t* spmCsrCounterpart = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	job_batch_t** subMtxArrPerThread = args->jobBatchArrPerThread;
	job_batch_t** jobBatchArrPerThread = NULL;
	int numThreads = args->numThreads;

	converter_extractPartialJDSArrMultiple(
			spmCsrCounterpart, numThreads, subMtxArrPerThread,
			jobBatchCountPerThread, permutation, &jobBatchArrPerThread);

	// local cleanup
	job_batch_deleteDenseArrMultiple(subMtxArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_JDS(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_JDS\n");

	spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_CSR(args);

	spm_cmp_t* spmCsrCounterpart = args->spmCsr;
	DECIMAL* permutation = args->permutation;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	job_batch_t** subMtxArrPerThread = args->jobBatchArrPerThread;
	job_batch_t** jobBatchArrPerThread = NULL;
	int numThreads = args->numThreads;

	converter_extractPartialJDSArrMultiple(
			spmCsrCounterpart, numThreads, subMtxArrPerThread,
			jobBatchCountPerThread, permutation, &jobBatchArrPerThread);

	// local cleanup
	job_batch_deleteDenseArrMultiple(subMtxArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

// ----------------------------------------------------------------------------------------------------------------

void spmxv_static_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR\n");

	spmxv_static_prepFromSubMtxTree_rowWise_JDS(args);

	job_batch_t** jobBatchArrPerThread = NULL;
	job_batch_t** partialJdsBatchArrPerThread = args->jobBatchArrPerThread;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	int numThreads = args->numThreads;

	converter_extractHybridJDSCSRArrMultiple(
			numThreads, partialJdsBatchArrPerThread,
			jobBatchCountPerThread, &jobBatchArrPerThread);

	// clean up
	job_batch_deleteDenseArrMultiple(
			partialJdsBatchArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromSubMtxTreeScatter_rowWise_hybrid_JDS_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromSubMtxTreeScatter_rowWise_hybrid_JDS_CSR\n");

	spmxv_static_prepFromSubMtxTreeScatter_rowWise_JDS(args);

	job_batch_t** jobBatchArrPerThread = NULL;
	job_batch_t** partialJdsBatchArrPerThread = args->jobBatchArrPerThread;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	int numThreads = args->numThreads;

	converter_extractHybridJDSCSRArrMultiple(
			numThreads, partialJdsBatchArrPerThread,
			jobBatchCountPerThread, &jobBatchArrPerThread);

	// clean up
	job_batch_deleteDenseArrMultiple(
			partialJdsBatchArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromSubMtxList_rowWise_hybrid_JDS_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromSubMtxList_rowWise_hybrid_JDS_CSR\n");

	spmxv_static_prepFromSubMtxList_rowWise_JDS(args);

	job_batch_t** jobBatchArrPerThread = NULL;
	job_batch_t** partialJdsBatchArrPerThread = args->jobBatchArrPerThread;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	int numThreads = args->numThreads;

	converter_extractHybridJDSCSRArrMultiple(
			numThreads, partialJdsBatchArrPerThread,
			jobBatchCountPerThread, &jobBatchArrPerThread);

	// clean up
	job_batch_deleteDenseArrMultiple(
			partialJdsBatchArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromSubMtxListScatter_rowWise_hybrid_JDS_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromSubMtxListScatter_rowWise_hybrid_JDS_CSR\n");

	spmxv_static_prepFromSubMtxListScatter_rowWise_JDS(args);

	job_batch_t** jobBatchArrPerThread = NULL;
	job_batch_t** partialJdsBatchArrPerThread = args->jobBatchArrPerThread;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	int numThreads = args->numThreads;

	converter_extractHybridJDSCSRArrMultiple(
			numThreads, partialJdsBatchArrPerThread,
			jobBatchCountPerThread, &jobBatchArrPerThread);

	// clean up
	job_batch_deleteDenseArrMultiple(
			partialJdsBatchArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR\n");

	spmxv_static_prepFromPartitionVector_rowWise_colNet_JDS(args);

	job_batch_t** jobBatchArrPerThread = NULL;
	job_batch_t** partialJdsBatchArrPerThread = args->jobBatchArrPerThread;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	int numThreads = args->numThreads;

	converter_extractHybridJDSCSRArrMultiple(
			numThreads, partialJdsBatchArrPerThread,
			jobBatchCountPerThread, &jobBatchArrPerThread);

	// clean up
	job_batch_deleteDenseArrMultiple(
			partialJdsBatchArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_hybrid_JDS_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_hybrid_JDS_CSR\n");

	spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_JDS(args);

	job_batch_t** jobBatchArrPerThread = NULL;
	job_batch_t** partialJdsBatchArrPerThread = args->jobBatchArrPerThread;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	int numThreads = args->numThreads;

	converter_extractHybridJDSCSRArrMultiple(
			numThreads, partialJdsBatchArrPerThread,
			jobBatchCountPerThread, &jobBatchArrPerThread);

	// clean up
	job_batch_deleteDenseArrMultiple(
			partialJdsBatchArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_hybrid_JDS_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_hybrid_JDS_CSR\n");

	spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_JDS(args);

	job_batch_t** jobBatchArrPerThread = NULL;
	job_batch_t** partialJdsBatchArrPerThread = args->jobBatchArrPerThread;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	int numThreads = args->numThreads;

	converter_extractHybridJDSCSRArrMultiple(
			numThreads, partialJdsBatchArrPerThread,
			jobBatchCountPerThread, &jobBatchArrPerThread);

	// clean up
	job_batch_deleteDenseArrMultiple(
			partialJdsBatchArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

void spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_hybrid_JDS_CSR(static_args_t* args)
{
	DEBUG("spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_hybrid_JDS_CSR\n");

	spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_JDS(args);

	job_batch_t** jobBatchArrPerThread = NULL;
	job_batch_t** partialJdsBatchArrPerThread = args->jobBatchArrPerThread;
	int* jobBatchCountPerThread = args->jobBatchCountPerThread;
	int numThreads = args->numThreads;

	converter_extractHybridJDSCSRArrMultiple(
			numThreads, partialJdsBatchArrPerThread,
			jobBatchCountPerThread, &jobBatchArrPerThread);

	// clean up
	job_batch_deleteDenseArrMultiple(
			partialJdsBatchArrPerThread, jobBatchCountPerThread, numThreads);

	// return values
	args->jobBatchArrPerThread = jobBatchArrPerThread;
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
