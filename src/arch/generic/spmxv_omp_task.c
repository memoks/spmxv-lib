
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "inc/spmxv_sequential.h"
#include "inc/spmxv_static.h"
#include "include/data_structure/sub_mtx.h"

#include "inc/spmxv_omp_task.h"

// Core
// ---------------------------------------------------------------------------------------------------------------------

// Row Parallel
// ---------------------------------------------------------------------------------------------------------------------

void spmxv_omp_task_multiply_CSR(
		spm_cmp_t* spmCsr, vector_real_t* x, vector_real_t* y_inout,
		job_batch_t** batchArrPerThread, int* batchCountPerThread, int numThreads,
		lg_t*** execHistoryPerThread_out)
{
	lg_t** execHistoryPerThread = (lg_t**) malloc(sizeof(lg_t*) * numThreads);

	#pragma omp parallel num_threads(numThreads)
	{
		int threadId = omp_get_thread_num();
		job_batch_t* threadBatchArr = batchArrPerThread[threadId];
		int threadBatchCount = batchCountPerThread[threadId];

		execHistoryPerThread[threadId] = lg_new();

		int i;
		for(i = 0; i < threadBatchCount; ++i)
		{
			#pragma omp task
			{
				threadId = omp_get_thread_num();

				job_batch_t* currBatch = &threadBatchArr[i];

				#if PROGRAM_MODE >= TRACE_MODE
					char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
					job_batch_toStringDetailed(currBatch, spmCsr, temp);
					PRINTF("Thread-%d executing %s.\n", threadId, temp);
				#endif

				lg_addTailData((void*) currBatch, execHistoryPerThread[threadId]);
				JOB_BATCH_EXECUTE_SUBMTX_CSR(currBatch, spmCsr, x, y_inout);
			}
		}
	}

	// return values
	*execHistoryPerThread_out = execHistoryPerThread;
}

void spmxv_omp_task_multiply_JDS(
		vector_real_t* x, vector_real_t* y_inout,
		job_batch_t** batchArrPerThread, int* batchCountPerThread, int numThreads,
		lg_t*** execHistoryPerThread_out)
{
	lg_t** execHistoryPerThread = (lg_t**) malloc(sizeof(lg_t*) * numThreads);

	#pragma omp parallel num_threads(numThreads)
	{
		int threadId = omp_get_thread_num();
		job_batch_t* threadBatchArr = batchArrPerThread[threadId];
		int threadBatchCount = batchCountPerThread[threadId];

		execHistoryPerThread[threadId] = lg_new();

		int i;
		for(i = 0; i < threadBatchCount; ++i)
		{
			#pragma omp task
			{
				threadId = omp_get_thread_num();

				job_batch_t* currBatch = &threadBatchArr[i];

				#if PROGRAM_MODE >= TRACE_MODE
					char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
					job_batch_toStringDetailed(currBatch, NULL, temp);
					PRINTF("Thread-%d executing %s.\n", threadId, temp);
				#endif

				lg_addTailData((void*) currBatch, execHistoryPerThread[threadId]);
				JOB_BATCH_EXECUTE_PARTIAL_JDS(currBatch, x, y_inout);
			}
		}
	}

	// return values
	*execHistoryPerThread_out = execHistoryPerThread;
}

// TODO not working for intel compiler
void spmxv_omp_task_multiply_hybrid_JDS_CSR(
		vector_real_t* x, vector_real_t* y_inout,
		job_batch_t** batchArrPerThread, int* batchCountPerThread,
		int numThreads, lg_t*** execHistoryPerThread_out)
{
	lg_t** execHistoryPerThread = (lg_t**) malloc(sizeof(lg_t*) * numThreads);

	#pragma omp parallel num_threads(numThreads)
	{
		int threadId = omp_get_thread_num();
		job_batch_t* threadBatchArr = batchArrPerThread[threadId];
		int threadBatchCount = batchCountPerThread[threadId];

		execHistoryPerThread[threadId] = lg_new();

		int i;
		for(i = 0; i < threadBatchCount; ++i)
		{
			#pragma omp task
			{
				threadId = omp_get_thread_num();

				job_batch_t* currBatch = &threadBatchArr[i];

				#if PROGRAM_MODE >= TRACE_MODE
					char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
					job_batch_toStringDetailed(currBatch, NULL, temp);
					PRINTF("Thread-%d executing %s.\n", threadId, temp);
				#endif

				lg_addTailData((void*) currBatch, execHistoryPerThread[threadId]);
				JOB_BATCH_EXECUTE_HYBRID_JDS_CSR(currBatch, x, y_inout);
			}
		}
	}

	// return values
	*execHistoryPerThread_out = execHistoryPerThread;
}

// Pre-multiplication routines / data preparations
// ---------------------------------------------------------------------------------------------------------------------

void spmxv_omp_task_prepFromSubMtxTree_rowWise_CSR(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromSubMtxTree_rowWise_CSRn");
	spmxv_static_prepFromSubMtxTree_rowWise_CSR((static_args_t*) args);
}

void spmxv_omp_task_prepFromSubMtxTreeScatter_rowWise_CSR(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromSubMtxTreeScatter_rowWise_CSR\n");
	spmxv_static_prepFromSubMtxTreeScatter_rowWise_CSR((static_args_t*) args);
}

void spmxv_omp_task_prepFromSubMtxList_rowWise_CSR(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromSubMtxList_rowWise_CSR\n");
	spmxv_static_prepFromSubMtxList_rowWise_CSR((static_args_t*) args);
}

void spmxv_omp_task_prepFromPartitionVector_rowWise_CSR(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromPartitionVector_rowWise_CSR\n");
	spmxv_static_prepFromPartitionVector_rowWise_colNet_CSR((static_args_t*) args);
}

void spmxv_omp_task_prepFromPartitionVectorScatter_rowWise_CSR(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromPartitionVectorScatter_rowWise_CSR\n");
	spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_CSR((static_args_t*) args);
}

// ---------------------------------------------------------------------------------------------------------------------

void spmxv_omp_task_prepFromSubMtxTree_rowWise_JDS(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromSubMtxTree_rowWise_JDS\n");
	spmxv_static_prepFromSubMtxTree_rowWise_JDS((static_args_t*) args);
}

void spmxv_omp_task_prepFromSubMtxTreeScatter_rowWise_JDS(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromSubMtxTreeScatter_rowWise_JDS\n");
	spmxv_static_prepFromSubMtxTreeScatter_rowWise_JDS((static_args_t*) args);
}

void spmxv_omp_task_prepFromSubMtxList_rowWise_JDS(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromSubMtxList_rowWise_JDS\n");
	spmxv_static_prepFromSubMtxList_rowWise_JDS((static_args_t*) args);
}

void spmxv_omp_task_prepFromPartitionVector_rowWise_JDS(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromPartitionVector_rowWise_JDS\n");
	spmxv_static_prepFromPartitionVector_rowWise_colNet_JDS((static_args_t*) args);
}

void spmxv_omp_task_prepFromPartitionVectorScatter_rowWise_JDS(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromPartitionVectorScatter_rowWise_JDS\n");
	spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_JDS((static_args_t*) args);
}

// ---------------------------------------------------------------------------------------------------------------------

void spmxv_omp_task_prepFromSubMtxTree_colNet_hybrid_JDS_CSR(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromSubMtxTree_colNet_hybrid_JDS_CSR\n");
	spmxv_static_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR((static_args_t*) args);
}

void spmxv_omp_task_prepFromSubMtxTreeScatter_colNet_hybrid_JDS_CSR(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromSubMtxTreeScatter_colNet_hybrid_JDS_CSR\n");
	spmxv_static_prepFromSubMtxTreeScatter_rowWise_hybrid_JDS_CSR((static_args_t*) args);
}

void spmxv_omp_task_prepFromSubMtxList_colNet_hybrid_JDS_CSR(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromSubMtxList_colNet_hybrid_JDS_CSR\n");
	spmxv_static_prepFromSubMtxList_rowWise_hybrid_JDS_CSR((static_args_t*) args);
}

void spmxv_omp_task_prepFromPartitionVector_colNet_hybrid_JDS_CSR(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromPartitionVector_colNet_hybrid_JDS_CSR\n");
	spmxv_static_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR((static_args_t*) args);
}

void spmxv_omp_task_prepFromPartitionVectorScatter_colNet_hybrid_JDS_CSR(omp_task_args_t* args)
{
	DEBUG("spmxv_omp_task_prepFromPartitionVectorScatter_colNet_hybrid_JDS_CSR\n");
	spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_hybrid_JDS_CSR((static_args_t*) args);
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
