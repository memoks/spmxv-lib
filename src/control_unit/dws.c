
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>

#include "include/config.h"
#include "include/control_unit/cu_options.h"
#include "include/control_unit/dws.h"

#include "arch/generic/inc/spmxv_sequential.h"
#include "arch/generic/inc/spmxv_omp_loop.h"
#include "arch/generic/inc/spmxv_omp_task.h"
#include "arch/generic/inc/spmxv_static.h"
#include "arch/generic/inc/spmxv_gws.h"
#include "arch/generic/inc/spmxv_dws.h"

dws_args_t* dws_args_new()
{
	dws_args_t* args = (dws_args_t*) malloc(sizeof(dws_args_t));
	dws_args_initDefault(args);
	return args;
}

void dws_args_initDefault(dws_args_t* args)
{
	static_args_initDefault(&args->staticArgs);

	args->algorithmType = -1;
	args->stealingScheme = -1;

	args->spmCsr = NULL;
	args->permutation = NULL;
	args->rowPartitioning = NULL;
	args->columnPartitioning = NULL;
	args->subMtxCount = -1;

	args->warps = NULL;
	args->numWarps = -1;
	args->threadsPerWarp = -1;
	args->stealTreshold = DEFAULT_STEAL_TRESHOLD;
	args->execHistoryPerWarp = NULL;
	args->execHistoryPerThread = NULL;

	// function pointers
	// ---------------------------------------------------------------

	args->prep = NULL;
	args->spmxv_dws_warmup = NULL;
	args->spmxv_dws = NULL;

	args->scheduler_init = NULL;
	args->scheduler_findVictim = NULL;
	args->scheduler_terminate = NULL;
}

static char* binaryRepresentation(int n)
{
	static char buff[35]; buff[34] = '\0';
	int rem = 0;
	int i = 33;

	while(n > 0)
	{
		rem = n % 2;
		n = n / 2;

		buff[i] = rem + '0';
		--i;
	}

	for(; 2 <= i; --i)
		buff[i] = '0';

	buff[0] = '0';
	buff[1] = 'b';

	return buff;
}

void dws_args_init(
		dws_args_t* args, fast_run_t* fr,
		int algorithmType, int stealingScheme, int stealTreshold)
{
	DEBUG("DWS_ARGS\n");
	DEBUG("algorithmType=%d stealingScheme=%d stealTreshold=%d\n",
			algorithmType, stealingScheme, stealTreshold);

	static_args_init(&args->staticArgs, fr);
	args->spmxv_dws = args->staticArgs.spmxv_static;

	args->algorithmType = algorithmType;
	args->stealingScheme = stealingScheme;

	args->spmCsr = fr->spmCsr;
	args->permutation = fr->permutation;
	args->rowPartitioning = fr->rowPartitioning;
	args->columnPartitioning = fr->columnPartitioning;

	int numBlocks = args->staticArgs.numBlocks;
	int threadsPerBlock = args->staticArgs.threadsPerBlock;
	int numThreads = args->staticArgs.numThreads;
	int storageFormat = args->staticArgs.storageFormat;
	int orderingType = args->staticArgs.orderingType;
	int partitionType = args->staticArgs.partitionType;
	int partitionMethod = args->staticArgs.partitionMethod;

	if(stealingScheme == STEALING_SCHEME_RING)
	{
		args->scheduler_init = ring_scheduler_init;
		args->scheduler_findVictim = ring_scheduler_findVictim;
		args->scheduler_terminate = ring_scheduler_terminate;
	}
	else if(stealingScheme == STEALING_SCHEME_TREE)
	{
		args->scheduler_init = tree_scheduler_init;
		args->scheduler_findVictim = tree_scheduler_findVictim;
		args->scheduler_terminate = tree_scheduler_terminate;
	}

	if(storageFormat == SPM_STORAGE_CSR)
	{
		args->spmxv_dws = spmxv_static_rowWise_CSR;

		if(algorithmType == DWS_ALGORITHM_TYPE_DEFAULT ||
			algorithmType == DWS_ALGORITHM_TYPE_SCATTER)
		{
			args->spmxv_dws_warmup = spmxv_dws_rowWise_CSR;
			args->numWarps = numThreads;
			args->threadsPerWarp = 1;
		}
		else if(algorithmType == DWS_ALGORITHM_TYPE_SHARED_BLOCK)
		{
			args->spmxv_dws_warmup = spmxv_dws_sharedBlock_rowWise_CSR;
			args->numWarps = numBlocks;
			args->threadsPerWarp = threadsPerBlock;
		}
		if(orderingType == ORDERING_TYPE_NONE &&
			partitionType == PARTITION_TYPE_1D_ROW_SLICE)
		{
			if(partitionMethod == PARTITION_METHOD_REGULAR)
				args->prep = spmxv_dws_prepFromSubMtxList_rowWise_CSR;
			else if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
				args->prep = spmxv_dws_prepFromSubMtxTree_rowWise_CSR;
		}
		else if(orderingType == ORDERING_TYPE_COLUMN_NET &&
			partitionType == PARTITION_TYPE_1D_ROW_SLICE)
		{
			if(partitionMethod == PARTITION_METHOD_REGULAR)
			{
				if(algorithmType == DWS_ALGORITHM_TYPE_SCATTER)
					args->prep = spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_CSR;
				else if(algorithmType == DWS_ALGORITHM_TYPE_DEFAULT || algorithmType == DWS_ALGORITHM_TYPE_SHARED_BLOCK)
					args->prep = spmxv_dws_prepFromPartitionVector_rowWise_colNet_CSR;
			}
			else if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
			{
				if(algorithmType == DWS_ALGORITHM_TYPE_SCATTER)
					args->prep = spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_CSR;
				else if(algorithmType == DWS_ALGORITHM_TYPE_DEFAULT || algorithmType == DWS_ALGORITHM_TYPE_SHARED_BLOCK)
					args->prep = spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_CSR;
			}
		}
	}
	else if(storageFormat == SPM_STORAGE_JDS)
	{
		args->spmxv_dws = spmxv_static_rowWise_JDS;

		if(algorithmType == DWS_ALGORITHM_TYPE_DEFAULT ||
			algorithmType == DWS_ALGORITHM_TYPE_SCATTER)
		{
			args->spmxv_dws_warmup = spmxv_dws_rowWise_JDS;
			args->numWarps = numThreads;
			args->threadsPerWarp = 1;
		}
		else if(algorithmType == DWS_ALGORITHM_TYPE_SHARED_BLOCK)
		{
			args->spmxv_dws_warmup = spmxv_dws_sharedBlock_rowWise_JDS;
			args->numWarps = numBlocks;
			args->threadsPerWarp = threadsPerBlock;
		}

		if(orderingType == ORDERING_TYPE_NONE &&
			partitionType == PARTITION_TYPE_1D_ROW_SLICE)
		{
			if(partitionMethod == PARTITION_METHOD_REGULAR)
				args->prep = spmxv_dws_prepFromSubMtxList_rowWise_JDS;
			else if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
				args->prep = spmxv_dws_prepFromSubMtxTree_rowWise_JDS;
		}
		else if(orderingType == ORDERING_TYPE_COLUMN_NET &&
			partitionType == PARTITION_TYPE_1D_ROW_SLICE)
		{
			if(partitionMethod == PARTITION_METHOD_REGULAR)
			{
				if(algorithmType == DWS_ALGORITHM_TYPE_SCATTER)
					args->prep = spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_JDS;
				else if(algorithmType == DWS_ALGORITHM_TYPE_DEFAULT || algorithmType == DWS_ALGORITHM_TYPE_SHARED_BLOCK)
					args->prep = spmxv_dws_prepFromPartitionVector_rowWise_colNet_JDS;
			}
			else if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
			{
				if(algorithmType == DWS_ALGORITHM_TYPE_SCATTER)
					args->prep = spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_JDS;
				else if(algorithmType == DWS_ALGORITHM_TYPE_DEFAULT || algorithmType == DWS_ALGORITHM_TYPE_SHARED_BLOCK)
					args->prep = spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_JDS;
			}
		}
	}
	else if(storageFormat == SPM_STORAGE_HYBRID_JDS_CSR)
	{
		args->spmxv_dws = spmxv_static_rowWise_hybrid_JDS_CSR;

		if(algorithmType == DWS_ALGORITHM_TYPE_DEFAULT ||
			algorithmType == DWS_ALGORITHM_TYPE_SCATTER)
		{
			args->spmxv_dws_warmup = spmxv_dws_rowWise_hybrid_JDS_CSR;
			args->numWarps = numThreads;
			args->threadsPerWarp = 1;
		}
		else if(algorithmType == DWS_ALGORITHM_TYPE_SHARED_BLOCK)
		{
			args->spmxv_dws_warmup = spmxv_dws_sharedBlock_rowWise_hybrid_JDS_CSR;
			args->numWarps = numBlocks;
			args->threadsPerWarp = threadsPerBlock;
		}

		if(orderingType == ORDERING_TYPE_NONE &&
			partitionType == PARTITION_TYPE_1D_ROW_SLICE)
		{
			if(partitionMethod == PARTITION_METHOD_REGULAR)
				args->prep = spmxv_dws_prepFromSubMtxList_rowWise_hybrid_JDS_CSR;
			else if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
				args->prep = spmxv_dws_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR;
		}
		else if(orderingType == ORDERING_TYPE_COLUMN_NET &&
			partitionType == PARTITION_TYPE_1D_ROW_SLICE)
		{
			if(partitionMethod == PARTITION_METHOD_REGULAR)
			{
				if(algorithmType == DWS_ALGORITHM_TYPE_SCATTER)
					args->prep = spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_hybrid_JDS_CSR;
				else if(algorithmType == DWS_ALGORITHM_TYPE_DEFAULT || algorithmType == DWS_ALGORITHM_TYPE_SHARED_BLOCK)
					args->prep = spmxv_dws_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR;
			}
			else if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
			{
				if(algorithmType == DWS_ALGORITHM_TYPE_SCATTER)
					args->prep = spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_hybrid_JDS_CSR;
				else if(algorithmType == DWS_ALGORITHM_TYPE_DEFAULT || algorithmType == DWS_ALGORITHM_TYPE_SHARED_BLOCK)
				{
					args->prep = spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_hybrid_JDS_CSR;
				}
			}
		}
	}
}

void dws_args_print(dws_args_t* args)
{
	PRINTF("DWS-ARGS\n");
	static_args_print(&args->staticArgs);
	PRINTF("algorithm type: %d\n", args->algorithmType);
	PRINTF("stealing scheme: %d\n", args->stealingScheme);
	PRINTF("steal treshold: %d\n", args->stealTreshold);
	PRINTF("num-warps: %d threads-per-warp: %d\n", args->numWarps, args->threadsPerWarp);
}

void dws_args_warmup(dws_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	lg_deleteShallowMultiple(args->execHistoryPerThread, args->staticArgs.numThreads);
	lg_deleteShallowMultiple(args->execHistoryPerWarp, args->numWarps);
	args->spmxv_dws_warmup(args, x, y_inout);
	block_reformMultiple(args->warps, args->execHistoryPerWarp, args->numWarps);
}

void dws_args_setup(dws_args_t* args)
{
	if(args->algorithmType == DWS_ALGORITHM_TYPE_SHARED_BLOCK)
	{
		lg_job_batch_toArrayMultiple(
				args->execHistoryPerThread, args->staticArgs.numThreads,
				&args->staticArgs.jobBatchArrPerThread, &args->staticArgs.jobBatchCountPerThread);
	}
	else // DWS_ALGORITHM_TYPE_DEFAULT DWS_ALGORITHM_TYPE_SCATTER
	{
		lg_job_batch_toArrayMultiple(
				args->execHistoryPerWarp, args->numWarps,
				&args->staticArgs.jobBatchArrPerThread, &args->staticArgs.jobBatchCountPerThread);
	}
}

void dws_args_reset(dws_args_t* args)
{
	block_resetMultiple(args->warps, args->numWarps);
}

void dws_args_terminate(dws_args_t* args)
{
	static_args_terminate(&args->staticArgs);
	block_deleteMultiple(args->warps, args->numWarps);

	if(args->execHistoryPerThread != NULL)
		lg_deleteShallowMultiple(args->execHistoryPerThread, args->staticArgs.numThreads);

	if(args->execHistoryPerWarp != NULL)
		lg_deleteShallowMultiple(args->execHistoryPerWarp, args->numWarps);
}

void dws_args_delete(dws_args_t* args)
{
	dws_args_terminate(args);
	free(args);
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
