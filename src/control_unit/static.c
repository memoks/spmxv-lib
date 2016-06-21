
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/control_unit/cu_options.h"
#include "arch/generic/inc/spmxv_static.h"

#include "include/control_unit/static.h"

static_args_t* static_args_new()
{
	static_args_t* args = (static_args_t*) malloc(sizeof(static_args_t));
	static_args_initDefault(args);

	return args;
}

void static_args_initDefault(static_args_t* args)
{
	args->quintets = NULL;
	args->nnz = -1;
	args->rowCount = -1;
	args->colCount = -1;

	// ASSUMPTION: This parameter should be left as it is.
	// If you want to use static-scatter algorithm then set this manually outside
	// before static_args_init function call.
	args->algorithmType = STATIC_ALGORITHM_TYPE_DEFAULT;

	args->storageFormat = -1;
	args->orderingType = -1;
	args->partitionType = -1;
	args->partitionMethod = -1;

	args->numBlocks = -1;
	args->threadsPerBlock = -1;
	args->numThreads = -1;

	args->simdLength = -1;
	args->targetedCacheSizeKB = -1;
	args->initialPartitioningData = NULL;
	args->subMtxHead = NULL;
	args->subMtxList = NULL;
	args->permutation = NULL;
	args->rowPartitioning = NULL;
	args->columnPartitioning = NULL;
	args->subMtxCount = -1;

	args->jobBatchArrPerThread = NULL;
	args->jobBatchCountPerThread = NULL;

	// function pointers
	// ----------------------------------------------------

	args->prep = NULL;
	args->spmxv_static = NULL;
}


void static_args_init(
		static_args_t* args, fast_run_t* fr)
{
	DEBUG("STATIC_ARGS\n");
	DEBUG("algorithmType=%d storageFormat=%d orderingType=%d partitionType=%d partitionMethod=%d\n",
		args->algorithmType, fr->storageFormat, fr->orderingType, fr->partitionType, fr->partitionMethod);

	args->quintets = fr->quintets;
	args->nnz = fr->nnz;
	args->rowCount = fr->rowCount;
	args->colCount = fr->colCount;

	args->storageFormat = fr->storageFormat;
	args->orderingType = fr->orderingType;
	args->partitionType = fr->partitionType;
	args->partitionMethod = fr->partitionMethod;

	args->numBlocks = fr->numBlocks;
	args->threadsPerBlock = fr->threadsPerBlock;
	args->numThreads = fr->numBlocks * fr->threadsPerBlock;

	args->simdLength = fr->simdLength;
	args->targetedCacheSizeKB = fr->targetedCacheSizeKB;
	args->initialPartitioningData = &fr->ipd;
	args->spmCsr = fr->spmCsr;
	args->permutation = fr->permutation;
	args->rowPartitioning = fr->rowPartitioning;
	args->columnPartitioning = fr->columnPartitioning;

	int algorithmType = args->algorithmType;
	int storageFormat = args->storageFormat;
	int orderingType = args->orderingType;
	int partitionMethod = args->partitionMethod;
	int partitionType = args->partitionType;

	if(orderingType == ORDERING_TYPE_NONE &&
		partitionMethod == PARTITION_METHOD_REGULAR)
	{
		args->subMtxList = args->initialPartitioningData->subMtxList;
		args->subMtxCount = lg_size(args->subMtxList);
	}
	else if(orderingType == ORDERING_TYPE_NONE &&
		partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
	{
		args->subMtxHead = args->initialPartitioningData->subMtxHead;
		args->subMtxCount = sub_mtx_tree_getLeafCount(args->subMtxHead);
	}
	else if(orderingType == ORDERING_TYPE_COLUMN_NET &&
		partitionType == PARTITION_TYPE_1D_ROW_SLICE)
	{
		args->subMtxCount = args->rowPartitioning->length;
	}

	// function pointers
	// -------------------------------------------------------------------------------------------

	if(storageFormat == SPM_STORAGE_CSR)
	{
		args->spmxv_static = spmxv_static_rowWise_CSR;

		if(orderingType == ORDERING_TYPE_NONE)
		{
			if(algorithmType == STATIC_ALGORITHM_TYPE_DEFAULT)
			{
				if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
					args->prep = spmxv_static_prepFromSubMtxTree_rowWise_CSR;
				else if(partitionMethod == PARTITION_METHOD_REGULAR)
					args->prep = spmxv_static_prepFromSubMtxList_rowWise_CSR;
			}
			else if(algorithmType == STATIC_ALGORITHM_TYPE_SCATTER)
			{
				if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
					args->prep = spmxv_static_prepFromSubMtxTreeScatter_rowWise_CSR;
				else if(partitionMethod == PARTITION_METHOD_REGULAR)
					args->prep = spmxv_static_prepFromSubMtxListScatter_rowWise_CSR;
			}
		}
		else if(orderingType == ORDERING_TYPE_COLUMN_NET)
		{
			if(algorithmType == STATIC_ALGORITHM_TYPE_DEFAULT)
			{
				if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
					args->prep = spmxv_static_prepFromPartitionVector_rowWise_colNet_CSR;
				else if(partitionMethod == PARTITION_METHOD_REGULAR)
					args->prep = spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_CSR;
			}
			else if(algorithmType == STATIC_ALGORITHM_TYPE_SCATTER)
			{
				if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
					args->prep = spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_CSR;
				else if(partitionMethod == PARTITION_METHOD_REGULAR)
					args->prep = spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_CSR;
			}
		}
	}
	else if(storageFormat == SPM_STORAGE_JDS)
	{
		args->spmxv_static = spmxv_static_rowWise_JDS;

		if(algorithmType == STATIC_ALGORITHM_TYPE_DEFAULT)
		{
			if(orderingType == ORDERING_TYPE_NONE)
			{
				if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
					args->prep = spmxv_static_prepFromSubMtxTree_rowWise_JDS;
				else if(partitionMethod == PARTITION_METHOD_REGULAR)
					args->prep = spmxv_static_prepFromSubMtxList_rowWise_JDS;
			}
			else if(orderingType == ORDERING_TYPE_COLUMN_NET)
			{
				if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
					args->prep = spmxv_static_prepFromPartitionVector_rowWise_colNet_JDS;
				else if(partitionMethod == PARTITION_METHOD_REGULAR)
					args->prep = spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_JDS;
			}
		}
		else if(algorithmType == STATIC_ALGORITHM_TYPE_SCATTER)
		{
			if(orderingType == ORDERING_TYPE_NONE)
			{
				if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
					args->prep = spmxv_static_prepFromSubMtxTreeScatter_rowWise_JDS;
				else if(partitionMethod == PARTITION_METHOD_REGULAR)
					args->prep = spmxv_static_prepFromSubMtxListScatter_rowWise_JDS;
			}
			else if(orderingType == ORDERING_TYPE_COLUMN_NET)
			{
				if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
					args->prep = spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_JDS;
				else if(partitionMethod == PARTITION_METHOD_REGULAR)
					args->prep = spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_JDS;
			}
		}
	}
	else if(storageFormat == SPM_STORAGE_HYBRID_JDS_CSR)
	{
		args->spmxv_static = spmxv_static_rowWise_hybrid_JDS_CSR;

		if(algorithmType == STATIC_ALGORITHM_TYPE_DEFAULT)
		{
			if(orderingType == ORDERING_TYPE_NONE)
			{
				if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
					args->prep = spmxv_static_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR;
				else if(partitionMethod == PARTITION_METHOD_REGULAR)
					args->prep = spmxv_static_prepFromSubMtxList_rowWise_hybrid_JDS_CSR;
			}
			else if(orderingType == ORDERING_TYPE_COLUMN_NET)
			{
				if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
					args->prep = spmxv_static_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR;
				else if(partitionMethod == PARTITION_METHOD_REGULAR)
					args->prep = spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_hybrid_JDS_CSR;
			}
		}
		else if(algorithmType == STATIC_ALGORITHM_TYPE_SCATTER)
		{
			if(orderingType == ORDERING_TYPE_NONE)
			{
				if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
					args->prep = spmxv_static_prepFromSubMtxTreeScatter_rowWise_hybrid_JDS_CSR;
				else if(partitionMethod == PARTITION_METHOD_REGULAR)
					args->prep = spmxv_static_prepFromSubMtxListScatter_rowWise_hybrid_JDS_CSR;
			}
			else if(orderingType == ORDERING_TYPE_COLUMN_NET)
			{
				if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
					args->prep = spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_hybrid_JDS_CSR;
				else if(partitionMethod == PARTITION_METHOD_REGULAR)
					args->prep = spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_hybrid_JDS_CSR;
			}
		}
	}
}

void static_args_print(static_args_t* args)
{
	PRINTF("STATIC-ARGS\n");
	PRINTF("NNZ: %d, row count: %d, column count: %d\n", args->nnz, args->rowCount, args->colCount);

	PRINTF("algorithm type: %d\n", args->algorithmType);
	PRINTF("storage format: %d\n", args->storageFormat);
	PRINTF("ordering type: %d\n", args->orderingType);
	PRINTF("partition type: %d\n", args->partitionType);
	PRINTF("partition method: %d\n", args->partitionMethod);

	PRINTF("num-blocks: %d, threads-per-block: %d, num-threads: %d\n",
			args->numBlocks, args->threadsPerBlock, args->numThreads);

	PRINTF("simd-length: %d, targeted cache size in KB: %f\n", args->simdLength, args->targetedCacheSizeKB);
	PRINTF("sub-mtx-count: %d\n", args->subMtxCount);

	if(args->jobBatchCountPerThread != NULL)
	{
		PRINTF("job-batch-count per thread:");
		int i;
		for(i = 0; i < args->numThreads; ++i)
			PRINTF(" %d", args->jobBatchCountPerThread[i]);
		PRINTF("\n");
	}
}

void static_args_terminate(static_args_t* args)
{
	job_batch_deleteDenseArrMultiple(args->jobBatchArrPerThread, args->jobBatchCountPerThread, args->numThreads);
	free(args->jobBatchCountPerThread);
}

void static_args_delete(static_args_t* args)
{
	static_args_terminate(args);
	free(args);
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
