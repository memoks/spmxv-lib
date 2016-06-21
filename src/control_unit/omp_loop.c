
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>

#include "include/config.h"
#include "include/control_unit/cu_options.h"
#include "arch/generic/inc/spmxv_omp_loop.h"

#include "include/control_unit/omp_loop.h"

omp_loop_args_t* omp_loop_args_new()
{
	omp_loop_args_t* args = (omp_loop_args_t*) malloc(sizeof(omp_loop_args_t));
	omp_loop_args_initDefault(args);

	return args;
}

void omp_loop_args_initDefault(omp_loop_args_t* args)
{
	static_args_initDefault((static_args_t*) args);

	args->chunkSize = -1;

	// function pointers
	// ------------------------------------------------------------------------------------------

	args->prep = NULL;
	args->spmxv_omp_loop_dynamic = NULL;
}

void omp_loop_args_init(
		omp_loop_args_t* args, fast_run_t* fr, int chunkSize)
{
	DEBUG("OMP_LOOP_ARGS\n");
	DEBUG("chunksSize=%d\n", chunkSize);

	// Keep these to restore these attributes later on fast_run data structure
	int numBlocks = fr->numBlocks;
	int threadsPerBlock = fr->threadsPerBlock;
	int numThreads = fr->numThreads;

	fr->numBlocks = 1;
	fr->threadsPerBlock = 1;
	fr->numThreads = 1;
	static_args_init((static_args_t*) args, fr);
	args->chunkSize = chunkSize;

	// Restore attributes on fast_run data structure
	fr->numBlocks = numBlocks;
	fr->threadsPerBlock = threadsPerBlock;
	fr->numThreads = numThreads;

	int storageFormat = args->staticArgs.storageFormat;
	int orderingType = args->staticArgs.orderingType;
	int partitionType = args->staticArgs.partitionType;
	int partitionMethod = args->staticArgs.partitionMethod;

	// function pointers
	// ------------------------------------------------------------------------------------------

	if(storageFormat == SPM_STORAGE_CSR &&
		partitionType == PARTITION_TYPE_1D_ROW_SLICE)
	{
		args->spmxv_omp_loop_dynamic = spmxv_omp_loop_dynamic_CSR;
		args->spmxv_omp_loop_guided = spmxv_omp_loop_guided_CSR;

		if(orderingType == ORDERING_TYPE_NONE)
		{
			if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
				args->prep = spmxv_omp_loop_prepFromSubMtxTree_rowWise_CSR;
			else if(partitionMethod == PARTITION_METHOD_REGULAR)
				args->prep = spmxv_omp_loop_prepFromSubMtxList_rowWise_CSR;
		}
		else if(orderingType == ORDERING_TYPE_COLUMN_NET)
		{
			// Partitioning is done outside. And ordering (initial task assignment)
			// doesn't make any difference for omp_loop algorithm (since only one
			// array is required - instead of many arrays where task assignment matters).
			// So, only one partitioning method is needed.
			// if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
			args->prep = spmxv_omp_loop_prepFromPartitionVector_rowWise_colNet_CSR;
			// else if(partitionMethod == PARTITION_METHOD_REGULAR)
			//	args->prep = spmxv_omp_loop_prepFromPartitionVector_rowWise_colNet_CSR;
		}
	}
	else if(storageFormat == SPM_STORAGE_JDS &&
		partitionType == PARTITION_TYPE_1D_ROW_SLICE)
	{
		args->spmxv_omp_loop_dynamic = spmxv_omp_loop_dynamic_JDS;
		args->spmxv_omp_loop_guided = spmxv_omp_loop_guided_JDS;

		if(orderingType == ORDERING_TYPE_NONE)
		{
			if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
				args->prep = spmxv_omp_loop_prepFromSubMtxTree_rowWise_JDS;
			else if(partitionMethod == PARTITION_METHOD_REGULAR)
				args->prep = spmxv_omp_loop_prepFromSubMtxList_rowWise_JDS;
		}
		else if(orderingType == ORDERING_TYPE_COLUMN_NET)
		{
			// Partitioning is done outside. And ordering (initial task assignment)
			// doesn't make any difference for omp_loop algorithm (since only one
			// array is required - instead of many arrays where task assignment matters).
			// So, only one partitioning method is needed.
			// if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
			args->prep = spmxv_omp_loop_prepFromPartitionVector_rowWise_colNet_JDS;
			// else if(partitionMethod == PARTITION_METHOD_REGULAR)
			//	args->prep = spmxv_omp_loop_prepFromSubMtxList_rowWise_JDS;
		}
	}
	else if(storageFormat == SPM_STORAGE_HYBRID_JDS_CSR &&
		partitionType == PARTITION_TYPE_1D_ROW_SLICE)
	{
		args->spmxv_omp_loop_dynamic = spmxv_omp_loop_dynamic_hybrid_JDS_CSR;
		args->spmxv_omp_loop_guided = spmxv_omp_loop_guided_hybrid_JDS_CSR;

		if(orderingType == ORDERING_TYPE_NONE)
		{
			if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
				args->prep = spmxv_omp_loop_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR;
			else if(partitionMethod == PARTITION_METHOD_REGULAR)
				args->prep = spmxv_omp_loop_prepFromSubMtxList_rowWise_hybrid_JDS_CSR;
		}
		else if(orderingType == ORDERING_TYPE_COLUMN_NET)
		{
			// Partitioning is done outside. And ordering (initial task assignment)
			// doesn't make any difference for omp_loop algorithm (since only one
			// array is required - instead of many arrays where task assignment matters).
			// So, only one partitioning method is needed.
			// if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
			args->prep = spmxv_omp_loop_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR;
			// else if(partitionMethod == PARTITION_METHOD_REGULAR)
			//	args->prep = spmxv_omp_loop_prepFromSubMtxList_rowWise_hybrid_JDS_CSR;
		}
	}
}

void omp_loop_args_print(omp_loop_args_t* args)
{
	PRINTF("OMP-LOOP-ARGS\n");
	static_args_print(&args->staticArgs);
	PRINTF("chunk size: %d\n", args->chunkSize);
}

void omp_loop_args_terminate(omp_loop_args_t* args)
{
	static_args_terminate(&args->staticArgs);
}

void omp_loop_args_delete(omp_loop_args_t* args)
{
	omp_loop_args_terminate(args);
	free(args);
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
