
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "arch/generic/inc/spmxv_gws.h"
#include "include/data_structure/list_generic.h"
#include "include/control_unit/cu_options.h"
#include "include/control_unit/fast_run.h"

#include "include/control_unit/gws.h"

gws_args_t* gws_args_new()
{
	gws_args_t* args = (gws_args_t*) malloc(sizeof(gws_args_t));
	gws_args_initDefault(args);
	return args;
}

void gws_args_initDefault(gws_args_t* args)
{
	static_args_initDefault(&args->staticArgs);

	args->globalQueue = NULL;
	args->execHistoryPerThread = NULL;

	// function pointers
	// -------------------------------------------------------------------------------

	args->prep = NULL;
	args->spmxv_gws_warmup = NULL;
	args->spmxv_gws = NULL;
}

void gws_args_init(gws_args_t* args, fast_run_t* fr)
{
	DEBUG("GWS_ARGS\n");
	static_args_init(&args->staticArgs, fr);
	args->spmxv_gws = args->staticArgs.spmxv_static;

	int storageFormat = args->staticArgs.storageFormat;
	int orderingType = args->staticArgs.orderingType;
	int partitionType = args->staticArgs.partitionType;
	int partitionMethod = args->staticArgs.partitionMethod;

	if(storageFormat == SPM_STORAGE_CSR && partitionType == PARTITION_TYPE_1D_ROW_SLICE)
	{
		args->spmxv_gws_warmup = spmxv_gws_multiply_CSR;

		if(orderingType == ORDERING_TYPE_NONE)
		{
			if(partitionMethod == PARTITION_METHOD_REGULAR)
				args->prep = spmxv_gws_prepFromSubMtxList_colNet_CSR;
			else if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
				args->prep = spmxv_gws_prepFromSubMtxTree_colNet_CSR;
		}
		else if(orderingType == ORDERING_TYPE_COLUMN_NET)
		{
			// since there is only one queue, no need to consider
			// initial task assignment (partition method)
			args->prep = spmxv_gws_prepFromPartitionVector_colNet_CSR;
		}
	}
	else if(storageFormat == SPM_STORAGE_JDS && partitionType == PARTITION_TYPE_1D_ROW_SLICE)
	{
		args->spmxv_gws_warmup = spmxv_gws_multiply_JDS;

		if(orderingType == ORDERING_TYPE_NONE)
		{
			if(partitionMethod == PARTITION_METHOD_REGULAR)
				args->prep = spmxv_gws_prepFromSubMtxList_colNet_JDS;
			else if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
				args->prep = spmxv_gws_prepFromSubMtxTree_colNet_JDS;
		}
		else if(orderingType == ORDERING_TYPE_COLUMN_NET)
		{
			// since there is only one queue, no need to consider
			// initial task assignment (partition method)
			args->prep = spmxv_gws_prepFromPartitionVector_colNet_JDS;
		}
	}
	else if(storageFormat == SPM_STORAGE_HYBRID_JDS_CSR && partitionType == PARTITION_TYPE_1D_ROW_SLICE)
	{
		args->spmxv_gws_warmup = spmxv_gws_multiply_hybrid_JDS_CSR;

		if(orderingType == ORDERING_TYPE_NONE)
		{
			if(partitionMethod == PARTITION_METHOD_REGULAR)
				args->prep = spmxv_gws_prepFromSubMtxList_colNet_hybrid_JDS_CSR;
			else if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
				args->prep = spmxv_gws_prepFromSubMtxTree_colNet_hybrid_JDS_CSR;
		}
		else if(orderingType == ORDERING_TYPE_COLUMN_NET)
		{
			// since there is only one queue, no need to consider
			// initial task assignment (partition method)
			args->prep = spmxv_gws_prepFromPartitionVector_colNet_hybrid_JDS_CSR;
		}
	}
}

void gws_args_warmup(gws_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	lg_deleteShallowMultiple(args->execHistoryPerThread, args->staticArgs.numThreads);
	job_queue_reset(args->globalQueue);
	args->spmxv_gws_warmup(args, x, y_inout);
}

void gws_args_print(gws_args_t* args)
{
	PRINTF("GWS-ARGS\n");
	static_args_print(&args->staticArgs);
}

void gws_args_setup(gws_args_t* args)
{
	lg_job_batch_toArrayMultiple(
			args->execHistoryPerThread, args->staticArgs.numThreads,
			&args->staticArgs.jobBatchArrPerThread, &args->staticArgs.jobBatchCountPerThread);
}

void gws_args_reset(gws_args_t* args)
{

}

void gws_args_terminate(gws_args_t* args)
{
	static_args_terminate(&args->staticArgs);
	job_queue_delete(args->globalQueue);
	lg_deleteShallowMultiple(args->execHistoryPerThread, args->staticArgs.numThreads);
}

void gws_args_delete(gws_args_t* args)
{
	gws_args_terminate(args);
	free(args);
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
