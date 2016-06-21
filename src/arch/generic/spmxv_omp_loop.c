
#include <omp.h>

#include "inc/kernel.h"
#include "inc/spmxv_omp_loop.h"
#include "include/io/converter.h"
#include "include/scheduler/job_batch.h"

#include "inc/spmxv_static.h"

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

// SpMxV routines
// ---------------------------------------------------------------------------------------------------------------------------------

inline void spmxv_omp_loop_dynamic_CSR(omp_loop_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	spm_cmp_t* spmCsr = args->staticArgs.spmCsr;
	DECIMAL chunkSize = args->chunkSize;

	DECIMAL i;
	DECIMAL j;
	#pragma omp parallel for private(i, j) schedule(dynamic, chunkSize)
	SPMXV_CSR_PARTIAL_FOR_LOOP(
			0, spmCsr->rowCount, spmCsr->ptr, spmCsr->ind, spmCsr->values, x->data, y_inout->data)
}

inline void spmxv_omp_loop_guided_CSR(omp_loop_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	spm_cmp_t* spmCsr = args->staticArgs.spmCsr;

	DECIMAL i;
	DECIMAL j;
	#pragma omp parallel for private(i, j) schedule(guided)
	SPMXV_CSR_PARTIAL_FOR_LOOP(
			0, spmCsr->rowCount, spmCsr->ptr, spmCsr->ind, spmCsr->values, x->data, y_inout->data)
}

inline void spmxv_omp_loop_dynamic_CSC(
		spm_cmp_t* spmCsc, vector_real_t* x, vector_real_t* y_inout, DECIMAL chunkSize)
{
	DECIMAL i;
	DECIMAL j;
	#pragma omp parallel for private(i, j) schedule(dynamic, chunkSize)
	SPMXV_CSC_PARTIAL_FOR_LOOP_ATOMIC(
			0, spmCsc->rowCount, spmCsc->ptr, spmCsc->ind, spmCsc->values, x->data, y_inout->data)
}

inline void spmxv_omp_loop_dynamic_JDS(
		omp_loop_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	job_batch_t* jobBatchArr = args->staticArgs.jobBatchArrPerThread[0];
	int length = args->staticArgs.jobBatchCountPerThread[0];
	DECIMAL chunkSize = args->chunkSize;

	DECIMAL i;
	#pragma omp parallel for private(i) schedule(dynamic, chunkSize)
	for(i = 0; i < length; ++i)
	{
		job_batch_t* currBatch = &jobBatchArr[i];
		JOB_BATCH_EXECUTE_PARTIAL_JDS(currBatch, x, y_inout);
	}
}

inline void spmxv_omp_loop_guided_JDS(
		omp_loop_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	job_batch_t* jobBatchArr = args->staticArgs.jobBatchArrPerThread[0];
	int length = args->staticArgs.jobBatchCountPerThread[0];

	DECIMAL i;
	#pragma omp parallel for private(i) schedule(guided)
	for(i = 0; i < length; ++i)
	{
		job_batch_t* currBatch = &jobBatchArr[i];
		JOB_BATCH_EXECUTE_PARTIAL_JDS(currBatch, x, y_inout);
	}
}


inline void spmxv_omp_loop_dynamic_hybrid_JDS_CSR(
		omp_loop_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	job_batch_t* jobBatchArr = args->staticArgs.jobBatchArrPerThread[0];
	int length = args->staticArgs.jobBatchCountPerThread[0];
	DECIMAL chunkSize = args->chunkSize;

	DECIMAL i;
	#pragma omp parallel for private(i) schedule(dynamic, chunkSize)
	for(i = 0; i < length; ++i)
	{
		job_batch_t* currBatch = &jobBatchArr[i];
		JOB_BATCH_EXECUTE_HYBRID_JDS_CSR(currBatch, x, y_inout);
	}
}

inline void spmxv_omp_loop_guided_hybrid_JDS_CSR(
		omp_loop_args_t* args, vector_real_t* x, vector_real_t* y_inout)
{
	job_batch_t* jobBatchArr = args->staticArgs.jobBatchArrPerThread[0];
	int length = args->staticArgs.jobBatchCountPerThread[0];

	DECIMAL i;
	#pragma omp parallel for private(i) schedule(guided)
	for(i = 0; i < length; ++i)
	{
		job_batch_t* currBatch = &jobBatchArr[i];
		JOB_BATCH_EXECUTE_HYBRID_JDS_CSR(currBatch, x, y_inout);
	}
}

// ---------------------------------------------------------------------------------------------------------------------------------

void spmxv_omp_loop_prepFromSubMtxTree_rowWise_CSR(omp_loop_args_t* args)
{
	DEBUG("spmxv_omp_loop_prepFromSubMtxTree_rowWise_CSR\n");
	spmxv_static_prepFromSubMtxTree_rowWise_CSR((static_args_t*) args);
}

void spmxv_omp_loop_prepFromSubMtxList_rowWise_CSR(omp_loop_args_t* args)
{
	DEBUG("spmxv_omp_loop_prepFromSubMtxList_rowWise_CSR\n");
	spmxv_static_prepFromSubMtxList_rowWise_CSR((static_args_t*) args);
}

void spmxv_omp_loop_prepFromPartitionVector_rowWise_colNet_CSR(omp_loop_args_t* args)
{
	DEBUG("spmxv_omp_loop_prepFromPartitionVector_rowWise_colNet_CSR\n");
	spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_CSR((static_args_t*) args);
}

// ---------------------------------------------------------------------------------------------------------------------------------

void spmxv_omp_loop_prepFromSubMtxTree_rowWise_JDS(omp_loop_args_t* args)
{
	DEBUG("spmxv_omp_loop_prepFromSubMtxTree_rowWise_JDS\n");
	spmxv_static_prepFromSubMtxTree_rowWise_JDS((static_args_t*) args);
}

void spmxv_omp_loop_prepFromSubMtxList_rowWise_JDS(omp_loop_args_t* args)
{
	DEBUG("spmxv_omp_loop_prepFromSubMtxList_rowWise_JDS\n");
	spmxv_static_prepFromSubMtxList_rowWise_JDS((static_args_t*) args);
}

void spmxv_omp_loop_prepFromPartitionVector_rowWise_colNet_JDS(omp_loop_args_t* args)
{
	DEBUG("spmxv_omp_loop_prepFromPartitionVector_rowWise_colNet_JDS\n");
	spmxv_static_prepFromPartitionVector_rowWise_colNet_JDS((static_args_t*) args);
}

// ---------------------------------------------------------------------------------------------------------------------------------

void spmxv_omp_loop_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR(omp_loop_args_t* args)
{
	DEBUG("spmxv_omp_loop_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR\n");

	int numThreads = args->staticArgs.numThreads;
	args->staticArgs.numThreads = 1;
	spmxv_static_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR((static_args_t*) args);
	args->staticArgs.numThreads = numThreads;
}

void spmxv_omp_loop_prepFromSubMtxList_rowWise_hybrid_JDS_CSR(omp_loop_args_t* args)
{
	DEBUG("spmxv_omp_loop_prepFromSubMtxList_rowWise_hybrid_JDS_CSR\n");

	int numThreads = args->staticArgs.numThreads;
	args->staticArgs.numThreads = 1;
	spmxv_static_prepFromSubMtxList_rowWise_hybrid_JDS_CSR((static_args_t*) args);
	args->staticArgs.numThreads = numThreads;
}

void spmxv_omp_loop_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR(omp_loop_args_t* args)
{
	DEBUG("spmxv_omp_loop_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR\n");
	spmxv_static_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR((static_args_t*) args);
}


#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

