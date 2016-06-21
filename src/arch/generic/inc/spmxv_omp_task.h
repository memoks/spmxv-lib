
#ifndef SPMXV_OMP_TASK_H_
#define SPMXV_OMP_TASK_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"
#include "include/scheduler/job_batch.h"
#include "include/data_structure/spm_storage.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/list_generic.h"
#include "include/control_unit/omp_task.h"

/**
 * Uses OpenMP task primitive to balance execution load through
 * run time. This is different from naive OpenMP loop
 * implementations in that it uses blocking. However load balancing
 * is not locality aware (and some kind of locking is utilized
 * by OpenMP runtime).
 */


// Core
// -----------------------------------------------------------------

extern void spmxv_omp_task_multiply_CSR(
		spm_cmp_t* spmCsr, vector_real_t* x, vector_real_t* y_inout,
		job_batch_t** batchArrPerThread, int* batchCountPerThread,
		int numThreads,
		lg_t*** executionHistoryPerThread_out);

extern void spmxv_omp_task_multiply_JDS(
		vector_real_t* x, vector_real_t* y_inout,
		job_batch_t** batchArrPerThread, int* batchCountPerThread,
		int numThreads,
		lg_t*** executionHistoryPerThread_out);

extern void spmxv_omp_task_multiply_hybrid_JDS_CSR(
		vector_real_t* x, vector_real_t* y_inout,
		job_batch_t** batchArrPerThread, int* batchCountPerThread,
		int numThreads,
		lg_t*** executionHistoryPerThread_out);

// Pre-multiplication routines / data preparations
// ------------------------------------------------------------------

// This implementation uses the same pre-multiplication routines
// with spmxv_static implementation. These are just wrapper functions.

extern void spmxv_omp_task_prepFromSubMtxTree_rowWise_CSR(omp_task_args_t* args);
extern void spmxv_omp_task_prepFromSubMtxTreeScatter_rowWise_CSR(omp_task_args_t* args);

extern void spmxv_omp_task_prepFromSubMtxList_rowWise_CSR(omp_task_args_t* args);

extern void spmxv_omp_task_prepFromPartitionVector_rowWise_CSR(omp_task_args_t* args);
extern void spmxv_omp_task_prepFromPartitionVectorScatter_rowWise_CSR(omp_task_args_t* args);

// ------------------------------------------------------------------

extern void spmxv_omp_task_prepFromSubMtxTree_rowWise_JDS(omp_task_args_t* args);
extern void spmxv_omp_task_prepFromSubMtxTreeScatter_rowWise_JDS(omp_task_args_t* args);

extern void spmxv_omp_task_prepFromSubMtxList_rowWise_JDS(omp_task_args_t* args);

extern void spmxv_omp_task_prepFromPartitionVector_rowWise_JDS(omp_task_args_t* args);
extern void spmxv_omp_task_prepFromPartitionVectorScatter_rowWise_JDS(omp_task_args_t* args);

// ------------------------------------------------------------------

extern void spmxv_omp_task_prepFromSubMtxTree_colNet_hybrid_JDS_CSR(omp_task_args_t* args);
extern void spmxv_omp_task_prepFromSubMtxTreeScatter_colNet_hybrid_JDS_CSR(omp_task_args_t* args);

extern void spmxv_omp_task_prepFromSubMtxList_colNet_hybrid_JDS_CSR(omp_task_args_t* args);

extern void spmxv_omp_task_prepFromPartitionVector_colNet_hybrid_JDS_CSR(omp_task_args_t* args);
extern void spmxv_omp_task_prepFromPartitionVectorScatter_colNet_hybrid_JDS_CSR(omp_task_args_t* args);


#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SPMXV_OMP_TASK_H_ */
