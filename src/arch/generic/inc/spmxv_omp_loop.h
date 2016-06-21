
#ifndef SPMXV_OMP_LOOP_H_
#define SPMXV_OMP_LOOP_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/data_structure/spm_storage.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/list_generic.h"
#include "include/scheduler/job_batch.h"
#include "include/control_unit/omp_loop.h"

/**
 * Very simple parallel implementations using OpenMP loop construct.
 */

// SpMxV routines
// -----------------------------------------------------------------------------------------------------

extern void spmxv_omp_loop_dynamic_CSR(
		omp_loop_args_t* args, vector_real_t* x, vector_real_t* y_inout);

extern void spmxv_omp_loop_guided_CSR(
		omp_loop_args_t* args, vector_real_t* x, vector_real_t* y_inout);

// -----------------------------------------------------------------------------------------------------

extern void spmxv_omp_loop_dynamic_CSC(
		spm_cmp_t* spmCsc, vector_real_t* x, vector_real_t* y_inout,
		DECIMAL chunkSize);

// -----------------------------------------------------------------------------------------------------

// TODO test
extern void spmxv_omp_loop_dynamic_JDS(
		omp_loop_args_t* args, vector_real_t* x, vector_real_t* y_inout);

// TODO test
extern void spmxv_omp_loop_guided_JDS(
		omp_loop_args_t* args, vector_real_t* x, vector_real_t* y_inout);

// -----------------------------------------------------------------------------------------------------

extern void spmxv_omp_loop_dynamic_hybrid_JDS_CSR(
		omp_loop_args_t* args, vector_real_t* x, vector_real_t* y_inout);

extern void spmxv_omp_loop_guided_hybrid_JDS_CSR(
		omp_loop_args_t* args, vector_real_t* x, vector_real_t* y_inout);


// SpMxV data preparation
// -----------------------------------------------------------------------------------------------------

extern void spmxv_omp_loop_prepFromSubMtxTree_rowWise_CSR(omp_loop_args_t* args);
extern void spmxv_omp_loop_prepFromSubMtxList_rowWise_CSR(omp_loop_args_t* args);
extern void spmxv_omp_loop_prepFromPartitionVector_rowWise_colNet_CSR(omp_loop_args_t* args);

// -----------------------------------------------------------------------------------------------------

extern void spmxv_omp_loop_prepFromSubMtxList_rowWise_JDS(omp_loop_args_t* args);
extern void spmxv_omp_loop_prepFromSubMtxTree_rowWise_JDS(omp_loop_args_t* args);
extern void spmxv_omp_loop_prepFromPartitionVector_rowWise_colNet_JDS(omp_loop_args_t* args);

// -----------------------------------------------------------------------------------------------------

extern void spmxv_omp_loop_prepFromSubMtxList_rowWise_hybrid_JDS_CSR(omp_loop_args_t* args);
extern void spmxv_omp_loop_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR(omp_loop_args_t* args);
extern void spmxv_omp_loop_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR(
		omp_loop_args_t* args);


#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SPMXV_OMP_LOOP_H_ */
