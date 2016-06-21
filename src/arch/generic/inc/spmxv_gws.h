
#ifndef SPMXV_GWS_H_
#define SPMXV_GWS_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"
#include "include/data_structure/spm_storage.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/list_generic.h"
#include "include/scheduler/job_queue.h"
#include "include/scheduler/job_batch.h"
#include "include/control_unit/gws.h"


/**
 * Global work stealing implementation for SpMxV in which there is
 * only one queue and all processes compete to detach work items from
 * that.
 *
 * Execution history per thread is recorded over a warm-up period to
 * discard locking overhead. After warm-up period static algorithms
 * are used with scheduling decisions fetched from execution history.
 */

extern void spmxv_gws_multiply_CSR(
		gws_args_t* args, vector_real_t* x, vector_real_t* y_inout);

extern void spmxv_gws_multiply_JDS(
		gws_args_t* args, vector_real_t* x, vector_real_t* y_inout);

extern void spmxv_gws_multiply_hybrid_JDS_CSR(
		gws_args_t* args, vector_real_t* x, vector_real_t* y_inout);


// Intermediate data preparation functions
// -----------------------------------------------------------------------

extern void spmxv_gws_prepFromSubMtxList_colNet_CSR(gws_args_t* args);
extern void spmxv_gws_prepFromSubMtxTree_colNet_CSR(gws_args_t* args);

extern void spmxv_gws_prepFromPartitionVector_colNet_CSR(gws_args_t* args);

// -----------------------------------------------------------------------

extern void spmxv_gws_prepFromSubMtxList_colNet_JDS(gws_args_t* args);
extern void spmxv_gws_prepFromSubMtxTree_colNet_JDS(gws_args_t* args);

extern void spmxv_gws_prepFromPartitionVector_colNet_JDS(gws_args_t* args);

// -----------------------------------------------------------------------

extern void spmxv_gws_prepFromSubMtxList_colNet_hybrid_JDS_CSR(gws_args_t* args);
extern void spmxv_gws_prepFromSubMtxTree_colNet_hybrid_JDS_CSR(gws_args_t* args);

extern void spmxv_gws_prepFromPartitionVector_colNet_hybrid_JDS_CSR(gws_args_t* args);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SPMXV_GWS_H_ */
