
#ifndef SPMXV_DWS_H_
#define SPMXV_DWS_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/data_structure/spm_storage.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/sub_mtx.h"
#include "include/data_structure/list_generic.h"
#include "include/scheduler/job_queue.h"
#include "include/scheduler/block_info.h"
#include "include/scheduler/ring_sched.h"
#include "include/scheduler/tree_sched.h"
#include "include/control_unit/dws.h"

/**
 * These algorithms utilize distributed work stealing scheme
 * where multiple queues are involved. After all the work items
 * which belong to one block (queue) finished, that block will
 * look to steal work items from other groups.
 *
 * A block can have a single or multiple execution contexts.
 */

/**
 * Regular distributed work stealing routine in which each block
 * has only one execution context.
 */
extern void spmxv_dws_rowWise_CSR(
		dws_args_t* args,
		vector_real_t* x, vector_real_t* y_out);

/**
 * Each block can have multiple execution contexts.
 */
extern void spmxv_dws_sharedBlock_rowWise_CSR(
		dws_args_t* args,
		vector_real_t* x, vector_real_t* y_inout);

extern void spmxv_dws_rowWise_JDS(
		dws_args_t* args,
		vector_real_t* x, vector_real_t* y_inout);

extern void spmxv_dws_sharedBlock_rowWise_JDS(
		dws_args_t* args,
		vector_real_t* x, vector_real_t* y_inout);

extern void spmxv_dws_rowWise_hybrid_JDS_CSR(
		dws_args_t* args,
		vector_real_t* x, vector_real_t* y_inout);

extern void spmxv_dws_sharedBlock_rowWise_hybrid_JDS_CSR(
		dws_args_t* args,
		vector_real_t* x, vector_real_t* y_inout);


// Pre-multiplication routines / data preparations
// ------------------------------------------------------------------

// TODO test all
extern void spmxv_dws_prepFromSubMtxTree_rowWise_CSR(dws_args_t* args);
extern void spmxv_dws_prepFromSubMtxTreeScatter_rowWise_CSR(dws_args_t* args);

extern void spmxv_dws_prepFromSubMtxList_rowWise_CSR(dws_args_t* args);

extern void spmxv_dws_prepFromPartitionVector_rowWise_colNet_CSR(dws_args_t* args);
extern void spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_CSR(dws_args_t* args);

extern void spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_CSR(dws_args_t* args);
extern void spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_CSR(dws_args_t* args);
// ------------------------------------------------------------------

extern void spmxv_dws_prepFromSubMtxTree_rowWise_JDS(dws_args_t* args);
extern void spmxv_dws_prepFromSubMtxTreeScatter_rowWise_JDS(dws_args_t* args);

extern void spmxv_dws_prepFromSubMtxList_rowWise_JDS(dws_args_t* args);

extern void spmxv_dws_prepFromPartitionVector_rowWise_colNet_JDS(dws_args_t* args);
extern void spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_JDS(dws_args_t* args);

extern void spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_JDS(dws_args_t* args);
extern void spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_JDS(dws_args_t* args);
// ------------------------------------------------------------------

extern void spmxv_dws_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR(dws_args_t* args);
extern void spmxv_dws_prepFromSubMtxTreeScatter_rowWise_hybrid_JDS_CSR(dws_args_t* args);

extern void spmxv_dws_prepFromSubMtxList_rowWise_hybrid_JDS_CSR(dws_args_t* args);

extern void spmxv_dws_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR(dws_args_t* args);
extern void spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_hybrid_JDS_CSR(dws_args_t* args);

extern void spmxv_dws_prepFromPartitionVectorPowerOf2_rowWise_colNet_hybrid_JDS_CSR(dws_args_t* args);
extern void spmxv_dws_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_hybrid_JDS_CSR(dws_args_t* args);

// Data structure clean-up
// ------------------------------------------------------------------

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif // SPMXV_DWS_H_

