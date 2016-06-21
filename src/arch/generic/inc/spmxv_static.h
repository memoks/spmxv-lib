/*
 * spmxv_static.h
 *
 * Static algorithms do not attempt to balance execution load.
 * Each execution context works on the sub-matrices they are
 * initially assigned to.
 */

#ifndef SPMXV_STATIC_H_
#define SPMXV_STATIC_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"
#include "include/io/input_parser.h"
#include "include/data_structure/spm_storage.h"
#include "include/data_structure/vector.h"
#include "include/scheduler/job_batch.h"
#include "include/control_unit/static.h"


/**
 * Static implementation of distributed work stealing algorithm in
 * which each thread has its own dense job-batch array (instead of
 * a queue). Sub-matrix execution order of each thread depends
 * solely on data preparation function used.
 */
extern void spmxv_static_rowWise_CSR(
		static_args_t* args, vector_real_t* x, vector_real_t* y_inout);

extern void spmxv_static_rowWise_JDS(
		static_args_t* args, vector_real_t* x, vector_real_t* y_inout);

extern void spmxv_static_rowWise_hybrid_JDS_CSR(
		static_args_t* args, vector_real_t* x, vector_real_t* y_inout);


// Pre-multiplication routines / data preparations
// ---------------------------------------------------------------

// now that we have 2 preparation methods for block and scatter
// there is no need for "int mergeJobBatches" as a parameter ??

extern void spmxv_static_prepFromSubMtxTree_rowWise_CSR(static_args_t* args);
extern void spmxv_static_prepFromSubMtxTreeScatter_rowWise_CSR(static_args_t* args);

extern void spmxv_static_prepFromSubMtxList_rowWise_CSR(static_args_t* args);
extern void spmxv_static_prepFromSubMtxListScatter_rowWise_CSR(static_args_t* args);

extern void spmxv_static_prepFromPartitionVector_rowWise_colNet_CSR(static_args_t* args);
extern void spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_CSR(static_args_t* args);

extern void spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_CSR(static_args_t* args);
extern void spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_CSR(static_args_t* args);

// ---------------------------------------------------------------

extern void spmxv_static_prepFromSubMtxTree_rowWise_JDS(static_args_t* args);
extern void spmxv_static_prepFromSubMtxTreeScatter_rowWise_JDS(static_args_t* args);

extern void spmxv_static_prepFromSubMtxList_rowWise_JDS(static_args_t* args);
extern void spmxv_static_prepFromSubMtxListScatter_rowWise_JDS(static_args_t* args);

extern void spmxv_static_prepFromPartitionVector_rowWise_colNet_JDS(static_args_t* args);
extern void spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_JDS(static_args_t* args);

extern void spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_JDS(static_args_t* args);
extern void spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_JDS(static_args_t* args);

// ---------------------------------------------------------------

extern void spmxv_static_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR(static_args_t* args);
extern void spmxv_static_prepFromSubMtxTreeScatter_rowWise_hybrid_JDS_CSR(static_args_t* args);

extern void spmxv_static_prepFromSubMtxList_rowWise_hybrid_JDS_CSR(static_args_t* args);
extern void spmxv_static_prepFromSubMtxListScatter_rowWise_hybrid_JDS_CSR(static_args_t* args);

extern void spmxv_static_prepFromPartitionVector_rowWise_colNet_hybrid_JDS_CSR(static_args_t* args);
extern void spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_hybrid_JDS_CSR(static_args_t* args);

extern void spmxv_static_prepFromPartitionVectorPowerOf2_rowWise_colNet_hybrid_JDS_CSR(static_args_t* args);
extern void spmxv_static_prepFromPartitionVectorScatterPowerOf2_rowWise_colNet_hybrid_JDS_CSR(static_args_t* args);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SPMXV_STATIC_H_ */
