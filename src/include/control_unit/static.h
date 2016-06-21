
#ifndef SRC_INCLUDE_CONTROL_UNIT_STATIC_H_
#define SRC_INCLUDE_CONTROL_UNIT_STATIC_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"
#include "include/data_structure/quintet.h"
#include "include/data_structure/sub_mtx.h"
#include "include/data_structure/spm_storage.h"
#include "include/data_structure/vector.h"
#include "include/scheduler/job_batch.h"
#include "include/task_decomposition/partitioning.h"
#include "include/control_unit/fast_run.h"

// TODO comment
struct static_args
{
	quintet_t* quintets;
	DECIMAL nnz;
	DECIMAL rowCount;
	DECIMAL colCount;

	int algorithmType;
	int storageFormat;
	int orderingType;
	int partitionType;
	int partitionMethod;

	int numBlocks;
	int threadsPerBlock;
	int numThreads;

	int simdLength;
	float targetedCacheSizeKB;
	ipd_t* initialPartitioningData;
	sub_mtx_tree_t* subMtxHead;
	lg_t* subMtxList;
	spm_cmp_t* spmCsr;
	DECIMAL* permutation;
	vector_int_t* rowPartitioning;
	vector_int_t* columnPartitioning;
	int subMtxCount;

	job_batch_t** jobBatchArrPerThread;
	int* jobBatchCountPerThread;

	// function pointers
	// -----------------------------------------------------------------------------------------

	// for spmxv multiply and data-preparation routines
	// (warm-up doesn't change for static routines - same as spmxv)
	void (*prep) (struct static_args* args);
	void (*spmxv_static) (struct static_args* args, vector_real_t* x, vector_real_t* y_out);
};

typedef struct static_args static_args_t;

extern static_args_t* static_args_new();
extern void static_args_initDefault(static_args_t* args);
extern void static_args_init(static_args_t* args, fast_run_t* fr);
extern void static_args_print(static_args_t* args);
extern void static_args_terminate(static_args_t* args);
extern void static_args_delete(static_args_t* args);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SRC_INCLUDE_CONTROL_UNIT_STATIC_STATIC_H_ */
