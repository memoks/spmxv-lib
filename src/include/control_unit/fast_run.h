
#ifndef SRC_INCLUDE_CONTROL_UNIT_FAST_RUN_H_
#define SRC_INCLUDE_CONTROL_UNIT_FAST_RUN_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif


#include "include/config.h"
#include "include/data_structure/spm_storage.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/sub_mtx.h"
#include "include/task_decomposition/partitioning.h"

// TODO comment
struct fast_run
{
	// these will not be freed upon termination
	// -------------------------------------------------------------------------------
	quintet_t* quintets;
	int rowCount;
	int colCount;
	int nnz;

	int numBlocks;
	int threadsPerBlock;
	int numThreads;

	ipd_t ipd;

	spm_cmp_t* spmCsr;

	spm_cmp_t* _spmCsr;
	spm_cmp_t* _spmCsrCounterpart;

	DECIMAL* permutation;
	vector_int_t* _rowOrderLookupJDS;
	vector_int_t* yVectorRowOrderLookup;

	int subMtxCount;

	// function parameters for control_unit
	// -------------------------------------------------------------------------------
	int storageFormat;
	int orderingType;
	int partitionType;
	int partitionMethod;

	int targetedCacheSizeKB;
	int simdLength;

	// order & partitioning data
	// -------------------------------------------------------------------------------
	vector_int_t* rowOrderLookup;
	partitioning_vector_t* rowPartitioning;
	vector_int_t* colOrderLookup;
	vector_int_t* columnPartitioning;
};

typedef struct fast_run fast_run_t;

extern fast_run_t* fast_run_new(void);
extern void fast_run_initDefault(fast_run_t* fr);
extern void fast_run_init(
		fast_run_t* fr,
		quintet_t* quintets, int rowCount, int colCount, int nnz,
		int numBlocks, int threadsPerBlock,
		int targetedCacheSizeKB, int simdLength,
		int storageFormat, int orderingType,
		int partitionType, int partitionMethod,
		vector_int_t* rowOrderLookup,
		vector_int_t* rowPartitioning,
		vector_int_t* colOrderLookup,
		vector_int_t* colPartitioning);
extern void fast_run_terminate(fast_run_t* fr);
extern void fast_run_delete(fast_run_t* fr);


#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SRC_INCLUDE_CONTROL_UNIT_FAST_RUN_H_ */
