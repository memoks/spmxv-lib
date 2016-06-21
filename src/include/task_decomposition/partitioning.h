
#ifndef PARTITIONING_H_
#define PARTITIONING_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"
#include "include/data_structure/sub_mtx.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/list_generic.h"
#include "include/data_structure/spm_storage.h"


/*
 * ASSUMPTION: Partitioning vector is an array of length [partition-count + 1]
 * (accumulated like CSR row-pointer), where:
 *
 * [partitionVector[i + 1] - partitionVector[i]] denotes row/column count of i.th partition
 *
 * Suppose we have a row partitioningVector={0, 3, 5, 7, 9}
 * <ul>
 *   <li>Then there are 4 partitions in total</li>
 *   <li>0th partition has 3 rows</li>
 *   <li>1st partition has 2 rows</li>
 *   <li>2nd partition has 2 rows</li>
 *   <li>3rd partition has 2 rows</li>
 * </ul>
 */

typedef vector_int_t partitioning_vector_t;

// TODO comment all
// ----------------------------------------------------------------------------------------------

// initial partitioning data
struct ipd
{
	lg_t* subMtxList;
	sub_mtx_tree_t* subMtxHead;
};

typedef struct ipd ipd_t;

extern ipd_t* ipd_new(void);
extern void ipd_initDefault(ipd_t* ipd);
extern void ipd_terminate(ipd_t* ipd);
extern void ipd_delete(ipd_t* ipd);

// ----------------------------------------------------------------------------------------------

extern void partitioning_1DRowWise(
		spm_cmp_t* spmCsr, int partitionMethod,
		ipd_t* ipd, sub_mtx_dim_t** subMtxArr_out, int* subMtxCount_out,
		REAL targetedCacheSizeInKB, DECIMAL rowIncrementSize, DECIMAL thresholdRowCount,
		REAL (*calculateSizeCallBack)(DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount));

extern void partitioning_1DRowSliceRecursiveBipartition_CSR(
		spm_cmp_t* spmCsr, REAL targetedCacheSizeInKB,
		sub_mtx_tree_t** subMtxTree_out,
		REAL (*calculateSizeCallBack)(DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount));

extern void partitioning_1DRowSliceLinear_CSR(
		spm_cmp_t* spmCsr, REAL targetedSizeInKB,
		DECIMAL rowIncrementSize, DECIMAL thresholdRowCount, lg_t** subMtxList_out,
		REAL (*calculateSizeCallBack)(DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount));

extern void partitioning_1DStaticRowSliceLinear_CSR(
		DECIMAL rowCount, DECIMAL colCount, int rowSliceSize,
		lg_t** subMtxList_out);

// ----------------------------------------------------------------------------------------------

extern REAL partitioning_calculateSizeInKBDefault(
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount);

// TODO test
extern DECIMAL* partitioning_calculateSubMtxCountPerBlockPowerOf2(
		int numBlocks, DECIMAL subMtxCount);

// TODO comment
extern partitioning_vector_t* partitioning_generatePartitioningFromPartVector(
		vector_int_t* partVector, int partCount);

// TODO test
extern vector_int_t* partitioning_generateOrderingFromPartVector(
		vector_int_t* partVector, int partCount);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* PARTITIONING_H_ */
