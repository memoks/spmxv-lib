
#ifndef JOBBATCH_H_
#define JOBBATCH_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <omp.h>

#include "include/config.h"
#include "include/data_structure/spm_storage.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/list.h"
#include "include/data_structure/list_generic.h"
#include "include/data_structure/sub_mtx.h"

struct partial_jds
{
	spm_jds_t* spmJds;
	sub_mtx_dim_t subMtx;
};

typedef struct partial_jds partial_jds_t;

struct hybrid_jds_csr
{
	spm_jds_t* jds;
	spm_cmp_t* csr;
	sub_mtx_dim_t subMtx;
};

typedef struct hybrid_jds_csr hybrid_jds_csr_t;

union job_batch_data
{
	partial_jds_t partialJds;
	hybrid_jds_csr_t hybrid_jds_csr;
	sub_mtx_dim_t subMtx;
};

typedef union job_batch_data job_batch_data_t;

enum job_batch_type
{
	JOB_SUB_MTX_CSR,
	JOB_PARTIAL_JDS,
	JOB_HYBRID_JDS_CSR,
	JOB_TRUE_HYBRID,
	JOB_NONE
};

typedef enum job_batch_type job_batch_type_t;


struct job_batch
{
	job_batch_data_t data;
	job_batch_type_t type;
	list_head_t head;
};

typedef struct job_batch job_batch_t;


#define JOB_BATCH_EXECUTE_SUBMTX_CSR(					\
		/* job_batch_t* */ jb,							\
		/* spm_csr_t* */ spm,							\
		/* vector_real_t* */ x,							\
		/* vector_real_t* */ y_inout)					\
														\
	do													\
	{													\
		sub_mtx_dim_t* subMtx = &jb->data.subMtx;		\
		SUB_MTX_EXECUTE_CSR(subMtx, spm, x, y_inout);	\
	}													\
	while(0)

#define JOB_BATCH_EXECUTE_PARTIAL_JDS(					\
		/* job_batch_t* */ jb,							\
		/* vector_real_t* */ x,							\
		/* vector_real_t* */ y_inout)					\
														\
	do													\
	{													\
		spm_jds_t* spmJds = jb->data.partialJds.spmJds;	\
		spmxv_jds_partial_optimized(					\
				jb->data.partialJds.subMtx.start.i,		\
				0,										\
				spmJds->idiagLength,					\
				spmJds->idiag,							\
				spmJds->jdiag,							\
				spmJds->dj,								\
				x->data,								\
				y_inout->data);							\
	}													\
	while(0)											\

#ifdef JDS_BASED_HYBRID_DATA_EXTRACTION
#define JOB_BATCH_EXECUTE_HYBRID_JDS_CSR(jb, x, y) JOB_BATCH_EXECUTE_HYBRIDJDSCSR_JDS(jb, x, y)
#else
#define JOB_BATCH_EXECUTE_HYBRID_JDS_CSR(jb, x, y) JOB_BATCH_EXECUTE_HYBRIDJDSCSR_CSR(jb, x, y)
// For debug mode
// #define JOB_BATCH_EXECUTE_HYBRID_JDS_CSR(jb, x, y) job_batch_execute_hybrid_jds_csr(jb, x, y);
#endif

#define JOB_BATCH_EXECUTE_HYBRIDJDSCSR_JDS(						\
		/* job_batch_t* */ jb,										\
		/* vector_real_t* */ x,										\
		/* vector_real_t* */ y_inout)								\
																	\
	do																\
	{																\
		sub_mtx_dim_t* subMtx = &jb->data.hybrid_jds_csr.subMtx;	\
																	\
		/* execute partial JDS */									\
		if(jb->data.hybrid_jds_csr.jds != NULL)						\
		{															\
			spm_jds_t* spmJds = jb->data.hybrid_jds_csr.jds;		\
			spmxv_jds_partial_optimized(							\
					subMtx->start.i,								\
					0,												\
					spmJds->idiagLength,							\
					spmJds->idiag,									\
					spmJds->jdiag,									\
					spmJds->dj,										\
					x->data,										\
					y_inout->data);									\
		}															\
																	\
		/* execute partial CSR */									\
		if(jb->data.hybrid_jds_csr.csr != NULL)						\
		{															\
			spm_cmp_t* spmCsr = jb->data.hybrid_jds_csr.csr;		\
			spmxv_csr_add(											\
					0,												\
					spmCsr->rowCount,								\
					spmCsr->ptr,									\
					spmCsr->ind,									\
					spmCsr->values,									\
					x->data,										\
					&y_inout->data[subMtx->start.i]);				\
		}															\
	}																\
	while(0)

#define JOB_BATCH_EXECUTE_HYBRIDJDSCSR_CSR(						\
		/* job_batch_t* */ jb,										\
		/* vector_real_t* */ x,										\
		/* vector_real_t* */ y_inout)								\
																	\
	do																\
	{																\
		sub_mtx_dim_t* subMtx = &jb->data.hybrid_jds_csr.subMtx;	\
		DECIMAL rowShift = 0;										\
		/* execute partial CSR */									\
		if(jb->data.hybrid_jds_csr.csr != NULL)						\
		{															\
			spm_cmp_t* spmCsr = jb->data.hybrid_jds_csr.csr;		\
			spmxv_csr_add(											\
					0,												\
					spmCsr->rowCount,								\
					spmCsr->ptr,									\
					spmCsr->ind,									\
					spmCsr->values,									\
					x->data,										\
					&y_inout->data[subMtx->start.i]);				\
			rowShift = spmCsr->rowCount;							\
		}															\
																	\
		/* execute partial JDS */									\
		if(jb->data.hybrid_jds_csr.jds != NULL)						\
		{															\
			spm_jds_t* spmJds = jb->data.hybrid_jds_csr.jds;		\
			spmxv_jds_partial_optimized(							\
					subMtx->start.i + rowShift,						\
					0,												\
					spmJds->idiagLength,							\
					spmJds->idiag,									\
					spmJds->jdiag,									\
					spmJds->dj,										\
					x->data,										\
					y_inout->data);									\
		}															\
	}																\
	while(0)

/*
extern void job_batch_execute_hybrid_jds_csr(
		job_batch_t* jb, vector_real_t* x, vector_real_t* y_inout);
*/

extern job_batch_t* job_batch_new(void);
extern job_batch_t* job_batch_newSubMtx(sub_mtx_dim_t* subMtx);
extern job_batch_t* job_batch_newPartialJds(
		sub_mtx_dim_t* subMtx, spm_jds_t* spmPartialJds);
extern job_batch_t* job_batch_newHybridJdsCsr(
		sub_mtx_dim_t* subMtx,
		spm_jds_t* spmJds, spm_cmp_t* spmCsr);

// TODO test
extern job_batch_t* job_batch_newTrueHybrid(
		sub_mtx_dim_t* subMtx,
		DECIMAL csrNNZ, DECIMAL csrRowCount, DECIMAL csrColCount,
		DECIMAL jdsNNZ, DECIMAL jdsRowCount,
		DECIMAL jdsColCount, DECIMAL jdsIdiagLength);

extern void job_batch_initDefault(job_batch_t* jb);
extern void job_batch_initSubMtx(job_batch_t* jb, sub_mtx_dim_t* subMtx);
extern void job_batch_initPartialJds(
		job_batch_t* jb, sub_mtx_dim_t* subMtx, spm_jds_t* spmPartialJds);


extern void job_batch_initHybridJdsCsrDefault(job_batch_t* jb);
/**
 * ASSUMTION: Function copies sub-matrix structure. But it does not
 * copy JDS and CSR structures. So do not just free them after
 * calling this function.
 */
extern void job_batch_initHybridJdsCsr(
		job_batch_t* jb, sub_mtx_dim_t* subMtx,
		spm_jds_t* spmJds, spm_cmp_t* spmCsr);

// TODO test
extern void job_batch_initTrueHybrid(
		job_batch_t* jb, sub_mtx_dim_t* subMtx,
		DECIMAL csrNNZ, DECIMAL csrRowCount, DECIMAL csrColCount,
		DECIMAL jdsNNZ, DECIMAL jdsRowCount,
		DECIMAL jdsColCount, DECIMAL jdsIdiagLength);
extern void job_batch_initTrueHybridDefault(job_batch_t* jb);

extern void job_batch_copy(job_batch_t* source, job_batch_t* destination);
extern job_batch_t* job_batch_copyToPtr(job_batch_t* source);
extern void job_batch_deleteNonPtr(job_batch_t* jb);
extern sub_mtx_dim_t* job_batch_getSubMtx(job_batch_t* jb);
extern void job_batch_extractStatisticsPartialJds(
		job_batch_t* jobBatch,
		REAL* maxr_out, REAL* avgr_out, REAL* minr_out,
		REAL* maxc_out, REAL* avgc_out, REAL* minc_out);
extern void job_batch_extractStatisticsArrPartialJds(
		job_batch_t* jobBatchArr, int length,
		REAL* maxr_out, REAL* avgr_out, REAL* minr_out,
		REAL* maxc_out, REAL* avgc_out, REAL* minc_out);
extern void job_batch_extractStatisticsArrMultiplePartialJds(
		job_batch_t** jobBatchArrPerBlock,
		int* jobBatchCountPerBlock, int numBlocks,
		REAL* maxr_out, REAL* avgr_out, REAL* minr_out,
		REAL* maxc_out, REAL* avgc_out, REAL* minc_out);

/**
 * Removes given job batch from list and completely deletes it.
 */
extern void job_batch_delete(job_batch_t* jobBatch);

/**
 * Completely deletes every job batch item in a given list,
 * except, original head which is populated by NULL data.
 */
extern void job_batch_deleteAll(job_batch_t* head);

/**
 * Different from job_batch_deleteAll function in that this only
 * removes job batch from the given list (again original head
 * remains untouched).
 */
extern void job_batch_deleteAllShallow(job_batch_t* head);

extern void job_batch_deleteDenseArr(
		job_batch_t* jobBatchArr, int length);
extern void job_batch_deleteDenseArrMultiple(
		job_batch_t** jobBatchArrPerBlock,
		int* jobBatchCountPerBlock, int numBlocks);
extern job_batch_t* job_batch_getBatch(list_head_t* listHead);
extern void job_batch_addBatchFront(job_batch_t* new, job_batch_t* head);
extern void job_batch_addBatchTail(job_batch_t* new, job_batch_t* head);
extern job_batch_t* job_batch_getNext(job_batch_t* batch);
extern job_batch_t* job_batch_getPrev(job_batch_t* batch);
extern DECIMAL job_batch_getRowCount(job_batch_t* jobBatch);
extern void job_batch_toString(
		job_batch_t* jobBatch, char* buff_inout);
extern void job_batch_print(job_batch_t* jobBatch);
extern void job_batch_printArr(
		char* message, job_batch_t* jobBatchArr, int length);
extern void job_batch_printArrMultiple(
		char* message, job_batch_t** jobBatchArrPerBlock,
		int* jobBatchCountPerBlock, int numBlocks);

// functions that have abstraction layer for the type of job-batch
// depending on storage type, implementations of these functions may vary.
// -----------------------------------------------------------------------

extern DECIMAL job_batch_getNNZ(job_batch_t* jobBatch, void* spm);
extern void job_batch_toStringDetailed(
		job_batch_t* jobBatch, void* spm, char* buff_inout);
extern void job_batch_printDetailed(job_batch_t* jobBatch, void* spm);
extern void job_batch_printArrDetailed(
		char* message, job_batch_t* jobBatchArr, int length, void* spm);
extern void job_batch_printArrMultipleDetailed(
		char* message, job_batch_t** jobBatchArrPerBlock,
		int* jobBatchCountPerBlock, int numBlocks, void* spm);

// SpMxV data preparation & manipulation functions
// ------------------------------------------------------------------------
/**
 * Creates and returns a list of job batch from given sub-matrix list.
 * Job batch list will use copy data (no reference).
 *
 * @subMtxList: source data used to create job batches.
 * @jobBatchList_out: (return value) generated job batch list.
 */
// TODO test
extern void job_batch_extractJobBatchListFromSubMtxList_rowWiseColNet_CSR(
		lg_t* subMtxList, lg_t** jobBatchList_out);


/**
 * Creates and fills a list from a given sub-matrix tree node.
 * Uses row wise distribution and column net partition model
 * in which only the tree leafs are considered jobs.
 *
 * @spmCsr: Complete sparse-matrix in CSR format
 * @node: tree node (partition) from whose child nodes
 *  (sub-partitions) new job batches will be created.
 * @jobBatchList_out: (return value) job batch list created by
 *  function itself. It holds the job batches created from the
 *  leafs under a given tree node.
 */
extern void job_batch_extractJobBatchListFromSubMtxTree_rowWiseColNet_CSR(
		spm_cmp_t* spmCsr, tree_node_t* node, lg_t** jobBatchList_out);

/**
 * <p>Creates and fills a list directly from a vector which provides
 * partitions generated by hypergraph tool using column net model
 * for row wise distribution scheme.</p>
 *
 * @spmCsr: full sparse-matrix in CSR format.
 * @rowPartitioningVector: a vector containing partition info in
 *  which each element is # of rows (sum of which is a partition).
 * @startIndex: starting from this index of partition vector to
 *  generate job batches.
 * @endIndex: until this index (endIndex itself is not included)
 * @jobBatchList_out: (return value) job batch list that contains
 *  all the job batches created using partition vector.
 *
 *
 * Suppose we have;
 * <ul>
 *  <li>rowPartitioning = {0, 20, 40, 60, 80, 100}</li>
 *  <li>startIndex = 1</li>
 *  <li>endIndex = 3</li>
 * </ul>
 *
 * Then we will have;
 * jobBatchList_out = {(20-40), (40,60)}
 *
 *
 * <p>See input_readPartitionVectorFromFile function description in
 * input_parser.h for more information regarding rowPartitioning
 * vector.</p>
 */
extern void job_batch_extractJobBatchListFromPartitionVector_rowWiseColNet_CSR(
		spm_cmp_t* spmCsr, vector_int_t* rowPartitioningVector,
		DECIMAL startIndex, DECIMAL endIndex, lg_t** jobBatchList_out);

/**
 * Distributes the elements in a given job-batch list to
 * multiple lists in a scatter fashion.
 *
 * @source: list whose elements will be distributed among multiple lists.
 * @scatterSize: number of lists that source list will be distributed to.
 * @jobBatchLists_out: (return value) array of lists whose union was
 * equal to given source array.
 *
 *
 * Assume;
 * source = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9} and scatterSize = 4
 *
 * After function call situation will be;
 * source = {}
 * jobBatchLists[0] = {0, 4, 8}
 * jobBatchLists[1] = {1, 5, 9}
 * jobBatchLists[2] = {2, 6}
 * jobBatchLists[3] = {3, 7}
 *
 * This functions is used to improve cache usage between
 * threads in the same block.
 */
extern void job_batch_scatterJobBatchList(
		lg_t* source, int scatterSize, lg_t*** jobBatchLists_out);

/**
 * Merges 2 adjacent job batches in given job batch list if;
 * <ul>
 * <li>sub-matrix end index of the first one</li>
 * <li>sub-matrix start index of the second one</li>
 * </ul>
 * are equal. Note that number of elements in give list will be
 * decreased by the number of merge operations performed.
 * This was done to reduce outher for loop overhead in SpMxV.
 *
 * @jobBatchList (lg_t*): job batch list whose elements will be merged.
 *
 *
 * Suppose we have;
 * jobBatchList = {(1, 3), (4, 6), (6, 8), (3, 4)}
 *
 * after function call we will have;
 * jobBatchList = {(1, 3), (4, 8), (3, 4)}
 */
// TODO test
extern void job_batch_mergeJobBatchList_CSR(lg_t* jobBatchList);

// TODO test
extern void lg_job_batch_extractSubMtx(
		lg_t* jobBatchList, lg_t* subMtxList_inout);

/**
 * Function prototypes for generic list.
 */
extern FUNC_PROTOTYPE_LG_TOARRAY(job_batch, job_batch_t);
extern FUNC_PROTOTYPE_LG_TOARRAY_MULTIPLE(job_batch, job_batch_t);
extern FUNC_PROTOTYPE_LG_COPYADD(job_batch, job_batch_t, Tail);
extern FUNC_PROTOTYPE_LG_COPYADD(job_batch, job_batch_t, Front);
extern FUNC_PROTOTYPE_LG_PRINT(job_batch, job_batch_t);
extern FUNC_PROTOTYPE_LG_PRINT_MULTIPLE(job_batch, job_batch_t);
extern FUNC_PROTOTYPE_LG_DELETEDEEP(job_batch, job_batch_t);
extern FUNC_PROTOTYPE_LG_DELETEDEEP_MULTIPLE(job_batch, job_batch_t);
extern FUNC_PROTOTYPE_LG_PARTITIONDEEP(job_batch, job_batch_t, decompose);
extern FUNC_PROTOTYPE_LG_PARTITIONDEEP(job_batch, job_batch_t, scatter);
extern FUNC_PROTOTYPE_LG_PARTITIONDEEP(job_batch, job_batch_t, split);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* JOBBATCH_H_ */
