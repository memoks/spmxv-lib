
#ifndef SUB_MTX_H_
#define SUB_MTX_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <omp.h>

#include "include/config.h"
#include "include/data_structure/list.h"
#include "include/data_structure/tree.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/spm_storage.h"
#include "include/data_structure/list_generic.h"

#include "arch/generic/inc/kernel.h"

// Index
// --------------------------------------------------------------------

struct index
{
	DECIMAL i; // row
	DECIMAL j; // column
};

typedef struct index index_t;

extern void index_initDefault(index_t* index);
extern void index_init(index_t* index, DECIMAL i, DECIMAL j);
extern index_t* index_new(DECIMAL i, DECIMAL j);
extern int index_equals(index_t* left, index_t* right);
extern void index_delete(index_t* index);


// Sub-matrix dimensions
// --------------------------------------------------------------------

struct sub_mtx_dim
{
	index_t start;
	index_t length;
};

typedef struct sub_mtx_dim sub_mtx_dim_t;

extern sub_mtx_dim_t* sub_mtx_new(void);
extern void sub_mtx_initDefault(sub_mtx_dim_t* subMtx);
extern void sub_mtx_init(sub_mtx_dim_t* subMtx,
		DECIMAL startRow, DECIMAL startCol,
		DECIMAL rowCount, DECIMAL colCount);
extern void sub_mtx_copy(
		sub_mtx_dim_t* source, sub_mtx_dim_t* dest);
extern sub_mtx_dim_t* sub_mtx_copyToPtr(sub_mtx_dim_t* subMtx);
extern void sub_mtx_delete(sub_mtx_dim_t* subMtx);
extern void sub_mtx_deleteNonPtr(sub_mtx_dim_t* subMtx);
extern void sub_mtx_toString(sub_mtx_dim_t* subMtx, char* buff_inout);
extern void sub_mtx_print(sub_mtx_dim_t* subMtx);
extern DECIMAL sub_mtx_equals(sub_mtx_dim_t* left, sub_mtx_dim_t* right);

/**
 * ASSUMPTION: assumes that overwritten sub-matrix end row is
 * the same as other sub-matrix's start row.
 *
 * Merges two give sequential sub-matrices by simply extending
 * the row count of the former.
 */
extern void sub_mtx_merge(
		sub_mtx_dim_t* overwritten, sub_mtx_dim_t* other);

/**
 * Function prototypes for generic list.
 */
extern FUNC_PROTOTYPE_LG_TOARRAY(sub_mtx, sub_mtx_dim_t);
extern FUNC_PROTOTYPE_LG_TOARRAY_MULTIPLE(sub_mtx, sub_mtx_dim_t);
extern FUNC_PROTOTYPE_LG_COPYADD(sub_mtx, sub_mtx_dim_t, Tail);
extern FUNC_PROTOTYPE_LG_COPYADD(sub_mtx, sub_mtx_dim_t, Front);
extern FUNC_PROTOTYPE_LG_PRINT(sub_mtx, sub_mtx_dim_t);
extern FUNC_PROTOTYPE_LG_PRINT_MULTIPLE(sub_mtx, sub_mtx_dim_t);
extern FUNC_PROTOTYPE_LG_DELETEDEEP(sub_mtx, sub_mtx_dim_t);
extern FUNC_PROTOTYPE_LG_DELETEDEEP_MULTIPLE(sub_mtx, sub_mtx_dim_t);
extern FUNC_PROTOTYPE_LG_PARTITIONDEEP(sub_mtx, sub_mtx_dim_t*, decompose);
extern FUNC_PROTOTYPE_LG_PARTITIONDEEP(sub_mtx, sub_mtx_dim_t*, scatter);
extern FUNC_PROTOTYPE_LG_PARTITIONDEEP(sub_mtx, sub_mtx_dim_t*, split);

// functions whose implementation varies depending on the storage format
// ---------------------------------------------------------------------

// CSR (Compressed Row Storage) Format
// ---------------------------------------------------------------------
extern void sub_mtx_print_CSR(sub_mtx_dim_t* subMtx, spm_cmp_t* spmCsr);
extern void sub_mtx_toString_CSR(
		sub_mtx_dim_t* subMtx, spm_cmp_t* spmCsr, char* buff_inout);
extern DECIMAL sub_mtx_getNNZ_CSR(sub_mtx_dim_t* subMtx, spm_cmp_t* spmCsr);
extern DECIMAL sub_mtx_isEmpty_CSR(sub_mtx_dim_t* subMtx, spm_cmp_t* spmCsr);

// JDS (Jagged Diagonal Storage) Format
// ---------------------------------------------------------------------
/**
 * ASSUMPTION: Functions which take both sub-matrix and SPM-JDS as
 * parameters assumes that sub-matrix is a chunk of SPM. Thus extracts
 * information about that chunk from global SPM.
 */
extern void sub_mtx_print_JDS(sub_mtx_dim_t* subMtx, spm_jds_t* spmJds);
extern void sub_mtx_toString_JDS(
		sub_mtx_dim_t* subMtx, spm_jds_t* spmJds, char* buff_inout);
extern DECIMAL sub_mtx_getNNZ_JDS(sub_mtx_dim_t* subMtx, spm_jds_t* spmJds);
extern DECIMAL sub_mtx_isEmpty_JDS(sub_mtx_dim_t* subMtx, spm_jds_t* spmJds);
extern void sub_mtx_extractRowStats_JDS(
		spm_jds_t* spmJds, REAL* maxr_out, REAL* avgr_out, REAL* minr_out);
extern void sub_mtx_extractColumnStats_JDS(
		spm_jds_t* spmJds, REAL* maxc_out, REAL* avgc_out, REAL* minc_out);

// Kernel wrapper interface for sub-matrix execution
// ---------------------------------------------------------------------

#define SUB_MTX_EXECUTE_CSR(						\
		/* sub_mtx_dim_t* */ subMtx,				\
		/* spm_cmp_t* */ spm, 						\
		/* vector_real_t* */ x, 					\
		/* vector_real_t* */ y_inout)				\
													\
	spmxv_csr_partial(								\
			subMtx->start.i, 						\
			subMtx->start.i + subMtx->length.i,		\
			spm->ptr,	 							\
			spm->ind, 								\
			spm->values, 							\
			x->data, 								\
			y_inout->data)

/**
 * ASSUMPTION: From CSR, group of row-slices are
 * implied. Therefore only y pointer is changed.
 */
#define SUB_MTX_EXECUTE_CSR_PARTIAL(				\
		/* sub_mtx_dim_t* */ subMtx,				\
		/* spm_cmp_t* */ spm,						\
		/* vector_real_t* */ x,						\
		/* vector_real_t* */ y_inout)				\
													\
	spmxv_csr_partial(								\
			0,										\
			spm->rowCount,							\
			spm->ptr,								\
			spm->ind,								\
			spm->values,							\
			x->data,								\
			&y_inout->data[subMtx->start.i])

// (current sub_matrix - JDS ordered sub_matrix) map for JDS format
// --------------------------------------------------------------------

struct sub_mtx_map
{
	sub_mtx_dim_t* initialSubMtx;
	DECIMAL maxNNZLength;
	sub_mtx_dim_t* afterOrderingSubMtx;
};

typedef struct sub_mtx_map sub_mtx_map_t;

extern sub_mtx_map_t* sub_mtx_map_new(
		sub_mtx_dim_t* subMtx, DECIMAL maxNNZLength);
extern void sub_mtx_map_init(
		sub_mtx_map_t* subMtxMap_inout,
		sub_mtx_dim_t* subMtx, DECIMAL maxNNZLength);
extern void sub_mtx_map_deleteNonPtr(sub_mtx_map_t* subMtxMap);
extern void sub_mtx_map_delete(sub_mtx_map_t* subMtxMap);
extern void sub_mtx_map_deleteNonPtrArray(
		sub_mtx_map_t* subMtxMapArr, int length);
extern void sub_mtx_map_print(sub_mtx_map_t* subMtxMap);
extern void sub_mtx_map_printArr(
		sub_mtx_map_t* subMtxMaps, int length);
extern int sub_mtx_map_cmpInverse(
		const void* left, const void* right);
extern void sub_mtx_map_sortDescending(
		sub_mtx_map_t* subMtxMapArr, int length);

/**
 * This function is used to successfully convert quintets
 * to JDS format without loosing partitioning information,
 * with minimum amount of zero-padding.
 *
 * ASSUMPTION: Before calling this function quintets must
 * be partially sorted by now (meaning NNZs are sorted in
 * a sub-matrix, but out of that scope (between different
 * sub-matrices), relative places of each NNZ remains the
 * same).
 *
 * Sorts sub-matrices by their longest non-zero element
 * count in descending order. It is aimed that, non-zeros
 * in the first sub-matrix in this array will be the first
 * non-zeros in JDS structure.
 *
 * After sub-matrix ordering new places of existing
 * sub-matrices are calculated and returned with sub-matrix
 * map array.
 */
extern void sub_mtx_map_extractOrderingForQuintets(
		sub_mtx_map_t* subMtxMapArr_inout, int length,
		vector_int_t* rowOrderLookup_inout);

// Sub-matrix tree
// --------------------------------------------------------------------

struct sub_mtx_tree
{
	tree_node_t node;

	int done;
	omp_lock_t writeLock;
	sub_mtx_dim_t* subMtx;

	omp_lock_t* borderLock;
};

typedef struct sub_mtx_tree sub_mtx_tree_t;

extern void sub_mtx_tree_init(sub_mtx_tree_t* subMtx);
extern sub_mtx_tree_t* sub_mtx_tree_new(
		DECIMAL startRow, DECIMAL startCol,
		DECIMAL rowCount, DECIMAL colCount);
extern void sub_mtx_tree_delete(tree_node_t* head);
extern void sub_mtx_tree_deleteSingle(sub_mtx_tree_t* subMtx);
extern void sub_mtx_tree_addLeft(
		sub_mtx_tree_t* parent, sub_mtx_tree_t* left);
extern void sub_mtx_tree_addRight(
		sub_mtx_tree_t* parent, sub_mtx_tree_t* left);
extern sub_mtx_tree_t* sub_mtx_tree_getLeft(sub_mtx_tree_t* parent);
extern sub_mtx_tree_t* sub_mtx_tree_getRight(sub_mtx_tree_t* parent);
extern sub_mtx_tree_t* sub_mtx_tree_getSubMtxTree(tree_node_t* node);
extern void sub_mtx_tree_printSingle(sub_mtx_tree_t* subMtxNode);
extern void sub_mtx_tree_toStringSingle(
		sub_mtx_tree_t* subMtxNode, char* buff_inout);
extern void sub_mtx_tree_printInOrder(tree_node_t* node);
extern void sub_mtx_tree_printPostOrder(tree_node_t* node);
extern void sub_mtx_tree_printLeafs(tree_node_t* node);
extern int sub_mtx_tree_isFull(sub_mtx_tree_t* subMtx);
extern int sub_mtx_tree_isLeaf(sub_mtx_tree_t* subMtx);
extern int sub_mtx_tree_getLeafCount(sub_mtx_tree_t* subMtx);

extern tree_node_t** sub_mtx_tree_getJobDistribution(
		sub_mtx_tree_t* head, int numBlocks);

extern void sub_mtx_tree_markTreeNodes(
		tree_node_t* node,
		tree_node_t*** nodesPerBlock_out, int start, int count);

extern void sub_mtx_tree_getLeafContentsShallow(
		sub_mtx_tree_t* subMtxTreeNode, lg_t** subMtxList_out);

extern void sub_mtx_tree_getLeafContentsDeep(
		sub_mtx_tree_t* subMtxTreeNode, lg_t** subMtxList_out);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* SUB_MTX_H_ */
