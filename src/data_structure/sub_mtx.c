
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "include/data_structure/sub_mtx.h"

// Index
// -------------------------------------------------------------------------------------------------------

void index_initDefault(index_t* index)
{
	index_init(index, -1, -1);
}

void index_init(index_t* index, DECIMAL i, DECIMAL j)
{
	index->i = i;
	index->j = j;
}

index_t* index_new(DECIMAL i, DECIMAL j)
{
	index_t* index = (index_t*) malloc(sizeof(index_t));
	index->i = i;
	index->j = j;
	return index;
}

int index_equals(index_t* left, index_t* right)
{
	if(left->i == right->i && left->j == right->j)
		return TRUE;

	return FALSE;
}

void index_delete(index_t* index)
{
	free(index);
}

// Sub-matrix Dimensions
// -------------------------------------------------------------------------------------------------------

sub_mtx_dim_t* sub_mtx_new(void)
{
	sub_mtx_dim_t* subMtx = (sub_mtx_dim_t*) malloc(sizeof(sub_mtx_dim_t));
	sub_mtx_initDefault(subMtx);
	return subMtx;
}

void sub_mtx_initDefault(sub_mtx_dim_t* subMtx)
{
	index_initDefault(&subMtx->start);
	index_initDefault(&subMtx->length);
}

void sub_mtx_init(sub_mtx_dim_t* subMtx,
		DECIMAL startRow, DECIMAL startCol,
		DECIMAL rowCount, DECIMAL colCount)
{
	index_init(&subMtx->start, startRow, startCol);
	index_init(&subMtx->length, rowCount, colCount);
}

void sub_mtx_copy(sub_mtx_dim_t* source, sub_mtx_dim_t* dest)
{
	if(dest == NULL || source == NULL)
		return;

	index_init(&dest->start, source->start.i, source->start.j);
	index_init(&dest->length, source->length.i, source->length.j);
}

sub_mtx_dim_t* sub_mtx_copyToPtr(sub_mtx_dim_t* subMtx)
{
	sub_mtx_dim_t* copy = sub_mtx_new();
	sub_mtx_copy(subMtx, copy);
	return copy;
}

void sub_mtx_delete(sub_mtx_dim_t* subMtx)
{
	sub_mtx_deleteNonPtr(subMtx);
	free(subMtx);
}

void sub_mtx_deleteNonPtr(sub_mtx_dim_t* subMtx)
{

}

void sub_mtx_toString(sub_mtx_dim_t* subMtx, char* buff_inout)
{
	char temp[DEFAULT_STR_BUFF_SIZE];
	sprintf(temp, "Start_ind: (%d, %d) End_ind: (%d, %d) " \
			"Rows: %d Columns: %d",
			subMtx->start.i, subMtx->start.j,
			subMtx->start.i + subMtx->length.i,
			subMtx->start.j + subMtx->length.j,
			subMtx->length.i,
			subMtx->length.j);
	strcat(buff_inout, temp);
}

void sub_mtx_print(sub_mtx_dim_t* subMtx)
{
	if(subMtx == NULL)
	{
		PRINTF("\n");
		return;
	}

	char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
	sub_mtx_toString(subMtx, temp);
	PRINTF("%s\n", temp);
}

DECIMAL sub_mtx_equals(sub_mtx_dim_t* left, sub_mtx_dim_t* right)
{
	if(index_equals(&left->length, &right->length) &&
			index_equals(&left->start, &right->start))
		return TRUE;

	return FALSE;
}

void sub_mtx_merge(sub_mtx_dim_t* overwritten, sub_mtx_dim_t* other)
{
	overwritten->length.i += other->length.i;
}


/**
 * Function definitions for generic list.
 */
FUNC_DEFINITION_LG_TOARRAY(sub_mtx, sub_mtx_dim_t);
FUNC_DEFINITION_LG_TOARRAY_MULTIPLE(sub_mtx, sub_mtx_dim_t);
FUNC_DEFINITION_LG_COPYADD(sub_mtx, sub_mtx_dim_t, Tail);
FUNC_DEFINITION_LG_COPYADD(sub_mtx, sub_mtx_dim_t, Front);
FUNC_DEFINITION_LG_PRINT(sub_mtx, sub_mtx_dim_t);
FUNC_DEFINITION_LG_PRINT_MULTIPLE(sub_mtx, sub_mtx_dim_t);
FUNC_DEFINITION_LG_DELETEDEEP(sub_mtx, sub_mtx_dim_t);
FUNC_DEFINITION_LG_DELETEDEEP_MULTIPLE(sub_mtx, sub_mtx_dim_t);
FUNC_DEFINITION_LG_PARTITIONDEEP(sub_mtx, sub_mtx_dim_t*, decompose);
FUNC_DEFINITION_LG_PARTITIONDEEP(sub_mtx, sub_mtx_dim_t*, scatter);
FUNC_DEFINITION_LG_PARTITIONDEEP(sub_mtx, sub_mtx_dim_t*, split);


// functions whose implementation varies depending on the storage format
// -------------------------------------------------------------------------------------------------------

// CSR (Compressed Row Storage) Format
// -------------------------------------------------------------------------------------------------------

void sub_mtx_print_CSR(sub_mtx_dim_t* subMtx, spm_cmp_t* spmCsr)
{
	if(subMtx == NULL)
	{
		PRINTF("NULL");
		return;
	}

	char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
	sub_mtx_toString_CSR(subMtx, spmCsr, temp);
	PRINTF("%s\n", temp);
}

void sub_mtx_toString_CSR(sub_mtx_dim_t* subMtx, spm_cmp_t* spmCsr, char* buff_inout)
{
	char temp[DEFAULT_STR_BUFF_SIZE];
	sprintf(temp, "[SUB_MTX_CSR] " \
			"Start_ind: (%d, %d) End_ind: (%d, %d) " \
			"Rows: %d Columns: %d NNZ: %d",
			subMtx->start.i, subMtx->start.j,
			subMtx->start.i + subMtx->length.i,
			subMtx->start.j + subMtx->length.j,
			subMtx->length.i,
			subMtx->length.j,
			sub_mtx_getNNZ_CSR(subMtx, spmCsr));
	strcat(buff_inout, temp);
}

DECIMAL sub_mtx_getNNZ_CSR(sub_mtx_dim_t* subMtx, spm_cmp_t* spmCsr)
{
	DECIMAL startRow = subMtx->start.i;
	DECIMAL endRow = subMtx->start.i + subMtx->length.i;
	DECIMAL startCol = subMtx->start.j;
	DECIMAL endCol = subMtx->start.j + subMtx->length.j;

	DECIMAL nnzCount = 0;
	DECIMAL i;
	DECIMAL j;
	for(i = startRow; i < endRow; ++i)
	{
		for(j = spmCsr->ptr[i]; j < spmCsr->ptr[i + 1]; ++j)
		{
			if(spmCsr->ind[j] >= startCol && spmCsr->ind[j] < endCol)
				++nnzCount;
		}
	}

	return nnzCount;
}

DECIMAL sub_mtx_isEmpty_CSR(sub_mtx_dim_t* subMtx, spm_cmp_t* spmCsr)
{
	if(subMtx->length.i <= 0 || subMtx->length.j <= 0)
		return TRUE;

	if(sub_mtx_getNNZ_CSR(subMtx, spmCsr) <= 0)
		return TRUE;

	return FALSE;
}

// JDS (Jagged Diagonal Storage) Format
// -------------------------------------------------------------------------------------------------------

void sub_mtx_print_JDS(sub_mtx_dim_t* subMtx, spm_jds_t* spmJds)
{
	if(subMtx == NULL)
	{
		PRINTF("NULL");
		return;
	}

	char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
	sub_mtx_toString_JDS(subMtx, spmJds, temp);
	PRINTF("%s\n", temp);
}

void sub_mtx_toString_JDS(sub_mtx_dim_t* subMtx, spm_jds_t* spmJds, char* buff_inout)
{
	char temp[DEFAULT_STR_BUFF_SIZE];
	sprintf(temp, "[SUB_MTX_JDS] " \
			"Start_ind: (%d, %d) End_ind: (%d, %d) " \
			"Rows: %d Columns: %d NNZ: %d",
			subMtx->start.i, subMtx->start.j,
			subMtx->start.i + subMtx->length.i,
			subMtx->start.j + subMtx->length.j,
			subMtx->length.i,
			subMtx->length.j,
			sub_mtx_getNNZ_JDS(subMtx, spmJds));
	strcat(buff_inout, temp);
}

DECIMAL sub_mtx_getNNZ_JDS(sub_mtx_dim_t* subMtx, spm_jds_t* spmJds)
{
	return sub_mtx_getNNZ_CSR(subMtx, spmJds->csrCounterpart);
}

DECIMAL sub_mtx_isEmpty_JDS(sub_mtx_dim_t* subMtx, spm_jds_t* spmJds)
{
	return sub_mtx_isEmpty_CSR(subMtx, spmJds->csrCounterpart);
}

void sub_mtx_extractRowStats_JDS(
		spm_jds_t* spmJds, REAL* maxr_out, REAL* avgr_out, REAL* minr_out)
{
	// count # of non-zeros per row
	// ------------------------------------------------------------------------------------
	DECIMAL i;
	DECIMAL j;
	DECIMAL* rowNNZArr = (DECIMAL*) malloc(sizeof(DECIMAL) * spmJds->rowCount);
	for(i = 0; i < spmJds->rowCount; ++i)
		rowNNZArr[i] = 0;

	for(i = 0; i < spmJds->idiagLength; ++i)
	{
		for(j = spmJds->idiag[i]; j < spmJds->idiag[i + 1]; ++j)
		{
			DECIMAL rowInd = j - spmJds->idiag[i];
			++rowNNZArr[rowInd];
		}
	}

	REAL sum = 0.0;
	for(i = 0; i < spmJds->rowCount; ++i)
	{
		sum += rowNNZArr[i];
	}

	REAL maxr = rowNNZArr[0];
	REAL minr = rowNNZArr[spmJds->rowCount - 1];
	REAL avgr = sum / (REAL) spmJds->rowCount;

	// clean up
	free(rowNNZArr);

	// return values
	*maxr_out = maxr;
	*minr_out = minr;
	*avgr_out = avgr;
}

void sub_mtx_extractColumnStats_JDS(
		spm_jds_t* spmJds, REAL* maxc_out, REAL* avgc_out, REAL* minc_out)
{
	REAL maxc = spmJds->idiag[1];
	REAL minc = spmJds->idiag[spmJds->idiagLength - 1] - spmJds->idiag[spmJds->idiagLength - 2];
	REAL avgc = ((REAL) spmJds->idiag[spmJds->idiagLength - 1]) / ((REAL) spmJds->idiagLength);

	// return values
	*maxc_out = maxc;
	*minc_out = minc;
	*avgc_out = avgc;
}


// (current sub_matrix - JDS ordered sub_matrix) map for JDS format
// -------------------------------------------------------------------------------------------------------

sub_mtx_map_t* sub_mtx_map_new(sub_mtx_dim_t* subMtx, DECIMAL maxNNZLength)
{
	sub_mtx_map_t* subMtxMap = (sub_mtx_map_t*) malloc(sizeof(sub_mtx_map_t));
	sub_mtx_map_init(subMtxMap, subMtx, maxNNZLength);
	return subMtxMap;
}

void sub_mtx_map_init(sub_mtx_map_t* subMtxMap, sub_mtx_dim_t* initialSubMtx, DECIMAL maxNNZLength)
{
	subMtxMap->maxNNZLength = maxNNZLength;
	subMtxMap->afterOrderingSubMtx = NULL;
	subMtxMap->initialSubMtx = sub_mtx_new();
	sub_mtx_copy(initialSubMtx, subMtxMap->initialSubMtx);
}

void sub_mtx_map_deleteNonPtr(sub_mtx_map_t* subMtxMap)
{
	if(subMtxMap->initialSubMtx != NULL)
		sub_mtx_delete(subMtxMap->initialSubMtx);
	if(subMtxMap->afterOrderingSubMtx != NULL)
		sub_mtx_delete(subMtxMap->afterOrderingSubMtx);
}

void sub_mtx_map_delete(sub_mtx_map_t* subMtxMap)
{
	sub_mtx_map_deleteNonPtr(subMtxMap);
	free(subMtxMap);
}

void sub_mtx_map_deleteNonPtrArray(sub_mtx_map_t* subMtxMapArr, int length)
{
	int i;
	for(i = 0; i < length; ++i)
		sub_mtx_map_deleteNonPtr(&subMtxMapArr[i]);

	free(subMtxMapArr);
}

void sub_mtx_map_print(sub_mtx_map_t* subMtxMap)
{
	PRINTF("Sub-matrix-map: NNZ Length=%d\n", subMtxMap->maxNNZLength);
	PRINTF("Initial Sub-matrix:\t\t");
	sub_mtx_print(subMtxMap->initialSubMtx);
	PRINTF("After Ordering sub-matrix:\t");
	sub_mtx_print(subMtxMap->afterOrderingSubMtx);
}

void sub_mtx_map_printArr(sub_mtx_map_t* subMtxMaps, int length)
{
	int i;
	for(i = 0; i < length; ++i)
	{
		sub_mtx_map_print(&subMtxMaps[i]);
	}
}

int sub_mtx_map_cmpInverse(const void* left, const void* right)
{
	sub_mtx_map_t* leftMap = (sub_mtx_map_t*) left;
	sub_mtx_map_t* rightMap = (sub_mtx_map_t*) right;

	if(leftMap->maxNNZLength > rightMap->maxNNZLength)
		return -1;
	else if(leftMap->maxNNZLength < rightMap->maxNNZLength)
		return 1;
	else
		return 0;
}

void sub_mtx_map_sortDescending(sub_mtx_map_t* subMtxMapArr, int length)
{
	qsort(subMtxMapArr, length, sizeof(sub_mtx_map_t), sub_mtx_map_cmpInverse);
}

void sub_mtx_map_extractOrderingForQuintets(
		sub_mtx_map_t* subMtxMapArr, int length, vector_int_t* rowOrderLookup)
{
	int i;
	DECIMAL j;
	DECIMAL currRowIndex = 0;

	for(i = 0; i < length; ++i)
	{
		sub_mtx_map_t* currSubMtxMap = &subMtxMapArr[i];
		sub_mtx_dim_t* initialSubMtx = currSubMtxMap->initialSubMtx;
		sub_mtx_dim_t* afterOrderingSubMtx = sub_mtx_new();

		DECIMAL startRow = initialSubMtx->start.i;
		DECIMAL endRow = initialSubMtx->start.i + initialSubMtx->length.i;

		// here "currRowIndex" is "startingRowIndex"
		// columns and the row length is same for column-net model
		sub_mtx_init(afterOrderingSubMtx,
				currRowIndex, initialSubMtx->start.j,
				initialSubMtx->length.i, initialSubMtx->length.j);

		currSubMtxMap->afterOrderingSubMtx = afterOrderingSubMtx;

		for(j = startRow; j < endRow; ++j)
		{
			// j.th row in current matrix will be the currIndexRow in newly ordered matrix
			rowOrderLookup->data[j] = currRowIndex;
			++currRowIndex;
		}

	}
}

// Sub-matrix Tree
// -------------------------------------------------------------------------------------------------------

// Some Helper functions
// -------------------------------------------------------------------------------------------------------

static void __sub_mtx_tree_getLeafContentsShallow(tree_node_t* node, lg_t* subMtxList);
static void __sub_mtx_tree_getLeafContentsDeep(tree_node_t* node, lg_t* subMtxList);

// -------------------------------------------------------------------------------------------------------

void sub_mtx_tree_init(sub_mtx_tree_t* subMtxNode)
{
	subMtxNode->done = FALSE;
	subMtxNode->subMtx = NULL;
	tree_node_init(&subMtxNode->node);
	omp_init_lock(&subMtxNode->writeLock);

	subMtxNode->borderLock = NULL;
}

sub_mtx_tree_t* sub_mtx_tree_new(DECIMAL startRow, DECIMAL startCol, DECIMAL rowCount, DECIMAL colCount)
{
	sub_mtx_tree_t* subMtxNode = (sub_mtx_tree_t*) malloc(sizeof(sub_mtx_tree_t));

	tree_node_init(&subMtxNode->node);
	subMtxNode->done = FALSE;
	omp_init_lock(&subMtxNode->writeLock);
	subMtxNode->borderLock = NULL;

	subMtxNode->subMtx = (sub_mtx_dim_t*) malloc(sizeof(sub_mtx_dim_t));
	sub_mtx_init(subMtxNode->subMtx, startRow, startCol, rowCount, colCount);

	return subMtxNode;
}

void sub_mtx_tree_delete(tree_node_t* head)
{
	if(head == NULL)
		return;

	sub_mtx_tree_delete(head->left);
	sub_mtx_tree_delete(head->right);

	sub_mtx_tree_t* subMtxNode = sub_mtx_tree_getSubMtxTree(head);
	sub_mtx_tree_deleteSingle(subMtxNode);
}

void sub_mtx_tree_deleteSingle(sub_mtx_tree_t* subMtxNode)
{
	omp_destroy_lock(&subMtxNode->writeLock);
	sub_mtx_delete(subMtxNode->subMtx);
	free(subMtxNode);
}

void sub_mtx_tree_addLeft(sub_mtx_tree_t* parent, sub_mtx_tree_t* left)
{
	if(left == NULL)
		return;

	parent->node.left = &left->node;
	left->node.parent = &parent->node;
}

void sub_mtx_tree_addRight(sub_mtx_tree_t* parent, sub_mtx_tree_t* right)
{
	if(right == NULL)
		return;

	parent->node.right = &right->node;
	right->node.parent = &parent->node;
}

sub_mtx_tree_t* sub_mtx_tree_getRight(sub_mtx_tree_t* parent)
{
	return sub_mtx_tree_getSubMtxTree(parent->node.right);
}

sub_mtx_tree_t* sub_mtx_tree_getLeft(sub_mtx_tree_t* parent)
{
	return sub_mtx_tree_getSubMtxTree(parent->node.left);
}

sub_mtx_tree_t* sub_mtx_tree_getSubMtxTree(tree_node_t* node)
{
	if(node == NULL)
		return NULL;

	return tree_entry(node, sub_mtx_tree_t, node);
}

void sub_mtx_tree_printSingle(sub_mtx_tree_t* subMtxNode)
{
	char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
	sub_mtx_tree_toStringSingle(subMtxNode, temp);
	PRINTF("%s\n", temp);
}

void sub_mtx_tree_toStringSingle(sub_mtx_tree_t* subMtxNode, char* buff_inout)
{
	sub_mtx_toString(subMtxNode->subMtx, buff_inout);

	char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
	sprintf(temp, " DONE: %d", subMtxNode->done);
	strcat(buff_inout, temp);
}

void sub_mtx_tree_printInOrder(tree_node_t* node)
{
	if(node == NULL)
		return;

	sub_mtx_tree_printInOrder(node->left);
	sub_mtx_tree_t* subMtxNode = sub_mtx_tree_getSubMtxTree(node);
	sub_mtx_tree_printSingle(subMtxNode);
	sub_mtx_tree_printInOrder(node->right);
}

void sub_mtx_tree_printPostOrder(tree_node_t* node)
{
	if(node == NULL)
		return;

	sub_mtx_tree_printPostOrder(node->left);
	sub_mtx_tree_printPostOrder(node->right);
	sub_mtx_tree_t* subMtxNode = sub_mtx_tree_getSubMtxTree(node);
	sub_mtx_tree_printSingle(subMtxNode);
}

void sub_mtx_tree_printLeafs(tree_node_t* node)
{
	if(node == NULL)
		return;

	if(tree_node_isLeaf(node))
	{
		sub_mtx_print(sub_mtx_tree_getSubMtxTree(node)->subMtx);
	}
	else
	{
		sub_mtx_tree_printLeafs(node->left);
		sub_mtx_tree_printLeafs(node->right);
	}
}

int sub_mtx_tree_isFull(sub_mtx_tree_t* subMtx)
{
	return tree_node_isFull(&subMtx->node);
}

int sub_mtx_tree_isLeaf(sub_mtx_tree_t* subMtx)
{
	return tree_node_isLeaf(&subMtx->node);
}

int sub_mtx_tree_getLeafCount(sub_mtx_tree_t* subMtx)
{
	if(subMtx == NULL)
		return 0;

	return tree_node_getLeafCount(&subMtx->node);
}

tree_node_t** sub_mtx_tree_getJobDistribution(sub_mtx_tree_t* head, int numBlocks)
{
#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("fetching job distribution numBlocks: %d, total-node-count: %d, leaf-count: %d ...\n",
			numBlocks, tree_node_getNodeCount(&head->node), sub_mtx_tree_getLeafCount(head));
#endif

	tree_node_t** nodesPerBlock = (tree_node_t**) malloc(sizeof(tree_node_t*) * numBlocks);
	int i;
	for(i = 0; i < numBlocks; ++i)
		nodesPerBlock[i] = NULL;

	sub_mtx_tree_markTreeNodes(&head->node, &nodesPerBlock, 0, numBlocks);

	return nodesPerBlock;
}

void sub_mtx_tree_markTreeNodes(tree_node_t* node, tree_node_t*** nodesPerBlock_out, int start, int count)
{
	if(node == NULL)
		return;

	if(tree_node_isLeaf(node) || count <= 1)
	{
		(*nodesPerBlock_out)[start] = node;
	}
	else
	{
		int toLeftCount = count / 2 + count % 2;
		int toRightCount = count / 2;
		sub_mtx_tree_markTreeNodes(node->left, nodesPerBlock_out, start, toLeftCount);
		sub_mtx_tree_markTreeNodes(node->right, nodesPerBlock_out, start + toLeftCount, toRightCount);
	}
}

void sub_mtx_tree_getLeafContentsShallow(sub_mtx_tree_t* subMtxTreeNode, lg_t** subMtxList_out)
{
	lg_t* subMtxList = lg_new();
	__sub_mtx_tree_getLeafContentsShallow(&subMtxTreeNode->node, subMtxList);
	*subMtxList_out = subMtxList;
}

void sub_mtx_tree_getLeafContentsDeep(sub_mtx_tree_t* subMtxTreeNode, lg_t** subMtxList_out)
{
	lg_t* subMtxList = lg_new();
	__sub_mtx_tree_getLeafContentsDeep(&subMtxTreeNode->node, subMtxList);
	*subMtxList_out = subMtxList;
}

// Helper functions
// -------------------------------------------------------------------------------------------------------

static void __sub_mtx_tree_getSubMtxPerBlockDenseInOrder(tree_node_t* head,
		sub_mtx_dim_t** subMtxHistoryPerBlock, int* subMtxCurrIndexPerBlock)
{
	if(head == NULL)
		return;

	__sub_mtx_tree_getSubMtxPerBlockDenseInOrder(head->left, subMtxHistoryPerBlock, subMtxCurrIndexPerBlock);

	sub_mtx_tree_t* subMtxTree = sub_mtx_tree_getSubMtxTree(head);
	int blockIndex = subMtxTree->done - 1;
	sub_mtx_dim_t* blockSubMtxHistory = subMtxHistoryPerBlock[blockIndex];
	sub_mtx_copy(subMtxTree->subMtx, &blockSubMtxHistory[subMtxCurrIndexPerBlock[blockIndex]]);
	++subMtxCurrIndexPerBlock[blockIndex];

	__sub_mtx_tree_getSubMtxPerBlockDenseInOrder(head->right, subMtxHistoryPerBlock, subMtxCurrIndexPerBlock);
}

static void __sub_mtx_tree_getSubMtxPerBlockDensePostOrder(tree_node_t* head,
		sub_mtx_dim_t** subMtxHistoryPerBlock, int* subMtxCurrIndexPerBlock)
{
	if(head == NULL)
		return;

	__sub_mtx_tree_getSubMtxPerBlockDensePostOrder(head->left, subMtxHistoryPerBlock, subMtxCurrIndexPerBlock);
	__sub_mtx_tree_getSubMtxPerBlockDensePostOrder(head->right, subMtxHistoryPerBlock, subMtxCurrIndexPerBlock);

	sub_mtx_tree_t* subMtxTree = sub_mtx_tree_getSubMtxTree(head);
	int blockIndex = subMtxTree->done - 1;
	sub_mtx_dim_t* blockSubMtxHistory = subMtxHistoryPerBlock[blockIndex];
	sub_mtx_copy(subMtxTree->subMtx, &blockSubMtxHistory[subMtxCurrIndexPerBlock[blockIndex]]);
	++subMtxCurrIndexPerBlock[blockIndex];
}

static void __sub_mtx_tree_getNodeCountPerBlock(tree_node_t* head, int* nodeCountPerBlock)
{
	if(head == NULL)
		return;

	sub_mtx_tree_t* subMtxNode = sub_mtx_tree_getSubMtxTree(head);
	++nodeCountPerBlock[subMtxNode->done - 1];

	__sub_mtx_tree_getNodeCountPerBlock(head->left, nodeCountPerBlock);
	__sub_mtx_tree_getNodeCountPerBlock(head->right, nodeCountPerBlock);
}

static void __sub_mtx_tree_getLeafContentsShallow(tree_node_t* node, lg_t* subMtxList)
{
	if(node == NULL)
		return;

	if(tree_node_isLeaf(node))
	{
		sub_mtx_tree_t* subMtxTreeNode = sub_mtx_tree_getSubMtxTree(node);
		lg_addTailData((void*) subMtxTreeNode->subMtx, subMtxList);
	}
	else
	{
		__sub_mtx_tree_getLeafContentsShallow(node->left, subMtxList);
		__sub_mtx_tree_getLeafContentsShallow(node->right, subMtxList);
	}
}

static void __sub_mtx_tree_getLeafContentsDeep(tree_node_t* node, lg_t* subMtxList)
{
	if(node == NULL)
		return;

	if(tree_node_isLeaf(node))
	{
		sub_mtx_dim_t* subMtx = sub_mtx_tree_getSubMtxTree(node)->subMtx;
		sub_mtx_dim_t* copy = sub_mtx_copyToPtr(subMtx);
		lg_addTailData((void*) copy, subMtxList);
	}
	else
	{
		__sub_mtx_tree_getLeafContentsDeep(node->left, subMtxList);
		__sub_mtx_tree_getLeafContentsDeep(node->right, subMtxList);
	}
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
