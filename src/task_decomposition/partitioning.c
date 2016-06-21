
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "include/control_unit/cu_options.h"

#include "include/task_decomposition/partitioning.h"

// -------------------------------------------------------------------------------------------------------

ipd_t* ipd_new(void)
{
	ipd_t* ipd = (ipd_t*) malloc(sizeof(ipd_t));
	ipd->subMtxList = NULL;
	ipd->subMtxHead = NULL;
	return ipd;
}

void ipd_initDefault(ipd_t* ipd)
{
	ipd->subMtxList = NULL;
	ipd->subMtxHead = NULL;
}

void ipd_terminate(ipd_t* ipd)
{
	if(ipd->subMtxList != NULL)
		lg_sub_mtx_deleteDeep(ipd->subMtxList);

	if(ipd->subMtxHead != NULL)
		sub_mtx_tree_delete(&ipd->subMtxHead->node);
}

void ipd_delete(ipd_t* ipd)
{
	ipd_terminate(ipd);
	free(ipd);
}

// Some Helper function declarations
// -------------------------------------------------------------------------------------------------------

static sub_mtx_tree_t* __partitioning_1DRowSliceRecursiveBipartition_CSR(
		spm_cmp_t* spmCsr, REAL targetedCacheSizeInKB, DECIMAL startRowInd, DECIMAL endRowInd,
		REAL (*calculateSizeCallBack)(DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount));

static void __partitioning_calculateSubMtxCountPerBlockPowerOf2(
		int currNumBlocks, DECIMAL currSubMtxCount, DECIMAL* startAddr);

// -------------------------------------------------------------------------------------------------------

void partitioning_1DRowWise(
		spm_cmp_t* spmCsr, int partitionMethod,
		ipd_t* ipd, sub_mtx_dim_t** subMtxArr_out, int* subMtxCount_out,
		REAL targetedCacheSizeInKB, DECIMAL rowIncrementSize, DECIMAL thresholdRowCount,
		REAL (*calculateSizeCallBack)(DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount))
{
	if(partitionMethod == PARTITION_METHOD_REGULAR)
	{
		lg_t* subMtxList = NULL;
		partitioning_1DRowSliceLinear_CSR(
				spmCsr, targetedCacheSizeInKB,
				rowIncrementSize, rowIncrementSize, &subMtxList, calculateSizeCallBack);

		lg_sub_mtx_toArray(subMtxList, subMtxArr_out, subMtxCount_out);
		ipd->subMtxList = subMtxList;
	}
	else if(partitionMethod == PARTITION_METHOD_RECURSIVE_BIPARTITION)
	{
		sub_mtx_tree_t* subMtxHead = NULL;
		lg_t* tempSubMtxList = NULL;

		partitioning_1DRowSliceRecursiveBipartition_CSR(
				spmCsr, targetedCacheSizeInKB, &subMtxHead, calculateSizeCallBack);

		sub_mtx_tree_getLeafContentsDeep(subMtxHead, &tempSubMtxList);
		lg_sub_mtx_toArray(tempSubMtxList, subMtxArr_out, subMtxCount_out);

		ipd->subMtxHead = subMtxHead;
		lg_sub_mtx_deleteDeep(tempSubMtxList);
	}
	else
	{
		PRINTF("Unrecognized partition method: \"%d\".\n", partitionMethod);
		PRINTF("Valid options: %d, %d.",
				PARTITION_METHOD_REGULAR, PARTITION_METHOD_RECURSIVE_BIPARTITION);
		PRINTF("Exiting Program.");
		exit(EXIT_FAILURE);
	}
}

inline void partitioning_1DRowSliceRecursiveBipartition_CSR(
		spm_cmp_t* spmCsr, REAL targetedCacheSizeInKB, sub_mtx_tree_t** subMtxTree_out,
		REAL (*calculateSizeCallBack)(DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount))
{
	*subMtxTree_out = __partitioning_1DRowSliceRecursiveBipartition_CSR(
					spmCsr, targetedCacheSizeInKB, 0, spmCsr->rowCount, calculateSizeCallBack);
}

void partitioning_1DRowSliceLinear_CSR(
		spm_cmp_t* spmCsr, REAL targetedSizeInKB,
		DECIMAL rowIncrementSize, DECIMAL tresholdRowCount, lg_t** subMtxList_out,
		REAL (*calculateSizeCallBack)(DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount))
{
	lg_t* subMtxList = lg_new();

	sub_mtx_dim_t* subMtx = NULL;
	int prevRowInd = 0;
	int currRowInd = 0;
	while(currRowInd < spmCsr->rowCount - rowIncrementSize)
	{
		REAL partSize = calculateSizeCallBack(
				spmCsr->ptr[currRowInd] - spmCsr->ptr[prevRowInd], currRowInd - prevRowInd, spmCsr->colCount);

		if(partSize > targetedSizeInKB)
		{
			// part size has exceeded targeted size so take one less row
			if(currRowInd - prevRowInd > tresholdRowCount)
				currRowInd -= rowIncrementSize;

			// create sub-matrix and append it to list
			subMtx = sub_mtx_new();
			sub_mtx_init(subMtx, prevRowInd, 0, currRowInd - prevRowInd, spmCsr->colCount);
			lg_addTailData((void*) subMtx, subMtxList);

			prevRowInd = currRowInd;
		}

		currRowInd += rowIncrementSize;
	}

	// for one last sub-matrix
	currRowInd -= rowIncrementSize;
	if(currRowInd < spmCsr->rowCount)
	{
		// Left-over sub-matrix is big enough to keep vector-unit full or
		// compared to targeted cache size matrix is really small that it produces only one sub-matrix
		// so, create new sub-matrix.
		if(currRowInd >= tresholdRowCount || subMtx == NULL)
		{
			subMtx = sub_mtx_new();
			sub_mtx_init(subMtx, prevRowInd, 0, spmCsr->rowCount - prevRowInd, spmCsr->colCount);
			lg_addTailData((void*) subMtx, subMtxList);
		}
		// creating new sub-matrix will probably result in inefficient use of vector unit,
		// and it will probably won't cause a huge load imbalance,
		// so, just append remaining to the last sub-matrix .
		else
		{
			subMtx->length.i = spmCsr->rowCount - subMtx->start.i;
		}
	}

	// return values
	*subMtxList_out = subMtxList;
}

void partitioning_1DStaticRowSliceLinear_CSR(
		DECIMAL rowCount, DECIMAL colCount, int rowSliceSize, lg_t** subMtxList_out)
{
	lg_t* subMtxList = lg_new();

	DECIMAL i;
	for(i = rowSliceSize; i < rowCount; i += rowSliceSize)
	{
		sub_mtx_dim_t* subMtx = sub_mtx_new();
		sub_mtx_init(subMtx, i - rowSliceSize, 0, rowSliceSize, colCount);

		lg_addTailData((void*) subMtx, subMtxList);
	}

	i = i - rowSliceSize;
	sub_mtx_dim_t* subMtx = sub_mtx_new();
	sub_mtx_init(subMtx, i, 0, rowCount - i, colCount);
	lg_addTailData((void*) subMtx, subMtxList);

	// return values
	*subMtxList_out = subMtxList;
}

REAL partitioning_calculateSizeInKBDefault(DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount)
{
	return (sizeof(DECIMAL) * nnz 						// size of colInd array in CSR format
			+ sizeof(REAL) * nnz 						// size of values array in CSR format
			+ sizeof(DECIMAL) * rowCount 				// size of rowPtr array in CSR format
			+ sizeof(REAL) * rowCount 					// size of output vector
			// + sizeof(REAL) * rowCount 				// size of input vector (estimated)
			) / pow(2, 10); 							// convert bytes to KB
}

// TODO test
DECIMAL* partitioning_calculateSubMtxCountPerBlockPowerOf2(int numBlocks, DECIMAL subMtxCount)
{
	DECIMAL* subMtxCountPerBlock = (DECIMAL*) malloc(sizeof(DECIMAL) * numBlocks);
	__partitioning_calculateSubMtxCountPerBlockPowerOf2(numBlocks, subMtxCount, subMtxCountPerBlock);

	return subMtxCountPerBlock;
}

partitioning_vector_t* partitioning_generatePartitioningFromPartVector(vector_int_t* partVector, int partCount)
{
	partitioning_vector_t* pv = vector_int_new(partCount + 1);
	pv->data[0] = 0;

	int i;
	for(i = 0; i < partVector->length; ++i)
	{
		++pv->data[partVector->data[i] + 1];
	}

	for(i = 0; i < partCount; ++i)
	{
		pv->data[i + 1] += pv->data[i];
	}

	return pv;
}

vector_int_t* partitioning_generateOrderingFromPartVector(vector_int_t* partVector, int partCount)
{
	int* partSizes = (int*) calloc(partCount, sizeof(int));
	int i;

	for(i = 0; i < partVector->length; ++i)
	{
		++partSizes[partVector->data[i]];
	}

	vector_int_t** partOrderings = (vector_int_t**) malloc(sizeof(vector_int_t*) * partCount);
	int* currIndices = (int*) malloc(sizeof(int) * partCount);
	for(i = 0; i < partCount; ++i)
	{
		partOrderings[i] = vector_int_new(partSizes[i]);
		currIndices[i] = 0;
	}

	for(i = 0; i < partVector->length; ++i)
	{
		int partIndex = partVector->data[i];
		partOrderings[partIndex]->data[currIndices[partIndex]] = i;
		++currIndices[partIndex];
	}

	vector_int_t* temp = vector_int_new(partVector->length);
	i = 0;
	while(i < temp->length)
	{
		int p;
		for(p = 0; p < partCount; ++p)
		{
			int j;
			for(j = 0; j < partOrderings[p]->length; ++j)
			{
				temp->data[i++] = partOrderings[p]->data[j];
			}
		}
	}

	vector_int_t* ordering = vector_int_new(partVector->length);
	i = 0;
	while(i < ordering->length)
	{
		ordering->data[temp->data[i]] = i;
		++i;
	}

	// clean up
	for(i = 0; i < partCount; ++i)
	{
		vector_int_delete(partOrderings[i]);
	}
	free(partOrderings);
	free(currIndices);
	free(partSizes);
	vector_int_delete(temp);

	return ordering;
}


// Helper functions
// -------------------------------------------------------------------------------------------------------

// TODO test
void __partitioning_calculateSubMtxCountPerBlockPowerOf2(
		int currNumBlocks, DECIMAL currSubMtxCount, DECIMAL* startAddr)
{
	if(currNumBlocks < 1)
		return;

	if(currSubMtxCount == 1)
	{
		startAddr[0] = currSubMtxCount;
		int i;
		for(i = 1; i < currNumBlocks; ++i)
			startAddr[i] = 0;
	}
	else if(currNumBlocks == 1)
	{
		startAddr[0] = currSubMtxCount;
	}
	else
	{
		DECIMAL nextSubMtxCount = currSubMtxCount / 2;

		int leftBlockCount = currNumBlocks / 2;
		__partitioning_calculateSubMtxCountPerBlockPowerOf2(leftBlockCount, nextSubMtxCount, startAddr);

		int rightBlockCount = currNumBlocks / 2 + currNumBlocks % 2;
		__partitioning_calculateSubMtxCountPerBlockPowerOf2(rightBlockCount, nextSubMtxCount, &startAddr[leftBlockCount]);
	}
}

sub_mtx_tree_t* __partitioning_1DRowSliceRecursiveBipartition_CSR(
		spm_cmp_t* spmCsr, REAL targetedCacheSizeInKB, DECIMAL startRowInd, DECIMAL endRowInd,
		REAL (*calculateSizeCallBack)(DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount))
{
	sub_mtx_tree_t* parent = sub_mtx_tree_new(startRowInd, 0, endRowInd - startRowInd, spmCsr->colCount);

	DECIMAL nnz = spmCsr->ptr[endRowInd] - spmCsr->ptr[startRowInd];

	// calculated cache size
	REAL calcCacheSizeInKB = calculateSizeCallBack(nnz, endRowInd - startRowInd, spmCsr->colCount);

	// if sub-matrix size is bigger than targeted cache size recursively divide matrix into smaller chunks
	// else just return current sub-matrix to be inserted into partitioning tree as a leaf node.
	// also stop partitioning if there are no more than 1 rows left
	if((calcCacheSizeInKB > targetedCacheSizeInKB) && ((endRowInd - startRowInd) > 1))
	// if((calcCacheSizeInKB > targetedCacheSizeInKB) && ((endRowInd - startRowInd) > SIMD_LENGTH))
	{
		double lowerSize = 0.0;
		double prevLowerSize = lowerSize;
		double upperSize = calcCacheSizeInKB;
		double prevUpperSize = upperSize;
		double lowerNNZ = 0;
		double upperNNZ = spmCsr->ptr[endRowInd] - spmCsr->ptr[startRowInd];
		DECIMAL midRow = startRowInd;

		while(upperSize > lowerSize &&
				midRow < endRowInd &&
				endRowInd - midRow > 1)
				// endRowInd - midRow > (SIMD_LENGTH))
		{
			prevLowerSize = lowerSize;
			prevUpperSize = upperSize;

			++midRow;
			lowerNNZ += spmCsr->ptr[midRow] - spmCsr->ptr[midRow - 1];
			upperNNZ -= spmCsr->ptr[midRow] - spmCsr->ptr[midRow - 1];
			// midRow += SIMD_LENGTH;
			// lowerNNZ += spmCsr->ptr[midRow] - spmCsr->ptr[midRow - SIMD_LENGTH];
			// upperNNZ -= spmCsr->ptr[midRow] - spmCsr->ptr[midRow - SIMD_LENGTH];

			lowerSize = calculateSizeCallBack(lowerNNZ, midRow - startRowInd, spmCsr->colCount);

			upperSize = calcCacheSizeInKB - lowerSize;
		}

		// Find the most balanced partitioning possible
		if(fabs(prevLowerSize - prevUpperSize) < fabs(lowerSize - upperSize))
			--midRow;

		// recursion until calculated sub-mtx size drops below targeted cache size
		sub_mtx_tree_t* leftChild = __partitioning_1DRowSliceRecursiveBipartition_CSR(
				spmCsr, targetedCacheSizeInKB, startRowInd, midRow, calculateSizeCallBack);
		sub_mtx_tree_t* rightChild = __partitioning_1DRowSliceRecursiveBipartition_CSR(
				spmCsr, targetedCacheSizeInKB, midRow, endRowInd, calculateSizeCallBack);

		// add childsTo subMtxTreeNode
		sub_mtx_tree_addLeft(parent, leftChild);
		sub_mtx_tree_addRight(parent, rightChild);
	}

	return parent;
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
