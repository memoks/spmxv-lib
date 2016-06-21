
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>

#include "include/algorithm/algorithm.h"
#include "include/util/utility.h"

#include "include/data_structure/quintet.h"


// jds_row functions
// ------------------------------------------------------------------------------------------

int jds_row_cmpInverse(const void* left, const void* right)
{
	jds_row_t* j1 = (jds_row_t*) left;
	jds_row_t* j2 = (jds_row_t*) right;

	// sort in descending order
	if(j1->nnz > j2->nnz)
		return -1;
	else if(j1->nnz < j2->nnz)
		return 1;
	else
		return 0;
}

void jds_row_sortDescending(jds_row_t* jdsRows, DECIMAL length)
{
	qsort(jdsRows, length, sizeof(jds_row_t), jds_row_cmpInverse);
}

void jds_row_accumulate(jds_row_t* jdsRows, DECIMAL length)
{
	int i;
	for(i = 1; i < length; ++i)
		jdsRows[i + 1].nnz += jdsRows[i].nnz;
}

void jds_row_deaccumulate(jds_row_t* jdsRows, DECIMAL length)
{
	int i;
	for(i = length - 1; i > 0; --i)
		jdsRows[i].nnz -= jdsRows[i - 1].nnz;
}

void jds_row_printArr(jds_row_t* jdsRows, DECIMAL length)
{
	DECIMAL i;
	for(i = 0; i < length; ++i)
		PRINTF("rowId: %d => nnz: %d\n", jdsRows[i].rowId, jdsRows[i].nnz);
}

// quintet functions
// ------------------------------------------------------------------------------------------

int quintet_cmpColIndex(const void* left, const void* right)
{
	quintet_t* t1 = (quintet_t*) left;
	quintet_t* t2 = (quintet_t*) right;

	if(t1->j > t2->j)
		return 1;
	else if(t1->j < t2->j)
		return -1;
	else
	{
		if(t1->i > t2->i)
			return 1;
		else if(t1->i < t2->i)
			return -1;

		return 0;
	}
}

int quintet_cmpColValue(const void* left, const void* right)
{
	quintet_t* t1 = (quintet_t*) left;
	quintet_t* t2 = (quintet_t*) right;

	if(t1->jVal > t2->jVal)
		return 1;
	else if(t1->jVal < t2->jVal)
		return -1;
	else
	{
		if(t1->iVal > t2->iVal)
			return 1;
		else if(t1->iVal < t2->iVal)
			return -1;

		return 0;
	}
}

int quintet_cmpRowIndex(const void* left, const void* right)
{
	quintet_t* t1 = (quintet_t*) left;
	quintet_t* t2 = (quintet_t*) right;

	if(t1->i > t2->i)
		return 1;
	else if(t1->i < t2->i)
		return -1;
	else
	{
		if(t1->j > t2->j)
			return 1;
		else if(t1->j < t2->j)
			return -1;

		return 0;
	}
}

int quintet_cmpRowValue(const void* left, const void* right)
{
	quintet_t* t1 = (quintet_t*) left;
	quintet_t* t2 = (quintet_t*) right;

	if(t1->iVal > t2->iVal)
		return 1;
	else if(t1->iVal < t2->iVal)
		return -1;
	else
	{
		if(t1->jVal > t2->jVal)
			return 1;
		else if(t1->jVal < t2->jVal)
			return -1;

		return 0;
	}
}

void quintet_sortColIndex(quintet_t* quintets, DECIMAL length)
{
	qsort(quintets, length, sizeof(quintet_t), quintet_cmpColIndex);
}

void quintet_sortColValue(quintet_t* quintets, DECIMAL length)
{
	qsort(quintets, length, sizeof(quintet_t), quintet_cmpColValue);
}

void quintet_sortRowIndex(quintet_t* quintets, DECIMAL length)
{
	qsort(quintets, length, sizeof(quintet_t), quintet_cmpRowIndex);
}

void quintet_sortRowValue(quintet_t* quintets, DECIMAL length)
{
	qsort(quintets, length, sizeof(quintet_t), quintet_cmpRowValue);
}

void quintet_sortCustom(
		quintet_t* quintets, DECIMAL length,
		int (*input_cmpQuintetFunc) (const void*, const void*))
{
	qsort(quintets, length, sizeof(quintet_t), input_cmpQuintetFunc);
}


void quintet_overwriteIndex(quintet_t* quintets_inout, DECIMAL nnz)
{
	DECIMAL i;
	for(i = 0; i < nnz; ++i)
	{
		quintets_inout[i].i = quintets_inout[i].iVal;
		quintets_inout[i].j = quintets_inout[i].jVal;
	}
}

void quintet_generateCSRRowPtr(quintet_t* quintets, DECIMAL length, DECIMAL** rowPtr_out, DECIMAL* rowPtrLength_out)
{
	DECIMAL rowCount = quintets[length - 1].iVal - quintets[0].iVal + 1;
	DECIMAL* rowPtr = (DECIMAL*) malloc(sizeof(DECIMAL) * (rowCount + 1));

	DECIMAL prevRowInd = -1;
	DECIMAL i = 0;

	while(i < length)
	{
		quintet_t* currQuintet = &quintets[i];
		while(currQuintet->iVal > prevRowInd)
		{
			++prevRowInd;
			rowPtr[prevRowInd] = i;
		}

		++i;
	}

	rowPtr[rowCount] = i;

	// return values
	*rowPtr_out = rowPtr;
	*rowPtrLength_out = rowCount + 1;
}

void quintet_countElementsPerRow(
		quintet_t* quintets, DECIMAL length, DECIMAL rowCount, DECIMAL** elemCountPerRowArr_out)
{
	DECIMAL* elemCountPerRowArr = (DECIMAL*) calloc(rowCount, sizeof(DECIMAL));

	DECIMAL i = 0;
	while(i < length)
	{
		quintet_t* currQuintet = &quintets[i];
		++elemCountPerRowArr[currQuintet->iVal];
		++i;
	}

	// return values
	*elemCountPerRowArr_out = elemCountPerRowArr;
}


void quintet_adjustBinaryMatrixCoefficients(quintet_t* quintets_inout, int numQuintets)
{
	int i;
	int j;
	int prevRow = 0;
	int currRow = 0;
	int nnzCountPrevRow = 0;
	for(i = 0; i < numQuintets; ++i)
	{
		// ASSUMPTION: when reading quintets for the first time
		// their original index must be copied to index values
		// event though there is no ordering scheme used.
		// As a result iVal and jVal must never be empty.
		currRow = quintets_inout[i].iVal;

		if(currRow != prevRow)
		{
			REAL rowValue = (float) 1 / (float) nnzCountPrevRow;
			for(j = 0; j < nnzCountPrevRow; ++j)
				quintets_inout[i - j - 1].val = rowValue;
			nnzCountPrevRow = 0;
			prevRow = currRow;
		}

		++nnzCountPrevRow;
	}

	float rowValue = (float) 1 / (float) nnzCountPrevRow;
	for(j = 0; j < nnzCountPrevRow; ++j)
		quintets_inout[i - j - 1].val = rowValue;
	nnzCountPrevRow = 0;
}

void quintet_transpose(quintet_t* quintets_inout, DECIMAL nnz,
		DECIMAL* rowCount_inout, DECIMAL* colCount_inout)
{
	DECIMAL i;
	DECIMAL temp;
	for(i = 0; i < nnz; ++i)
	{
		quintet_t* currTrip = &quintets_inout[i];
		temp = currTrip->i;
		currTrip->i = currTrip->j;
		currTrip->j = temp;
	}

	temp = *rowCount_inout;
	*rowCount_inout = *colCount_inout;
	*colCount_inout = temp;
}


void quintet_addOrderingInfoByIndex(
		quintet_t* quintets_inout, DECIMAL nnz, vector_int_t* rowOrderLookup, vector_int_t* columnOrderLookup)
{
	DECIMAL i;

	if(rowOrderLookup != NULL)
	{
		for(i = 0; i < nnz; ++i)
		{
			quintets_inout[i].iVal = rowOrderLookup->data[quintets_inout[i].i];
		}
	}

	if(columnOrderLookup != NULL)
	{
		for(i = 0; i < nnz; ++i)
			quintets_inout[i].jVal = columnOrderLookup->data[quintets_inout[i].j];
	}

}

void quintet_addRowOrderingInfoByValue(quintet_t* quintets_inout, DECIMAL nnz, vector_int_t* rowOrderLookup)
{
	DECIMAL i;
	for(i = 0; i < nnz; ++i)
	{
		quintets_inout[i].iVal = rowOrderLookup->data[quintets_inout[i].iVal];
	}
}

void quintet_copyDefaultOrderingInfoToValue(quintet_t* quintets_inout, DECIMAL nnz)
{
	DECIMAL i;
	for(i = 0; i < nnz; ++i)
	{
		quintets_inout[i].iVal = quintets_inout[i].i;
		quintets_inout[i].jVal = quintets_inout[i].j;
	}
}

void quintet_printQuintets(char* m, quintet_t* quintets, int length)
{
	PRINTF("%s\n", m);
	int i;
	for(i = 0; i < length; ++i)
		PRINTF("%d => (%d, %d) (%d, %d) > %3.2f\n",
				i, quintets[i].i, quintets[i].j, quintets[i].iVal, quintets[i].jVal, quintets[i].val);
	PRINTF("\n");
}

void quintet_copyQuintet(quintet_t* source, quintet_t* destination)
{
	destination->i = source->i;
	destination->j = source->j;
	destination->val = source->val;
	destination->iVal = source->iVal;
	destination->jVal = source->jVal;
}

void quintet_copyQuintets(quintet_t* source, quintet_t* destination, DECIMAL length)
{
	DECIMAL i;
	for(i = 0; i < length; ++i)
	{
		quintet_copyQuintet(&source[i], &destination[i]);
	}
}

void quintet_orderForJDSCounterpart(
		quintet_t** quintets_inout, DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount,
		sub_mtx_dim_t* subMtxArr, int subMtxArrLength,
		vector_int_t** rowOrderLookup_out, DECIMAL** elemCountPerRow_out)
{
	quintet_t* quintets = *quintets_inout;
	DECIMAL* elemCountPerRow = (DECIMAL*) malloc(sizeof(DECIMAL) * rowCount);
	DECIMAL i;

	// Count NNZ per row to order rows later
	// ---------------------------------------------------------------------------------------

	// NNZ count per row
	// BEWARE of the first element of this array (will be used for accumulation / deaccumulation)
	jds_row_t* jdsRowAcc = (jds_row_t*) malloc(sizeof(jds_row_t) * (rowCount + 1));
	jdsRowAcc[0].nnz = 0;
	jdsRowAcc[0].rowId = -1;
	jds_row_t* jdsRows = &jdsRowAcc[1];

	// ordering info for SPM
	vector_int_t* rowOrderLookup = vector_int_new(rowCount);

	for(i = 0; i < rowCount; ++i)
	{
		jdsRows[i].nnz = 0;
		jdsRows[i].rowId = i;
		rowOrderLookup->data[i] = i;
	}

	// count NNZ per row
	for(i = 0; i < nnz; ++i)
		++jdsRows[quintets[i].i].nnz;

	// Sort rows by their nnz count
	int numThreads = omp_get_num_threads();
	int* subMtxCountAcc = (int*) malloc(sizeof(int) * (numThreads + 1));
	subMtxCountAcc[0] = 0;

	// sort quintets for CSR first
	algorithm_parallelMergeQuickSort(quintets, nnz, numThreads, quintet_cmpRowValue);

	#pragma omp parallel num_threads(numThreads)
	{
		int i;
		int j;

		// handle sub-matrix distribution
		// -----------------------------------------------------------------------------------------
		int threadId = omp_get_thread_num();
		int threadSubMtxCount = subMtxArrLength / numThreads;
		if(threadId < subMtxArrLength % numThreads)
			++threadSubMtxCount;

		subMtxCountAcc[threadId + 1] = threadSubMtxCount;

		// Sort rows (in a sub-matrix) by their nnz count so that the one with most nnz is first row
		// -----------------------------------------------------------------------------------------
		for(i = subMtxCountAcc[threadId]; i < subMtxCountAcc[threadId + 1]; ++i)
		{
			sub_mtx_dim_t* currSubMtx = &subMtxArr[i];

			jds_row_sortDescending(&jdsRows[currSubMtx->start.i], currSubMtx->length.i);

			// create permutation array for JDS format
			for(j = currSubMtx->start.i; j < currSubMtx->start.i + currSubMtx->length.i; ++j)
				rowOrderLookup->data[jdsRows[j].rowId] = j;
		}

		#pragma omp barrier
		#pragma omp single
		{
			accumulate(subMtxCountAcc, numThreads + 1);
			jds_row_accumulate(jdsRowAcc, rowCount + 1);
		}
		#pragma omp barrier

		// "Partially" sort quintets in sub-matrix
		// -----------------------------------------------------------------------------------------
		for(i = subMtxCountAcc[threadId]; i < subMtxCountAcc[threadId + 1]; ++i)
		{
			sub_mtx_dim_t* currSubMtx = &subMtxArr[i];
			DECIMAL subMtxQuintetStart = jdsRowAcc[currSubMtx->start.i].nnz;
			DECIMAL subMtxQuintetLength =
					jdsRowAcc[currSubMtx->start.i + currSubMtx->length.i].nnz - subMtxQuintetStart;

			quintet_sortRowValue(&quintets[subMtxQuintetStart], subMtxQuintetLength);

			for(j = currSubMtx->start.i; j < currSubMtx->start.i + currSubMtx->length.i; ++i)
				elemCountPerRow[j] = jdsRowAcc[j + 1].nnz - jdsRowAcc[j].nnz;
		}
	}

	// clean up
	free(subMtxCountAcc);
	free(jdsRowAcc);

	// return values
	*rowOrderLookup_out = rowOrderLookup;
	*elemCountPerRow_out = elemCountPerRow;
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
