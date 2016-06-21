
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "include/algorithm/algorithm.h"
#include "include/io/converter.h"

// Helper functions
// ---------------------------------------------------------------------------------------------------------------

extern DECIMAL __converter_CSR_to_ICSR_calculatePtrLength(DECIMAL rowCount, DECIMAL* ptr);
extern void __converter_CSR_to_ICSR_fillPtrArr(DECIMAL rowCount, DECIMAL* ptr, DECIMAL* incPtr);
extern void __converter_copyReal(REAL* source, REAL* destination, DECIMAL length);
extern void __converter_CSR_to_ICSR_fillIndArr(DECIMAL rowCount, DECIMAL colCount, DECIMAL* ptr, DECIMAL* ind, DECIMAL* incInd);
extern void __converter_smoothCSR(spm_cmp_t* csr);


// ---------------------------------------------------------------------------------------------------------------

void converter_CSR_to_ICSR(spm_cmp_t* spmCsr, spm_inc_t** spmIcsr_out)
{
	if(spmCsr == NULL)
		return;

	DECIMAL ptrLength = __converter_CSR_to_ICSR_calculatePtrLength(spmCsr->rowCount, spmCsr->ptr);

	spm_inc_t* spmIcsr = spm_inc_new(spmCsr->nnz, spmCsr->rowCount, spmCsr->colCount, ptrLength);

	__converter_CSR_to_ICSR_fillPtrArr(spmCsr->rowCount, spmCsr->ptr, spmIcsr->ptr);
	__converter_CSR_to_ICSR_fillIndArr(spmCsr->rowCount, spmCsr->colCount, spmCsr->ptr, spmCsr->ind, spmIcsr->ind);
	__converter_copyReal(spmCsr->values, spmIcsr->values, spmIcsr->nnz);

	// return values
	*spmIcsr_out = spmIcsr;
}

void converter_CSR_to_ICSR_multiple(spm_cmp_t** spmCsrs, spm_inc_t*** spmIcsrs_out, int numSpms)
{
	spm_inc_t** spmIcsrs = (spm_inc_t**) malloc(sizeof(spm_inc_t*) * numSpms);

	#pragma omp parallel num_threads(numSpms)
	{
		int blockId = omp_get_thread_num();
		spm_inc_t* spmIcsr = NULL;
		converter_CSR_to_ICSR(spmCsrs[blockId], &spmIcsr);

		spmIcsrs[blockId] = spmIcsr;
	}

	// return values
	*spmIcsrs_out = spmIcsrs;
}

// TODO not used, test
void converter_ICSR_to_CSR(spm_inc_t* spmIcsr, spm_cmp_t** spmCsr_out)
{
	if(spmIcsr->nnz == 0)
		return;

	spm_cmp_t* spmCsr = spm_cmp_new(spmIcsr->rowCount, spmIcsr->colCount, spmIcsr->nnz);
	spmCsr->ptr[0] = 0;

	DECIMAL i;
	DECIMAL j;
	DECIMAL acc = 0;
	DECIMAL buffIndex = 0;
	DECIMAL* buffer = (DECIMAL*) malloc(sizeof(DECIMAL) * spmIcsr->ptrLength);
	for(i = 0; i < spmIcsr->nnz; ++i)
	{
		acc += spmIcsr->ind[i];
		if(acc >= spmIcsr->colCount)
		{
			buffer[buffIndex] = i;
			++buffIndex;
			acc = acc % spmIcsr->colCount;
		}

		spmCsr->ind[i] = acc;
	}

	DECIMAL ptrInd = 1;
	for(i = 1; i < spmIcsr->ptrLength; ++i)
	{
		for(j = 0; j < spmIcsr->ptr[i]; ++j)
		{
			spmCsr->ptr[ptrInd] = buffer[i - 1];
			++ptrInd;
		}
	}

	__converter_copyReal(spmIcsr->values, spmCsr->values, spmCsr->nnz);

	// local clean up
	free(buffer);

	// return values
	*spmCsr_out = spmCsr;
}

// quintet conversion functions
// ---------------------------------------------------------------------------------------------------------

void converter_quintetToCSR_alloc(
		quintet_t* quintets, spm_cmp_t** spmCsr_out,
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount)
{
	spm_cmp_t* spmCsr = spm_cmp_new(rowCount, colCount, nnz);
	converter_quintetToCSR(quintets, spmCsr);
	*spmCsr_out = spmCsr;
}

void converter_quintetToCSR(quintet_t* quintets, spm_cmp_t* spmCsr_inout)
{
	DECIMAL i;
	DECIMAL j;
	DECIMAL row = 0;
	DECIMAL prevRow = -1;
	DECIMAL rowIndex = 0;
	for(i = 0; i < spmCsr_inout->nnz; ++i)
	{
		// ASSUMPTION: when reading quintets for the first time
		// their original index must be copied to index values
		// even though there is no ordering scheme used.
		// As a result iVal and jVal must never be empty.
		row = quintets[i].iVal;

		if(row != prevRow)
		{
			// in case there are empty rows
			for(j = prevRow; j < row; ++j)
			{
				spmCsr_inout->ptr[rowIndex] = i;
				++rowIndex;
			}

			prevRow = row;
		}

		// ASSUMPTION: when reading quintets for the first time
		// their original index must be copied to index values
		// even though there is no ordering scheme used.
		// As a result iVal and jVal must never be empty.
		spmCsr_inout->ind[i] = quintets[i].jVal;

		spmCsr_inout->values[i] = quintets[i].val;
	}

	spmCsr_inout->ptr[rowIndex] = i;
}

void converter_extractCsrNnzPerRow(quintet_t* quintets, DECIMAL length, DECIMAL* nnzPerRow_inout)
{
	if(length < 0)
		return;

	DECIMAL i = 0;
	DECIMAL ptrIndex = 0;
	DECIMAL currRow = quintets[0].iVal;
	while(i < length)
	{
		if(quintets[i].iVal != currRow)
		{
			++currRow;
			++ptrIndex;
			continue;
		}

		++nnzPerRow_inout[ptrIndex];
		++i;
	}
}

void converter_extractJdsNnzPerColumn(quintet_t* quintets, DECIMAL length, DECIMAL* nnzPerColumn_inout)
{
	if(length < 0)
		return;

	DECIMAL i = 0;
	DECIMAL colIndex = 0;
	DECIMAL prevRow = 0;
	while(i < length)
	{
		if(prevRow >= quintets[i].iVal)
		{
			++colIndex;
		}

		++nnzPerColumn_inout[colIndex];
		++i;
	}
}

extern void converter_quintetToJDS_alloc(
		quintet_t* quintets, spm_jds_t** spmJds_out,
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount, DECIMAL idiagLength)
{
	spm_jds_t* spmJds = spm_jds_new(nnz, rowCount, colCount, idiagLength);
	converter_quintetToJDS(quintets, spmJds);
	*spmJds_out = spmJds;
}

extern void converter_quintetToJDS(quintet_t* quintets, spm_jds_t* spmJds_inout)
{
	DECIMAL i = 0;
	DECIMAL j;
	DECIMAL prevRowIndex = 0;
	DECIMAL idiagIndex = 0;

	spmJds_inout->idiag[idiagIndex] = 0;
	while(i < spmJds_inout->nnz)
	{
		quintet_t* currQuintet = &quintets[i];
		while(prevRowIndex <= currQuintet->iVal)
		{
			spmJds_inout->jdiag[i] = currQuintet->jVal;
			spmJds_inout->dj[i] = currQuintet->val;

			++i;
			if(i >= spmJds_inout->nnz)
				break;

			currQuintet = &quintets[i];
		}

		++idiagIndex;
		spmJds_inout->idiag[idiagIndex] = i;
	}
}

void converter_accumulate(DECIMAL* arr_inout, DECIMAL length)
{
	DECIMAL i;
	for(i = 1; i < length; ++i)
		arr_inout[i] += arr_inout[i - 1];

}

void converter_quintetToCSC(
		quintet_t* quintets, spm_cmp_t** spmCsc_out,
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount)
{
	spm_cmp_t* spmCsc = spm_cmp_new(colCount, rowCount, nnz);

	int i;
	int j;
	int col = 0;
	int prevCol = -1;
	int colIndex = 0;
	for(i = 0; i < nnz; ++i)
	{
		// ASSUMPTION: when reading quintets for the first time
		// their original index must be copied to index values
		// even though there is no ordering scheme used.
		// As a result iVal and jVal must never be empty.
		if(quintets[i].jVal < 0)
			col = quintets[i].j;
		else
			col = quintets[i].jVal;

		if(col != prevCol)
		{
			// in case there are empty columns
			for(j = prevCol; j < col; ++j)
			{
				spmCsc->ptr[colIndex] = i;
				++colIndex;
			}

			prevCol = col;
		}

		// ASSUMPTION: when reading quintets for the first time
		// their original index must be copied to index values
		// even though there is no ordering scheme used.
		// As a result iVal and jVal must never be empty.
		spmCsc->ind[i] = quintets[i].iVal;

		spmCsc->values[i] = quintets[i].val;
	}

	spmCsc->ptr[colIndex] = i;

	// return values
	*spmCsc_out = spmCsc;
}

void converter_quintetToJDSOrderedQuintet(
		quintet_t* quintets,DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount,
		sub_mtx_dim_t* subMtxArr, int subMtxCount, vector_int_t** rowOrderLookup_out)
{
	DECIMAL i;

	// Count NNZ per row to order rows later
	// ---------------------------------------------------------------------------------------

	// NNZ count per row
	jds_row_t* jdsRows = (jds_row_t*) malloc(sizeof(jds_row_t) * rowCount);

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("Allocated jds_row_t array.\n");
#endif

	// ordering info for SPM
	vector_int_t* rowOrderLookup = vector_int_new(rowCount);

	for(i = 0; i < rowCount; ++i)
	{
		jdsRows[i].nnz = 0;
		jdsRows[i].rowId = i;
		rowOrderLookup->data[i] = i;
	}

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("Initialized jds_row_t array.\n");
#endif

	// count NNZ per row
	for(i = 0; i < nnz; ++i)
		jdsRows[quintets[i].i].nnz += 1;

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("Counted NNZ per row.\n");
	PRINTF("\nTotal NNZ: %d\n", nnz);
#endif

	// Calculate zero padded NNZ length
	// partially sort JDS rows (do not sort quintets yet)
	// ---------------------------------------------------------------------------------------

	// non-padded usual CSR
	spm_cmp_t* spmCsr = NULL;
	converter_quintetToCSR_alloc(quintets, &spmCsr, nnz, rowCount, colCount);

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("Created CSR.\n");
#endif

	for(i = 0; i < subMtxCount; ++i)
	{
		sub_mtx_dim_t* currSubMtx = &subMtxArr[i];

		// partially sort JDS rows (by NNZ elements per row in descending order)
		jds_row_sortDescending(&jdsRows[currSubMtx->start.i], currSubMtx->length.i);
	}

	// extract row ordering information from JDS-rows
	// create permutation array for JDS format
	// sort quintets
	// create new CSR counterpart in which quintets are ordered for JDS format
	// ---------------------------------------------------------------------------------------

	for(i = 0; i < rowCount; ++i)
		rowOrderLookup->data[jdsRows[i].rowId] = i;

	quintet_addOrderingInfoByIndex(quintets, nnz, rowOrderLookup, NULL);
	// TODO using self-implemented hybrid sort
	// input_sortQuintetRowValue(quintets, nnz);
	// TODO changed to parallel partial sort since quick sort of global array (which is almost sorted) is inefficient
	// algorithm_parallelMergeQuickSort(quintets, nnz, omp_get_num_threads(), input_cmpQuintetRowValue);
	// Changed to parallel partial quick sort at last
	#pragma omp parallel for private(i)
	for(i = 0; i < subMtxCount; ++i)
	{
		sub_mtx_dim_t* currSubMtx = &subMtxArr[i];
		DECIMAL jobBatchNNZStartIndex = spmCsr->ptr[currSubMtx->start.i];
		DECIMAL jobBatchNNZLength = spmCsr->ptr[currSubMtx->start.i + currSubMtx->length.i] - jobBatchNNZStartIndex;

		quintet_sortRowValue(&quintets[jobBatchNNZStartIndex], jobBatchNNZLength);
	}


	spm_cmp_delete(spmCsr);
	converter_quintetToCSR_alloc(quintets, &spmCsr, nnz, rowCount, colCount);

	// clean up
	spm_cmp_delete(spmCsr);
	free(jdsRows);

	// return values
	*rowOrderLookup_out = rowOrderLookup;
}

void converter_quintetToJDSCounterpartCSR(
		quintet_t* quintets, DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount,
		sub_mtx_dim_t* subMtxArr, int subMtxCount,
		spm_cmp_t** spmCsrCounterpart_out, DECIMAL** permutation_out,
		vector_int_t** rowOrderLookup_out)
{
#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("Converting quintets to JDS...\n");
#endif

	// Create a new quintet array which is zero-padded version of the
	// given quintet array. Quintets are padded according to given job-batches.
	// Also paddedQuintets array is partially ordered.
	// --------------------------------------------------------------------------------------
	vector_int_t* rowOrderLookupJDS = NULL;
	sub_mtx_map_t* subMtxMaps = NULL;

	converter_quintetToJDSOrderedQuintet(
			quintets, nnz, rowCount, colCount, subMtxArr, subMtxCount, &rowOrderLookupJDS);

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("Quintets are partially ordered for JDS format.\n");
#endif

	// create final version of CSR which will be converted to JDS directly.
	// --------------------------------------------------------------------------------------
	spm_cmp_t* spmCsrCounterpart = NULL;
	converter_quintetToCSR_alloc(quintets, &spmCsrCounterpart, nnz, rowCount, colCount);

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("CSR counterpart of JDS is created.\n");
#endif

	// generate permutation array from partial ordering (for JDS format)
	// --------------------------------------------------------------------------------------

	DECIMAL* permutation = (DECIMAL*) malloc(sizeof(DECIMAL) * rowCount);
	DECIMAL i;
	for(i = 0; i < rowCount; ++i)
	{
		permutation[rowOrderLookupJDS->data[i]] = i;
	}

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("\"permutation\" array is created.\n");
#endif

#if PROGRAM_MODE >= TRACE_MODE
	PRINTF("Local cleanup done.\n");
#endif

	// return values
	*permutation_out = permutation;
	*spmCsrCounterpart_out = spmCsrCounterpart;
	*rowOrderLookup_out = rowOrderLookupJDS;
}

void converter_CSRToJDS(spm_cmp_t* spmCsr, DECIMAL* permutation, spm_jds_t** spmJds_out)
{
	sub_mtx_dim_t subMtx;
	sub_mtx_init(&subMtx, 0, 0, spmCsr->rowCount, spmCsr->colCount);
	converter_CSRPartToJDS(spmCsr, &subMtx, permutation, spmJds_out);
}

extern void converter_CSRPartToJDS(
		spm_cmp_t* spmCsr, sub_mtx_dim_t* subMtx, DECIMAL* permutation, spm_jds_t** spmJds_out)
{
	DECIMAL i;
	DECIMAL j;
	DECIMAL maxNNZ = spm_cmp_findMaxPtrGap(
			spmCsr, subMtx->start.i, subMtx->start.i + subMtx->length.i);
	DECIMAL subMtxTotalNNZ = sub_mtx_getNNZ_CSR(subMtx, spmCsr);

	// Fill JDS from CSR
	DECIMAL currIndex = 0;
	DECIMAL currColInd = 0;
	spm_jds_t* spmJds = spm_jds_new(subMtxTotalNNZ, subMtx->length.i, spmCsr->colCount, maxNNZ);
	for(currColInd = 0; currColInd < maxNNZ; ++currColInd)
	{
		spmJds->idiag[currColInd] = currIndex;
		for(i = subMtx->start.i; i < subMtx->start.i + subMtx->length.i; ++i)
		{
			int rowNNZ = spmCsr->ptr[i + 1] - spmCsr->ptr[i];
			if(currColInd < rowNNZ)
			{
				spmJds->dj[currIndex] = spmCsr->values[spmCsr->ptr[i] + currColInd];
				spmJds->jdiag[currIndex] = spmCsr->ind[spmCsr->ptr[i] + currColInd];
				++currIndex;
			}
		}
	}
	spmJds->idiag[currColInd] = currIndex;

	// Create partial permutation array for that specific sub-matrix
	spm_jds_addPermutation(spmJds, &permutation[subMtx->start.i]);


	// fill CSR counterpart
	spm_cmp_t* spmCsrCounterpart = spm_cmp_new(subMtx->length.i, subMtx->length.j, subMtxTotalNNZ);

	DECIMAL csrCounterpartInd = 0;
	DECIMAL subMtxColStartInd = subMtx->start.j;
	DECIMAL subMtxColEndInd = subMtx->start.j + subMtx->length.j;
	for(i = subMtx->start.i; i < subMtx->start.i + subMtx->length.i; ++i)
	{
		spmCsrCounterpart->ptr[i - subMtx->start.i] = csrCounterpartInd;
		for(j = spmCsr->ptr[i]; j < spmCsr->ptr[i + 1]; ++j)
		{
			if(spmCsr->ind[j] >= subMtxColStartInd && spmCsr->ind[j] < subMtxColEndInd)
			{
				spmCsrCounterpart->ind[csrCounterpartInd] = spmCsr->ind[j];
				spmCsrCounterpart->values[csrCounterpartInd] = spmCsr->values[j];
				++csrCounterpartInd;
			}
		}
	}
	spmCsrCounterpart->ptr[i - subMtx->start.i] = csrCounterpartInd;

	spmJds->csrCounterpart = spmCsrCounterpart;

	// return values
	*spmJds_out = spmJds;
}

void converter_extractPartialJDSArr(
		spm_cmp_t* spmCsrCounterpart, job_batch_t* subMtxArr, int subMtxCount,
		DECIMAL* permutation, job_batch_t** partialJdsBatchArr_out)
{
	job_batch_t* partialJDSBatchArr = (job_batch_t*) malloc(sizeof(job_batch_t) * subMtxCount);

	int i;
	for(i = 0; i < subMtxCount; ++i)
	{
		spm_jds_t* partialJds = NULL;
		sub_mtx_dim_t* subMtx = &subMtxArr[i].data.subMtx;
		converter_CSRPartToJDS(spmCsrCounterpart, subMtx, permutation, &partialJds);

		job_batch_initPartialJds(&partialJDSBatchArr[i], subMtx, partialJds);
	}

	// return values
	*partialJdsBatchArr_out = partialJDSBatchArr;
}

void converter_extractPartialJDSArrMultiple(
		spm_cmp_t* spmCsrCounterpart, int numBlocks,
		job_batch_t** subMtxArrPerBlock, int* subMtxCountPerBlock,
		DECIMAL* permutation, job_batch_t*** partialJdsBatchArrPerBlock_out)
{
	job_batch_t** partialJdsBatchArrPerBlock = (job_batch_t**) malloc(sizeof(job_batch_t*) * numBlocks);

	#pragma omp parallel num_threads(numBlocks)
	{
		int blockId = omp_get_thread_num();
		job_batch_t* partialJDSBatchArr = NULL;
		converter_extractPartialJDSArr(
				spmCsrCounterpart, subMtxArrPerBlock[blockId], subMtxCountPerBlock[blockId],
				permutation, &partialJDSBatchArr);

		partialJdsBatchArrPerBlock[blockId] = partialJDSBatchArr;
	}

	// return values
	*partialJdsBatchArrPerBlock_out = partialJdsBatchArrPerBlock;
}

void converter_extractHybridJDSCSR(
		job_batch_t* partialJdsBatch, job_batch_t* hybridJdsCsr_inout)
{
#ifdef JDS_BASED_HYBRID_DATA_EXTRACTION
	converter_extractHybridJDSCSR_JDS(partialJdsBatch, hybridJdsCsr_inout);
#else
	converter_extractHybridJDSCSR_CSR(partialJdsBatch, hybridJdsCsr_inout);
#endif
}

// TODO add support for true hybrid form in this
void converter_extractHybridJDSCSR_JDS(
		job_batch_t* partialJdsBatch, job_batch_t* hybridJdsCsr_inout)
{
	DECIMAL optColIndex = 0;
	DECIMAL optRowIndex = 0;
	spm_jds_t* spmJds = partialJdsBatch->data.partialJds.spmJds;
	spm_jds_findOptimumCut_JDS(spmJds, &optColIndex, &optRowIndex);

	spm_jds_t* jds = NULL;
	spm_cmp_t* csr = NULL;
	DECIMAL jdsNNZ = spmJds->idiag[optColIndex];
	DECIMAL csrNNZ = spmJds->nnz - jdsNNZ;

	DECIMAL i;
	DECIMAL j;
	DECIMAL k;

	// create JDS partition
	if(jdsNNZ > 0)
	{
		jds = spm_jds_new(jdsNNZ, spmJds->rowCount, optColIndex, optColIndex);
		for(i = 0; i < jdsNNZ; ++i)
		{
			jds->dj[i] = spmJds->dj[i];
			jds->jdiag[i] = spmJds->jdiag[i];
		}

		for(i = 0; i < jds->idiagLength + 1; ++i)
			jds->idiag[i] = spmJds->idiag[i];

		spm_jds_addPermutation(jds, spmJds->permutation);
	}

	// create CSR partition
	if(csrNNZ > 0)
	{
		csr = spm_cmp_new(optRowIndex, spmJds->colCount - optColIndex, csrNNZ);

		DECIMAL csrNNZIndex = 0;

		for(i = 0; i < csr->rowCount; ++i)
		{
			csr->ptr[i] = csrNNZIndex;
			for(j = spmJds->csrCounterpart->ptr[i] + optColIndex; j < spmJds->csrCounterpart->ptr[i + 1]; ++j)
			{
				csr->values[csrNNZIndex] = spmJds->csrCounterpart->values[j];
				csr->ind[csrNNZIndex] = spmJds->csrCounterpart->ind[j];
				++csrNNZIndex;
			}
		}
		csr->ptr[i] = csrNNZIndex;
	}

	// return values
	job_batch_initHybridJdsCsr(
			hybridJdsCsr_inout, &partialJdsBatch->data.partialJds.subMtx, jds, csr);
}

// TODO this has support for true hybrid form but it is done through a pragma, which has to be changed
void converter_extractHybridJDSCSR_CSR(
		job_batch_t* partialJdsBatch, job_batch_t* hybridJdsCsr_inout)
{
	DECIMAL optColIndex = 0;
	DECIMAL optRowIndex = 0;
	spm_jds_t* spmJds = partialJdsBatch->data.partialJds.spmJds;
	spm_jds_findOptimumCut_CSR(spmJds->csrCounterpart, &optColIndex, &optRowIndex);

	spm_jds_t* jds = NULL;
	spm_cmp_t* csr = NULL;
	DECIMAL csrNNZ = spmJds->csrCounterpart->ptr[optRowIndex];
	DECIMAL jdsNNZ = spmJds->csrCounterpart->nnz - csrNNZ;

	DECIMAL i;
	DECIMAL j;
	DECIMAL k;

#ifdef TRUE_HYBRID_FORM
	job_batch_initTrueHybrid(
			hybridJdsCsr_inout, &partialJdsBatch->data.partialJds.subMtx,
			csrNNZ, optRowIndex, spmJds->colCount,
			jdsNNZ, spmJds->rowCount - optRowIndex, optColIndex, optColIndex);

	// create CSR partition
	if(csrNNZ > 0)
	{
		csr = hybridJdsCsr_inout->data.hybrid_jds_csr.csr;

		for(i = 0; i < csr->rowCount; ++i)
		{
			csr->ptr[i] = spmJds->csrCounterpart->ptr[i];
			for(j = spmJds->csrCounterpart->ptr[i]; j < spmJds->csrCounterpart->ptr[i + 1]; ++j)
			{
				csr->values[j] = spmJds->csrCounterpart->values[j];
				csr->ind[j] = spmJds->csrCounterpart->ind[j];
			}
		}
		csr->ptr[i] = spmJds->csrCounterpart->ptr[i];
	}

	// create JDS partition
	if(jdsNNZ > 0)
	{
		jds = hybridJdsCsr_inout->data.hybrid_jds_csr.jds;

		DECIMAL jdsNNZIndex = 0;

		for(i = 0; i < jds->idiagLength; ++i)
		{
			jds->idiag[i] = jdsNNZIndex;
			for(j = spmJds->idiag[i] + optRowIndex;
				(j < spmJds->idiag[i + 1]);
				++j)
			{
				jds->dj[jdsNNZIndex] = spmJds->dj[j];
				jds->jdiag[jdsNNZIndex] = spmJds->jdiag[j];
				++jdsNNZIndex;
			}
		}
		jds->idiag[i] = jdsNNZIndex;

		spm_jds_addPermutation(jds, &spmJds->permutation[optRowIndex]);
	}
#else
	// create CSR partition
	if(csrNNZ > 0)
	{
		csr = spm_cmp_new(optRowIndex, spmJds->colCount, csrNNZ);

		for(i = 0; i < csr->rowCount; ++i)
		{
			csr->ptr[i] = spmJds->csrCounterpart->ptr[i];
			for(j = spmJds->csrCounterpart->ptr[i]; j < spmJds->csrCounterpart->ptr[i + 1]; ++j)
			{
				csr->values[j] = spmJds->csrCounterpart->values[j];
				csr->ind[j] = spmJds->csrCounterpart->ind[j];
			}
		}
		csr->ptr[i] = spmJds->csrCounterpart->ptr[i];
	}

	// create JDS partition
	if(jdsNNZ > 0)
	{
		jds = spm_jds_new(jdsNNZ, spmJds->rowCount - optRowIndex, optColIndex, optColIndex);

		DECIMAL jdsNNZIndex = 0;

		for(i = 0; i < jds->idiagLength; ++i)
		{
			jds->idiag[i] = jdsNNZIndex;
			for(j = spmJds->idiag[i] + optRowIndex;
				(j < spmJds->idiag[i + 1]);
				++j)
			{
				jds->dj[jdsNNZIndex] = spmJds->dj[j];
				jds->jdiag[jdsNNZIndex] = spmJds->jdiag[j];
				++jdsNNZIndex;
			}
		}
		jds->idiag[i] = jdsNNZIndex;

		spm_jds_addPermutation(jds, &spmJds->permutation[optRowIndex]);
	}

	job_batch_initHybridJdsCsr(
			hybridJdsCsr_inout, &partialJdsBatch->data.partialJds.subMtx, jds, csr);
#endif

}

void converter_extractHybridJDSCSRArr(
		job_batch_t* partialJdsBatchArr, int length, job_batch_t** hybridJdsCsrArr_out)
{
	job_batch_t* hybridJdsCsrArr = (job_batch_t*) malloc(sizeof(job_batch_t) * length);
	int i;
	for(i = 0; i < length; ++i)
	{
		job_batch_t* currPartialJdsBatch = &partialJdsBatchArr[i];
		job_batch_t* currHybridJdsCsrBatch = &hybridJdsCsrArr[i];

		converter_extractHybridJDSCSR(currPartialJdsBatch, currHybridJdsCsrBatch);
	}

	// return values
	*hybridJdsCsrArr_out = hybridJdsCsrArr;
}

void converter_extractHybridJDSCSRArrMultiple(
		int numBlocks, job_batch_t** partialJdsBatchArrPerBlock,
		int* partialJdsBatchCountPerBlock, job_batch_t*** hybridJdsCsrArrPerBlock_out)
{
	job_batch_t** hybridJdsCsrArrPerBlock =
			(job_batch_t**) malloc(sizeof(job_batch_t*) * numBlocks);

	#pragma omp parallel num_threads(numBlocks)
	{
		int blockId = omp_get_thread_num();
		int blockBatchCount = partialJdsBatchCountPerBlock[blockId];
		job_batch_t* blockPartialJdsBatchArr = partialJdsBatchArrPerBlock[blockId];
		job_batch_t* blockHybridBatchArr = NULL;

		converter_extractHybridJDSCSRArr(
				blockPartialJdsBatchArr, blockBatchCount, &blockHybridBatchArr);

		// thread output
		hybridJdsCsrArrPerBlock[blockId] = blockHybridBatchArr;
	}

	// return values
	*hybridJdsCsrArrPerBlock_out = hybridJdsCsrArrPerBlock;
}



// Helper functions
// ---------------------------------------------------------------------------------------------------------

DECIMAL __converter_CSR_to_ICSR_calculatePtrLength(DECIMAL rowCount, DECIMAL* ptr)
{
	DECIMAL ptrLength = 0;
	DECIMAL i = 0;
	while(i < rowCount)
	{
		if(ptr[i + 1] - ptr[i] > 0)
			++ptrLength;
		++i;
	}

	return ptrLength;
}

void __converter_CSR_to_ICSR_fillPtrArr(DECIMAL rowCount, DECIMAL* ptr, DECIMAL* incPtr)
{
	DECIMAL currInc = 0;
	DECIMAL incIndex = 0;

	DECIMAL i = 0;
	while(i < rowCount)
	{
		if(ptr[i + 1] - ptr[i] > 0)
		{
			incPtr[incIndex] = currInc;
			currInc = 1;
			++incIndex;
		}
		else
			++currInc;

		++i;
	}
}

void __converter_CSR_to_ICSR_fillIndArr(
		DECIMAL rowCount, DECIMAL colCount, DECIMAL* ptr, DECIMAL* ind, DECIMAL* incInd_inout)
{
	DECIMAL lastColInd = colCount;
	DECIMAL i;
	DECIMAL j;
	for(i = 0; i < rowCount; ++i)
	{
		if(ptr[i + 1] - ptr[i] > 0)
		{
			DECIMAL start = ptr[i];
			DECIMAL end = ptr[i + 1];
			incInd_inout[start] = colCount - lastColInd + ind[start];
			++start;

			while(start < end)
			{
				incInd_inout[start] = ind[start] - ind[start - 1];
				++start;
			}

			lastColInd = ind[end - 1];
		}
	}
}

void __converter_copyReal(REAL* source, REAL* destination, DECIMAL length)
{
	DECIMAL i;
	for(i = 0; i < length; ++i)
		destination[i] = source[i];
}

void __converter_smoothCSR(spm_cmp_t* csr)
{
	DECIMAL i;
	for(i = csr->rowCount; i > 0; --i)
	{
		if(csr->ptr[i] < csr->ptr[i - 1])
			csr->ptr[i - 1] = csr->ptr[i];
	}
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
