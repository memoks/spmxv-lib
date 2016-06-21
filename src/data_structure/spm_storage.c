
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdio.h>
#include <stdlib.h>
#include <alloca.h>
#include <string.h>

#include "include/data_structure/spm_storage.h"

#ifdef __ICC
#include <mkl.h>
#endif

// functions for traditional CSR / CSC storage format
// -------------------------------------------------------------------------

spm_cmp_t* spm_cmp_new(DECIMAL rowCount, DECIMAL colCount, DECIMAL nnz)
{
	spm_cmp_t* spm = NULL;
#ifdef __ICC
	spm = (spm_cmp_t*) mkl_malloc(sizeof(spm_cmp_t), ALIGNMENT);
#else
	posix_memalign(&spm, ALIGNMENT, sizeof(spm_cmp_t));
#endif
	spm_cmp_initDefault(spm);


#ifdef __ICC
	spm->values = (REAL*) mkl_malloc(nnz * sizeof(REAL), ALIGNMENT);
	spm->ind = (DECIMAL*) mkl_malloc(nnz * sizeof(DECIMAL), ALIGNMENT);
	spm->ptr = (DECIMAL*) mkl_malloc((rowCount + 1) * sizeof(DECIMAL), ALIGNMENT);
#else
	posix_memalign(&spm->values, ALIGNMENT, sizeof(REAL) * nnz);
	posix_memalign(&spm->ind, ALIGNMENT, sizeof(DECIMAL) * nnz);
	posix_memalign(&spm->ptr, ALIGNMENT, sizeof(DECIMAL) * (rowCount + 1));
#endif

	spm->nnz = nnz;
	spm->rowCount = rowCount;
	spm->colCount = colCount;

	DECIMAL i;
	for(i = 0; i < nnz; ++i)
	{
		spm->values[i] = 0.0;
		spm->ind[i] = 0;
	}

	for(i = 0; i < rowCount + 1; ++i)
	{
		spm->ptr[i] = nnz;
	}

	return spm;
}

void spm_cmp_initDefault(spm_cmp_t* spm)
{
	spm->nnz = 0;
	spm->rowCount = 0;
	spm->colCount = 0;
	spm->values = NULL;
	spm->ind = NULL;
	spm->ptr = NULL;
}

void spm_cmp_copy(spm_cmp_t* source, spm_cmp_t* destination)
{
	destination->nnz = source->nnz;
	destination->rowCount = source->rowCount;
	destination->colCount = source->colCount;

	DECIMAL i;
	for(i = 0; i < destination->nnz; ++i)
	{
		destination->values[i] = source->values[i];
		destination->ind[i] = source->ind[i];
	}

	for(i = 0; i < destination->rowCount + 1; ++i)
		destination->ptr[i] = source->ptr[i];
}

REAL spm_cmp_get(spm_cmp_t* spm, DECIMAL ptrIndex, DECIMAL indIndex)
{
	int ind;
	for(ind = spm->ptr[ptrIndex]; ind < spm->ptr[ptrIndex + 1]; ++ind)
	{
		if(indIndex == spm->ind[ind])
			return spm->values[ind];
	}

	return 0.0;
}

void spm_cmp_print(spm_cmp_t* spm)
{
	if(spm == NULL)
		return;

	spm_cmp_printAttributes(spm);
	int i;

	printf("values array:");
	for(i = 0; i < spm->nnz; ++i)
		printf(" %3.2f", spm->values[i]);
	printf("\n");

	printf("indexArr array:");
	for(i = 0; i < spm->nnz; ++i)
		printf(" %d", spm->ind[i]);
	printf("\n");

	printf("ptrArray array:");
	for(i = 0; i < spm->rowCount + 1; ++i)
		printf(" %d", spm->ptr[i]);
	printf("\n");
}

void spm_cmp_printMultiple(char* m, spm_cmp_t** spms, int numSpms)
{
	printf("%s\n", m);
	int i;
	for(i = 0; i < numSpms; ++i)
		spm_cmp_print(spms[i]);
	printf("\n");
}

void spm_cmp_print2DFormat(char *m, spm_cmp_t* spm)
{
	printf("%s\n", m);
	spm_cmp_printAttributes(spm);

	int i;
	int j;
	for(i = 0; i < spm->rowCount; ++i)
	{
		int colInd = 0;
		for(j = spm->ptr[i]; j < spm->ptr[i + 1]; ++j)
		{
			while(colInd != spm->ind[j])
			{
				printf("%3.2f  ", 0.0);
				++colInd;
			}

			printf("%3.2f  ", spm->values[j]);
			++colInd;
		}

		while(colInd < spm->colCount)
		{
			printf("%3.2f  ", 0.0);
			++colInd;
		}

		printf("\n");
	}

	printf("\n");
}

void spm_cmp_print2DFormatMultiple(char* m, spm_cmp_t** spms, int numSpms)
{
	printf("%s\n", m);
	printf("---------------------");
	int i;
	for(i = 0; i < numSpms; ++i)
		spm_cmp_print2DFormat("", spms[i]);
	printf("---------------------\n\n");
}

void spm_cmp_printToFile(char* filePath, char* m, spm_cmp_t* spm)
{
	FILE* f = fopen(filePath, "w");
	fprintf(f, "%s\n", m);

	int i;
	int j;
	for(i = 0; i < spm->rowCount; ++i)
		for(j = 0; j < spm->nnz; ++j)
			fprintf(f, "%d %d\n", i, spm->ind[j]);

	fclose(f);
}

void spm_cmp_delete(spm_cmp_t* spm)
{
	spm_cmp_deleteNonPtr(spm);

#ifdef __ICC
	mkl_free(spm);
#else
	free(spm);
#endif
}

void spm_cmp_deleteNonPtr(spm_cmp_t* spm)
{
	if(spm == NULL)
		return;

#ifdef __ICC
	mkl_free(spm->values);
	mkl_free(spm->ind);
	mkl_free(spm->ptr);
#else
	free(spm->values);
	free(spm->ind);
	free(spm->ptr);
#endif
}

void spm_cmp_extractNonZeroStatistics_CSR(
		spm_cmp_t* spmCsr,
		REAL* rowMinNNZ_out, REAL* rowAvgNNZ_out, REAL* rowMaxNNZ_out,
		REAL* columnMinNNZ_out, REAL* columnAvgNNZ_out, REAL* columnMaxNNZ_out)
{
	DECIMAL i;

	// ------------------------------------------------------------------------

	REAL avgr = ((REAL) spmCsr->nnz) / ((REAL) spmCsr->rowCount);
	REAL avgc = ((REAL) spmCsr->nnz) / ((REAL) spmCsr->colCount);

	// ------------------------------------------------------------------------

	REAL maxr = spmCsr->ptr[1] - spmCsr->ptr[0];
	REAL minr = spmCsr->ptr[1] - spmCsr->ptr[0];

	for(i = 0; i < spmCsr->rowCount; ++i)
	{
		if(maxr < spmCsr->ptr[i + 1] - spmCsr->ptr[i])
			maxr = spmCsr->ptr[i + 1] - spmCsr->ptr[i];
		else if(minr > spmCsr->ptr[i + 1] - spmCsr->ptr[i])
			minr = spmCsr->ptr[i + 1] - spmCsr->ptr[i];
	}

	// ------------------------------------------------------------------------

	DECIMAL* nonZeroCountPerCol = (DECIMAL*) malloc(sizeof(DECIMAL) * spmCsr->colCount);
	for(i = 0; i < spmCsr->colCount; ++i)
		nonZeroCountPerCol[i] = 0;

	for(i = 0; i < spmCsr->nnz; ++i)
		++nonZeroCountPerCol[spmCsr->ind[i]];

	REAL maxc = nonZeroCountPerCol[0];
	REAL minc = nonZeroCountPerCol[0];

	for(i = 0; i < spmCsr->colCount; ++i)
	{
		if(maxc < nonZeroCountPerCol[i])
			maxc = nonZeroCountPerCol[i];
		else if(minc > nonZeroCountPerCol[i])
			minc = nonZeroCountPerCol[i];
	}

	// ------------------------------------------------------------------------

	// clean up
	free(nonZeroCountPerCol);

	// return values
	*rowMinNNZ_out = minr;
	*rowAvgNNZ_out = avgr;
	*rowMaxNNZ_out = maxr;
	*columnMinNNZ_out = minc;
	*columnAvgNNZ_out = avgc;
	*columnMaxNNZ_out = maxc;
}

void spm_cmp_plotShovedLeft(spm_cmp_t* spm)
{
	if(spm == NULL)
		return;

	char buff[DEFAULT_STR_BUFF_SIZE] = "\0";
	spm_cmp_plotToBufferShovedLeft(spm, buff);
	PRINTF("%s\n", buff);
}

void spm_cmp_plotToBufferShovedLeft(spm_cmp_t* spm, char* buff_inout)
{
	if(spm == NULL)
		return;

	DECIMAL i;
	DECIMAL j;
	for(i = 0; i < spm->rowCount; ++i)
	{
		for(j = spm->ptr[i]; j < spm->ptr[i + 1]; ++j)
		{
			strcat(buff_inout, "*");
		}

		strcat(buff_inout, "\n");
	}
}

void spm_cmp_printAttributes(spm_cmp_t* spm)
{
	char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
	spm_cmp_toStringAttributes(spm, temp);
	PRINTF("%s\n", temp);
}

void spm_cmp_toStringAttributes(spm_cmp_t* spm, char* buff_inout)
{
	if(spm == NULL)
		return;

	char temp[DEFAULT_STR_BUFF_SIZE];
	sprintf(temp, "SPM_CMP: %d rows, %d cols, and %d elements",
			spm->rowCount, spm->colCount, spm->nnz);
	strcat(buff_inout, temp);
}

DECIMAL spm_cmp_findMaxPtrGap(spm_cmp_t* spm, int startPtrIndex, int endPtrIndex)
{
	DECIMAL maxGap = 0;
	DECIMAL i;
	DECIMAL j;
	for(i = startPtrIndex; i < endPtrIndex; ++i)
	{
		DECIMAL ithRowNNZ = spm->ptr[i + 1] - spm->ptr[i];
		if(maxGap < ithRowNNZ)
			maxGap = ithRowNNZ;
	}

	return maxGap;
}

// functions for ICSR / ICSC storage format
// -------------------------------------------------------------------------------------------------------------

spm_inc_t* spm_inc_new(DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount, DECIMAL ptrLength)
{
	spm_inc_t* spmInc = NULL;
#ifdef __ICC
	spmInc = (spm_inc_t*) mkl_malloc(sizeof(spm_inc_t), ALIGNMENT);
#else
	posix_memalign(&spmInc, ALIGNMENT, sizeof(spm_inc_t));
#endif

	spm_inc_init(spmInc, nnz, rowCount, colCount, ptrLength);
	return spmInc;
}

void spm_inc_init(
		spm_inc_t* spmInc, DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount, DECIMAL ptrLength)
{
	spmInc->nnz = nnz;
	spmInc->rowCount = rowCount;
	spmInc->colCount = colCount;
	spmInc->ptrLength = ptrLength;

#ifdef __ICC
	spmInc->values = (REAL*) mkl_malloc((spmInc->nnz + 1) * sizeof(REAL), ALIGNMENT);
	spmInc->ind = (DECIMAL*) mkl_malloc((spmInc->nnz + 1) * sizeof(DECIMAL), ALIGNMENT);
	spmInc->ptr = (DECIMAL*) mkl_malloc((spmInc->ptrLength) * sizeof(DECIMAL), ALIGNMENT);
#else
	posix_memalign(&spmInc->values, ALIGNMENT, sizeof(REAL) * (spmInc->nnz + 1));
	posix_memalign(&spmInc->ind, ALIGNMENT, sizeof(DECIMAL) * (spmInc->nnz + 1));
	posix_memalign(&spmInc->ptr, ALIGNMENT, sizeof(DECIMAL) * (spmInc->ptrLength));
#endif

	// this is not secure for partial spm_inc implementation
	// spmInc->ind[nnz] = nnz;
	spmInc->ind[nnz] = colCount + rowCount + nnz;
	spmInc->values[nnz] = 0.0;
}

void spm_inc_initDefault(spm_inc_t* spmInc)
{
	spmInc->nnz = 0;
	spmInc->rowCount = 0;
	spmInc->colCount = 0;
	spmInc->ptrLength = 0;
	spmInc->values = NULL;
	spmInc->ind = NULL;
	spmInc->ptr = NULL;
}

void spm_inc_copy(spm_inc_t* source, spm_inc_t* destination)
{
	spm_inc_init(destination, source->nnz, source->rowCount, source->colCount, source->ptrLength);

	DECIMAL i;
	for(i = 0; i < destination->nnz + 1; ++i)
	{
		destination->values[i] = source->values[i];
		destination->ind[i] = source->ind[i];
	}

	for(i = 0; i < destination->ptrLength; ++i)
		destination->ptr[i] = source->ptr[i];
}

void spm_inc_delete(spm_inc_t* spm)
{
	if(spm == NULL)
		return;

	spm_inc_deleteNonPtr(spm);

#ifdef __ICC
	mkl_free(spm);
#else
	free(spm);
#endif

}

void spm_inc_deleteNonPtr(spm_inc_t* spm)
{
	spm->nnz = 0;
	spm->rowCount = 0;
	spm->colCount = 0;
	spm->ptrLength = 0;

#ifdef __ICC
	mkl_free(spm->values);
	mkl_free(spm->ind);
	mkl_free(spm->ptr);
#else
	free(spm->values);
	free(spm->ind);
	free(spm->ptr);
#endif
}

void spm_inc_deleteMultiple(spm_inc_t** spms, int numSpms)
{
	int i;
	for(i = 0; i < numSpms; ++i)
		spm_inc_delete(spms[i]);

	free(spms);
}

void spm_inc_print(spm_inc_t* spmInc)
{
	DECIMAL i;
	printf("ICSR %d rows, %d columns, %d nnz.\n",
			spmInc->rowCount, spmInc->colCount, spmInc->nnz);

	printf("Values:");
	for(i = 0; i < spmInc->nnz; ++i)
		printf(" %3.2f", spmInc->values[i]);

	printf("\nIndex-Increment:");
	for(i = 0; i < spmInc->nnz + 1; ++i)
		printf(" %d", spmInc->ind[i]);

	printf("\nPtr-increment:");
	for(i = 0; i < spmInc->ptrLength; ++i)
		printf(" %d", spmInc->ptr[i]);

	printf("\n");
}

void spm_inc_printMultiple(char* m, spm_inc_t** spmIncs, int numSpms)
{
	printf("%s\n", m);
	int i;
	for(i = 0; i < numSpms; ++i)
	{
		spm_inc_print(spmIncs[i]);
		printf("\n");
	}
}

void spm_inc_printAttributes(spm_inc_t* spm)
{
	char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
	spm_inc_toStringAttributes(spm, temp);
	PRINTF("%s\n", temp);
}

void spm_inc_toStringAttributes(spm_inc_t* spm, char* buff_inout)
{
	char temp[DEFAULT_STR_BUFF_SIZE];
	sprintf(temp, "SPM_INC: %d rows, %d columns, and %d elements...",
			spm->rowCount, spm->colCount, spm->nnz);
	strcat(buff_inout, temp);
}

// functions JDS format
// -------------------------------------------------------------------------------------------------------------

spm_jds_t* spm_jds_new(DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount, DECIMAL idiagLength)
{
	spm_jds_t* spmJds = NULL;

#ifdef __ICC
	spmJds = (spm_jds_t*) mkl_malloc(sizeof(spm_jds_t), ALIGNMENT);
	spmJds->dj = (REAL*) mkl_malloc(nnz * sizeof(REAL), ALIGNMENT);
	spmJds->jdiag = (DECIMAL*) mkl_malloc(nnz * sizeof(DECIMAL), ALIGNMENT);
	spmJds->idiag = (DECIMAL*) mkl_malloc((idiagLength + 1) * sizeof(DECIMAL), ALIGNMENT);
#else
	posix_memalign(&spmJds, ALIGNMENT, sizeof(spm_jds_t));
	posix_memalign(&spmJds->dj, ALIGNMENT, sizeof(REAL) * nnz);
	posix_memalign(&spmJds->jdiag, ALIGNMENT, sizeof(DECIMAL) * nnz);
	posix_memalign(&spmJds->idiag, ALIGNMENT, sizeof(DECIMAL) * (idiagLength + 1));
#endif

	spmJds->permutation = NULL;
	spmJds->colCount = colCount;
	spmJds->rowCount = rowCount;
	spmJds->nnz = nnz;
	spmJds->idiagLength = idiagLength;
	spmJds->csrCounterpart = NULL;

	spm_jds_initDefault(spmJds);

	return spmJds;
}

void spm_jds_initDefault(spm_jds_t* spmJds)
{
	DECIMAL i;
	for(i = 0; i < spmJds->nnz; ++i)
	{
		spmJds->dj[i] = 0;
		spmJds->jdiag[i] = 0;
	}

	for(i = 0; i < spmJds->idiagLength + 1; ++i)
		spmJds->idiag[i] = 0;
}

void spm_jds_copy(spm_jds_t* source, spm_jds_t* destination)
{
	DECIMAL i;
	for(i = 0; i < source->nnz; ++i)
	{
		destination->dj[i] = source->dj[i];
		destination->jdiag[i] = source->jdiag[i];
	}

	for(i = 0; i < source->idiagLength + 1; ++i)
		destination->idiag[i] = source->idiag[i];

	spm_jds_addPermutation(destination, source->permutation);

	if(source->csrCounterpart != NULL)
	{
		destination->csrCounterpart = spm_cmp_new(
				source->csrCounterpart->rowCount,
				source->csrCounterpart->colCount,
				source->csrCounterpart->nnz);

		spm_cmp_copy(source->csrCounterpart, destination->csrCounterpart);
	}
}

void spm_jds_printAttributes(spm_jds_t* spmJds)
{
	char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
	spm_jds_toStringAttributes(spmJds, temp);
	PRINTF("%s\n", temp);
}

void spm_jds_toStringAttributes(spm_jds_t* spmJds, char* buff_inout)
{
	if(spmJds == NULL)
		return;

	char temp[DEFAULT_STR_BUFF_SIZE];
	sprintf(temp, "SPM_JDS: %d rows, %d cols, %d idiagLength, and %d elements",
			spmJds->rowCount, spmJds->colCount, spmJds->idiagLength, spmJds->nnz);
	strcat(buff_inout, temp);
}

void spm_jds_print(spm_jds_t* spmJds)
{
	if(spmJds == NULL)
		return;

	spm_jds_printAttributes(spmJds);

	DECIMAL i;

	PRINTF("DJ:");
	for(i = 0; i < spmJds->nnz; ++i)
		PRINTF(" %0.2f", spmJds->dj[i]);

	PRINTF("\nJDIAG:");
	for(i = 0; i < spmJds->nnz; ++i)
		PRINTF(" %d", spmJds->jdiag[i]);

	PRINTF("\nIDIAG:");
	for(i = 0; i < spmJds->idiagLength + 1; ++i)
		PRINTF(" %d", spmJds->idiag[i]);

	PRINTF("\nPermutation:");
	for(i = 0; i < spmJds->rowCount; ++i)
		PRINTF(" %d", spmJds->permutation[i]);

	PRINTF("\n");
}

void spm_jds_deleteNonPtr(spm_jds_t* spmJds)
{
	if(spmJds == NULL)
		return;

#ifdef __ICC
	mkl_free(spmJds->dj);
	mkl_free(spmJds->jdiag);
	mkl_free(spmJds->idiag);

	if(spmJds->permutation != NULL)
		mkl_free(spmJds->permutation);
#else

	free(spmJds->dj);
	free(spmJds->jdiag);
	free(spmJds->idiag);

	if(spmJds->permutation != NULL)
		free(spmJds->permutation);
#endif

	if(spmJds->csrCounterpart != NULL)
		spm_cmp_delete(spmJds->csrCounterpart);
}

void spm_jds_delete(spm_jds_t* spmJds)
{
	if(spmJds == NULL)
		return;

	spm_jds_deleteNonPtr(spmJds);

#ifdef __ICC
	mkl_free(spmJds);
#else
	free(spmJds);
#endif

}

void spm_jds_addPermutation(spm_jds_t* spmJds, DECIMAL* permutation)
{
#ifdef __ICC
	spmJds->permutation = (DECIMAL*) mkl_malloc(spmJds->rowCount * sizeof(DECIMAL), ALIGNMENT);
#else
	posix_memalign(&spmJds->permutation, ALIGNMENT, sizeof(DECIMAL) * spmJds->rowCount);
#endif

	int i;
	for(i = 0; i < spmJds->rowCount; ++i)
		spmJds->permutation[i] = permutation[i];
}

void spm_jds_extractNonZeroStatistics_JDS(
		spm_jds_t* spmJds,
		REAL* rowMinNNZ_out, REAL* rowAvgNNZ_out, REAL* rowMaxNNZ_out,
		REAL* columnMinNNZ_out, REAL* columnAvgNNZ_out, REAL* columnMaxNNZ_out)
{
	spm_cmp_extractNonZeroStatistics_CSR(
			spmJds->csrCounterpart,
			rowMinNNZ_out, rowAvgNNZ_out, rowMaxNNZ_out,
			columnMinNNZ_out, columnAvgNNZ_out, columnMaxNNZ_out);
}

void spm_jds_plot_JDS(spm_jds_t* spmJds)
{
	if(spmJds == NULL)
		return;

	char buff[DEFAULT_STR_BUFF_SIZE] = "\0";
	spm_jds_plotToBuffer_JDS(spmJds, buff);
	PRINTF("%s\n", buff);
}

void spm_jds_plotToBuffer_JDS(spm_jds_t* spmJds, char* buff_inout)
{
	if(spmJds == NULL)
		return;

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

	for(i = 0; i < spmJds->rowCount; ++i)
	{
		for(j = 0; j < rowNNZArr[i]; ++j)
			strcat(buff_inout, "*");
		strcat(buff_inout, "\n");
	}

	// clean up
	free(rowNNZArr);
}

void spm_jds_findOptimumCut_JDS(
		spm_jds_t* spmJds, DECIMAL* optColInd_out, DECIMAL* optRowInd_out)
{
	DECIMAL optColInd = 0;
	DECIMAL optRowInd = 0;
	DECIMAL i;
	DECIMAL j;
	DECIMAL maxArea = 0;
	DECIMAL y = 0;
	DECIMAL x = 0;
	DECIMAL effectiveRowCount = spmJds->rowCount;
	DECIMAL effectiveColCount = spmJds->colCount;
	// INTEGER effectiveColCount = spmJds->idiagLength;
	for(i = 0; i < spmJds->idiagLength; ++i)
	{
		y = effectiveRowCount - (spmJds->idiag[i + 1] - spmJds->idiag[i]);
		x = effectiveColCount - i;
		// INTEGER area = x * y;
		DECIMAL area = x + y;
		if(maxArea <= area)
		{
			maxArea = area;
			optColInd = i;
			optRowInd = effectiveRowCount - y;
		}
	}

	y = spmJds->rowCount;
	x = spmJds->colCount - i;
	// INTEGER area = x * y;
	DECIMAL area = x + y;
	if(maxArea <= area)
	{
		maxArea = area;
		optColInd = i;
		optRowInd = 0;
	}

	// return values
	*optColInd_out = optColInd;
	*optRowInd_out = optRowInd;
}

void spm_jds_findOptimumCut_CSR(
		spm_cmp_t* spmCsr, DECIMAL* optColInd_out, DECIMAL* optRowInd_out)
{
	DECIMAL optColInd = 0;
	DECIMAL optRowInd = 0;
	DECIMAL i;
	DECIMAL j;
	DECIMAL currIdiagStart = 0;
	DECIMAL maxPerimeter = 0;
	DECIMAL effectiveRowCount = spmCsr->rowCount;
	DECIMAL effectiveColCount = spmCsr->colCount;
	DECIMAL y = 0;
	DECIMAL x = 0;
	for(i = 0; i < spmCsr->rowCount; ++i)
	{
		y = effectiveRowCount - i;
		x = effectiveColCount - (spmCsr->ptr[i + 1] - spmCsr->ptr[i]);

		DECIMAL perimeter = x + y;
		if(maxPerimeter <= perimeter)
		{
			maxPerimeter = perimeter;
			optColInd = effectiveColCount - x;
			optRowInd = i;
		}
	}

	y = spmCsr->rowCount - i;
	x = spmCsr->colCount;

	DECIMAL perimeter = x + y;
	if(maxPerimeter <= perimeter)
	{
		maxPerimeter = perimeter;
		optColInd = 0;
		optRowInd = i;
	}

	// return values
	*optColInd_out = optColInd;
	*optRowInd_out = optRowInd;
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
