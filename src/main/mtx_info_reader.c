
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "include/data_structure/spm_storage.h"
#include "include/io/input_parser.h"
#include "include/io/converter.h"
#include "include/io/cli.h"

int main(int argc, char* argsv[])
{
	if(argc < 3)
	{
		printf("Usage: ./%s mtx_file_path (partitionSize_KB)+\n", argsv[0]);
		printf("Example: ./%s atmosmodd 1024", argsv[0]);
		printf("Example: ./%s Freescale1 1024 512 256 128", argsv[0]);
		printf("Don't forget the results are in double precision!!");
		return EXIT_SUCCESS;
	}

	int nnz;
	int rowCount;
	int colCount;

	quintet_t* triplets;
	mm_info_t* mmInfo;
	input_readQuintets(argsv[1], &nnz, &rowCount, &colCount, &triplets, &mmInfo);
	quintet_sortRowIndex(triplets, nnz);

	spm_cmp_t* spm = spm_cmp_new(rowCount, colCount, nnz);
	converter_quintetToCSR(triplets, spm);

	float avgNonZeroPerRow = (float) nnz / (float) rowCount;
	float avgNonZeroPerCol = (float) nnz / (float) colCount;
	int maxNonZeroPerRow = spm->ptr[1] - spm->ptr[0];
	int minNonZeroPerRow = spm->ptr[1] - spm->ptr[0];

	int i;
	for(i = 0; i < spm->rowCount; ++i)
	{
		if(maxNonZeroPerRow < spm->ptr[i + 1] - spm->ptr[i])
			maxNonZeroPerRow = spm->ptr[i + 1] - spm->ptr[i];
		else if(minNonZeroPerRow > spm->ptr[i + 1] - spm->ptr[i])
			minNonZeroPerRow = spm->ptr[i + 1] - spm->ptr[i];
	}

	int* nonZeroCountPerCol = malloc(sizeof(int) * spm->colCount);
	for(i = 0; i < spm->colCount; ++i)
		nonZeroCountPerCol[i] = 0;

	for(i = 0; i < spm->nnz; ++i)
		++nonZeroCountPerCol[spm->ind[i]];

	int maxNonZeroPerCol = nonZeroCountPerCol[0];
	int minNonZeroPerCol = nonZeroCountPerCol[0];

	for(i = 0; i < spm->colCount; ++i)
	{
		if(maxNonZeroPerCol < nonZeroCountPerCol[i])
			maxNonZeroPerCol = nonZeroCountPerCol[i];
		else if(minNonZeroPerCol > nonZeroCountPerCol[i])
			minNonZeroPerCol = nonZeroCountPerCol[i];
	}

//	float mtxSizeKBSinglePrecision =
//			(spm->nnz * sizeof(int) + spm->nnz * sizeof(float) + (spm->rowCount + 1) * sizeof(int)) / pow(2, 10);
//	float operationSizeKBSinglePrecision = mtxSizeKBSinglePrecision +
//			(spm->rowCount + spm->colCount) * sizeof(float) / pow(2, 10);

	float mtxSizeKBDoublePrecision =
			(spm->nnz * sizeof(int) + spm->nnz * sizeof(double) + (spm->rowCount + 1) * sizeof(int)) / pow(2, 10);
	float operationSizeKBDoublePrecision = mtxSizeKBDoublePrecision +
			(spm->rowCount + spm->colCount) * sizeof(double) / pow(2, 10);

	/*
	printf("MTX: %s\n", argsv[1]);
	printf("nnz: %d\n", nnz);
	printf("row-count: %d\n", rowCount);
	printf("col-count: %d\n", colCount);
	printf("average non-zeros per row: %0.2f\n", avgNonZeroPerRow);
	printf("average non-zeros per column: %0.2f\n", avgNonZeroPerCol);
	printf("max non-zeros per row: %d\n", maxNonZeroPerRow);
	printf("max non-zeros per column: %d\n", maxNonZeroPerCol);
	printf("min non-zeros per row: %d\n", minNonZeroPerRow);
	printf("min non-zeros per column: %d\n", minNonZeroPerCol);
	printf("single-precision mtx-size (KB): %3.2f, operation-size (KB): %3.2f\n", mtxSizeKBSinglePrecision, operationSizeKBSinglePrecision);
	printf("double-precision mtx-size (KB): %3.2f, operation-size (KB): %3.2f\n", mtxSizeKBDoublePrecision, operationSizeKBDoublePrecision);

	for(i = 2; i < argc; ++i)
	{
		float cacheSize = atof(argsv[i]);
		printf("sub-matrix count for cache-size: %3.2f => single-precision: %3.2f  double-precision: %3.2f.\n",
				cacheSize, mtxSizeKBDoublePrecision / cacheSize, mtxSizeKBDoublePrecision / cacheSize);
	}

	*/

	char *mtxName = NULL;
	char *mtxPath = NULL;
	printf("matrixList.append(Matrix(\"%s\", %d, %d, \"colnet\", {", argsv[1], mmInfo->isSymmetric, !mmInfo->isReal);
	printf("\"%.0fK\": %d", atof(argsv[2]), (int) (mtxSizeKBDoublePrecision / atof(argsv[2])) + 1);
	for(i = 3; i < argc; ++i)
	{
		float cacheSize = atof(argsv[i]);
		printf(", \"%.0fK\": %d", cacheSize ,(int) (mtxSizeKBDoublePrecision / cacheSize) + 1);
	}
	printf("}))\n");

	// printf("\n");

	free(triplets);
	free(nonZeroCountPerCol);
	spm_cmp_delete(spm);
	mm_info_delete(mmInfo);

	return EXIT_SUCCESS;
}
