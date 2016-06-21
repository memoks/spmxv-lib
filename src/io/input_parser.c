
#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>

#include "include/config.h"
#include "include/task_decomposition/partitioning.h"
#include "include/algorithm/algorithm.h"
#include "include/io/converter.h"


#include "include/io/input_parser.h"

// Some macros to save same space and hide code blocks
// -----------------------------------------------------------------------------

// Opens a file stream.
#define OPEN_FILE_STREAM(fileName, mode)			\
	FILE *f = fopen(fileName, mode);				\
	if(f == NULL)									\
	{												\
		PRINTF("Cannot find file: %s\n", fileName);	\
		exit(EXIT_FAILURE);							\
	}

// Closes a file stream and a quick clean-up service.
#define CLOSE_FILE_STREAM(line, f)					\
	if(line != NULL) 								\
		free(line);									\
	if(f != NULL)									\
		fclose(f);


void input_printQuintetsToFileInMMF(
		FILE *f, quintet_t* quintets, int rowCount, int colCount, int nnz, int isBinary)
{
	fprintf(f, "%c MTX header\n", '%');
	fprintf(f, "%d %d %d\n", rowCount, colCount, nnz);

	int i;
	if(isBinary)
	{
		for(i = 0; i < nnz; ++i)
		{
			fprintf(f, "%d %d\n", quintets[i].i + 1, quintets[i].j + 1);
		}
	}
	else
	{
		for(i = 0; i < nnz; ++i)
		{
			fprintf(f, "%d %d %f\n", quintets[i].i + 1, quintets[i].j + 1, quintets[i].val);
		}
	}
}

void input_readQuintets(
		char* mmfPath,
		DECIMAL *nnz_out, DECIMAL* rowCount_out, DECIMAL* colCount_out,
		quintet_t** quintets_out, mm_info_t** mmInfo_out)
{
	OPEN_FILE_STREAM(mmfPath, "r")
	DECIMAL* rowIndices = NULL;
	DECIMAL* colIndices = NULL;
	REAL* values = NULL;
	mm_info_t* mmInfo = NULL;
	DECIMAL nnz;

	mm_read(f, rowCount_out, colCount_out, &nnz, &values, &rowIndices, &colIndices, &mmInfo);
	quintet_t* quintets = (quintet_t*) malloc(sizeof(quintet_t) * (nnz));

	DECIMAL i;
	for(i = 0; i < nnz; ++i)
	{
		quintets[i].i = rowIndices[i];
		quintets[i].j = colIndices[i];
		quintets[i].val = values[i];
		quintets[i].iVal = quintets[i].i;
		quintets[i].jVal = quintets[i].j;
	}

	if(mmInfo->isPattern)
		quintet_adjustBinaryMatrixCoefficients(quintets, nnz);

	// local clean up
	CLOSE_FILE_STREAM(NULL, f)
	free(rowIndices);
	free(colIndices);
	free(values);

	// return values
	*quintets_out = quintets;
	*mmInfo_out = mmInfo;
	*nnz_out = nnz;
}


void input_readOrderLookup(char* fileName, quintet_t* quintets_in_out,
		DECIMAL nnz, DECIMAL rowCount, DECIMAL columnCount,
		vector_int_t** rowOrderLookup_out, vector_int_t** colOrderLookup_out)
{
	*rowOrderLookup_out = vector_int_new(rowCount);
	*colOrderLookup_out = vector_int_new(columnCount);

	OPEN_FILE_STREAM(fileName, "r")

	char* line = NULL;
	size_t len = 0;

	int elem;
	int i;

	// Read row ordering
	for(i = 0; i < (*rowOrderLookup_out)->length; ++i)
	{
		getline(&line, &len, f);
		sscanf(line, "%d", &elem);
		(*rowOrderLookup_out)->data[i] = elem;

	}

	// Read column ordering
	for(i = 0; i < (*colOrderLookup_out)->length; ++i)
	{
		getline(&line, &len, f);
		sscanf(line, "%d", &elem);
		(*colOrderLookup_out)->data[i] = elem;
	}

	// add ordering information to quintets
	for(i = 0; i < nnz; ++i)
	{
		quintets_in_out[i].iVal = (*rowOrderLookup_out)->data[quintets_in_out[i].i];
		quintets_in_out[i].jVal = (*colOrderLookup_out)->data[quintets_in_out[i].j];
	}

	CLOSE_FILE_STREAM(line, f)
}

int input_getLineCount(char* fileName)
{
	OPEN_FILE_STREAM(fileName, "r")

	char* line = NULL;
	size_t len = 0;
	size_t read = 0;
	int lineCount = 0;

	while((read = getline(&line, &len, f)) != -1)
	{
		++lineCount;
	}

	CLOSE_FILE_STREAM(line, f)

	return lineCount;
}

void input_readVector(char* fileName, vector_int_t* vector_inout)
{
	OPEN_FILE_STREAM(fileName, "r")

	char* line = NULL;
	size_t len = 0;

	int i;
	for(i = 0; i < vector_inout->length; ++i)
	{
		getline(&line, &len, f);
		vector_inout->data[i] = atoi(line);
	}

	CLOSE_FILE_STREAM(line, f)

//	#if PROGRAM_MODE >= TRACE_MODE
//		vector_int_print("\nPrinting ordering... ", vector_inout);
//	#endif

}

void input_kPatohOrder(
		char* mmfPath, char* mmfDirPath, char* mmfName, int cacheSizeKB, char* partitioningTypeStr,
		quintet_t* quintets_inout, DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount,
		vector_int_t** rowOrderLookup_out, vector_int_t** rowPartitioning_out,
		vector_int_t** colOrderLookup_out, vector_int_t** colPartitioning_out)
{
	// read order lookups
	// ------------------------------------------------------------------------------------------------------------
	char rowOrderLookupPath[DEFAULT_STR_BUFF_SIZE];
	sprintf(rowOrderLookupPath, "%s/%s/%s_vertexOrder_%dK", mmfDirPath, partitioningTypeStr, mmfName, cacheSizeKB);
	vector_int_t* rowOrderLookup = input_readVectorIntFromFile(rowOrderLookupPath);

	char colOrderLookupPath[DEFAULT_STR_BUFF_SIZE];
	sprintf(colOrderLookupPath, "%s/%s/%s_netOrder_%dK", mmfDirPath, partitioningTypeStr, mmfName, cacheSizeKB);
	vector_int_t* colOrderLookup = input_readVectorIntFromFile(colOrderLookupPath);

	// inject ordering info into quintets and then quintets
	// ------------------------------------------------------------------------------------------------------------
	quintet_addOrderingInfoByIndex(quintets_inout, nnz, rowOrderLookup, colOrderLookup);
	// TODO Using self-implemented hybrid sort
	// input_sortQuintetRowValue(quintets_inout, nnz);
	algorithm_parallelMergeQuickSort(quintets_inout, nnz, omp_get_num_threads(), quintet_cmpRowValue);
	quintet_overwriteIndex(quintets_inout, nnz);

	// read partitioning information
	// ------------------------------------------------------------------------------------------------------------
	char rowPartitioningPath[DEFAULT_STR_BUFF_SIZE];
	sprintf(rowPartitioningPath, "%s/%s/%s_vertexDimension_%dK", mmfDirPath, partitioningTypeStr, mmfName, cacheSizeKB);
	vector_int_t* rowPartitioning = input_readPartitionVectorFromFile(rowPartitioningPath);

	// TODO why no col partitioning path??
	sprintf(rowPartitioningPath, "%s/%s/%s_netDimension_%dK", mmfDirPath, partitioningTypeStr, mmfName, cacheSizeKB);
	char colPartitioningPath[DEFAULT_STR_BUFF_SIZE];
	vector_int_t* colPartitioning = input_readPartitionVectorFromFile(rowPartitioningPath);

	// return values
	*rowOrderLookup_out = rowOrderLookup;
	*colOrderLookup_out = colOrderLookup;
	*rowPartitioning_out = rowPartitioning;
	*colPartitioning_out = colPartitioning;
}

void input_powerOf2Order(
		char* mmfDirPath, char* mmfName, int cacheSizeKB,
		quintet_t* quintets_inout, DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount,
		vector_int_t** rowOrderLookup_out, vector_int_t** rowPartitioning_out,
		vector_int_t** colOrderLookup_out, vector_int_t** colPartitioning_out)
{
	// read partitioning stats
	// ------------------------------------------------------------------------------------------------------------
	char partStatsPath[DEFAULT_STR_BUFF_SIZE];
	int partCount = 0;
	sprintf(partStatsPath, "%s/%s_%dKB.stats", mmfDirPath, mmfName, cacheSizeKB);
	input_readPartitionStatFile(partStatsPath, &partCount);

	// read partitioning information
	// ------------------------------------------------------------------------------------------------------------
	char partVectorPath[DEFAULT_STR_BUFF_SIZE];
	sprintf(partVectorPath, "%s/%s_%dKB.partVector", mmfDirPath, mmfName, cacheSizeKB);
	vector_int_t* partVector = input_readPartVector(partVectorPath);

	// generate partitioning info (how many rows a partition has?)
	// ------------------------------------------------------------------------------------------------------------
	partitioning_vector_t* rowPartitioning = partitioning_generatePartitioningFromPartVector(partVector, partCount);

	// generate ordering info (from resulting vector, ith row will be data[i]th row after ordering
	// ------------------------------------------------------------------------------------------------------------
	vector_int_t* rowOrderLookup = partitioning_generateOrderingFromPartVector(partVector, partCount);

	// inject ordering info into quintets and then quintets
	// ------------------------------------------------------------------------------------------------------------
	quintet_addOrderingInfoByIndex(quintets_inout, nnz, rowOrderLookup, NULL);
	algorithm_parallelMergeQuickSort(quintets_inout, nnz, omp_get_num_threads(), quintet_cmpRowValue);
	quintet_overwriteIndex(quintets_inout, nnz);

	// clean up
	vector_int_delete(partVector);

	// return values
	*rowOrderLookup_out = rowOrderLookup;
	*rowPartitioning_out = rowPartitioning;
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

void input_orderVector(vector_real_t** inputAddress_inout, vector_int_t* orderLookup)
{
	if(orderLookup == NULL)
		return;

	vector_real_t* input = *inputAddress_inout;

	if(input->length != orderLookup->length)
		return;

	vector_real_t* inputOrdered = vector_real_new(input->length);
	int i;
	for(i = 0; i < inputOrdered->length; ++i)
		inputOrdered->data[orderLookup->data[i]] = input->data[i];

	vector_real_delete(*inputAddress_inout);
	*inputAddress_inout = inputOrdered;
}

void input_reorderVector(vector_real_t** inputAddress_inout, vector_int_t* orderLookup)
{
	if(orderLookup == NULL)
		return;

	vector_real_t* input = *inputAddress_inout;

	if(input->length != orderLookup->length)
		return;

	vector_real_t* inputOrdered = vector_real_new(input->length);
	int i;
	for(i = 0; i < inputOrdered->length; ++i)
		inputOrdered->data[i] = input->data[orderLookup->data[i]];

	// swap ordered and unordered vector
	vector_real_delete(*inputAddress_inout);
	*inputAddress_inout = inputOrdered;
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif


vector_real_t* input_readVectorRealFromFile(char* filePath)
{
	OPEN_FILE_STREAM(filePath, "r")

	char* line = NULL;
	size_t len = 0;

	getline(&line, &len, f);
	int vecLength = 0;
	sscanf(line, "%d", &vecLength);

	vector_real_t* v = vector_real_new(vecLength);
	int i;
	for(i = 0; i < vecLength; ++i)
	{
		getline(&line, &len, f);

#ifdef SINGLE_PRECISION
		sscanf(line, "%f", &(v->data[i]));
#else
		sscanf(line, "%lf", &(v->data[i]));
#endif
	}

	CLOSE_FILE_STREAM(line, f)

	return v;
}

vector_int_t* input_readVectorIntFromFile(char* filePath)
{
	OPEN_FILE_STREAM(filePath, "r")

	char* line = NULL;
	size_t len = 0;

	getline(&line, &len, f);
	int vecLength = 0;
	sscanf(line, "%d", &vecLength);

	vector_int_t* v = vector_int_new(vecLength);
	int i;
	for(i = 0; i < vecLength; ++i)
	{
		getline(&line, &len, f);
		sscanf(line, "%d", &(v->data[i]));
	}

	CLOSE_FILE_STREAM(line, f)

	return v;
}

vector_int_t* input_readRowPartitionVectorFromPowerOf2File(char* filePath, int rowCount)
{
	OPEN_FILE_STREAM(filePath, "r")

	char* line = NULL;
	size_t len = 0;

	vector_int_t* v = vector_int_new(rowCount);
	int i;
	for(i = 0; i < rowCount; ++i)
	{
		getline(&line, &len, f);
		sscanf(line, "%d", &(v->data[i]));
	}

	CLOSE_FILE_STREAM(line, f)

	return v;
}

// ASSUMPTION: partitioningInfos are accumulated just like CSR rowPtr
// partition-count = partitioningInfo->length - 1;
vector_int_t* input_readPartitionVectorFromFile(char* filePath)
{
	OPEN_FILE_STREAM(filePath, "r")

	char* line = NULL;
	size_t len = 0;

	getline(&line, &len, f);
	int vecLength = 0;
	sscanf(line, "%d", &vecLength);

	vector_int_t* v = vector_int_new(vecLength + 1);
	v->length = vecLength;
	v->data[0] = 0;

	int i;
	for(i = 1; i <= vecLength; ++i)
	{
		getline(&line, &len, f);
		sscanf(line, "%d", &(v->data[i]));
		v->data[i] += v->data[i - 1];
	}

	CLOSE_FILE_STREAM(line, f)

	return v;
}

void input_printVectorRealToFile(char* filePath, vector_real_t* v)
{
	OPEN_FILE_STREAM(filePath, "w")

	fprintf(f, "%d\n", v->length);

	int i;
	for(i = 0; i < v->length; ++i)
		fprintf(f, "%3.2lf\n", v->data[i]);

	CLOSE_FILE_STREAM(NULL, f)
}

void input_printVectorIntToFile(char* filePath, vector_int_t* v)
{
	OPEN_FILE_STREAM(filePath, "w")

	fprintf(f, "%d\n", v->length);

	int i;
	for(i = 0; i < v->length; ++i)
		fprintf(f, "%d\n", v->data[i]);

	CLOSE_FILE_STREAM(NULL, f)
}

void input_readPartitionStatFile(char* filePath, int* partCount_out)
{
	OPEN_FILE_STREAM(filePath, "r")

	char* line = NULL;
	size_t len = 0;

	getline(&line, &len, f);
	sscanf(line, "%d", partCount_out);

	CLOSE_FILE_STREAM(line, f)
}

vector_int_t* input_readPartVector(char* filePath)
{
	int length = input_getLineCount(filePath);
	vector_int_t* partVector = vector_int_new(length);

	OPEN_FILE_STREAM(filePath, "r")

	char* line = NULL;
	size_t len = 0;

	int i;
	for(i = 0; i < length; ++i)
	{
		getline(&line, &len, f);
		sscanf(line, "%d", &partVector->data[i]);
	}

	CLOSE_FILE_STREAM(line, f)

	return partVector;
}
