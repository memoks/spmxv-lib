
/*
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "mkl.h"


struct triplet
{
	MKL_INT i; 		// row index when ordered
	MKL_INT j;		// column index when ordered
	double val;		// entry value
};

typedef struct triplet triplet_t;

struct spm_csr
{
	double* values;
	int* colInd;
	int* rowPtr;

	int nnz;
	int rowCount;
	int colCount;

	int startRow;
	int startCol;
};

typedef struct spm_csr spm_csr_t;


// Some macros to save same space and hide code blocks
// -----------------------------------------------------------------------------

// Opens a file stream.
#define OPEN_FILE_STREAM(fileName, mode)			\
	FILE *f = fopen(fileName, mode);				\
	if(f == NULL)									\
	{												\
		printf("Cannot find file: %s\n", fileName);	\
		exit(EXIT_FAILURE);							\
	}

// Closes a file stream and a quick clean-up service.
#define CLOSE_FILE_STREAM(line, f)					\
	if(line) free(line);							\
	fclose(f);

// -----------------------------------------------------------------------------

extern triplet_t* input_readTriplets(char* fileName, int *nnz_out, int* rowCount_out, int* colCount_out);
extern void input_sortTripletIndex(triplet_t* triplets, int length);
extern int input_cmpTripletIndex(const void* left, const void* right);
extern spm_csr_t* spm_csr_new(int rowCount, int colCount, int nnz, int startRow, int startCol);


int main(int argsc, char* argsv[])
{
	if(argsc < 1)
	{
		printf("Usage: ./bin <mtx_file_name>\n");
	}

	int nnz;
	int rowCount;
	int colCount;
	triplet_t* triplets = input_readTriplets(argsv[1], &nnz, &rowCount, &colCount);
	spm_csr_t* spm = input_tripletToCSR(triplets, nnz, rowCount, colCount);

	mkl_mic_enable();

	// offload all the work to device-0
	mkl_mic_set_workdivision(MKL_TARGET_MIC, 0, 1);


	mkl_cspblas_dcsrgemv('N', spm->rowCount, spm->values, spm->rowPtr, , vecin, vecout);


	return EXIT_SUCCESS;
}


triplet_t* input_readTriplets(char* fileName, int *nnz_out, int* rowCount_out, int* colCount_out)
{
	OPEN_FILE_STREAM(fileName, "r")

	char* line = NULL;
	size_t len = 0;
	size_t read = 0;

	// read first line
	getline(&line, &len, f);
	sscanf(line, "%d %d %d", rowCount_out, colCount_out, nnz_out);

	triplet_t* tripletArr = (triplet_t*) malloc(sizeof(triplet_t) * *nnz_out);

	int i = 0;
	int j = 0;
	int index = 0;
	while((read = getline(&line, &len, f)) != -1)
	{
		sscanf(line, "%d %d", &i, &j);
		tripletArr[index].i = i - 1;
		tripletArr[index].j = j - 1;
		tripletArr[index].val = 1.0;
		++index;
	}

	CLOSE_FILE_STREAM(line, f)

	input_sortTripletIndex(tripletArr, *nnz_out);

	return tripletArr;
}

int input_cmpTripletIndex(const void* left, const void* right)
{
	triplet_t* t1 = (triplet_t*) left;
	triplet_t* t2 = (triplet_t*) right;

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

void input_sortTripletIndex(triplet_t* triplets, int length)
{
	qsort(triplets, length, sizeof(triplet_t), input_cmpTripletIndex);
}

spm_csr_t* input_tripletToCSR(triplet_t* triplets, int nnz, int rowCount, int colCount)
{
	return input_tripletToCSRPartial(triplets, 0, rowCount, nnz, rowCount, colCount);
}

spm_csr_t* input_tripletToCSRPartial(triplet_t* triplets,
		int rowStart, int rowEnd, int nnz, int totalRowCount, int totalColCount)
{
	MKL_INT startIndex = 0;
	if(rowStart != 0)
		startIndex = input_binarySearchRow(triplets, nnz, rowStart);

	MKL_INT endIndex = nnz;
	if(rowEnd != totalRowCount)
		endIndex = input_binarySearchRow(triplets, nnz, rowEnd);

	spm_csr_t* spm = spm_csr_new(rowEnd - rowStart, totalColCount, endIndex - startIndex, rowStart, 0);

	int i;
	int row = 0;
	int prevRow = -1;
	int rowIndex = 0;
	for(i = 0; i < endIndex - startIndex; ++i)
	{
		row = triplets[startIndex + i].i;

		if(row != prevRow)
		{
			prevRow = row;
			spm->rowPtr[rowIndex] = i;
			++rowIndex;
		}

		spm->colInd[i] = triplets[i + startIndex].j;
		spm->values[i] = triplets[i + startIndex].val;
	}

	spm->rowPtr[rowIndex] = endIndex - startIndex;

	return spm;
}

spm_csr_t* spm_csr_new(int rowCount, int colCount, int nnz, int startRow, int startCol)
{
	spm_csr_t* spm = (spm_csr_t*) malloc(sizeof(spm_csr_t));

	spm->values = (double*) mkl_malloc(sizeof(double) * nnz, 64);
	spm->colInd = (MKL_INT*) mkl_malloc(sizeof(MKL_INT) * nnz, 64);
	spm->rowPtr = (MKL_INT*) mkl_malloc(sizeof(MKL_INT) * (rowCount + 1), 64);

	spm->nnz = nnz;
	spm->rowCount = rowCount;
	spm->colCount = colCount;

	spm->startRow = startRow;
	spm->startCol = startCol;

	return spm;
}

*/
