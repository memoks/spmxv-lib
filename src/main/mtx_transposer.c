
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "include/config.h"
#include "include/io/converter.h"
#include "include/io/input_parser.h"

extern int parseMtxPaths(char* matrices, int* mtxCount_out, char*** mtxPaths_out);
extern int parseHasValues(char* hasValueStr, int* length_out, int** hasValuesArr_out);

int main(int argsc, char* argsv[])
{
	if(argsc < 3)
	{
		printf("USAGE ./bin mtx_path hasValues\n");
		printf("EX> ./bin input/demo.mtx 1\n");
		return EXIT_FAILURE;
	}

	char* mtxPath = argsv[1];
	int hasValues = atoi(argsv[2]);
	printf("%s hasValues: %d\n", mtxPath, hasValues);

	DECIMAL nnz = 0;
	DECIMAL rowCount = 0;
	DECIMAL colCount = 0;
	char* desc = NULL;

	quintet_t* quintets = input_readQuintets(mtxPath, &nnz, &rowCount, &colCount, &desc, 0, hasValues);
	quintet_transpose(quintets, nnz, &rowCount, &colCount);
	quintet_sortRowIndex(quintets, nnz);

	char mtxPathTransposed[4096];
	mtxPathTransposed[0] = '\0';
	strcat(mtxPathTransposed, mtxPath);
	strcat(mtxPathTransposed, "_t");

	input_writeQuintetsInMtxFormat(quintets, nnz, rowCount, colCount, hasValues, desc, mtxPathTransposed);

	// Clean up
	free(quintets);
	free(desc);


	return EXIT_SUCCESS;
}

int parseMtxPaths(char* matrices, int* mtxCount_out, char*** mtxPaths_out)
{
	int mtxCount = 1;
	int i;
	for(i = 0; i < strlen(matrices); ++i)
	{
		if(matrices[i] == ',')
			++mtxCount;
	}

	char** mtxPaths = (char**) malloc(sizeof(char*) * mtxCount);

	int mtxIndex = 0;
	char* str = strtok(matrices, ",");
	while(str != NULL)
	{
		mtxPaths[mtxIndex] = (char*) malloc(sizeof(char) * (strlen(str) + 1));
		mtxPaths[mtxIndex] = str;

		++mtxIndex;
		str = strtok(NULL, ",");
	}

	*mtxPaths_out = mtxPaths;
	*mtxCount_out = mtxCount;

	return EXIT_SUCCESS;
}

int parseHasValues(char* hasValueStr, int* length_out, int** hasValuesArr_out)
{
	int mtxCount = 1;
	int i;
	for(i = 0; i < strlen(hasValueStr); ++i)
	{
		if(hasValueStr[i] == ',')
			++mtxCount;
	}

	int* hasValuesArr = (int*) malloc(sizeof(int) * mtxCount);

	int mtxIndex = 0;
	char* str = strtok(hasValueStr, ",");
	while(str != NULL)
	{
		hasValuesArr[mtxIndex] = atoi(str);

		++mtxIndex;
		str = strtok(NULL, ",");
	}

	*hasValuesArr_out = hasValuesArr;
	*length_out = mtxCount;

	return EXIT_SUCCESS;
}
