/*

#include <stdlib.h>
#include <stdio.h>

#include "include/input_parser.h"
#include "include/data_structure/block_tree.h"
#include "include/scheduler/scheduler.h"


int main(void)
{
	int nnz = 0;
	int rowCount = 0;
	int colCount = 0;
	triplet_t* triplets = input_readOrderedTriplets(
			"input/test/spm_ptr_test_small", "input/test/testa",
			&nnz, &rowCount, &colCount);


	int numBlocks = 2;
	vector_int_t* blockData = input_readBlockFile("input/test/testb");
	vector_int_t* rowBlockArr = block_createRowBlockArr(blockData, numBlocks);
	spm_csr_t** partialSpms = scheduler_distributeMatrix(triplets, nnz, rowCount, colCount, rowBlockArr, NULL, numBlocks);

	printf("Created.\n");

	int i;
	for(i = 0; i < numBlocks; ++i)
	{
		spm_csr_t* spm = partialSpms[i];
		spm_print(spm);
	}



	for(i = 0; i < numBlocks; ++i)
		spm_delete(partialSpms[i]);

	free(partialSpms);
	vector_int_delete(blockData);
	vector_int_delete(rowBlockArr);
	free(triplets);



	return EXIT_SUCCESS;
}
*/
