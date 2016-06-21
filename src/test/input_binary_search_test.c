/*

#include <stdlib.h>
#include <stdio.h>

#include "include/data_structure/vector.h"
#include "include/data_structure/spm.h"
#include "include/input_parser.h"


int main(void)
{
	int nnz;
	int rowCount;
	int colCount;

	triplet_t* triplets = input_readTriplets("input/test/spm_ptr_test_small",
			&nnz, &rowCount, &colCount);
	input_qsortRowColIndex(triplets, nnz);

	input_readPartVectorLookup("input/test/part_vec_test_small", triplets, nnz, rowCount, colCount);
	input_qsortRowColValue(triplets, nnz);


	printf("NNZ: %d Row-count: %d Col-count: %d\n", nnz, rowCount, colCount);
	int i;
	for(i = 0; i < nnz; ++i)
	{
		printf("%d - (%d, %d) (%d, %d) => %3.2f\n", i, triplets[i].i, triplets[i].j,
				triplets[i].iVal, triplets[i].jVal, triplets[i].val);
	}

	int searchIndex = 5;
	printf("Testing binary search for rowIndex=%d... ", searchIndex);
	printf("done\nResult=%d\n", input_binarySearchRow(triplets, nnz, searchIndex));

	searchIndex = 0;
	printf("Testing binary search for rowIndex=%d... ", searchIndex);
	printf("done\nResult=%d\n", input_binarySearchRow(triplets, nnz, searchIndex));

	searchIndex = 6;
	printf("Testing binary search for rowIndex=%d... ", searchIndex);
	printf("done\nResult=%d\n", input_binarySearchRow(triplets, nnz, searchIndex));

	free(triplets);

	return EXIT_SUCCESS;
}

*/
