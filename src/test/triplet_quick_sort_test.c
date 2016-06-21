/*

#include <stdlib.h>
#include <stdio.h>

#include "include/input_parser.h"

int main(void)
{
	int nnz = 0;
	int rowCount = 0;
	int colCount = 0;
	triplet_t* triplets = input_readTriplets("input/test/spm_ptr_test_small", &nnz, &rowCount, &colCount);

	input_qsortRowCol(triplets, nnz);

	int i;
	for(i = 0; i < nnz; ++i)
	{
		printf("(%d, %d) => %d\n", triplets[i].i, triplets[i].j, triplets[i].val);
	}

	return EXIT_SUCCESS;
}

*/
