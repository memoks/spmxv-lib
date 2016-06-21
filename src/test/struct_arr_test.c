/*

#include <stdlib.h>
#include <stdio.h>

#include "include/input_parser.h"


triplet_t* createTriplet(int length)
{
	triplet_t* arr = malloc(sizeof(triplet_t) * length);

	int i;
	for(i = 0; i < length; ++i)
	{
		arr[i].i = i;
		arr[i].j = i;
		arr[i].val = i;
	}

	return arr;
}

int main(void)
{
	int length = 20;
	triplet_t* arr = createTriplet(length);

	int i;
	for(i = 0; i < length; ++i)
	{
		printf("(%d, %d) => %d\n", arr[i].i, arr[i].j, arr[i].val);
	}

	free(arr);

	return EXIT_SUCCESS;
}

*/
