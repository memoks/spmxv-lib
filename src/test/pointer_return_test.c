
#include <stdlib.h>
#include <stdio.h>

/*

void allocArr(int length, int** arr)
{
	*arr = (int*) malloc(sizeof(int) * length);
	int i;
	for(i = 0; i < length; ++i)
		(*arr)[i] = i;
}

int main(void)
{
	int* arr = NULL;
	printf("Address => %p\n", &arr);
	allocArr(10, &arr);

	int i;
	for(i = 0; i < 10; ++i)
		printf("%d\n", arr[i]);

	free(arr);

	return EXIT_SUCCESS;
}

*/
