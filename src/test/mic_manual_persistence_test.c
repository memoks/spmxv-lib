
#include <stdlib.h>
#include <stdio.h>

#define ALLOC alloc_if(1)
#define REUSE alloc_if(0)
#define FREE free_if(1)
#define RETAIN free_if(0)

__attribute__((target(mic))) void printData(char* m, int* data, int length);

int main(void)
{
	__attribute__((target(mic))) int* data;
	__attribute__((target(mic))) int length;
	int i;

	length = 10;
	data = (int*) malloc(sizeof(int) * length);
	for(i = 0; i < length; ++i)
		data[i] = i;

	printData("On CPU", data, length);

	#pragma offload target(mic) \
		in(length) \
		nocopy(data)
	{
		printf("(OFF-1) data: %p\t&data: %p\n", data, &data);
		data = (int*) malloc(sizeof(int) * length);

		int i;
		for(i = 0; i < length; ++i)
			data[i] = length - i;
	}

	#pragma offload target(mic) \
		nocopy(length) \
		nocopy(data)
	{
		printf("(OFF-2) data: %p\t&data: %p\n", data, &data);
		printData("On Coprocessor", data, length);

	}

	#pragma offload target(mic) \
		nocopy(length) \
		nocopy(data)
	{
		printf("(OFF-3) data: %p\t&data: %p\n", data, &data);
	}

	#pragma offload target(mic) \
		nocopy(data)
	{
		free(data);
		printf("(OFF-4) data: %p\t&data: %p\n", data, &data);
	}
}

__attribute__((target(mic))) void printData(char* m, int* data, int length)
{
	int i;
	printf("%s (%d):", m, length);
	for(i = 0; i < length; ++i)
		printf(" %d", data[i]);
	printf("\n");
}
