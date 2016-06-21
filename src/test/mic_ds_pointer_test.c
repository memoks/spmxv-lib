
#include <stdlib.h>
#include <stdio.h>

#define ALLOC alloc_if(1)
#define REUSE alloc_if(0)
#define FREE free_if(1)
#define RETAIN free_if(1)

struct test
{
	int length;
	int* data;
};

typedef struct test test_t;

__attribute__((target(mic))) test_t* test_new(int length, int* data)
{
	test_t* test = (test_t*) malloc(sizeof(test_t));
	test->length = length;
	test->data = data;

	return test;
}

__attribute__((target(mic))) void test_delete(test_t* test)
{
	free(test->data);
	// free(test);
}

int main(int argsc, char* argsv[])
{
	// Allocate memory on the host
	__attribute__((target(mic))) test_t* testPtr;
	testPtr = (test_t*) malloc(sizeof(test_t));
	testPtr->length = 10;
	testPtr->data = (int*) malloc(sizeof(int) * 10);

	// fill in some values for array
	int i;
	for(i = 0; i < testPtr->length; ++i)
		testPtr->data[i] = i;

	__attribute__((target(mic))) int length = testPtr->length;
	__attribute__((target(mic))) int* data = testPtr->data;

	#pragma offload target(mic) \
		nocopy(testPtr) \
		in(length) \
		in(data:length(length) ALLOC RETAIN)
	{
		testPtr = test_new(length, data);
		printf("(METHOD) addresses on coprocessor &test %p\t test->data %p\n", &testPtr, &testPtr->data);
	}

	#pragma offload target(mic) nocopy(testPtr)
	{
		printf("(MAIN) addresses on coprocessor &test %p\t test->data %p\n", &testPtr, &testPtr->data);

		int i;
		printf("test: ");
		for(i = 0; i < 10; ++i)
			printf(" %d", testPtr->data[i]);
		printf("\n");

		free(testPtr->data);
		free(testPtr);

		if(testPtr == NULL)
			printf("testPtr is NULL\n");
		else
			printf("testPtr is NOT null\n");
	}

	return EXIT_SUCCESS;
}
