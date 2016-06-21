#include <stdio.h>
#include <stdlib.h>

#define SIZE 10

#define ALLOC alloc_if(1) free_if(0)
#define REUSE alloc_if(0) free_if(0)
#define FREE  alloc_if(0) free_if(1)

// Example of Non-Bitwise Object Transfer, All Data Elements Needed
typedef struct
{
    int m1;
    int *m2;
} nbwcs;

void send_inputs(nbwcs* struct1)
{
	int m1;
	int *m2;

	// Initialize the struct
	struct1 = (nbwcs*) malloc(sizeof(nbwcs));
	struct1->m1 = 10;
	struct1->m2 = (int*) malloc(SIZE * sizeof(int));
	for (int i=0; i<SIZE; i++)
	{
		struct1->m2[i] = i;
	}

	// In this offload data is transferred
	m1 = struct1->m1;
	m2 = struct1->m2;
	#pragma offload target(mic:0) in(m1) in(m2[0:SIZE] : ALLOC) nocopy(struct1:length(1) ALLOC)
	{
		struct1->m1 = m1;
		struct1->m2 = m2;
		printf("MIC offload1: struct1.m2[0] = %d, struct1.m2[SIZE-1] = %d\n", struct1->m2[0], struct1->m2[SIZE-1]);
		fflush(0);
	}
}

void use_the_data(nbwcs* struct1)
{
	// In this offload data is used and updated
	#pragma offload target(mic:0) nocopy(struct1)
	{
		for (int i=0; i<SIZE; i++)
		{
			struct1->m2[i] += i;
		}

		printf("MIC offload2: struct1.m2[0] = %d, struct1.m2[SIZE-1] = %d\n", struct1->m2[0], struct1->m2[SIZE-1]);
		fflush(0);
	}
}

void receive_results(nbwcs* struct1)
{
	int *m2;
	// In this offload data is used,, updated, freed on MIC and brought back to the CPU
	m2 = struct1->m2;
	#pragma offload target(mic:0) out(m2[0:SIZE] : FREE) nocopy(struct1:length(1) FREE)
	{
		for (int i=0; i<SIZE; i++)
		{
			struct1->m2[i] += i;
		}
		printf("MIC offload3: struct1.m2[0] = %d, struct1.m2[SIZE-1] = %d\n", struct1->m2[0], struct1->m2[SIZE-1]);
		fflush(0);
	}

	printf("CPU: struct1.m2[0] = %d, struct1.m2[SIZE-1] = %d\n", struct1->m2[0], struct1->m2[SIZE-1]);
}

int main()
{
	__declspec(target(mic)) nbwcs* struct1;

	send_inputs(struct1);
	use_the_data(struct1);
	receive_results(struct1);
	return 0;
}
