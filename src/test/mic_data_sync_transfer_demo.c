
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(void)
{
	printf("demo\n");

	// Arrays allocated here will reside in coprocessor
	#pragma offload_attribute(push, target(mic))

	int i;
	int offLen = 240;
	float* x = NULL;// = malloc(sizeof(float) * offLen);
	float* y = NULL;// = malloc(sizeof(float) * offLen);
	float* z = malloc(sizeof(float) * offLen);

	#pragma offload_attribute(pop)

	printf("Allocation...\n");

	// when nothing is specified,
	// default behaviour will be: alloc_if(0) free_if(0)
	// no need to specify free_if here
	#pragma offload target(mic) \
		nocopy(x:length(offLen) alloc_if(1)) \
		nocopy(y:length(offLen) alloc_if(1))
	{
		#pragma omp parallel for
		for(i = 0; i < offLen; ++i)
		{
			x[i] = i;
			y[i] = offLen - i;
		}
	}

	printf("Offloaded code...\n");

	// no need to specify free_if for x and y because of the default behaviour
	#pragma offload target(mic) \
		nocopy(x:length(offLen) free_if(1)) \
		nocopy(y:length(offLen) free_if(1)) \
		out(z:length(offLen) alloc_if(1) free_if(1))
	{
		#pragma omp parallel for
		for(i = 0; i < offLen; ++i)
		{
			z[i] = x[i] + y[i];
		}
	}

	// no need to specify free_if for x and y because of the default behaviour
	#pragma offload target(mic) \
		nocopy(x:length(offLen) free_if(1)) \
		nocopy(y:length(offLen) free_if(1)) \
		out(z:length(offLen) alloc_if(1) free_if(1))
	{
		#pragma omp parallel for
		for(i = 0; i < offLen; ++i)
		{
			z[i] = x[i] + y[i];
		}
	}

	printf("Printing result...\n");
	for(i = 0; i < offLen; ++i)
		printf("%3.2f ", z[i]);
	printf("\n");

	if(x != NULL) free(x);
	if(y != NULL) free(y);
	if(z != NULL) free(z);

	return EXIT_SUCCESS;
}
