
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

__attribute__((target(mic))) calculation(float* res, float* x, float* y, int length);

int main(void)
{
	// Arrays allocated here will reside
	// in processor not coprocessor !!
	#pragma offload_attribute(push, target(mic))

	int i;
	int offLen = 240;
	float* x = malloc(sizeof(float) * offLen);
	float* y = malloc(sizeof(float) * offLen);
	float* z = malloc(sizeof(float) * offLen);
	float* zProc = malloc(sizeof(float) * offLen);

	#pragma offload_attribute(pop)

	int signalVal;
	#pragma offload target(mic:0) signal(&signalVal) \	cd
		nocopy(x:length(offLen) alloc_if(1) free_if(1)) \
		nocopy(y:length(offLen) alloc_if(1) free_if(1)) \
		out(z:length(offLen) alloc_if(1) free_if(1))
	{
		calculation(z, x, y, offLen);
	}

	// concurrent calculation on the processor
	calculation(zProc, x, y, offLen);
	printf("Processor done, waiting for the coprocessor...\n");

	// unlike "offload", "offload_*" clauses are "0" code clauses
	#pragma offload_wait target(mic:0) wait(&signalVal)

	printf("Results...\n");
	for(i = 0; i < offLen; ++i)
		printf("%3.2f %3.2f\n", z[i], zProc[i]);
	printf("\n");

	if(x != NULL) free(x);
	if(y != NULL) free(y);
	if(z != NULL) free(z);

	printf("Done!!\n");

	return EXIT_SUCCESS;
}

__attribute__((target(mic))) calculation(float* res, float* x, float* y, int length)
{
	int i;
	#pragma omp parallel for
	for(i = 0; i < length; ++i)
	{
		x[i] = i;
		y[i] = length - i;
	}

	#pragma omp parallel for
	for(i = 0; i < length; ++i)
	{
		res[i] = x[i] + y[i];
	}
}

