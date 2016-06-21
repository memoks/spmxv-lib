
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(void)
{
	// all the variables here will be defined in coprocessor as well
	#pragma offload_attribute(push, target(mic))

	float alpha = 1.0f;
	float beta = 1.0f;

	int i;
	int length = 100;
	float* x;
	float* y;
	float* out;

	#pragma offload_attribute(pop)

	// --------------------------------------------------------------

	x = malloc(sizeof(int) * length);
	y = malloc(sizeof(int) * length);
	out = malloc(sizeof(int) * length);

	for(i = 0; i < length; ++i)
	{
		x[i] = i;
		y[i] = length - i - 1;
	}

	#pragma offload_transfer target(mic:0) \
		in(x:length(length) alloc_if(1) free_if(0)) \
		in(y:length(length) alloc_if(1) free_if(0))


	#pragma offload target(mic:0) \
		nocopy(x:length(length)) \
		nocopy(y:length(length)) \
		out(out:length(length) alloc_if(1) free_if(1))
	{
		#pragma omp parallel for
		for(i = 0; i < length; ++i)
			out[i] = alpha * x[i] + beta * y[i];
	}

	for(i = 0; i < length; ++i)
		printf("%3.2f ", out[i]);

	printf("\n");
}
