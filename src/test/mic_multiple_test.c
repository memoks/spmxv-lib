
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(void)
{
	int numCoprocessors = _Offload_number_of_devices();
	printf("There are %d coprocessors, offloading code...\n", numCoprocessors);

	int i;
	#pragma omp parallel for num_threads(numCoprocessors) schedule(static, 1)
	for(i = 0; i < numCoprocessors; ++i)
	{
		#pragma offload target(mic:i)
		{
			#pragma omp parallel
			#pragma omp master
				printf("MIC-%d ready with %d threads\n", i, omp_get_num_threads());
		}
	}

	return EXIT_SUCCESS;
}
