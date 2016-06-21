
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(int argsc, char* argsv[])
{
	#ifdef __INTEL_OFFLOAD

		printf("Using offload... \n");
		fflush(0);

		#pragma offload target(mic)
		{
			#pragma omp parallel
			#pragma omp master
			{
				printf("Offloaded code is running on coprocessor with %d threads...\n", omp_get_num_threads());
				fflush(0);
			}
		}

		#pragma omp parallel
		#pragma omp master
		{
			printf("Meanwhile on the processor with %d threads...\n", omp_get_num_threads());
			fflush(0);
		}
	#else

		#ifdef __MIC__

			#pragma omp parallel
			#pragma omp master
			{
				printf("Using coprocessor native execution with %d threads... \n", omp_get_num_threads());
				fflush(0);
			}
		#else

			#pragma omp parallel
			#pragma omp master
			{
				printf("Using only processor with %d threads... \n", omp_get_num_threads());
				fflush(0);
			}

		#endif

	#endif

	return EXIT_SUCCESS;
}

