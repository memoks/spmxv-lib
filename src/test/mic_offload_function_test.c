
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

__attribute__((target(mic))) offloadFunc(char* message);

/*
This code doesn't output as expected. However by looking at
the execution time it seems to be working. Must be something
to do with interrupt handling mechanism.
*/

int main(void)
{
        int signalVar;

        // coprocessor activity
        #pragma offload target(mic:0) signal(&signalVar)
        {
                offloadFunc("COPROCESSOR");
        }

        // concurrent processor activity
        offloadFunc("PROCESSOR_ASYNC");

        // wait for coprocessor to finish
        #pragma offload_wait target(mic:0) wait(&signalVar)

        offloadFunc("PROCESSOR_SYNC");
}


__attribute__((target(mic))) offloadFunc(char* message)
{
        #pragma omp parallel
        #pragma omp master
                printf("%s: %d available threads\n", message, omp_get_num_threads());
}
