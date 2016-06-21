/*

#include <stdlib.h>
#include <stdio.h>

#include "include/scheduler/job_batch.h"
#include "include/scheduler/job_queue.h"


int main(void)
{
	job_queue_t* jq = job_queue_new(0, 10);

	int i;
	for(i = 0; i < 10; ++i)
		job_queue_addBatch(jq, job_batch_new(i + 1));

	job_batch_t* removedBatch = job_queue_removeFirstBatch(jq);
	printf("Removed Batch: ");
	job_batch_print(removedBatch);

	job_queue_print(jq);

	job_batch_delete(removedBatch);
	job_queue_delete(jq);

	return EXIT_SUCCESS;
}
*/
