/*

#include <stdlib.h>
#include <stdio.h>

#include "include/scheduler/job_batch.h"
#include "include/scheduler/job_queue.h"

int main(void)
{
	job_queue_t* jq = job_queue_new(0);

	int i;
	for(i = 0; i < 10; ++i)
		job_queue_addBatch(jq, job_batch_new(i + 1));

	job_batch_t* startBatch;
	int count;
	job_queue_getLastBatchHalf(jq, &startBatch, &count);
	printf("Printing stolen chunk: %d\n", count);
	list_head_t* head = &startBatch->head;
	for(i = 0; i < count; ++i)
	{
		job_batch_print(job_batch_getBatch(head));
		head = head->next;
	}

	printf("Printing job_queue:\n");
	job_queue_print(jq);

	job_batch_t* batch = job_queue_getFirstBatch(jq);
	printf("Executing first ");
	job_batch_print(batch);

	batch = job_queue_getLastBatch(jq);
	printf("Executing last ");
	job_batch_print(batch);

	job_queue_getLastBatchHalf(jq, &startBatch, &count);


	job_batch_t* stolen1 = job_batch_new(i + 1);
	job_batch_t* stolen2 = job_batch_new(i + 2);
	job_queue_addStolenBatch(jq, stolen1);
	job_queue_addStolenBatch(jq, stolen2);
	job_queue_removeStolenBatch(jq);
	job_queue_print(jq);

	job_queue_delete(jq);
	job_batch_delete(stolen1);
	job_batch_delete(stolen2);

	return EXIT_SUCCESS;
}

*/
