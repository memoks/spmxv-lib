/*

#include <stdlib.h>
#include <stdio.h>

#include "include/config.h"
#include "include/scheduler/job_batch.h"
#include "include/scheduler/job_queue.h"
#include "include/scheduler/scheduler.h"


int main(void)
{
	int i;
	int j;
	int numQueues = 4;
	job_queue_t** queues = (job_queue_t**) malloc(sizeof(job_queue_t*) * numQueues);

	for(i = 0; i < numQueues; ++i)
	{
		queues[i] = job_queue_new(i);
		for(j = 0; j < 10; ++j)
			job_queue_addBatch(queues[i], job_batch_new(i));
	}

	job_queue_printMultiple("Printing initial states...", queues, numQueues);

	queues[1]->victimId = 2;
	scheduler_localityAwareSteal_multi(1, numQueues, queues);

	job_queue_printMultiple("Printing states after steal from queue-2...", queues, numQueues);


	job_queue_deleteMultiple(queues, numQueues);

	return EXIT_SUCCESS;
}
*/
