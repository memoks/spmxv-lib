
/*

#include <stdlib.h>
#include <stdio.h>

#include "include/scheduler/job_batch.h"
#include "include/data_structure/list.h"


int main(void)
{
	job_batch_t* batch = job_batch_new(10);

	printf("Address of batch: %p\n", batch);
	printf("Address of batch->node: %p \t offset: %lu\n", batch->head, (unsigned long) &((job_batch_t*) 0)->head);
	printf("Address of batch->bytes: %p \t offset: %lu\n", &batch->bytes, (unsigned long) &((job_batch_t*) 0)->bytes);

	// sizeof(index_t*), (unsigned long) &((struct block *)0)->start

	job_batch_delete(batch);

	return EXIT_SUCCESS;
}
*/
