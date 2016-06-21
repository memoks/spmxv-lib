
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "include/config.h"
#include "include/io/converter.h"
#include "include/data_structure/list.h"

#include "include/scheduler/job_queue.h"

job_queue_t* job_queue_new(int stealTreshold)
{
	job_queue_t* jq = (job_queue_t*) malloc(sizeof(job_queue_t));
	job_queue_init(jq, stealTreshold);

	return jq;
}

void job_queue_init(job_queue_t* queue, int stealTreshold)
{
	queue->batchCount = 0;
	omp_init_lock(&queue->writeLock);
	queue->batchHead = job_batch_new();

	queue->executableBatchCount = 0;
	queue->execPtrStart = queue->batchHead;
	queue->execPtrEnd = queue->batchHead;

	queue->stealTreshold = stealTreshold;
	queue->stolenBatchCount = 0;
	queue->stolenBatchList = lg_new();

	queue->state = JOB_QUEUE_STATE_INITIAL;
}

void job_queue_deleteNonPtr(job_queue_t* queue)
{
	omp_destroy_lock(&queue->writeLock);
	job_batch_deleteAll(queue->batchHead);
	job_batch_delete(queue->batchHead);
	lg_deleteShallow(queue->stolenBatchList);
}

void job_queue_delete(job_queue_t* queue)
{
	job_queue_deleteNonPtr(queue);
	free(queue);
}

void job_queue_printDetailed(job_queue_t* queue, void* spm)
{
	job_queue_printAttributes(queue);

	DECIMAL totalNNZ = 0;
	list_head_t* curr;
	list_for_each(curr, &queue->batchHead->head)
	{
		job_batch_t* currBatch = job_batch_getBatch(curr);
		job_batch_printDetailed(currBatch, spm);

		totalNNZ += job_batch_getNNZ(currBatch, spm);
	}

	PRINTF("TOTAL NNZ: %d\n", totalNNZ);
}

void job_queue_printAttributes(job_queue_t* queue)
{
	PRINTF("[JOB_QUEUE]\n");
	PRINTF("Local batch count: %d, Steal Treshold: %d.\n",
			queue->batchCount, queue->stealTreshold);
}

void job_queue_print(job_queue_t* queue)
{
	job_queue_printAttributes(queue);
	list_head_t* curr;
	list_for_each(curr, &queue->batchHead->head)
	{
		job_batch_t* currBatch = job_batch_getBatch(curr);
		job_batch_print(currBatch);
	}
}

void job_queue_printCurrentState(job_queue_t* queue, void* spm)
{
	job_queue_printDetailed(queue, spm);

	PRINTF("Executable batch count: %d\n", queue->executableBatchCount);
	job_batch_t* currExecBatch = queue->execPtrStart;
	int i = 0;
	while(currExecBatch->head.next != &queue->execPtrEnd->head)
	{
		currExecBatch = job_batch_getBatch(currExecBatch->head.next);
		printf("%d. ", i++);
		job_batch_printDetailed(currExecBatch, spm);
	}

	lg_job_batch_print("Stolen batches", queue->stolenBatchList);
	PRINTF("\n");
}

void job_queue_copy(job_queue_t* source, job_queue_t* destination)
{
	if(source == NULL)
		return;

	destination->batchCount = source->batchCount;
	omp_init_lock(&destination->writeLock);
	destination->batchHead = source->batchHead;

	destination->executableBatchCount = source->executableBatchCount;
	destination->execPtrStart = source->execPtrStart;
	destination->execPtrEnd = source->execPtrEnd;

	destination->stealTreshold = source->stealTreshold;
	destination->stolenBatchCount = source->stolenBatchCount;
	destination->stolenBatchList = source->stolenBatchList;

	destination->state = source->state;
}

inline void job_queue_reset(job_queue_t* queue)
{
	queue->execPtrStart = queue->batchHead;
	queue->execPtrEnd = queue->batchHead;
	queue->executableBatchCount = queue->batchCount;
	queue->state = JOB_QUEUE_STATE_INITIAL;

	queue->stolenBatchCount = 0;
	lg_removeAll(queue->stolenBatchList);
}

// SpMxV data preparation functions
// -------------------------------------------------------------------------------------------------------------------

void job_queue_fillFromSubMtxTree_rowWiseColNet_CSR(spm_cmp_t* spmCsr, job_queue_t* jq, tree_node_t* head)
{
	lg_t* jobBatchList = NULL;
	job_batch_extractJobBatchListFromSubMtxTree_rowWiseColNet_CSR(spmCsr, head, &jobBatchList);
	job_queue_fillFromJobBatchList(jq, jobBatchList);

	// clean up
	lg_job_batch_deleteDeep(jobBatchList);
}

void job_queue_fillFromSubMtxTreeScatter_rowWiseColNet_CSR(
		spm_cmp_t* spmCsr, job_queue_t** jqPtrArr, int jqCount, tree_node_t* head)
{
	lg_t* jobBatchList = NULL;
	job_batch_extractJobBatchListFromSubMtxTree_rowWiseColNet_CSR(spmCsr, head, &jobBatchList);
	job_queue_fillFromJobBatchListScatter(jqPtrArr, jqCount, jobBatchList);

	// clean up
	lg_job_batch_deleteDeep(jobBatchList);
}

void job_queue_fillFromJobBatchList(job_queue_t* jq, lg_t* jobBatchList)
{
	if(jobBatchList == NULL)
		return;

	list_head_t* currListHead = NULL;
	list_head_t* tempListHead = NULL;
	list_for_each_safe(currListHead, tempListHead, &jobBatchList->headNode->listHead)
	{
		lng_t* currListNode = lng_getListGeneric(currListHead);
		job_batch_t* jobBatch = (job_batch_t*) lng_remove(currListNode);
		job_queue_addBatchWithoutLocking(jq, jobBatch);
	}
}

void job_queue_fillFromJobBatchListScatter(
		job_queue_t** jqPtrArr, int jqCount, lg_t* jobBatchList)
{
	if(jobBatchList == NULL)
		return;

	int currIndex = 0;
	list_head_t* currListHead = NULL;
	list_head_t* temp = NULL;
	list_for_each_safe(currListHead, temp, &jobBatchList->headNode->listHead)
	{
		lng_t* currListNode = lng_getListGeneric(currListHead);
		job_batch_t* jobBatch = lng_remove(currListNode);

		job_queue_addBatch(jqPtrArr[currIndex], jobBatch);

		++currIndex;
		currIndex = currIndex % jqCount;
	}
}

void job_queue_fillFromParitioningVector_rowWiseColNet_CSR(
		spm_cmp_t* spmCsr, job_queue_t* jq, vector_int_t* rowPartitioning,
		int startIndex, int endIndex)
{
	lg_t* jobBatchList = NULL;
	job_batch_extractJobBatchListFromPartitionVector_rowWiseColNet_CSR(
			spmCsr, rowPartitioning, startIndex, endIndex, &jobBatchList);

	job_queue_addJobBatchListWithoutLocking(jq, jobBatchList);

	// clean up
	lg_deleteShallow(jobBatchList);
}

void job_queue_fillFromParitioningVectorScatter_rowWiseColNet_CSR(
		spm_cmp_t* spmCsr, job_queue_t** jobQueuePtrArr, int jobQueueCount,
		vector_int_t* rowPartitioning, int startIndex, int endIndex)
{
	lg_t* jobBatchList = NULL;
	job_batch_extractJobBatchListFromPartitionVector_rowWiseColNet_CSR(
			spmCsr, rowPartitioning, startIndex, endIndex, &jobBatchList);

	job_queue_fillFromJobBatchListScatter(jobQueuePtrArr, jobQueueCount, jobBatchList);
	lg_deleteShallow(jobBatchList);
}

void job_queue_convertToJDSPartial(
		spm_cmp_t* spmCsrCounterpart, DECIMAL* permutation, job_queue_t* jobQueue_inout)
{
	list_head_t* head = &jobQueue_inout->batchHead->head;
	list_head_t* curr = NULL;
	list_for_each(curr, head)
	{
		job_batch_t* currBatch = job_batch_getBatch(curr);

		// since we will overwrite job_batch, let's copy the sub-matrix so that it will be healthier
		sub_mtx_dim_t currSubMtx;
		sub_mtx_copy(&currBatch->data.subMtx, &currSubMtx);

		// create partial JDS of current sub-matrix
		spm_jds_t* spmJds = NULL;
		converter_CSRPartToJDS(spmCsrCounterpart, &currSubMtx, permutation, &spmJds);

		// convert (or overwrite in this case) but do not touch list structure
		currBatch->type = JOB_PARTIAL_JDS;
		sub_mtx_copy(&currSubMtx, &currBatch->data.partialJds.subMtx);
		currBatch->data.partialJds.spmJds = spmJds;
	}
}

void job_queue_convertFromPartialJDSToHybridJDSCSR(
		job_queue_t* partialJdsJobQueue,
		job_queue_t* hybridJdsCsrJobQueue_inout)
{
	list_head_t* head = &partialJdsJobQueue->batchHead->head;
	list_head_t* curr = NULL;
	list_for_each(curr, head)
	{
		job_batch_t* partialJdsBatch = job_batch_getBatch(curr);

		// extract hybrid job batch from partial JDS batch
		job_batch_t* hybridJdsCsrBatch = job_batch_new();
		converter_extractHybridJDSCSR(partialJdsBatch, hybridJdsCsrBatch);

		// add extracted batch to new queue
		job_queue_addBatchWithoutLocking(hybridJdsCsrJobQueue_inout, hybridJdsCsrBatch);
	}
}

// Job Queue Runtime Management Functions
// -------------------------------------------------------------------------------------------------------------------

inline int job_queue_isEmpty(job_queue_t* queue)
{
	return &queue->batchHead->head == queue->batchHead->head.next;
}

inline int job_queue_hasExecutableBatch(job_queue_t* queue)
{
	return queue->execPtrStart->head.next != &queue->execPtrEnd->head;
}

inline int job_queue_hasStolenBatch(job_queue_t* queue)
{
	return !lg_isEmpty(queue->stolenBatchList);
}

inline void job_queue_emptyQueue(job_queue_t* queue)
{
	job_batch_deleteAllShallow(queue->batchHead);
	queue->batchCount = 0;
	queue->executableBatchCount = 0;
}


void job_queue_addBatchWithoutLocking(job_queue_t* jq, job_batch_t* batch)
{
	job_batch_addBatchTail(batch, jq->batchHead);
	++jq->batchCount;
	++jq->executableBatchCount;
}

void job_queue_addBatch(job_queue_t* jq, job_batch_t* batch)
{
	omp_set_lock(&jq->writeLock);
	{
		job_queue_addBatchWithoutLocking(jq, batch);
	}
	omp_unset_lock(&jq->writeLock);
}

void job_queue_addJobBatchListWithoutLocking(job_queue_t* queue, lg_t* jobBatchList)
{
	list_head_t* currListHead = NULL;
	list_for_each(currListHead, &jobBatchList->headNode->listHead)
	{
		lng_t* currListNode = lng_getListGeneric(currListHead);
		job_queue_addBatchWithoutLocking(queue, (job_batch_t*) currListNode->dataPtr);
	}
}

job_batch_t* job_queue_getFirstBatchWithoutLocking(job_queue_t* jq)
{
	job_batch_t* executable = NULL;

	if(job_queue_hasExecutableBatch(jq))
	{
		executable = job_batch_getNext(jq->execPtrStart);
		jq->execPtrStart = executable;
		--jq->executableBatchCount;
	}

	return executable;
}

job_batch_t* job_queue_getFirstBatch(job_queue_t* jq)
{
	job_batch_t* executable = NULL;
	omp_set_lock(&jq->writeLock);
	{
		if(job_queue_hasExecutableBatch(jq))
		{
			executable = job_batch_getNext(jq->execPtrStart);
			jq->execPtrStart = executable;
			--jq->executableBatchCount;
		}
	}
	omp_unset_lock(&jq->writeLock);

	return executable;
}

job_batch_t* job_queue_removeFirstBatch(job_queue_t* jq)
{
	job_batch_t* batch = NULL;
	omp_set_lock(&jq->writeLock);
	{
		if(!job_queue_isEmpty(jq))
		{
			batch = job_batch_getNext(jq->batchHead);
			list_del(&batch->head);
			--jq->batchCount;
		}
	}
	omp_unset_lock(&jq->writeLock);

	return batch;
}

job_batch_t* job_queue_getLastBatch(job_queue_t* jq)
{
	job_batch_t* executable = NULL;
	omp_set_lock(&jq->writeLock);
	{
		if(job_queue_hasExecutableBatch(jq))
		{
			executable = job_batch_getPrev(jq->execPtrEnd);
			jq->execPtrEnd = executable;
			--jq->executableBatchCount;
		}
	}
	omp_unset_lock(&jq->writeLock);

	return executable;
}

void job_queue_getLastBatchHalf(job_queue_t* jq, job_batch_t** start_out, int* count_out)
{
	omp_set_lock(&jq->writeLock);
	{
		if(jq->executableBatchCount > jq->stealTreshold)
		{
			*count_out = jq->executableBatchCount / 2;
			jq->executableBatchCount -= *count_out;

			list_head_t* head = &jq->execPtrEnd->head;
			int i;
			for(i = 0; i < *count_out; ++i)
				head = head->prev;

			*start_out = job_batch_getBatch(head);
			jq->execPtrEnd = *start_out;
		}
		else
		{
			*start_out = NULL;
			*count_out = 0;
		}
	}
	omp_unset_lock(&jq->writeLock);
}

void job_queue_addStolenBatchWithoutLocking(job_queue_t* jq, job_batch_t* batch)
{
	// omp_set_lock(&jq->writeLock);
	{
		lg_addTailData((void*) batch, jq->stolenBatchList);
		++jq->stolenBatchCount;
	}
	// omp_unset_lock(&jq->writeLock);
}

job_batch_t* job_queue_removeStolenBatchWithoutLocking(job_queue_t* jq)
{
	job_batch_t* batch = NULL;

	// omp_set_lock(&jq->writeLock);
	{
		if(job_queue_hasStolenBatch(jq))
		{
			batch = (job_batch_t*) lg_removeTail(jq->stolenBatchList);
		}
	}
	// omp_unset_lock(&jq->writeLock);

	return batch;
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
