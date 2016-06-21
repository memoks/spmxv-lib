/*
 * WorkQueue.h
 *
 *  Created on: Jul 31, 2013
 *      Author: matara
 */

#ifndef JOBQUEUE_H_
#define JOBQUEUE_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <omp.h>

#include "include/config.h"
#include "include/scheduler/job_batch.h"
#include "include/data_structure/sub_mtx.h"
#include "include/data_structure/tree.h"
#include "include/data_structure/spm_storage.h"
#include "include/data_structure/list_generic.h"

/**
 * Each node/die/core has a job queue in which job-batches are stored as
 * doubly linked list. Has a lock to support simultaneous accesses.
 */
struct job_queue
{
	int batchCount;
	omp_lock_t writeLock;
	job_batch_t* batchHead;

	int executableBatchCount;
	job_batch_t* execPtrStart;
	job_batch_t* execPtrEnd;

	int stealTreshold;
	int stolenBatchCount;

	lg_t* stolenBatchList;

	int state;
};

typedef struct job_queue job_queue_t;

extern job_queue_t* job_queue_new(int stealTreshold);
extern void job_queue_init(job_queue_t* jq, int stealTreshold);
extern void job_queue_deleteNonPtr(job_queue_t* queue);
extern void job_queue_delete(job_queue_t* queue);
extern void job_queue_printDetailed(job_queue_t* queue, void* spm);
extern void job_queue_printAttributes(job_queue_t* queue);
extern void job_queue_print(job_queue_t* queue);
extern void job_queue_printCurrentState(
		job_queue_t* queue, void* spm);
extern void job_queue_copy(
		job_queue_t* source, job_queue_t* destination);

/**
 * Resets the execution state of a given job queue to its
 * default value.
 */
extern void job_queue_reset(job_queue_t* queue);

/**
 * Removes all the job batches in a given queue without deleting
 * any job-batch. Keep a reference to job batches in queue to
 * avoid memory leak.
 */
extern void job_queue_emptyQueue(job_queue_t* queue);



// SpMxV data preparation functions
// ------------------------------------------------------------------------

/**
 * Fills a job queue from a given sub-matrix tree node.
 * Uses row wise distribution and column net partition model
 * in which only the tree leafs are considered jobs.
 *
 * @spmCsr: Complete sparse-matrix in CSR format
 * @jq: job queue to which created job batches from tree node
 *  will be appended.
 * @head: tree node (partition) from whose child nodes
 * (sub-partitions) new job batches will be created.
 *
 */
extern void job_queue_fillFromSubMtxTree_rowWiseColNet_CSR(
		spm_cmp_t* spmCsr, job_queue_t* jq, tree_node_t* head);

/**
 * Fills multiple job queues from a given sub-matrix tree
 * node. Uses row wise distribution and column net partition
 * model in which only the tree leafs are considered jobs.
 *
 * @spmCsr: Complete sparse-matrix in CSR format
 * @jq: job queue to which created job batches from tree node
 *  will be appended.
 * @head: tree node (partition) from whose child nodes
 * (sub-partitions) new job batches will be created.
 *
 *
 * Assuming;
 * <ul>
 * <li>sparse matrix has 300 rows, </li>
 * <li>every sub-partition is composed of 20 rows, </li>
 * <li>jqCount is 3 </li>
 *
 * After method call job-queue contents will look like;
 * <ul>
 * <li>jq1: 0-20, 60-80, 120-140, 180-200, 240-260 </li>
 * <li>jq2: 20-40, 80-100, 140-160, 200-220, 260-280 </li>
 * <li>jq3: 40-60, 100-120, 160-180, 220-240, 280-300 </li>
 * </ul>
 *
 * Function is designed to use caches more efficiently by
 * multiple threads in a single block.
 */
extern void job_queue_fillFromSubMtxTreeScatter_rowWiseColNet_CSR(
		spm_cmp_t* spmCsr, job_queue_t** jqPtrArr, int jqCount,
		tree_node_t* head);

/**
 * Fills a given job queue from a given job batch list. Items in
 * job batch list will be removed and appended to job queue.
 * After function call job batch list will be empty.
 *
 * @jq: job queue to which items will be appended.
 * @jobBatchList: job batch list whose items will be removed.
 */
extern void job_queue_fillFromJobBatchList(
		job_queue_t* jq, lg_t* jobBatchList);

/**
 * Fills multiple job queues from a given list using
 * scatter distribution method. Job-batch elements in list
 * will only be removed and added to job queues.
 *
 * @jqPtrArr: Array of job-queue pointers to be filled.
 * @jqCount: Length of job-queue array.
 * @jobBatchList: List initially containing job batches.
 *
 *
 * Assuming;
 * <ul>
 * <li>
 *   list = {(0-20), (20-40), (40-60), (60-80), (80-100)
 *   		(100-120), (120-140), (140-160), (160-180),
 *   		(180, 200)}
 * </li>
 * <li>jqCount = 3 </li>
 *
 * After method call job-queue contents will look like;
 * <ul>
 * <li>list = {}</li>
 * <li>jqPtrArr[0]: (0-20), (60-80), (120-140), (180-200)</li>
 * <li>jqPtrArr[1]: (20-40), (80-100), (140-160)</li>
 * <li>jqPtrArr[2]: (40-60), (100-120), (160-180)</li>
 * </ul>
 *
 * Function is designed to use caches more efficiently by
 * multiple threads in a single block.
 */
extern void job_queue_fillFromJobBatchListScatter(
		job_queue_t** jqPtrArr, int jqCount, lg_t* jobBatchList);

/**
 * <p>Fills given job queue directly from a vector which provides
 * partitions generated by hypergraph tool using column net model
 * for row wise distribution scheme.</p>
 *
 * @spmCsr: full sparse-matrix in CSR format.
 * @jq: job-queue to be appended with newly created job-batches.
 * @rowPartitioning: a vector containing partition info in which
 *  each element is # of rows (sum of which is a partition)
 * @startIndex: starting from this index of partition vector to
 *  generate job batches.
 * @endIndex: until this index (endIndex itself is not included)
 *
 *
 * Suppose we have;
 * <ul>
 *  <li>rowPartitioning = {0, 20, 40, 60, 80, 100}</li>
 *  <li>startIndex = 1</li>
 *  <li>endIndex = 3</li>
 * </ul>
 *
 * Then, (assuming jq is initialy empty) we will have;
 * jq = {(20-40), (40,60)}
 *
 *
 * <p>See input_readPartitionVectorFromFile function description in
 * input_parser.h for more information regarding rowPartitioning
 * vector.</p>
 */
extern void job_queue_fillFromParitioningVector_rowWiseColNet_CSR(
		spm_cmp_t* spmCsr, job_queue_t* jq, vector_int_t* rowPartitioning,
		int startIndex, int endIndex);

/**
 * <p>Fills given job queues directly from a vector which provides
 * partitions generated by hypergraph tool using column net model
 * for row wise distribution scheme. Job queues are filled using
 * scatter fashion.</p>
 *
 * @spmCsr: full sparse-matrix in CSR format.
 * @jobQueuePtrArr: job-queues to be appended with newly
 *  created job-batches.
 * @jobQueueCount: length of jobQueuePtrArr.
 * @rowPartitioning: a vector containing partition info in which
 *  each element is # of rows (sum of which is a partition)
 * @startIndex: starting from this index of partition vector to
 *  generate job batches.
 * @endIndex: until this index (endIndex itself is not included)
 *
 *
 * Suppose we have;
 * <ul>
 *  <li>
 *      rowPartitioning =
 *      {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200}
 *  </li>
 *  <li>jobQueueCount = 3</li>
 *  <li>startIndex = 2</li>
 *  <li>endIndex = 8</li>
 * </ul>
 *
 * Job queue contents will look like this after function call;
 * <ul>
 * <li>jqPtrArr[0]: (40-60), (100-120)</li>
 * <li>jqPtrArr[1]: (60-80), (120-140)</li>
 * <li>jqPtrArr[2]: (80-100), (140-160)/li>
 * </ul>
 *
 * <p>See input_readPartitionVectorFromFile function description in
 * input_parser.h for more information regarding rowPartitioning
 * vector. </p>
 */
extern void job_queue_fillFromParitioningVectorScatter_rowWiseColNet_CSR(
		spm_cmp_t* spmCsr, job_queue_t** jobQueuePtrArr, int jobQueueCount,
		vector_int_t* rowPartitioning, int startIndex, int endIndex);

// TODO this should be deleted completely
extern void job_queue_convertToJDSPartial(
		spm_cmp_t* spmCsrCounterpart, DECIMAL* permutation,
		job_queue_t* jobQueue_inout);

// TODO this should be deleted completely
extern void job_queue_convertFromPartialJDSToHybridJDSCSR(
		job_queue_t* partialJdsJobQueue,
		job_queue_t* hybridJdsCsrJobQueue_inout);

// Job Queue Runtime Management Functions
// Functions below are used in runtime to ensure queue state management
// without causing data-corruption.
// ------------------------------------------------------------------------

/**
 * If there are no local job-batches in a queue than returns 1,
 * return 0 otherwise.
 */
extern int job_queue_isEmpty(job_queue_t* queue);

/**
 * Returns 1 if there are any executable batch residing in execution
 * queue, returns 0 otherwise.
 */
extern int job_queue_hasExecutableBatch(job_queue_t* queue);

/**
 * Returns 1 if there are any stolen batch references residing in
 * steal queue, returns 0 otherwise.
 */
extern int job_queue_hasStolenBatch(job_queue_t* jq);

/**
 * Adds a job-batch at the end of queue. Doesn't use any locks to support
 * simultaneous accesses by itself. Therefore shouldn't be used in
 * multi-threaded algorithms.
 */
extern void job_queue_addBatchWithoutLocking(
		job_queue_t* queue, job_batch_t* batch);

/**
 * Does exactly the same thing job_queue_addBatchWithoutLocking function
 * does. But locks the queue while doing so.
 */
extern void job_queue_addBatch(
		job_queue_t* queue, job_batch_t* batch);

/**
 * Does exactly the same thing job_queue_addBatchWithoutLocking function
 * however adds a whole batch list instead of just one job-batch. And it
 * doesn't add a copy of each batch from batchList. Adds a pointer
 * reference instead.
 *
 * ASSUMPTION: Do not delete the contents of batchList if you want to use
 * batches in job_queue later.
 */
extern void job_queue_addJobBatchListWithoutLocking(
		job_queue_t* queue, lg_t* jobBatchList);

/**
 * Returns a JOB_BATCH reference pointed by the executionPtrStart->next
 * (front of the queue) and increases the executionPtrStart.
 * After calling this function a second time, the returned JOB_BATCH
 * reference will be the one pointed by previously returned JOB_BATCH
 * reference (by next pointer). Doesn't use any locking schemes to
 * support simultaneous accesses.
 */
extern job_batch_t* job_queue_getFirstBatchWithoutLocking(
		job_queue_t* queue);

/**
 * Does exactly the same thing job_queue_getFirstBatchWithoutLocking
 * function does. But locks the queue while doing so.
 */
extern job_batch_t* job_queue_getFirstBatch(
		job_queue_t* queue);

/**
 * Removes and returns the current (top most) job batch in a given
 * job queue. Uses locking schemes.
 */
extern job_batch_t* job_queue_removeFirstBatch(
		job_queue_t* queue);

/**
 * Returns a JOB_BATCH reference pointed by the executionPtrEnd
 * (back of the queue) and decreases the executionPtrEnd. Uses a
 * lock routine to support simultaneous accesses.
 */
extern job_batch_t* job_queue_getLastBatch(
		job_queue_t* queue);

/**
 * Does exactly the same thing that job_queue_getLastBatch function
 * does except returns half of the queue instead of just one work
 * item.
 */
extern void job_queue_getLastBatchHalf(
		job_queue_t* queue, job_batch_t** start_out, int* count_out);

/**
 * Adds a JOB_BATCH reference" to steal queue in given job queue.
 */
extern void job_queue_addStolenBatchWithoutLocking(
		job_queue_t* queue, job_batch_t* batch);

/**
 * Removes a JOB_BATCH reference from the steal queue in jq. Reference
 * should not be freed afterwards.
 */
extern job_batch_t* job_queue_removeStolenBatchWithoutLocking(
		job_queue_t* queue);


#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* JOBQUEUE_H_ */
