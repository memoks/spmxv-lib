
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#ifdef __ICC
#include <mkl.h>
#endif

#include "include/task_decomposition/partitioning.h"

#include "include/scheduler/job_batch.h"

// Job Batch Helper functions
// --------------------------------------------------------------------------------------------------------------------------

static void __job_batch_extractJobBatchListFromSubMtxTree_rowWiseColNet_CSR(
		spm_cmp_t* spmCsr, lg_t* jobBatchList, tree_node_t* head);

// --------------------------------------------------------------------------------------------------------------------------

void job_batch_execute_hybrid_jds_csr(
		job_batch_t* jb, vector_real_t* x, vector_real_t* y_inout)
{
	sub_mtx_dim_t* subMtx = &jb->data.hybrid_jds_csr.subMtx;

	job_batch_printDetailed(jb, NULL);
	spm_cmp_print(jb->data.hybrid_jds_csr.csr);

	// execute partial CSR
	DECIMAL rowShift = 0;
	if(jb->data.hybrid_jds_csr.csr != NULL)
	{
		spm_cmp_t* spmCsr = jb->data.hybrid_jds_csr.csr;
		spmxv_csr_add(
				0,
				spmCsr->rowCount,
				spmCsr->ptr,
				spmCsr->ind,
				spmCsr->values,
				x->data,
				&y_inout->data[subMtx->start.i]);
		rowShift = spmCsr->rowCount;
	}

	DECIMAL i;
	for(i = subMtx->start.i; i < subMtx->start.i + subMtx->length.i; ++i)
		printf("%d. %f\n", i, y_inout->data[i]);
	printf("\n");

	spm_jds_print(jb->data.hybrid_jds_csr.jds);
	if(jb->data.hybrid_jds_csr.jds != NULL)
	{
		spm_jds_t* spmJds = jb->data.hybrid_jds_csr.jds;
		spmxv_jds_partial_optimized(
				subMtx->start.i,
				0,
				spmJds->idiagLength,
				spmJds->idiag,
				spmJds->jdiag,
				spmJds->dj,
				x->data,
				y_inout->data);
	}

	for(i = subMtx->start.i; i < subMtx->start.i + subMtx->length.i; ++i)
		printf("%d. %f\n", i, y_inout->data[i]);
	printf("\n");
}


// Job Batch functions
// --------------------------------------------------------------------------------------------------------------------------

job_batch_t* job_batch_new(void)
{
	job_batch_t* jb = (job_batch_t*) malloc(sizeof(job_batch_t));
	job_batch_initDefault(jb);

	return jb;
}

job_batch_t* job_batch_newSubMtx(sub_mtx_dim_t* subMtx)
{
	job_batch_t* jb = (job_batch_t*) malloc(sizeof(job_batch_t));
	job_batch_initSubMtx(jb, subMtx);

	return jb;
}

job_batch_t* job_batch_newPartialJds(sub_mtx_dim_t* subMtx, spm_jds_t* spmPartialJds)
{
	job_batch_t* jb = (job_batch_t*) malloc(sizeof(job_batch_t));
	job_batch_initPartialJds(jb, subMtx, spmPartialJds);

	return jb;
}

job_batch_t* job_batch_newHybridJdsCsr(
		sub_mtx_dim_t* subMtx, spm_jds_t* spmJds, spm_cmp_t* spmCsr)
{
	job_batch_t* jb = (job_batch_t*) malloc(sizeof(job_batch_t));
	job_batch_initHybridJdsCsr(jb, subMtx, spmJds, spmCsr);

	return jb;
}

job_batch_t* job_batch_newTrueHybrid(
		sub_mtx_dim_t* subMtx,
		DECIMAL csrNNZ, DECIMAL csrRowCount, DECIMAL csrColCount,
		DECIMAL jdsNNZ, DECIMAL jdsRowCount, DECIMAL jdsColCount, DECIMAL jdsIdiagLength)
{
	job_batch_t* jb = (job_batch_t*) malloc(sizeof(job_batch_t));
	job_batch_initTrueHybrid(
			jb, subMtx, csrNNZ, csrRowCount, csrColCount, jdsNNZ, jdsRowCount, jdsColCount, jdsIdiagLength);

	return jb;
}

void job_batch_initDefault(job_batch_t* jb)
{
	INIT_LIST_HEAD(&jb->head);
	jb->type = JOB_NONE;
}

void job_batch_initSubMtx(job_batch_t* jb, sub_mtx_dim_t* subMtx)
{
	INIT_LIST_HEAD(&jb->head);
	jb->type = JOB_SUB_MTX_CSR;

	sub_mtx_initDefault(&jb->data.subMtx);
	sub_mtx_copy(subMtx, &jb->data.subMtx);
}

void job_batch_initPartialJds(
		job_batch_t* jb, sub_mtx_dim_t* subMtx, spm_jds_t* spmPartialJds)
{
	INIT_LIST_HEAD(&jb->head);
	jb->type = JOB_PARTIAL_JDS;

	sub_mtx_initDefault(&jb->data.partialJds.subMtx);
	sub_mtx_copy(subMtx, &jb->data.partialJds.subMtx);
	jb->data.partialJds.spmJds = spmPartialJds;
}

void job_batch_initHybridJdsCsrDefault(job_batch_t* jb)
{
	INIT_LIST_HEAD(&jb->head);
	jb->type = JOB_HYBRID_JDS_CSR;
	jb->data.hybrid_jds_csr.jds = NULL;
	jb->data.hybrid_jds_csr.csr = NULL;
	sub_mtx_initDefault(&jb->data.hybrid_jds_csr.subMtx);
}

void job_batch_initHybridJdsCsr(
		job_batch_t* jb, sub_mtx_dim_t* subMtx, spm_jds_t* spmJds, spm_cmp_t* spmCsr)
{
	job_batch_initHybridJdsCsrDefault(jb);

	sub_mtx_copy(subMtx, &jb->data.hybrid_jds_csr.subMtx);

	if(spmJds != NULL)
		jb->data.hybrid_jds_csr.jds = spmJds;

	if(spmCsr != NULL)
		jb->data.hybrid_jds_csr.csr = spmCsr;
}

void job_batch_initTrueHybrid(
		job_batch_t* jb, sub_mtx_dim_t* subMtx,
		DECIMAL csrNNZ, DECIMAL csrRowCount, DECIMAL csrColCount,
		DECIMAL jdsNNZ, DECIMAL jdsRowCount,
		DECIMAL jdsColCount, DECIMAL jdsIdiagLength)
{
	job_batch_initTrueHybridDefault(jb);
	sub_mtx_copy(subMtx, &jb->data.hybrid_jds_csr.subMtx);

	DECIMAL nnz = csrNNZ + jdsNNZ;
	spm_jds_t* spmJds = NULL;
	spm_cmp_t* spmCsr = NULL;

	DECIMAL totalPtrLength = 0;
	if(csrNNZ > 0)
		totalPtrLength += csrRowCount + 1;
	if(jdsNNZ > 0)
		totalPtrLength += jdsIdiagLength + 1;

	REAL* values = NULL;
	DECIMAL* colInd = NULL;
	DECIMAL* ptr = NULL;

#ifdef __ICC
	if(csrNNZ > 0)
		spmCsr = (spm_cmp_t*) mkl_malloc(sizeof(spm_cmp_t), ALIGNMENT);
	if(jdsNNZ > 0)
		spmJds = (spm_jds_t*) mkl_malloc(sizeof(spm_jds_t), ALIGNMENT);

	values = (REAL*) mkl_malloc(nnz * sizeof(REAL), ALIGNMENT);
	colInd = (DECIMAL*) mkl_malloc(nnz * sizeof(DECIMAL), ALIGNMENT);
	ptr = (DECIMAL*) mkl_malloc(totalPtrLength * sizeof(DECIMAL), ALIGNMENT);
#else
	if(csrNNZ > 0)
		posix_memalign(&spmCsr, ALIGNMENT, sizeof(spm_cmp_t));
	if(jdsNNZ > 0)
		posix_memalign(&spmJds, ALIGNMENT, sizeof(spm_jds_t));

	posix_memalign(&values, ALIGNMENT, sizeof(REAL) * nnz);
	posix_memalign(&colInd, ALIGNMENT, sizeof(DECIMAL) * nnz);
	posix_memalign(&ptr, ALIGNMENT, sizeof(DECIMAL) * totalPtrLength);
#endif

	DECIMAL carryOverPtrIndex = 0;
	if(csrNNZ > 0)
	{
		spmCsr->colCount = csrColCount;
		spmCsr->rowCount = csrRowCount;
		spmCsr->nnz = csrNNZ;
		spmCsr->values = values;
		spmCsr->ind = colInd;
		spmCsr->ptr = ptr;
		jb->data.hybrid_jds_csr.csr = spmCsr;
		carryOverPtrIndex = csrRowCount + 1;
	}

	if(jdsNNZ > 0)
	{
		spmJds->permutation = NULL;
		spmJds->colCount = jdsColCount;
		spmJds->rowCount = jdsRowCount;
		spmJds->nnz = jdsNNZ;
		spmJds->idiagLength = jdsIdiagLength;
		spmJds->csrCounterpart = NULL;
		spmJds->dj = &values[csrNNZ];
		spmJds->jdiag = &colInd[csrNNZ];
		spmJds->idiag = &ptr[carryOverPtrIndex];
		jb->data.hybrid_jds_csr.jds = spmJds;
	}

	// Initialize
	DECIMAL i;
	for(i = 0; i < nnz; ++i)
	{
		values[i] = 0.0;
		colInd[i] = 0;
	}
	for(i = 0; i < totalPtrLength; ++i)
		ptr[i] = 0;
}

void job_batch_initTrueHybridDefault(job_batch_t* jb)
{
	INIT_LIST_HEAD(&jb->head);
	jb->type = JOB_TRUE_HYBRID;
	jb->data.hybrid_jds_csr.jds = NULL;
	jb->data.hybrid_jds_csr.csr = NULL;
	sub_mtx_initDefault(&jb->data.hybrid_jds_csr.subMtx);
}

void job_batch_copy(job_batch_t* source, job_batch_t* destination)
{
	destination->type = source->type;
	if(source->type == JOB_SUB_MTX_CSR)
	{
		sub_mtx_copy(&source->data.subMtx, &destination->data.subMtx);
	}
	else if(source->type == JOB_PARTIAL_JDS)
	{
		// copy sub-matrix
		sub_mtx_copy(&source->data.partialJds.subMtx, &destination->data.partialJds.subMtx);
		spm_jds_t* sourceJds = source->data.partialJds.spmJds;

		// copy partial JDS
		destination->data.partialJds.spmJds = spm_jds_new(
				sourceJds->nnz, sourceJds->rowCount, sourceJds->colCount, sourceJds->idiagLength);
		spm_jds_copy(sourceJds, destination->data.partialJds.spmJds);
	}
	else if(source->type == JOB_HYBRID_JDS_CSR)
	{
		job_batch_initHybridJdsCsrDefault(destination);

		// copy sub-matrix
		sub_mtx_copy(&source->data.hybrid_jds_csr.subMtx, &destination->data.hybrid_jds_csr.subMtx);

		// copy partial JDS
		if(source->data.hybrid_jds_csr.jds != NULL)
		{
			spm_jds_t* sourceJds = source->data.hybrid_jds_csr.jds;
			destination->data.hybrid_jds_csr.jds = spm_jds_new(
					sourceJds->nnz, sourceJds->rowCount, sourceJds->colCount,
					sourceJds->idiagLength);

			spm_jds_copy(sourceJds, destination->data.hybrid_jds_csr.jds);
		}

		// copy partial CSR
		if(source->data.hybrid_jds_csr.csr != NULL)
		{
			spm_cmp_t* sourceCsr = source->data.hybrid_jds_csr.csr;
			destination->data.hybrid_jds_csr.csr = spm_cmp_new(
					sourceCsr->rowCount, sourceCsr->colCount, sourceCsr->nnz);

			spm_cmp_copy(sourceCsr, destination->data.hybrid_jds_csr.csr);
		}
	}
	else if(source->type == JOB_TRUE_HYBRID)
	{
		DECIMAL i;

		REAL* sourceValues = NULL;
		DECIMAL* sourceColInd = NULL;
		DECIMAL* sourcePtr = NULL;

		REAL** destinationValues = NULL;
		DECIMAL** destinationColInd = NULL;
		DECIMAL** destinationPtr = NULL;

		DECIMAL csrNNZ = 0;
		DECIMAL csrRowCount = 0;
		DECIMAL csrColCount = 0;

		DECIMAL jdsNNZ = 0;
		DECIMAL jdsRowCount = 0;
		DECIMAL jdsColCount = 0;
		DECIMAL jdsIdiagLength = 0;

		if(source->data.hybrid_jds_csr.csr != NULL)
		{
			csrNNZ = source->data.hybrid_jds_csr.csr->nnz;
			csrRowCount = source->data.hybrid_jds_csr.csr->rowCount;
			csrColCount = source->data.hybrid_jds_csr.csr->colCount;

			sourceValues = source->data.hybrid_jds_csr.csr->values;
			sourceColInd = source->data.hybrid_jds_csr.csr->ind;
			sourcePtr = source->data.hybrid_jds_csr.csr->ptr;
		}
		else
		{
			sourceValues = source->data.hybrid_jds_csr.jds->dj;
			sourceColInd = source->data.hybrid_jds_csr.jds->jdiag;
			sourcePtr = source->data.hybrid_jds_csr.jds->idiag;
		}

		if(source->data.hybrid_jds_csr.jds != NULL)
		{
			jdsNNZ = source->data.hybrid_jds_csr.jds->nnz;
			jdsRowCount = source->data.hybrid_jds_csr.jds->rowCount;
			jdsColCount = source->data.hybrid_jds_csr.jds->colCount;
			jdsIdiagLength = source->data.hybrid_jds_csr.jds->idiagLength;
		}

		job_batch_initTrueHybrid(destination, &source->data.hybrid_jds_csr.subMtx,
				csrNNZ, csrRowCount, csrColCount,
				jdsNNZ, jdsRowCount, jdsColCount, jdsIdiagLength);

		if(source->data.hybrid_jds_csr.csr != NULL)
		{
			destinationValues = &destination->data.hybrid_jds_csr.csr->values;
			destinationColInd = &destination->data.hybrid_jds_csr.csr->ind;
			destinationPtr = &destination->data.hybrid_jds_csr.csr->ptr;
		}
		else
		{
			destinationValues = &destination->data.hybrid_jds_csr.jds->dj;
			destinationColInd = &destination->data.hybrid_jds_csr.jds->jdiag;
			destinationPtr = &destination->data.hybrid_jds_csr.jds->idiag;
		}

		// copy
		DECIMAL totalPtr = 0;
		DECIMAL nnzIndex = 0;
		DECIMAL ptrIndex = 0;

		// copy csr
		if(source->data.hybrid_jds_csr.csr != NULL)
		{
			while(nnzIndex < csrNNZ)
			{
				(*destinationValues)[nnzIndex] = sourceValues[nnzIndex];
				(*destinationColInd)[nnzIndex] = sourceColInd[nnzIndex];

				++nnzIndex;
			}

			totalPtr = csrRowCount + 1;
			while(ptrIndex < totalPtr)
			{
				(*destinationPtr)[ptrIndex] = sourcePtr[ptrIndex];
				++ptrIndex;
			}
		}

		// copy jds
		if(source->data.hybrid_jds_csr.jds != NULL)
		{
			while(nnzIndex < (csrNNZ + jdsNNZ))
			{
				(*destinationValues)[nnzIndex] = sourceValues[nnzIndex];
				(*destinationColInd)[nnzIndex] = sourceColInd[nnzIndex];

				++nnzIndex;
			}

			totalPtr += jdsIdiagLength + 1;
			while(ptrIndex < totalPtr)
			{
				(*destinationPtr)[ptrIndex] = sourcePtr[ptrIndex];
				++ptrIndex;
			}
		}
	}
	else // if(source->type == JOB_NONE)
	{
		PRINTF("job_batch_copy: " \
				"Something is wrong with the prepared data. Exiting program...\n");
		exit(EXIT_FAILURE);
	}
}

job_batch_t* job_batch_copyToPtr(job_batch_t* source)
{
	job_batch_t* copy = job_batch_new();
	job_batch_copy(source, copy);
	return copy;
}

void job_batch_deleteNonPtr(job_batch_t* jb)
{
	if(jb == NULL)
		return;

	if(jb->type == JOB_SUB_MTX_CSR)
	{
		// there are no pointers for this type so do nothing
	}
	else if(jb->type == JOB_PARTIAL_JDS)
	{
		spm_jds_delete(jb->data.partialJds.spmJds);
	}
	else if(jb->type == JOB_HYBRID_JDS_CSR)
	{
		if(jb->data.hybrid_jds_csr.jds != NULL)
			spm_jds_delete(jb->data.hybrid_jds_csr.jds);

		if(jb->data.hybrid_jds_csr.csr != NULL)
			spm_cmp_delete(jb->data.hybrid_jds_csr.csr);
	}
	else if(jb->type == JOB_TRUE_HYBRID)
	{
		if(jb->data.hybrid_jds_csr.csr != NULL)
		{
			spm_cmp_delete(jb->data.hybrid_jds_csr.csr);
			if(jb->data.hybrid_jds_csr.jds != NULL)
			{
				free(jb->data.hybrid_jds_csr.jds->permutation);
				free(jb->data.hybrid_jds_csr.jds);
			}
		}
		else
		{
			spm_jds_delete(jb->data.hybrid_jds_csr.jds);
		}
	}
	else // if(jb->type == JOB_NONE)
	{
		// do nothing
	}
}

sub_mtx_dim_t* job_batch_getSubMtx(job_batch_t* jb)
{
	sub_mtx_dim_t* subMtx = NULL;

	if(jb->type == JOB_SUB_MTX_CSR)
	{
		subMtx = &jb->data.subMtx;
	}
	else if(jb->type == JOB_PARTIAL_JDS)
	{
		subMtx = &jb->data.partialJds.subMtx;
	}
	else if(jb->type == JOB_HYBRID_JDS_CSR)
	{
		subMtx = &jb->data.hybrid_jds_csr.subMtx;
	}
	else // if(jb->type == JOB_NONE)
	{
		PRINTF("job_batch_getSubMtx: "
				"Something is wrong with the prepared data. Exiting program...\n");
		exit(EXIT_FAILURE);
	}

	return subMtx;
}

void job_batch_extractStatisticsPartialJds(
		job_batch_t* jobBatch,
		REAL* maxr_out, REAL* avgr_out, REAL* minr_out,
		REAL* maxc_out, REAL* avgc_out, REAL* minc_out)
{
	sub_mtx_extractRowStats_JDS(jobBatch->data.partialJds.spmJds, maxr_out, avgr_out, minr_out);
	sub_mtx_extractColumnStats_JDS(
			jobBatch->data.partialJds.spmJds, maxc_out, avgc_out, minc_out);
}

void job_batch_extractStatisticsArrPartialJds(
		job_batch_t* jobBatchArr, int length,
		REAL* maxr_out, REAL* avgr_out, REAL* minr_out,
		REAL* maxc_out, REAL* avgc_out, REAL* minc_out)
{
	REAL maxr = 0.0, avgr = 0.0, minr = 0.0, maxc = 0.0, avgc = 0.0, minc = 0.0;
	DECIMAL i;
	for(i = 0; i < length; ++i)
	{
		REAL tempMaxr = 0.0, tempAvgr = 0.0, tempMinr = 0.0;
		REAL tempMaxc = 0.0, tempAvgc = 0.0, tempMinc = 0.0;
		job_batch_t* currBatch = &jobBatchArr[i];

		job_batch_extractStatisticsPartialJds(
				currBatch,
				&tempMaxr, &tempAvgr, &tempMinr,
				&tempMaxc, &tempAvgc, &tempMinc);

		maxr += tempMaxr;
		avgr += tempAvgr;
		minr += tempMinr;
		maxc += tempMaxc;
		avgc += tempAvgc;
		minc += tempMinc;
	}

	// return values
	*maxr_out = maxr / ((REAL) length);
	*avgr_out = avgr / ((REAL) length);
	*minr_out = minr / ((REAL) length);
	*maxc_out = maxc / ((REAL) length);
	*avgc_out = avgc / ((REAL) length);
	*minc_out = minc / ((REAL) length);
}

void job_batch_extractStatisticsArrMultiplePartialJds(
		job_batch_t** jobBatchArrPerBlock, int* jobBatchCountPerBlock, int numBlocks,
		REAL* maxr_out, REAL* avgr_out, REAL* minr_out,
		REAL* maxc_out, REAL* avgc_out, REAL* minc_out)
{
	REAL maxr = 0.0, avgr = 0.0, minr = 0.0, maxc = 0.0, avgc = 0.0, minc = 0.0;

	#pragma omp parallel num_threads(numBlocks) reduction(+:maxr, avgr, minr, maxc, avgc, minc)
	{
		int blockId = omp_get_thread_num();

		REAL tempMaxr = 0.0, tempAvgr = 0.0, tempMinr = 0.0;
		REAL tempMaxc = 0.0, tempAvgc = 0.0, tempMinc = 0.0;
		job_batch_t* currBatchArr = jobBatchArrPerBlock[blockId];
		int currBatchArrLength = jobBatchCountPerBlock[blockId];

		job_batch_extractStatisticsArrPartialJds(
				currBatchArr, currBatchArrLength,
				&tempMaxr, &tempAvgr, &tempMinr,
				&tempMaxc, &tempAvgc, &tempMinc);


		maxr += tempMaxr;
		avgr += tempAvgr;
		minr += tempMinr;
		maxc += tempMaxc;
		avgc += tempAvgc;
		minc += tempMinc;
	}

	// return values
	*maxr_out = maxr / ((REAL) numBlocks);
	*avgr_out = avgr / ((REAL) numBlocks);
	*minr_out = minr / ((REAL) numBlocks);
	*maxc_out = maxc / ((REAL) numBlocks);
	*avgc_out = avgc / ((REAL) numBlocks);
	*minc_out = minc / ((REAL) numBlocks);
}


void job_batch_delete(job_batch_t* jobBatch)
{
	// Check if node is already removed from the list. If not remove it.
	if(!(jobBatch->head.next == NULL && jobBatch->head.prev == NULL))
		list_del(&jobBatch->head);

	job_batch_deleteNonPtr(jobBatch);

	free(jobBatch);
}

void job_batch_deleteAll(job_batch_t* head)
{
	// Only when there are more than 1 batch (base) in the list start deleting
	if(head != job_batch_getNext(head))
	{
		list_head_t* curr;
		list_head_t* temp;
		// start from first element, not base element
		list_for_each_safe(curr, temp, &head->head)
		{
			job_batch_t* temp = job_batch_getBatch(curr);
			job_batch_delete(temp);
		}
	}
}

void job_batch_deleteAllShallow(job_batch_t* head)
{
	// Only when there are more than 1 batch (base) in the list start deleting
	if(head != job_batch_getNext(head))
	{
		list_head_t* curr;
		list_head_t* temp;
		// start from first element, not base element
		list_for_each_safe(curr, temp, &head->head)
		{
			job_batch_t* temp = job_batch_getBatch(curr);

			// Instead of deleting, only remove from queue
			if(!(temp->head.next == NULL && temp->head.prev == NULL))
				list_del(&temp->head);
		}
	}
}

void job_batch_deleteDenseArr(job_batch_t* jobBatchArr, int length)
{
	if(jobBatchArr == NULL)
		return;

	int i;
	for(i = 0; i < length; ++i)
		job_batch_deleteNonPtr(&jobBatchArr[i]);

	free(jobBatchArr);
}

void job_batch_deleteDenseArrMultiple(
		job_batch_t** jobBatchArrPerBlock, int* jobBatchCountPerBlock, int numBlocks)
{
	if(jobBatchArrPerBlock == NULL)
		return;

	int i;
	for(i = 0; i < numBlocks; ++i)
	{
		job_batch_deleteDenseArr(jobBatchArrPerBlock[i], jobBatchCountPerBlock[i]);
	}

	free(jobBatchArrPerBlock);
}

job_batch_t* job_batch_getBatch(list_head_t* listHead)
{
	return list_entry(listHead, job_batch_t, head);
}

void job_batch_addBatchFront(job_batch_t* new, job_batch_t* head)
{
	list_add(&new->head, &head->head);
}

void job_batch_addBatchTail(job_batch_t* new, job_batch_t* head)
{
	list_add_tail(&new->head, &head->head);
}

job_batch_t* job_batch_getNext(job_batch_t* batch)
{
	return job_batch_getBatch(batch->head.next);
}

job_batch_t* job_batch_getPrev(job_batch_t* batch)
{
	return job_batch_getBatch(batch->head.prev);
}

DECIMAL job_batch_getRowCount(job_batch_t* jobBatch)
{
	// currently we have only one job-batch type
	// if(jobBatch->type == JOB_SUB_MTX_CSR)
	return jobBatch->data.subMtx.length.i;
}

void job_batch_toString(job_batch_t* jobBatch, char* buff_inout)
{
	if(jobBatch->type == JOB_SUB_MTX_CSR)
	{
		strcat(buff_inout, "JOB_SUB_MTX_CSR: ");
		sub_mtx_toString(&jobBatch->data.subMtx, buff_inout);
	}
	else if(jobBatch->type == JOB_PARTIAL_JDS)
	{
		strcat(buff_inout, "JOB_PARTIAL_JDS: ");
		sub_mtx_toString(&jobBatch->data.partialJds.subMtx, buff_inout);
	}
	else if(jobBatch->type == JOB_HYBRID_JDS_CSR || jobBatch->type == JOB_TRUE_HYBRID)
	{
		strcat(buff_inout, "JOB_HYBRID_JDS_CSR: ");
		sub_mtx_toString(&jobBatch->data.hybrid_jds_csr.subMtx, buff_inout);
	}
	else // if(jobBatch->type == JOB_NONE)
	{
		// DO NOTHING
	}
}

void job_batch_print(job_batch_t* jobBatch)
{
	char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
	job_batch_toString(jobBatch, temp);
	PRINTF("%s\n", temp);
}

void job_batch_printArr(char* message, job_batch_t* jobBatchArr, int length)
{
	PRINTF("%s\n", message);

	int i;
	for(i = 0; i < length; ++i)
	{
		PRINTF("%d) ", i + 1);
		job_batch_print(&jobBatchArr[i]);
	}
}

void job_batch_printArrMultiple(
		char* message, job_batch_t** jobBatchArrPerBlock, int* jobBatchCountPerBlock, int numBlocks)
{
	PRINTF("%s\n", message);
	int i;
	for(i = 0; i < numBlocks; ++i)
	{
		char mess[128];
		sprintf(mess, "block-%d => %d job batches", i, jobBatchCountPerBlock[i]);
		job_batch_printArr(mess, jobBatchArrPerBlock[i], jobBatchCountPerBlock[i]);
	}
}

// --------------------------------------------------------------------------------------------------------------------------

DECIMAL job_batch_getNNZ(job_batch_t* jobBatch, void* spm)
{
	DECIMAL nnz = 0;

	if(jobBatch->type == JOB_SUB_MTX_CSR)
	{
		nnz = sub_mtx_getNNZ_CSR(&jobBatch->data.subMtx, (spm_cmp_t*) spm);
	}
	else if(jobBatch->type == JOB_PARTIAL_JDS)
	{
		nnz = jobBatch->data.partialJds.spmJds->nnz;
	}
	else if(jobBatch->type == JOB_HYBRID_JDS_CSR || jobBatch->type == JOB_TRUE_HYBRID)
	{
		if(jobBatch->data.hybrid_jds_csr.jds != NULL)
			nnz += jobBatch->data.hybrid_jds_csr.jds->nnz;
		if(jobBatch->data.hybrid_jds_csr.csr != NULL)
			nnz += jobBatch->data.hybrid_jds_csr.csr->nnz;
	}
	else // if(jobBatch->type == JOB_NONE)
	{
		// DO NOTHING
	}

	return nnz;
}

void job_batch_printDetailed(job_batch_t* batch, void* spm)
{
	if(batch == NULL)
		return;

	char temp[DEFAULT_STR_BUFF_SIZE] = "\0";
	job_batch_toStringDetailed(batch, spm, temp);
	PRINTF("%s\n", temp);
}

void job_batch_printArrDetailed(char* message, job_batch_t* jobBatchArr, int length, void* spm)
{
	if(message != NULL)
		PRINTF("%s\n", message);

	int i;
	for(i = 0; i < length; ++i)
		job_batch_printDetailed(&jobBatchArr[i], spm);
}

void job_batch_printArrMultipleDetailed(
		char* message, job_batch_t** jobBatchArrPerBlock, int* jobBatchCountPerBlock, int numBlocks, void* spm)
{
	PRINTF("%s\n", message);
	int i;
	int j;
	for(i = 0; i < numBlocks; ++i)
	{
		PRINTF("batch-count: %d\n", jobBatchCountPerBlock[i]);
		job_batch_printArrDetailed(NULL, jobBatchArrPerBlock[i], jobBatchCountPerBlock[i], spm);
		PRINTF("\n");
	}
}

void job_batch_toStringDetailed(job_batch_t* jobBatch, void* spm, char* buff_inout)
{
	if(jobBatch->type == JOB_SUB_MTX_CSR)
	{
		strcat(buff_inout, "JOB_SUB_MTX_CSR: ");
		sub_mtx_toString(&jobBatch->data.subMtx, buff_inout);
		char temp[DEFAULT_STR_BUFF_SIZE];
		sprintf(temp, " NNZ: %d", sub_mtx_getNNZ_CSR(&jobBatch->data.subMtx, (spm_cmp_t*) spm));
		strcat(buff_inout, temp);
	}
	else if(jobBatch->type == JOB_PARTIAL_JDS)
	{
		strcat(buff_inout, "JOB_PARTIAL_JDS: ");
		sub_mtx_toString(&jobBatch->data.partialJds.subMtx, buff_inout);
		char temp[DEFAULT_STR_BUFF_SIZE];
		sprintf(temp, " NNZ: %d\n", jobBatch->data.partialJds.spmJds->nnz);
		strcat(buff_inout, temp);

		spm_jds_toStringAttributes(jobBatch->data.partialJds.spmJds, buff_inout);
		// strcat(buff_inout, "\n");
		// spm_jds_plotToBuffer_JDS(jobBatch->data.partialJds.spmJds, buff_inout);
	}
	else if(jobBatch->type == JOB_HYBRID_JDS_CSR)
	{
		strcat(buff_inout, "JOB_HYBRID_JDS_CSR: ");
		sub_mtx_toString(&jobBatch->data.hybrid_jds_csr.subMtx, buff_inout);

		strcat(buff_inout, "\n");
		spm_jds_toStringAttributes(jobBatch->data.hybrid_jds_csr.jds, buff_inout);
		// strcat(buff_inout, "\n");
		// spm_jds_plotToBuffer_JDS(jobBatch->data.hybrid_jds_csr.jds, buff_inout);

		strcat(buff_inout, "\n");
		spm_cmp_toStringAttributes(jobBatch->data.hybrid_jds_csr.csr, buff_inout);
		// strcat(buff_inout, "\n");
		// spm_cmp_plotToBufferShovedLeft(jobBatch->data.hybrid_jds_csr.csr, buff_inout);
	}
	else if(jobBatch->type == JOB_TRUE_HYBRID)
	{
		strcat(buff_inout, "JOB_HYBRID_JDS_CSR: ");
		sub_mtx_toString(&jobBatch->data.hybrid_jds_csr.subMtx, buff_inout);

		strcat(buff_inout, "\n");
		spm_jds_toStringAttributes(jobBatch->data.hybrid_jds_csr.jds, buff_inout);
		// strcat(buff_inout, "\n");
		// spm_jds_plotToBuffer_JDS(jobBatch->data.hybrid_jds_csr.jds, buff_inout);

		strcat(buff_inout, "\n");
		spm_cmp_toStringAttributes(jobBatch->data.hybrid_jds_csr.csr, buff_inout);
		// strcat(buff_inout, "\n");
		// spm_cmp_plotToBuffer_JDS(jobBatch->data.hybrid_jds_csr.csr, buff_inout);
	}
	else // if(jobBatch->type == JOB_NONE)
	{
		// DO NOTHING
	}
}

void job_batch_extractJobBatchListFromSubMtxList_rowWiseColNet_CSR(
		lg_t* subMtxList, lg_t** jobBatchList_out)
{
	if(subMtxList == NULL)
		return;

	lg_t* jobBatchList = lg_new();

	list_head_t* currListHead = NULL;
	list_for_each(currListHead, &subMtxList->headNode->listHead)
	{
		sub_mtx_dim_t* currSubMtx = (sub_mtx_dim_t*) lng_getListGeneric(currListHead)->dataPtr;
		job_batch_t* new = job_batch_new();
		job_batch_initSubMtx(new, currSubMtx);

		lg_addTailData(new, jobBatchList);
	}

	*jobBatchList_out = jobBatchList;
}

void job_batch_extractJobBatchListFromSubMtxTree_rowWiseColNet_CSR(
		spm_cmp_t* spmCsr, tree_node_t* node, lg_t** jobBatchList_out)
{
	lg_t* jobBatchList = lg_new();

	__job_batch_extractJobBatchListFromSubMtxTree_rowWiseColNet_CSR(spmCsr, jobBatchList, node);

	*jobBatchList_out = jobBatchList;
}

void job_batch_extractJobBatchListFromPartitionVector_rowWiseColNet_CSR(
		spm_cmp_t* spmCsr, partitioning_vector_t* rowPartitioningVector,
		DECIMAL startIndex, DECIMAL endIndex, lg_t** jobBatchList_out)
{
	lg_t* jobBatchList = lg_new();

	int i;
	for(i = startIndex; i < endIndex; ++i)
	{
		sub_mtx_dim_t subMtx;
		sub_mtx_init(&subMtx,
				rowPartitioningVector->data[i],
				0,
// TODO come back here
				rowPartitioningVector->data[i + 1] - rowPartitioningVector->data[i],
//				rowCount,
				spmCsr->colCount);

		job_batch_t* batch = job_batch_newSubMtx(&subMtx);
		lg_addTailData((void*) batch, jobBatchList);
	}

	// return values
	*jobBatchList_out = jobBatchList;
}

void job_batch_scatterJobBatchList(
		lg_t* source, int scatterSize, lg_t*** jobBatchLists_out)
{
	if(source == NULL)
		return;

	lg_t** jobBatchLists = (lg_t**) malloc(sizeof(lg_t*) * scatterSize);
	int i;
	for(i = 0; i < scatterSize; ++i)
	{
		jobBatchLists[i] = lg_new();
	}

	int currListIndex = 0;
	list_head_t* currListHead = NULL;
	list_head_t* temp = NULL;
	list_for_each_safe(currListHead, temp, &source->headNode->listHead)
	{
		lng_t* currListNode = lng_getListGeneric(currListHead);
		void* data = lng_remove(currListNode);

		lg_addTailData(data, jobBatchLists[currListIndex]);

		++currListIndex;
		currListIndex = currListIndex % scatterSize;
	}

	// return values
	*jobBatchLists_out = jobBatchLists;
}

// TODO must be tested lg_remove method changed
void job_batch_mergeJobBatchList_CSR(lg_t* jobBatchList)
{
	if(jobBatchList == NULL)
		return;

	lng_t* currListItem = lng_next(jobBatchList->headNode);
	lng_t* nextListItem = NULL;
	job_batch_t* currJobBatch = (job_batch_t*) currListItem->dataPtr;
	job_batch_t* nextJobBatch = NULL;

	// merging is done by one pass window from left to right
	for(nextListItem = lng_next(currListItem);
		&currListItem->listHead != &jobBatchList->headNode->listHead &&
		&nextListItem->listHead != &jobBatchList->headNode->listHead;
		nextListItem = lng_next(currListItem), nextJobBatch = (job_batch_t*) currListItem->dataPtr)
	{
		// If two job batches are of same time, only then we consider merging
		if(currJobBatch->type == nextJobBatch->type)
		{
			sub_mtx_dim_t* curr = job_batch_getSubMtx(currJobBatch);
			sub_mtx_dim_t* next = job_batch_getSubMtx(nextJobBatch);

			// if left.end == right.start then there is an overlap
			// In other words, curr and next are mergable
			if((curr->start.i + curr->length.i) == curr->start.i)
			{
				sub_mtx_merge(curr, next);

				// delete next since it is already a part of curr-sub-matrix
				lg_remove(jobBatchList, nextListItem);
				job_batch_delete(nextJobBatch);
			}
			// Current job batch is not mergable with any other. So get the next one.
			else
			{
				currListItem = nextListItem;
			}
		}
	}
}

void lg_job_batch_extractSubMtx(lg_t* jobBatchList, lg_t* subMtxList_inout)
{
	if(jobBatchList == NULL)
		return;

	list_head_t* curr;

	list_for_each(curr, &jobBatchList->headNode->listHead)
	{
		lng_t* currListNode = lng_getListGeneric(curr);
		job_batch_t* currJobBatch = (job_batch_t*) currListNode->dataPtr;
		sub_mtx_dim_t* currSubMtx = job_batch_getSubMtx(currJobBatch);

		sub_mtx_dim_t* cpy = sub_mtx_copyToPtr(currSubMtx);
		lg_addTailData((void*) cpy, subMtxList_inout);
	}
}

/**
 * Function definitions for generic list.
 */
FUNC_DEFINITION_LG_TOARRAY(job_batch, job_batch_t);
FUNC_DEFINITION_LG_TOARRAY_MULTIPLE(job_batch, job_batch_t);
FUNC_DEFINITION_LG_COPYADD(job_batch, job_batch_t, Tail);
FUNC_DEFINITION_LG_COPYADD(job_batch, job_batch_t, Front);
FUNC_DEFINITION_LG_PRINT(job_batch, job_batch_t);
FUNC_DEFINITION_LG_PRINT_MULTIPLE(job_batch, job_batch_t);
FUNC_DEFINITION_LG_DELETEDEEP(job_batch, job_batch_t);
FUNC_DEFINITION_LG_DELETEDEEP_MULTIPLE(job_batch, job_batch_t);
FUNC_DEFINITION_LG_PARTITIONDEEP(job_batch, job_batch_t, decompose);
FUNC_DEFINITION_LG_PARTITIONDEEP(job_batch, job_batch_t, scatter);
FUNC_DEFINITION_LG_PARTITIONDEEP(job_batch, job_batch_t, split);


// Job Batch Helper functions
// --------------------------------------------------------------------------------------------------------------------------

static void __job_batch_extractJobBatchListFromSubMtxTree_rowWiseColNet_CSR(
		spm_cmp_t* spmCsr, lg_t* jobBatchList, tree_node_t* head)
{
	if(head == NULL)
		return;

	if(tree_node_isLeaf(head))
	{
		sub_mtx_tree_t* subMtxNode = sub_mtx_tree_getSubMtxTree(head);
		DECIMAL rowStart = subMtxNode->subMtx->start.i;
		DECIMAL rowLength = subMtxNode->subMtx->length.i;

		sub_mtx_dim_t subMtxRowPar;
		sub_mtx_init(&subMtxRowPar, rowStart, 0, rowLength, spmCsr->colCount);

		// If sub matrix is empty (either dimensions are 0 or NNZ count)
		// then don't bother adding it as a job
		if(!sub_mtx_isEmpty_CSR(&subMtxRowPar, spmCsr))
		{
			job_batch_t* batch = job_batch_newSubMtx(&subMtxRowPar);
			lg_addTailData(batch, jobBatchList);
		}
	}
	else
	{
		__job_batch_extractJobBatchListFromSubMtxTree_rowWiseColNet_CSR(spmCsr, jobBatchList, head->left);
		__job_batch_extractJobBatchListFromSubMtxTree_rowWiseColNet_CSR(spmCsr, jobBatchList, head->right);
	}
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
