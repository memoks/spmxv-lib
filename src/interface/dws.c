
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "include/scheduler/ring_scheduler.h"
#include "include/scheduler/tree_scheduler.h"
#include "arch/generic/inc/spmxv_dws.h"
#include "arch/generic/inc/spmxv_static.h"

#include "include/interface/dws.h"


// values to function pointers in dws_routine struct
// -----------------------------------------------------------------------------------------------------------

extern void dws_routine_initialize(dws_routine* self);
extern void dws_routine_warmUp(dws_routine* self, vector_real_t* x);
extern void dws_routine_multiply(dws_routine* self, vector_real_t* x, vector_real_t* y_out);
extern void dws_routine_scalarMult(dws_routine* self, vector_real_t* x, REAL scalar, vector_real_t* y_out);
extern REAL dws_routine_innerProd(dws_routine* self, vector_real_t* r1, vector_real_t* r2);
extern void dws_routine_terminate(dws_routine* self);

// -----------------------------------------------------------------------------------------------------------

dws_routine* dws_routine_new(cli_options_t* options)
{
	dws_routine* dr = (dws_routine*) malloc(sizeof(dws_routine));

	dr->cr = common_routine_new(options);
	dr->spm = dr->cr->spm;

	dr->blockInfos = NULL;
	dr->subMtxHead = NULL;
	dr->execHistoryPerBlock = NULL;
	dr->subMtxArrayPerBlock = NULL;
	dr->subMtxCountPerBlock = NULL;

	dr->initialize = dws_routine_initialize;
	dr->warmUp = dws_routine_warmUp;
	dr->mult = dws_routine_multiply;
	dr->scalarMult = dws_routine_scalarMult;
	dr->innerProd = dws_routine_innerProd;
	dr->terminate = dws_routine_terminate;

	return dr;
}

void dws_routine_cleanUp(dws_routine* dr)
{
	dr->terminate(dr);
	free(dr);
}

// dws_routine struct function definitions
// -----------------------------------------------------------------------------------------------------------

void dws_routine_initialize(dws_routine* self)
{
	/*
	if(self->cr->options->useOrdering)
		self->subMtxHead = sub_mtx_tree_createSubMtxTree(self->cr->subMtxData, self->cr->options->netType);
	else
		self->subMtxHead = sub_mtx_tree_createUnorderedSubMtxTreeColNet(self->spm, DEFAULT_CACHE_SIZE);

	if(self->cr->options->scheduler == SCHEDULER_RING)
	{
		dws_scheduler_init = ring_scheduler_init;
		dws_scheduler_localityAwareStealHalf = ring_scheduler_localityAwareStealHalf;
		dws_scheduler_terminate = ring_scheduler_terminate;
	}
	else if(self->cr->options->scheduler == SCHEDULER_TREE)
	{
		dws_scheduler_init = tree_scheduler_init;
		dws_scheduler_localityAwareStealHalf = tree_scheduler_localityAwareStealHalf;
		dws_scheduler_terminate = tree_scheduler_terminate;
	}
	*/

	// spmxv_dws_generateBlockInfosSubMtxTree_CSR(&self->blockInfos, self->cr->options->numBlocks, self->spm, self->subMtxHead);
}

void dws_routine_warmUp(dws_routine* self, vector_real_t* x)
{
	vector_real_t* yTemp = vector_real_new(self->spm->rowCount);
	vector_real_t* xTemp = vector_real_new(self->spm->colCount);
	vector_real_copy(x, xTemp);

	int i;
	for(i = 0; i < self->cr->options->runCountWarmUp; ++i)
	{
		batch_list_deleteAllShallowMultiple(self->execHistoryPerBlock, self->cr->options->numBlocks);
		block_resetMultiple(self->blockInfos, self->cr->options->numBlocks);
		// spmxv_dws_colnet_CSR(self->spm, xTemp, yTemp, self->blockInfos, self->cr->options->numBlocks, &self->execHistoryPerBlock);
	}

	vector_real_delete(xTemp);
	vector_real_delete(yTemp);

//	batch_list_convertToSubMtxDenseRowParallelMultiple(
//			self->execHistoryPerBlock, self->cr->options->numBlocks,
//			&self->subMtxArrayPerBlock, &self->subMtxCountPerBlock);

	sub_mtx_tree_delete(&self->subMtxHead->node);
}

void dws_routine_multiply(dws_routine* self, vector_real_t* x, vector_real_t* y_out)
{
	// spmxv_static_colNet_CSR(self->spm, x, y_out, self->subMtxArrayPerBlock, self->subMtxCountPerBlock, self->cr->options->numBlocks);
}

void dws_routine_scalarMult(dws_routine* self, vector_real_t* x, REAL scalar, vector_real_t* y_inout)
{
	#pragma omp parallel num_threads(self->cr->options->numBlocks)
	{
		int threadId = omp_get_thread_num();
		int subMtxCount = self->subMtxCountPerBlock[threadId];
		sub_mtx_dim_t* subMtxs = self->subMtxArrayPerBlock[threadId];

		int i;
		int j;
		for(i = 0; i < subMtxCount; ++i)
		{
			for(j = subMtxs[i].start.i;
				j < subMtxs[i].start.i + subMtxs[i].length.i;
				++j)
			{
				y_inout->data[j] = x->data[j] * scalar;
			}
		}
	}
}

REAL dws_routine_innerProd(dws_routine* self, vector_real_t* r1, vector_real_t* r2)
{
	REAL* results = (REAL*) malloc(sizeof(REAL) * self->cr->options->numBlocks);
	#pragma omp parallel num_threads(self->cr->options->numBlocks)
	{
		int threadId = omp_get_thread_num();
		int subMtxCount = self->subMtxCountPerBlock[threadId];
		sub_mtx_dim_t* subMtxs = self->subMtxArrayPerBlock[threadId];

		results[threadId] = 0.0;
		int i;
		int j;
		for(i = 0; i < subMtxCount; ++i)
		{
			for(j = subMtxs[i].start.i; j < subMtxs[i].start.i + subMtxs[i].length.i; ++j)
			{
				results[threadId] += r1->data[j] * r2->data[j];
			}
		}
	}

	REAL accRes = 0.0;

	int i;
	for(i = 0; i< self->cr->options->numBlocks; ++i)
		accRes += results[i];

	return accRes;
}

void dws_routine_terminate(dws_routine* self)
{
	spmxv_dws_cleanUp(self->blockInfos, self->cr->options->numBlocks);

	int i;
	for(i = 0; i < self->cr->options->numBlocks; ++i)
	{
		batch_list_deleteAllShallow(self->execHistoryPerBlock[i]);
		batch_list_deleteSingleShallow(self->execHistoryPerBlock[i]);
		if(self->subMtxArrayPerBlock[i] != NULL)
			free(self->subMtxArrayPerBlock[i]);
	}

	free(self->subMtxArrayPerBlock);
	free(self->subMtxCountPerBlock);
	free(self->execHistoryPerBlock);
}
