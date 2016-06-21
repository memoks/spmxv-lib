
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/io/converter.h"
#include "include/control_unit/cu_options.h"

#include "include/control_unit/fast_run.h"

fast_run_t* fast_run_new(void)
{
	fast_run_t* fr = (fast_run_t*) malloc(sizeof(fast_run_t));
	fast_run_initDefault(fr);
	return fr;
}

void fast_run_initDefault(fast_run_t* fr)
{
	fr->quintets = NULL;
	fr->rowCount = -1;
	fr->colCount = -1;
	fr->nnz = -1;

	fr->numBlocks = -1;
	fr->threadsPerBlock = -1;
	fr->numThreads = -1;

	ipd_initDefault(&fr->ipd);

	fr->spmCsr = NULL;
	fr->_spmCsr = NULL;
	fr->_spmCsrCounterpart = NULL;

	fr->permutation = NULL;
	fr->_rowOrderLookupJDS = NULL;
	fr->yVectorRowOrderLookup = NULL;

	fr->subMtxCount = -1;

	fr->storageFormat = -1;
	fr->orderingType = -1;
	fr->partitionType = -1;
	fr->partitionMethod = -1;

	fr->simdLength = -1;
	fr->targetedCacheSizeKB = -1;

	fr->rowOrderLookup = NULL;
	fr->rowPartitioning = NULL;
	fr->colOrderLookup = NULL;
	fr->columnPartitioning = NULL;
}

void fast_run_init(
		fast_run_t* fr,
		quintet_t* quintets, int rowCount, int colCount, int nnz,
		int numBlocks, int threadsPerBlock,
		int targetedCacheSizeKB, int simdLength,
		int storageFormat, int orderingType, int partitionType, int partitionMethod,
		vector_int_t* rowOrderLookup, vector_int_t* rowPartitioning,
		vector_int_t* colOrderLookup, vector_int_t* colPartitioning)
{
	fr->quintets = quintets;
	fr->rowCount = rowCount;
	fr->colCount = colCount;
	fr->nnz = nnz;

	fr->numBlocks = numBlocks;
	fr->threadsPerBlock = threadsPerBlock;
	fr->numThreads = numBlocks * threadsPerBlock;

	ipd_initDefault(&fr->ipd);

	fr->storageFormat = storageFormat;
	fr->orderingType = orderingType;
	fr->partitionType = partitionType;
	fr->partitionMethod = partitionMethod;

	fr->simdLength = simdLength;
	fr->targetedCacheSizeKB = targetedCacheSizeKB;

	fr->rowOrderLookup = rowOrderLookup;
	fr->rowPartitioning = rowPartitioning;
	fr->colOrderLookup = colOrderLookup;
	fr->columnPartitioning = colPartitioning;

	if(partitionType == PARTITION_TYPE_1D_ROW_SLICE)
	{
		if(orderingType == ORDERING_TYPE_NONE)
		{
			converter_quintetToCSR_alloc(fr->quintets, &fr->_spmCsr, nnz, rowCount, colCount);

			sub_mtx_dim_t* subMtxArr = NULL;
			int subMtxCount = 0;

			partitioning_1DRowWise(fr->_spmCsr, partitionMethod, &fr->ipd, &subMtxArr, &subMtxCount,
					targetedCacheSizeKB, simdLength, simdLength, partitioning_calculateSizeInKBDefault);

			if(storageFormat == SPM_STORAGE_CSR)
			{
				fr->spmCsr = fr->_spmCsr;
			}
			else if(storageFormat == SPM_STORAGE_JDS ||
					storageFormat == SPM_STORAGE_HYBRID_JDS_CSR)
			{
				// convert quintets to JDS
				converter_quintetToJDSCounterpartCSR(
						fr->quintets, nnz, rowCount, colCount, subMtxArr, subMtxCount,
						&fr->_spmCsrCounterpart, &fr->permutation, &fr->_rowOrderLookupJDS);

				fr->spmCsr = fr->_spmCsrCounterpart;
			}

			fr->subMtxCount = subMtxCount;

			free(subMtxArr);
		}
		else if(orderingType == ORDERING_TYPE_COLUMN_NET)
		{
			// In ordered SpMxV partition count is the length-1 of one of two partitioning vectors (row or column)
			fr->subMtxCount = fr->rowPartitioning->length - 1;

			converter_quintetToCSR_alloc(fr->quintets, &fr->_spmCsr, nnz, rowCount, colCount);

			if(storageFormat == SPM_STORAGE_CSR)
			{
				fr->spmCsr = fr->_spmCsr;
			}
			else if(storageFormat == SPM_STORAGE_JDS ||
					storageFormat == SPM_STORAGE_HYBRID_JDS_CSR)
			{
				lg_t* jobBatchList = NULL;
				job_batch_extractJobBatchListFromPartitionVector_rowWiseColNet_CSR(
						fr->_spmCsr, fr->rowPartitioning, 0, fr->rowPartitioning->length - 1, &jobBatchList);

				lg_t* subMtxList = lg_new();
				lg_job_batch_extractSubMtx(jobBatchList, subMtxList);
				lg_job_batch_deleteDeep(jobBatchList);

				int subMtxCount = 0;
				sub_mtx_dim_t* subMtxArr = NULL;
				lg_sub_mtx_toArray(subMtxList, &subMtxArr, &subMtxCount);

				converter_quintetToJDSCounterpartCSR(
						quintets, nnz, rowCount, colCount, subMtxArr, subMtxCount,
						&fr->_spmCsrCounterpart, &fr->permutation, &fr->_rowOrderLookupJDS);

				fr->yVectorRowOrderLookup = fr->_rowOrderLookupJDS;
				fr->spmCsr = fr->_spmCsrCounterpart;

				free(subMtxArr);
				lg_sub_mtx_deleteDeep(subMtxList);
			}
		}
	}
}

void fast_run_terminate(fast_run_t* fr)
{
	ipd_terminate(&fr->ipd);

	if(fr->_spmCsr != NULL) spm_cmp_delete(fr->_spmCsr);
	if(fr->_spmCsrCounterpart != NULL) spm_cmp_delete(fr->_spmCsrCounterpart);
	if(fr->permutation != NULL) free(fr->permutation);
	if(fr->_rowOrderLookupJDS != NULL) vector_int_delete(fr->_rowOrderLookupJDS);
	if(fr->rowOrderLookup != NULL) free(fr->rowOrderLookup);
	if(fr->rowPartitioning != NULL) free(fr->rowPartitioning);
	if(fr->colOrderLookup != NULL) free(fr->colOrderLookup);
	if(fr->columnPartitioning != NULL) free(fr->columnPartitioning);
}

void fast_run_delete(fast_run_t* fr)
{
	fast_run_terminate(fr);
	free(fr);
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

