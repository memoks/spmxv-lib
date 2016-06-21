
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "include/config.h"
#include "include/io/converter.h"
#include "include/io/input_parser.h"
#include "include/io/cli.h"
#include "include/timer/custom_timer.h"
#include "include/data_structure/sub_mtx.h"
#include "include/data_structure/comm_tree.h"
#include "include/data_structure/list_generic.h"
#include "include/scheduler/job_queue.h"
#include "include/scheduler/job_batch.h"
#include "include/scheduler/ring_scheduler.h"
#include "include/scheduler/tree_scheduler.h"
#include "include/scheduler/block_info.h"
#include "include/algorithm/algorithm.h"
#include "include/task_decomposition/partitioning.h"

#include "arch/generic/inc/spmxv_sequential.h"
#include "arch/generic/inc/spmxv_omp_loop.h"
#include "arch/generic/inc/spmxv_omp_task.h"
#include "arch/generic/inc/spmxv_static.h"
#include "arch/generic/inc/spmxv_gws.h"
#include "arch/generic/inc/spmxv_dws.h"

#ifdef __ICC
#include "mkl.h"
#include "mkl_types.h"
#include "mkl_cblas.h"
#include "arch/icc/inc/spmxv_mkl.h"
#include "offload.h"
#endif


#define WARM_UP(runCount, func_calls)					\
	do{													\
		int i;											\
		for(i = 0; i < runCount; ++i)					\
		{												\
			func_calls									\
		}												\
	} while(0)

#define MEASURE(runCount, func_calls, timeElapsed)		\
	do{													\
		startTimer();									\
		for(i = 0; i < runCount; ++i)					\
		{												\
			func_calls									\
		}												\
		stopTimer(&timeElapsed);						\
	} while(0)

#define RUN(runCount, func_calls, buffer)				\
	do{													\
		int i;											\
		double totalTime = 0;							\
		double halfTime = 0;							\
														\
		MEASURE(runCount, func_calls, totalTime);		\
		MEASURE(runCount, func_calls, halfTime);		\
														\
		double diff = 0.5 -totalTime / halfTime;		\
		char temp[128];									\
		if(diff <= 0.1 || diff >= -0.1)					\
		{												\
			sprintf(temp, "%lf;", totalTime);			\
			strcat(buffer, temp);						\
		}												\
		else											\
		{												\
			sprintf(temp, "%lf;", diff);				\
			strcat(buffer, temp);						\
		}												\
	} while(0)

#define CLEAN_UP_JOB_BATCHES(jobBatchArrays, jobBatchCountPerArray, numArrays)					\
	do{																							\
		job_batch_deleteDenseArrMultiple(jobBatchArrays, jobBatchCountPerArray, numArrays);		\
		free(jobBatchCountPerArray);															\
	} while(0)


extern int calculateSum(int* arr, int length);

int main(int argsc, char* argsv[])
{
	// read command line parameters
	// ---------------------------------------------------------------

	cli_options_t options;
	cli_options_init(&options);
	cli_parseInput(argsc, argsv, &options);

	int numBlocks = options.numBlocks;
	int threadsPerBlock = options.threadsPerBlock;
	int numThreads = options.numBlocks * options.threadsPerBlock;
	float targetedCacheSize = options.targetedCacheSizeKB;

	// read matrix
	quintet_t* quintets = NULL;
	mm_info_t* mmInfo = NULL;
	DECIMAL nnz = 0;
	DECIMAL rowCount = 0;
	DECIMAL colCount = 0;

	input_readQuintets(options.mmfPath, &nnz, &rowCount, &colCount, &quintets, &mmInfo);
	// Using self-implemented hybrid sort instead of "input_sortQuintetRowValue(quintets, nnz);" for speed
	algorithm_parallelMergeQuickSort(quintets, nnz, numThreads, quintet_cmpRowValue);

	PRINTF("%s;", options.mmfName);

	// Unordered runs
	// --------------------------------------------------------------------------------------------------

#ifdef ARCH_MIC_OFFLOAD
	#pragma offload target(mic:3) in(quintets[0:nnz], nnz, rowCount, colCount, numBlocks, numThreads, threadsPerBlock, targetedCacheSize)
#endif
	{
		int block_schedType = DEFAULT_SCHEDULER;

		// Output buffers for CSV format
		// ---------------------------------------------------------------

		// char statisticsHeaders[16384]; statisticsHeaders[0] = '\0';
		// char algorithmHeaders[16384]; algorithmHeaders[0] = '\0';
		// char orderingHeaders[16384]; orderingHeaders[0] = '\0';

		char sequentialResults[512]; sequentialResults[0] = '\0';
		#ifdef __ICC
			char mklResults[512]; mklResults[0] = '\0';
		#endif
		char staticResults[512]; staticResults[0] = '\0';
		char ompDynamicResults[512]; ompDynamicResults[0] = '\0';
		char ompGuidedResults[512]; ompGuidedResults[0] = '\0';
		// char ompTaskResults[512]; ompTaskResults[0] = '\0';
		char gwsResults[512]; gwsResults[0] = '\0';
		char dwsRingResults[512]; dwsRingResults[0] = '\0';
		char dwsRingScatterSharedBlockResults[512]; dwsRingScatterSharedBlockResults[0] = '\0';

		REAL rowMinNNZ = 0;
		REAL rowAvgNNZ = 0;
		REAL rowMaxNNZ = 0;
		REAL columnMinNNZ = 0;
		REAL columnAvgNNZ = 0;
		REAL columnMaxNNZ = 0;
		int unorderedPartitionCount = 0;

		// headers in CSV format in order
		// ---------------------------------------------------------------

//		strcat(algorithmHeaders, "Overall_Statistics;;;;;;;;;;;;");
//		strcat(algorithmHeaders, ALGORITHM_SEQUENTIAL_STR); strcat(algorithmHeaders, ";");
//		#ifdef __ICC
//			strcat(algorithmHeaders, ALGORITHM_MKL_STR); strcat(algorithmHeaders, ";");
//			strcat(orderingHeaders, "un;");
//		#endif
//		strcat(algorithmHeaders, ALGORITHM_STATIC_STR); strcat(algorithmHeaders, ";");
//		// strcat(algorithmHeaders, ALGORITHM_OMP_TASK_STR); strcat(algorithmHeaders, ";");
//		strcat(algorithmHeaders, ALGORITHM_DWS_STR); strcat(algorithmHeaders, ";");
//		strcat(algorithmHeaders, ALGORITHM_DWS_SCATTER_SHARED_BLOCK_STR); strcat(algorithmHeaders, ";");

		// Overall statistics headers
//		strcat(statisticsHeaders, "rows;columns;nnz;minr;avgr;maxr;minc;avgc;maxc;cache_size_KB;k_un;k_ord;");


//		int i;
//		for(i = 0; i < 5; ++i)
//			strcat(orderingHeaders, "un;");

		// ---------------------------------------------------------------


		// generate x vector randomly
		vector_real_t* x = vector_real_random(colCount);
		vector_real_t* ySequentialCSR = vector_real_new(rowCount);

		// ---------------------------------------------------------------

		// vector_real_t* yOmpTask = vector_real_new(rowCount);
		// vector_real_t* yOmpTaskWarmUp = vector_real_new(rowCount);

		REAL targetedCacheSizeInKB = options.targetedCacheSizeKB;
		spm_cmp_t* spmCsr = NULL;
		spm_cmp_t* spmCsrCounterpart = NULL;
		DECIMAL* permutation = NULL;
		vector_int_t* rowOrderLookupJDS = NULL;

		converter_quintetToCSR_alloc(quintets, &spmCsr, nnz, rowCount, colCount);
		sub_mtx_tree_t* subMtxHead = NULL;
		partitioning_1DRowSliceRecursiveBipartition_CSR(
				spmCsr, targetedCacheSizeInKB, &subMtxHead, partitioning_calculateSizeInKBDefault);

		spmxv_sequential_CSR(spmCsr, x, ySequentialCSR);

		// extract statistics
		spm_cmp_extractNonZeroStatistics_CSR(
				spmCsr, &rowMinNNZ, &rowAvgNNZ, &rowMaxNNZ, &columnMinNNZ, &columnAvgNNZ, &columnMaxNNZ);

		// convert quintets to JDS
		sub_mtx_dim_t* subMtxArr = NULL;
		int subMtxCount = 0;
		lg_t* subMtxList = NULL;

		sub_mtx_tree_getLeafContentsDeep(subMtxHead, &subMtxList);
		lg_sub_mtx_toArray(subMtxList, &subMtxArr, &subMtxCount);

		converter_quintetToJDSCounterpartCSR(
				quintets, nnz, rowCount, colCount, subMtxArr, subMtxCount,
				&spmCsrCounterpart, &permutation, &rowOrderLookupJDS);

		free(subMtxArr);

		// Sequential CSR Implementation
		// ----------------------------------------------------------------------------------------------------
		{
			vector_real_t* ySequential = vector_real_new(rowCount);

			WARM_UP(
				options.runCountWarmUp,
				spmxv_sequential_CSR(spmCsr, x, ySequential);
			);
			RUN(
				options.runCountMeasure,
				spmxv_sequential_CSR(spmCsr, x, ySequential);,
				sequentialResults
			);

			vector_real_compare("SEQUENTIAL_CSR", ySequentialCSR, ySequential);
			vector_real_delete(ySequential);
		}

	#ifdef __ICC
		// MKL Implementation
		// ----------------------------------------------------------------------------------------------------
		{
			vector_real_t* yMkl = vector_real_new(rowCount);

			WARM_UP(
				options.runCountWarmUp,
				spmxv_mkl(spmCsr, x, yMkl);
			);
			RUN(
				options.runCountMeasure,
				spmxv_mkl(spmCsr, x, yMkl);,
				mklResults
			);

			vector_real_compare("MKL", ySequentialCSR, yMkl);
			vector_real_delete(yMkl);
		}
		// ----------------------------------------------------------------------------------------------------
#endif

		// OpenMP default scheduling
		{
			job_batch_t* jobBatchArr = NULL;
			int jobBatchCount = 0;
			spmxv_omp_loop_prepFromSubMtxList_rowWise_hybrid_JDS_CSR(spmCsrCounterpart, subMtxList, permutation, &jobBatchArr, &jobBatchCount);

			// OpenMP dynamic scheduling
			{
				vector_real_t* yOmpDynamic = vector_real_new(rowCount);

				int chunkSize;
				for(chunkSize = 1; chunkSize < 5; ++chunkSize)
				{
					WARM_UP(
						options.runCountWarmUp,
						spmxv_omp_loop_dynamic_hybrid_JDS_CSR(jobBatchArr, jobBatchCount, x, yOmpDynamic, chunkSize);
					);
					RUN(
						options.runCountMeasure,
						spmxv_omp_loop_dynamic_hybrid_JDS_CSR(jobBatchArr, jobBatchCount, x, yOmpDynamic, chunkSize);,
						ompDynamicResults
					);

					vector_real_reset(yOmpDynamic);
					spmxv_omp_loop_dynamic_hybrid_JDS_CSR(jobBatchArr, jobBatchCount, x, yOmpDynamic, chunkSize);
					input_reorderVector(&yOmpDynamic, rowOrderLookupJDS);
				}

				vector_real_compare("OMP_Dynamic", ySequentialCSR, yOmpDynamic);
				vector_real_delete(yOmpDynamic);
			}

			// OpenMP guided scheduling
			{
				vector_real_t* yOmpGuided = vector_real_new(rowCount);

				WARM_UP(
					options.runCountWarmUp,
					spmxv_omp_loop_guided_hybrid_JDS_CSR(jobBatchArr, jobBatchCount, x, yOmpGuided);
				);
				RUN(
					options.runCountMeasure,
					spmxv_omp_loop_guided_hybrid_JDS_CSR(jobBatchArr, jobBatchCount, x, yOmpGuided);,
					ompGuidedResults
				);

				vector_real_reset(yOmpGuided);
				spmxv_omp_loop_guided_hybrid_JDS_CSR(jobBatchArr, jobBatchCount, x, yOmpGuided);
				input_reorderVector(&yOmpGuided, rowOrderLookupJDS);

				vector_real_compare("OMP_Guided", ySequentialCSR, yOmpGuided);
				vector_real_delete(yOmpGuided);
			}


			job_batch_deleteDenseArr(jobBatchArr, jobBatchCount);
		}

		// ----------------------------------------------------------------------------------------------------
		{
			job_batch_t** jobBatchArrPerThread = NULL;
			int* jobBatchCountPerThread = NULL;
			job_batch_t** jobBatchArrPerBlock = NULL;
			int* jobBatchCountPerBlock = NULL;

			lg_t** execHistoryPerThread = NULL;
			lg_t** execHistoryPerBlock = NULL;

			unorderedPartitionCount = tree_node_getLeafCount(&subMtxHead->node);

			// Static Implementation
			// ------------------------------------------------------------------------------------------------
			{
				vector_real_t* yStatic = vector_real_new(rowCount);

				spmxv_static_prepFromSubMtxTree_rowWise_hybrid_JDS_CSR(
						spmCsrCounterpart, subMtxHead, numThreads, permutation, &jobBatchArrPerThread, &jobBatchCountPerThread);

				WARM_UP(
					options.runCountWarmUp,
					spmxv_static_rowWise_hybrid_JDS_CSR(x, yStatic, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_hybrid_JDS_CSR(x, yStatic, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					staticResults
				);

				vector_real_reset(yStatic);
				spmxv_static_rowWise_hybrid_JDS_CSR(x, yStatic, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				input_reorderVector(&yStatic, rowOrderLookupJDS);
				vector_real_compare("STATIC", ySequentialCSR, yStatic);

				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				vector_real_delete(yStatic);
			}

/*
			// OMP victim any
			// ------------------------------------------------------------------------------------------------
			{
				job_batch_t** denseJobBatchArrPerThread = NULL;
				int* denseJobBatchCountPerThread = NULL;
				execHistoryPerThread = NULL;

				spmxv_omp_task_prepFromSubMtxTree_colNet_hybrid_JDS_CSR(
						spmCsrCounterpart, subMtxHead, numThreads, permutation,
						&denseJobBatchArrPerThread, &denseJobBatchCountPerThread);

				WARM_UP(
					options.runCountWarmUp,
					lg_deleteShallowMultiple(execHistoryPerThread, numThreads);
					spmxv_omp_task_multiply_hybrid_JDS_CSR(
							x, yOmpTaskWarmUp, denseJobBatchArrPerThread, denseJobBatchCountPerThread,
							numThreads, &execHistoryPerThread);
				);

				lg_job_batch_toArrayMultiple(execHistoryPerThread, numThreads, &jobBatchArrPerThread, &jobBatchCountPerThread);

				vector_real_reset(yOmpTaskWarmUp);
				lg_deleteShallowMultiple(execHistoryPerThread, numThreads);
				spmxv_omp_task_multiply_hybrid_JDS_CSR(
						x, yOmpTaskWarmUp, denseJobBatchArrPerThread, denseJobBatchCountPerThread,
						numThreads, &execHistoryPerThread);
				input_reorderVector(&yOmpTaskWarmUp, rowOrderLookupJDS);

				// Clean up dynamic scheduling data structures
				lg_deleteShallowMultiple(execHistoryPerThread, numThreads);


				RUN(
					options.runCountMeasure,
					spmxv_static_multiply_colNet_hybrid_JDS_CSR(
							x, yOmpTask, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					ompTaskResults
				);


				vector_real_reset(yOmpTask);
				spmxv_static_multiply_colNet_hybrid_JDS_CSR(x, yOmpTask, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				input_reorderVector(&yOmpTask, rowOrderLookupJDS);


				// Clean up
				CLEAN_UP_JOB_BATCHES(denseJobBatchArrPerThread, denseJobBatchCountPerThread, numThreads);
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}
*/



			// Distributed Work Stealing Algorithm (Ring Scheduler)
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				block_schedType = SCHEDULER_RING;

				execHistoryPerThread = NULL;
				block_info_t** blockInfos = NULL;
				spmxv_dws_prepFromSubMtxTree_colNet_hybrid_JDS_CSR(
						spmCsrCounterpart, permutation, subMtxHead, numThreads, options.stealTreshold, block_schedType, &blockInfos);

				vector_real_t* yDwsRingWarmUp = vector_real_new(rowCount);
				WARM_UP(
					options.runCountWarmUp,
					lg_deleteShallowMultiple(execHistoryPerThread, numThreads);
					spmxv_dws_rowWise_hybrid_JDS_CSR(x, yDwsRingWarmUp, blockInfos, numThreads, &execHistoryPerThread);
					block_reformMultiple(blockInfos, execHistoryPerThread, numThreads);
				);


				vector_real_reset(yDwsRingWarmUp);
				block_resetMultiple(blockInfos, numThreads);
				spmxv_dws_rowWise_hybrid_JDS_CSR(x, yDwsRingWarmUp, blockInfos, numThreads, &execHistoryPerBlock);
				lg_deleteShallowMultiple(execHistoryPerBlock, numThreads);
				input_reorderVector(&yDwsRingWarmUp, rowOrderLookupJDS);
				vector_real_compare("DWS_RING_WARMUP", ySequentialCSR, yDwsRingWarmUp);


				lg_job_batch_toArrayMultiple(execHistoryPerThread, numThreads, &jobBatchArrPerThread, &jobBatchCountPerThread);


				// Clean up dynamic scheduling data structures
				spmxv_dws_cleanUp(blockInfos, numThreads, block_schedType);
				lg_deleteShallowMultiple(execHistoryPerThread, numThreads);
				vector_real_delete(yDwsRingWarmUp);

				vector_real_t* yDwsRing = vector_real_new(rowCount);
				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_hybrid_JDS_CSR(
							x, yDwsRing, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					dwsRingResults
				);


				vector_real_reset(yDwsRing);
				spmxv_static_rowWise_hybrid_JDS_CSR(
						x, yDwsRing, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				input_reorderVector(&yDwsRing, rowOrderLookupJDS);
				vector_real_compare("DWS_RING", ySequentialCSR, yDwsRing);


				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				vector_real_delete(yDwsRing);
			}


			// Distributed Work Stealing Algorithm (Ring Scheduler) Shared Queue
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				block_schedType = SCHEDULER_RING;

				execHistoryPerThread = NULL;
				execHistoryPerBlock = NULL;
				block_info_t** blockInfos = NULL;
				spmxv_dws_prepFromSubMtxTree_colNet_hybrid_JDS_CSR(
						spmCsrCounterpart, permutation, subMtxHead, numBlocks, options.stealTreshold, block_schedType, &blockInfos);

				vector_real_t* yDwsRingScatterSharedBlockWarmUp = vector_real_new(rowCount);
				WARM_UP(
					options.runCountWarmUp,
					lg_deleteShallowMultiple(execHistoryPerThread, numThreads);
					lg_deleteShallowMultiple(execHistoryPerBlock, numBlocks);
					spmxv_dws_sharedBlock_rowWise_hybrid_JDS_CSR(
							x, yDwsRingScatterSharedBlockWarmUp, blockInfos, numBlocks, threadsPerBlock,
							&execHistoryPerThread, &execHistoryPerBlock);
					block_reformMultiple(blockInfos, execHistoryPerBlock, numBlocks);
				);

				lg_job_batch_toArrayMultiple(execHistoryPerThread, numThreads, &jobBatchArrPerThread, &jobBatchCountPerThread);


				vector_real_reset(yDwsRingScatterSharedBlockWarmUp);
				lg_deleteShallowMultiple(execHistoryPerThread, numThreads);
				lg_deleteShallowMultiple(execHistoryPerBlock, numBlocks);
				block_resetMultiple(blockInfos, numBlocks);
				spmxv_dws_sharedBlock_rowWise_hybrid_JDS_CSR(
						x, yDwsRingScatterSharedBlockWarmUp, blockInfos, numBlocks, threadsPerBlock,
						&execHistoryPerThread, &execHistoryPerBlock);
				input_reorderVector(&yDwsRingScatterSharedBlockWarmUp, rowOrderLookupJDS);
				vector_real_compare("DWS_RING_SHARED_BLOCK_WARMUP", ySequentialCSR, yDwsRingScatterSharedBlockWarmUp);


				// Clean up dynamic scheduling data structures
				spmxv_dws_cleanUp(blockInfos, numBlocks, block_schedType);
				lg_deleteShallowMultiple(execHistoryPerThread, numThreads);
				lg_deleteShallowMultiple(execHistoryPerBlock, numBlocks);
				vector_real_delete(yDwsRingScatterSharedBlockWarmUp);

				vector_real_t* yDwsRingScatterSharedBlock = vector_real_new(rowCount);
				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_hybrid_JDS_CSR(x, yDwsRingScatterSharedBlock,
							jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					dwsRingScatterSharedBlockResults
				);

				vector_real_reset(yDwsRingScatterSharedBlock);
				spmxv_static_rowWise_hybrid_JDS_CSR(
						x, yDwsRingScatterSharedBlock, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				input_reorderVector(&yDwsRingScatterSharedBlock, rowOrderLookupJDS);
				vector_real_compare("DWS_RING_SHARED_BLOCK", ySequentialCSR, yDwsRingScatterSharedBlock);

				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				vector_real_delete(yDwsRingScatterSharedBlock);
			}


			// Correctness test
			// ------------------------------------------------------------------------------------------------

			// vector_real_compare("OMP_TASK_WARMUP", comparisonVector, yOmpTaskWarmUp);
			// vector_real_compare("OMP_TASK", comparisonVector, yOmpTask);

			// Clean up
			// ------------------------------------------------------------------------------------------------

			// vector_real_delete(yOmpTask);
			// vector_real_delete(yOmpTaskWarmUp);

			free(permutation);
			vector_int_delete(rowOrderLookupJDS);
			sub_mtx_tree_delete(&subMtxHead->node);
			spm_cmp_delete(spmCsr);
			spm_cmp_delete(spmCsrCounterpart);
		}

		PRINTF("%d;%d;%d;%0.0f;%0.1f;%0.0f;%0.0f;%0.1f;%0.1f;%0.1f;%d;",
				rowCount, colCount, nnz, rowMinNNZ, rowAvgNNZ, rowMaxNNZ,
				columnMinNNZ, columnAvgNNZ, columnMaxNNZ, options.targetedCacheSizeKB,
				unorderedPartitionCount);

		PRINTF("%s", sequentialResults);
	#ifdef __ICC
		PRINTF("%s", mklResults);
	#endif
		PRINTF("%s", staticResults);
		// strcat(results, ompTaskResults);
		PRINTF("%s", ompDynamicResults);
		PRINTF("%s", ompGuidedResults);
		PRINTF("%s", gwsResults);
		PRINTF("%s", dwsRingResults);
		PRINTF("%s", dwsRingScatterSharedBlockResults);
		PRINTF("\n");


		vector_real_delete(ySequentialCSR);
		vector_real_delete(x);

		lg_sub_mtx_deleteDeep(subMtxList);
	}

	free(quintets);
	cli_options_deleteNonPtr(&options);
	mm_info_delete(mmInfo);

	return EXIT_SUCCESS;
}

int calculateSum(int* arr, int length)
{
	int sum = 0;

	int i;
	for(i = 0; i < length; ++i)
		sum += arr[i];

	return sum;
}
