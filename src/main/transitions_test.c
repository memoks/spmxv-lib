
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

int block_schedType = DEFAULT_SCHEDULER;

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

	// read matrix
	quintet_t* triplets = NULL;
	mm_info_t* mmInfo = NULL;
	DECIMAL nnz = 0;
	DECIMAL rowCount = 0;
	DECIMAL colCount = 0;

	input_readQuintets(options.mmfPath, &nnz, &rowCount, &colCount, &triplets, &mmInfo);
	// Using self-implemented hybrid sort instead of "input_sortTripletRowValue(triplets, nnz);" for speed
	algorithm_parallelMergeQuickSort(triplets, nnz, numThreads, quintet_cmpRowValue);

	// Don't show extra information ...
	// cli_print(&options);
	// PRINTF("\n");
	// mm_info_print(mmInfo);

	// Output buffers for CSV format
	// ---------------------------------------------------------------

	char statisticsHeaders[16384]; statisticsHeaders[0] = '\0';
	char algorithmHeaders[16384]; algorithmHeaders[0] = '\0';
	char orderingHeaders[16384]; orderingHeaders[0] = '\0';

	char sequentialResults[512]; sequentialResults[0] = '\0';
#ifdef __ICC
	char mklResults[512]; mklResults[0] = '\0';
#endif
	char staticResults[512]; staticResults[0] = '\0';
	char staticScatterResults[512]; staticScatterResults[0] = '\0';
	char ompTaskResults[512]; ompTaskResults[0] = '\0';
	char ompTaskScatterResults[512]; ompTaskScatterResults[0] = '\0';
	char gwsResults[512]; gwsResults[0] = '\0';
	char dwsRingResults[512]; dwsRingResults[0] = '\0';
	char dwsRingScatterResults[512]; dwsRingScatterResults[0] = '\0';
	char dwsRingScatterSharedBlockResults[512]; dwsRingScatterSharedBlockResults[0] = '\0';
	char dwsTreeResults[512]; dwsTreeResults[0] = '\0';
	char dwsTreeScatterResults[512]; dwsTreeScatterResults[0] = '\0';
	char dwsTreeScatterSharedBlockResults[512]; dwsTreeScatterSharedBlockResults[0] = '\0';

	char results[32868]; results[0] = '\0';

	REAL rowMinNNZ = 0;
	REAL rowAvgNNZ = 0;
	REAL rowMaxNNZ = 0;
	REAL columnMinNNZ = 0;
	REAL columnAvgNNZ = 0;
	REAL columnMaxNNZ = 0;
	int unorderedPartitionCount = 0;
	int orderedPartitionCount = 0;

	// headers in CSV format in order
	// ---------------------------------------------------------------

	strcat(algorithmHeaders, "Overall_Statistics;;;;;;;;;;;;");
	strcat(algorithmHeaders, ALGORITHM_SEQUENTIAL_STR); strcat(algorithmHeaders, ";;");
#ifdef __ICC
	strcat(algorithmHeaders, ALGORITHM_MKL_STR); strcat(algorithmHeaders, ";;");
	strcat(orderingHeaders, "un;ord;");
#endif
	strcat(algorithmHeaders, ALGORITHM_STATIC_STR); strcat(algorithmHeaders, ";;");
	strcat(algorithmHeaders, ALGORITHM_STATIC_SCATTER_STR); strcat(algorithmHeaders, ";;");
	strcat(algorithmHeaders, ALGORITHM_OMP_TASK_STR); strcat(algorithmHeaders, ";;");
	strcat(algorithmHeaders, ALGORITHM_OMP_TASK_SCATTER_STR); strcat(algorithmHeaders, ";;");
	strcat(algorithmHeaders, ALGORITHM_GWS_STR); strcat(algorithmHeaders, ";;");
	strcat(algorithmHeaders, ALGORITHM_DWS_STR); strcat(algorithmHeaders, ";;");
	strcat(algorithmHeaders, ALGORITHM_DWS_SCATTER_STR); strcat(algorithmHeaders, ";;");
	strcat(algorithmHeaders, ALGORITHM_DWS_SCATTER_SHARED_BLOCK_STR); strcat(algorithmHeaders, ";;");
	strcat(algorithmHeaders, "dws_tree;;");
	strcat(algorithmHeaders, "dws_tree_scatter;;");
	strcat(algorithmHeaders, "dws_tree_scatter_shared_block;;");

	// Overall statistics headers
	strcat(statisticsHeaders, "rows;columns;nnz;minr;avgr;maxr;minc;avgc;maxc;cache_size_KB;k_un;k_ord;");


	int i;
	for(i = 0; i < 12; ++i)
		strcat(orderingHeaders, "un;ord;");


	// ---------------------------------------------------------------

	// generate x vector randomly
	vector_real_t* x = vector_real_random(colCount);

	vector_real_t* ySequentialCSR = vector_real_new(rowCount);
	vector_real_t* ySequential = vector_real_new(rowCount);


	// Unordered runs
	// --------------------------------------------------------------------------------------------------

	{
		int i;
		vector_real_t* yMkl = vector_real_new(rowCount);
		vector_real_t* yStatic = vector_real_new(rowCount);
		vector_real_t* yStaticScatter = vector_real_new(rowCount);
		vector_real_t* yOmpTask = vector_real_new(rowCount);
		vector_real_t* yOmpTaskWarmUp = vector_real_new(rowCount);
		vector_real_t* yOmpTaskScatter = vector_real_new(rowCount);
		vector_real_t* yOmpTaskScatterWarmUp = vector_real_new(rowCount);
		vector_real_t* yGws = vector_real_new(rowCount);
		vector_real_t* yGwsWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsRing = vector_real_new(rowCount);
		vector_real_t* yDwsRingWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsRingScatter = vector_real_new(rowCount);
		vector_real_t* yDwsRingScatterWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsRingScatterSharedBlock = vector_real_new(rowCount);
		vector_real_t* yDwsRingScatterSharedBlockWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsTree = vector_real_new(rowCount);
		vector_real_t* yDwsTreeWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsTreeScatter = vector_real_new(rowCount);
		vector_real_t* yDwsTreeScatterWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsTreeScatterSharedBlock = vector_real_new(rowCount);
		vector_real_t* yDwsTreeScatterSharedBlockWarmUp = vector_real_new(rowCount);


		REAL targetedCacheSizeInKB = options.targetedCacheSizeKB;
		spm_cmp_t* spmCsr = NULL;
		spm_cmp_t* spmCsrCounterpart = NULL;
		DECIMAL* permutation = NULL;
		vector_int_t* rowOrderLookupJDS = NULL;

		converter_quintetToCSR_alloc(triplets, &spmCsr, nnz, rowCount, colCount);
		sub_mtx_tree_t* subMtxHead = NULL;
		partitioning_1DRowSliceRecursiveBipartition_CSR(
				spmCsr, targetedCacheSizeInKB, &subMtxHead, partitioning_calculateSizeInKBDefault);

		spmxv_sequential_CSR(spmCsr, x, ySequentialCSR);

		// extract statistics
		spm_cmp_extractNonZeroStatistics_CSR(
				spmCsr, &rowMinNNZ, &rowAvgNNZ, &rowMaxNNZ, &columnMinNNZ, &columnAvgNNZ, &columnMaxNNZ);

		// convert triplets to JDS
		job_batch_t* jobBatchArr = NULL;
		int jobBatchCount = 0;


		batch_list_t* batchList = NULL;

		batch_list_extractBatchListFromSubMtxTree_colNet_CSR(spmCsr, &subMtxHead->node, &batchList);

		batch_list_convertToJobBatchDense(batchList, &jobBatchArr, &jobBatchCount);



		batch_list_deleteAll(batchList);
		batch_list_deleteSingleShallow(batchList);

		converter_quintetToJDSCounterpartCSR(
				triplets, nnz, rowCount, colCount, jobBatchArr, jobBatchCount,
				&spmCsrCounterpart, &permutation, &rowOrderLookupJDS);
		free(jobBatchArr);


		// Sequential CSR Implementation
		// ----------------------------------------------------------------------------------------------------
		{
			WARM_UP(
				options.runCountWarmUp,
				spmxv_sequential_CSR(spmCsr, x, ySequential);
			);
			RUN(
				options.runCountMeasure,
				spmxv_sequential_CSR(spmCsr, x, ySequential);,
				sequentialResults
			);
		}

#ifdef __ICC
		// MKL Implementation
		// ----------------------------------------------------------------------------------------------------
		{
			WARM_UP(
				options.runCountWarmUp,
				spmxv_mkl(spmCsr, x, yMkl);
			);
			RUN(
				options.runCountMeasure,
				spmxv_mkl(spmCsr, x, yMkl);,
				mklResults
			);
		}
		// ----------------------------------------------------------------------------------------------------
#endif

		// Unordered algorithms with cache blocking
		// ----------------------------------------------------------------------------------------------------
		{
			job_batch_t** jobBatchArrPerThread = NULL;
			int* jobBatchCountPerThread = NULL;
			job_batch_t** jobBatchArrPerBlock = NULL;
			int* jobBatchCountPerBlock = NULL;

			batch_list_t** execHistoryPerThread = NULL;
			batch_list_t** execHistoryPerBlock = NULL;

			unorderedPartitionCount = tree_node_getLeafCount(&subMtxHead->node);

			// Static Implementation
			// ------------------------------------------------------------------------------------------------
			{
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


				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}

			// Static Scatter Implementation
			// ------------------------------------------------------------------------------------------------
			{
				spmxv_static_prepFromSubMtxTreeScatter_rowWise_hybrid_JDS_CSR(
						spmCsrCounterpart, subMtxHead, numBlocks, threadsPerBlock, permutation,
						&jobBatchArrPerThread, &jobBatchCountPerThread);

				WARM_UP(
					options.runCountWarmUp,
					spmxv_static_rowWise_hybrid_JDS_CSR(
							x, yStaticScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				);
				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_hybrid_JDS_CSR(
							x, yStaticScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					staticScatterResults
				);

				vector_real_reset(yStaticScatter);
				spmxv_static_rowWise_hybrid_JDS_CSR(x, yStaticScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				input_reorderVector(&yStaticScatter, rowOrderLookupJDS);

				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}

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
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					spmxv_omp_task_multiply_hybrid_JDS_CSR(
							x, yOmpTaskWarmUp, denseJobBatchArrPerThread, denseJobBatchCountPerThread,
							numThreads, &execHistoryPerThread);
				);

				batch_list_convertToJobBatchDenseMultiple(
						execHistoryPerThread, numThreads, FALSE, &jobBatchArrPerThread, &jobBatchCountPerThread);


				vector_real_reset(yOmpTaskWarmUp);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				spmxv_omp_task_multiply_hybrid_JDS_CSR(
						x, yOmpTaskWarmUp, denseJobBatchArrPerThread, denseJobBatchCountPerThread,
						numThreads, &execHistoryPerThread);
				input_reorderVector(&yOmpTaskWarmUp, rowOrderLookupJDS);


				// Clean up dynamic scheduling data structures
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_hybrid_JDS_CSR(
							x, yOmpTask, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					ompTaskResults
				);


				vector_real_reset(yOmpTask);
				spmxv_static_rowWise_hybrid_JDS_CSR(x, yOmpTask, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				input_reorderVector(&yOmpTask, rowOrderLookupJDS);


				// Clean up
				CLEAN_UP_JOB_BATCHES(denseJobBatchArrPerThread, denseJobBatchCountPerThread, numThreads);
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}

			// OMP victim any scatter
			// ------------------------------------------------------------------------------------------------
			{
				job_batch_t** denseJobBatchArrPerThread = NULL;
				int* denseJobBatchCountPerThread = NULL;
				execHistoryPerThread = NULL;

				spmxv_omp_task_prepFromSubMtxTreeScatter_colNet_hybrid_JDS_CSR(
						spmCsrCounterpart, subMtxHead, numBlocks, threadsPerBlock, permutation,
						&denseJobBatchArrPerThread, &denseJobBatchCountPerThread);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					spmxv_omp_task_multiply_hybrid_JDS_CSR(x, yOmpTaskScatterWarmUp,
							denseJobBatchArrPerThread, denseJobBatchCountPerThread, numThreads, &execHistoryPerThread);
				);

				batch_list_convertToJobBatchDenseMultiple(
						execHistoryPerThread, numThreads, FALSE, &jobBatchArrPerThread, &jobBatchCountPerThread);


				vector_real_reset(yOmpTaskScatterWarmUp);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				spmxv_omp_task_multiply_hybrid_JDS_CSR(x, yOmpTaskScatterWarmUp,
						denseJobBatchArrPerThread, denseJobBatchCountPerThread, numThreads, &execHistoryPerThread);
				input_reorderVector(&yOmpTaskScatterWarmUp, rowOrderLookupJDS);


				// Clean up dynamic scheduling data structures
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_hybrid_JDS_CSR(
							x, yOmpTaskScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					ompTaskScatterResults
				);

				vector_real_reset(yOmpTaskScatter);
				spmxv_static_rowWise_hybrid_JDS_CSR(
						x, yOmpTaskScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				input_reorderVector(&yOmpTaskScatter, rowOrderLookupJDS);

				// Clean up
				CLEAN_UP_JOB_BATCHES(denseJobBatchArrPerThread, denseJobBatchCountPerThread, numThreads);
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}


			// Global work stealing algorithm
			// ------------------------------------------------------------------------------------------------
			{
				execHistoryPerThread = NULL;
				job_queue_t* globalQueue = NULL;
				spmxv_gws_prepFromSubMtxTree_colNet_hybrid_JDS_CSR(
						spmCsrCounterpart, permutation, subMtxHead, options.stealTreshold, &globalQueue);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					job_queue_reset(globalQueue);
					spmxv_gws_multiply_hybrid_JDS_CSR(x, yGwsWarmUp, globalQueue, numThreads, &execHistoryPerThread);
				);

				batch_list_convertToJobBatchDenseMultiple(
						execHistoryPerThread, numThreads, FALSE, &jobBatchArrPerThread, &jobBatchCountPerThread);


				vector_real_reset(yGwsWarmUp);
				job_queue_reset(globalQueue);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				spmxv_gws_multiply_hybrid_JDS_CSR(x, yGwsWarmUp, globalQueue, numThreads, &execHistoryPerThread);
				input_reorderVector(&yGwsWarmUp, rowOrderLookupJDS);


				// Clean up dynamic scheduling data structures
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				job_queue_delete(globalQueue);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_hybrid_JDS_CSR(
							x, yGws, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					gwsResults
				);

				vector_real_reset(yGws);
				spmxv_static_rowWise_hybrid_JDS_CSR(
						x, yGws, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				input_reorderVector(&yGws, rowOrderLookupJDS);

				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}

			// Distributed Work Stealing Algorithm (Ring Scheduler)
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				block_schedType = SCHEDULER_RING;
				dws_scheduler_init = ring_scheduler_init;
				dws_scheduler_localityAwareStealHalf = ring_scheduler_localityAwareStealHalf;
				dws_scheduler_terminate = ring_scheduler_terminate;

				execHistoryPerThread = NULL;
				block_info_t** blockInfos = NULL;
				spmxv_dws_prepFromSubMtxTree_colNet_hybrid_JDS_CSR(
						spmCsrCounterpart, permutation, subMtxHead, numThreads, options.stealTreshold, &blockInfos);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					spmxv_dws_rowWise_hybrid_JDS_CSR(x, yDwsRingWarmUp, blockInfos, numThreads, &execHistoryPerThread);
					block_reformMultiple(blockInfos, execHistoryPerThread, numThreads);
				);


				vector_real_reset(yDwsRingWarmUp);
				block_resetMultiple(blockInfos, numThreads);
				spmxv_dws_rowWise_hybrid_JDS_CSR(x, yDwsRingWarmUp, blockInfos, numThreads, &execHistoryPerBlock);
				batch_list_deleteAllShallowMultiple(execHistoryPerBlock, numThreads);
				input_reorderVector(&yDwsRingWarmUp, rowOrderLookupJDS);


				batch_list_convertToJobBatchDenseMultiple(
						execHistoryPerThread, numThreads, FALSE, &jobBatchArrPerThread, &jobBatchCountPerThread);

				// Clean up dynamic scheduling data structures
				spmxv_dws_cleanUp(blockInfos, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

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


				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}

			// Distributed Work Stealing Algorithm Scatter (Ring Scheduler)
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				block_schedType = SCHEDULER_RING;
				dws_scheduler_init = ring_scheduler_init;
				dws_scheduler_localityAwareStealHalf = ring_scheduler_localityAwareStealHalf;
				dws_scheduler_terminate = ring_scheduler_terminate;

				execHistoryPerThread = NULL;
				block_info_t** blockInfos = NULL;
				spmxv_dws_prepFromSubMtxTreeScatter_colNet_hybrid_JDS_CSR(
						spmCsrCounterpart, permutation, subMtxHead, numBlocks, threadsPerBlock, options.stealTreshold, &blockInfos);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					spmxv_dws_rowWise_hybrid_JDS_CSR(x, yDwsRingScatterWarmUp, blockInfos, numThreads, &execHistoryPerThread);
					block_reformMultiple(blockInfos, execHistoryPerThread, numThreads);
				);

				batch_list_convertToJobBatchDenseMultiple(
						execHistoryPerThread, numThreads, FALSE, &jobBatchArrPerThread, &jobBatchCountPerThread);


				vector_real_reset(yDwsRingScatterWarmUp);
				block_resetMultiple(blockInfos, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				spmxv_dws_rowWise_hybrid_JDS_CSR(x, yDwsRingScatterWarmUp, blockInfos, numThreads, &execHistoryPerThread);
				input_reorderVector(&yDwsRingScatterWarmUp, rowOrderLookupJDS);


				// Clean up dynamic scheduling data structures
				spmxv_dws_cleanUp(blockInfos, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_hybrid_JDS_CSR(
							x, yDwsRingScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					dwsRingScatterResults
				);


				vector_real_reset(yDwsRingScatter);
				spmxv_static_rowWise_hybrid_JDS_CSR(
						x, yDwsRingScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				input_reorderVector(&yDwsRingScatter, rowOrderLookupJDS);


				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}


			// Distributed Work Stealing Algorithm (Ring Scheduler) Shared Queue
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				block_schedType = SCHEDULER_RING;
				dws_scheduler_init = ring_scheduler_init;
				dws_scheduler_localityAwareStealHalf = ring_scheduler_localityAwareStealHalf;
				dws_scheduler_terminate = ring_scheduler_terminate;

				execHistoryPerThread = NULL;
				execHistoryPerBlock = NULL;
				block_info_t** blockInfos = NULL;
				spmxv_dws_prepFromSubMtxTree_colNet_hybrid_JDS_CSR(
						spmCsrCounterpart, permutation, subMtxHead, numBlocks, options.stealTreshold, &blockInfos);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					batch_list_deleteAllShallowMultiple(execHistoryPerBlock, numBlocks);
					spmxv_dws_sharedBlock_rowWise_hybrid_JDS_CSR(
							x, yDwsRingScatterSharedBlockWarmUp, blockInfos, numBlocks, threadsPerBlock,
							&execHistoryPerThread, &execHistoryPerBlock);
					block_reformMultiple(blockInfos, execHistoryPerBlock, numBlocks);
				);

				batch_list_convertToJobBatchDenseMultiple(
						execHistoryPerThread, numThreads, FALSE, &jobBatchArrPerThread, &jobBatchCountPerThread);


				vector_real_reset(yDwsRingScatterSharedBlockWarmUp);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerBlock, numBlocks);
				block_resetMultiple(blockInfos, numBlocks);
				spmxv_dws_sharedBlock_rowWise_hybrid_JDS_CSR(
						x, yDwsRingScatterSharedBlockWarmUp, blockInfos, numBlocks, threadsPerBlock,
						&execHistoryPerThread, &execHistoryPerBlock);
				input_reorderVector(&yDwsRingScatterSharedBlockWarmUp, rowOrderLookupJDS);


				// Clean up dynamic scheduling data structures
				spmxv_dws_cleanUp(blockInfos, numBlocks);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerBlock, numBlocks);

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

				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}


			// Distributed Work Stealing Algorithm (Tree Scheduler)
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				block_schedType = SCHEDULER_TREE;
				dws_scheduler_init = tree_scheduler_init;
				dws_scheduler_localityAwareStealHalf = tree_scheduler_localityAwareStealHalf;
				dws_scheduler_terminate = tree_scheduler_terminate;

				execHistoryPerThread = NULL;
				block_info_t** blockInfos = NULL;
				spmxv_dws_prepFromSubMtxTree_colNet_hybrid_JDS_CSR(
						spmCsrCounterpart, permutation, subMtxHead, numThreads, options.stealTreshold, &blockInfos);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					spmxv_dws_rowWise_hybrid_JDS_CSR(x, yDwsTreeWarmUp, blockInfos, numThreads, &execHistoryPerThread);
					block_reformMultiple(blockInfos, execHistoryPerThread, numThreads);
				);

				batch_list_convertToJobBatchDenseMultiple(
						execHistoryPerThread, numThreads, FALSE, &jobBatchArrPerThread, &jobBatchCountPerThread);


				vector_real_reset(yDwsTreeWarmUp);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				block_resetMultiple(blockInfos, numThreads);
				spmxv_dws_rowWise_hybrid_JDS_CSR(x, yDwsTreeWarmUp, blockInfos, numThreads, &execHistoryPerThread);
				input_reorderVector(&yDwsTreeWarmUp, rowOrderLookupJDS);

				// Clean up dynamic scheduling data structures
				spmxv_dws_cleanUp(blockInfos, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_hybrid_JDS_CSR(
							x, yDwsTree, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					dwsTreeResults
				);

				vector_real_reset(yDwsTree);
				spmxv_static_rowWise_hybrid_JDS_CSR(
						x, yDwsTree, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				input_reorderVector(&yDwsTree, rowOrderLookupJDS);

				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}

			// Distributed Work Stealing Algorithm Scatter (Tree Scheduler)
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				block_schedType = SCHEDULER_TREE;
				dws_scheduler_init = tree_scheduler_init;
				dws_scheduler_localityAwareStealHalf = tree_scheduler_localityAwareStealHalf;
				dws_scheduler_terminate = tree_scheduler_terminate;

				execHistoryPerThread = NULL;
				block_info_t** blockInfos = NULL;
				spmxv_dws_prepFromSubMtxTreeScatter_colNet_hybrid_JDS_CSR(
						spmCsrCounterpart, permutation, subMtxHead, numBlocks, threadsPerBlock, options.stealTreshold, &blockInfos);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					spmxv_dws_rowWise_hybrid_JDS_CSR(x, yDwsTreeScatterWarmUp, blockInfos, numThreads, &execHistoryPerThread);
					block_reformMultiple(blockInfos, execHistoryPerThread, numThreads);
				);

				batch_list_convertToJobBatchDenseMultiple(
						execHistoryPerThread, numThreads, FALSE, &jobBatchArrPerThread, &jobBatchCountPerThread);


				vector_real_reset(yDwsTreeScatterWarmUp);
				block_resetMultiple(blockInfos, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				spmxv_dws_rowWise_hybrid_JDS_CSR(x, yDwsTreeScatterWarmUp, blockInfos, numThreads, &execHistoryPerThread);
				input_reorderVector(&yDwsTreeScatterWarmUp, rowOrderLookupJDS);


				// Clean up dynamic scheduling data structures
				spmxv_dws_cleanUp(blockInfos, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_hybrid_JDS_CSR(
							x, yDwsTreeScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					dwsTreeScatterResults
				);


				vector_real_reset(yDwsTreeScatter);
				spmxv_static_rowWise_hybrid_JDS_CSR(
						x, yDwsTreeScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				input_reorderVector(&yDwsTreeScatter, rowOrderLookupJDS);


				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}


			// Distributed Work Stealing Algorithm (Tree Scheduler) Shared Queue
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				block_schedType = SCHEDULER_TREE;
				dws_scheduler_init = tree_scheduler_init;
				dws_scheduler_localityAwareStealHalf = tree_scheduler_localityAwareStealHalf;
				dws_scheduler_terminate = tree_scheduler_terminate;

				execHistoryPerThread = NULL;
				execHistoryPerBlock = NULL;
				block_info_t** blockInfos = NULL;
				spmxv_dws_prepFromSubMtxTree_colNet_hybrid_JDS_CSR(
						spmCsrCounterpart, permutation, subMtxHead, numBlocks, options.stealTreshold, &blockInfos);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					batch_list_deleteAllShallowMultiple(execHistoryPerBlock, numBlocks);
					spmxv_dws_sharedBlock_rowWise_hybrid_JDS_CSR(
							x, yDwsTreeScatterSharedBlockWarmUp, blockInfos, numBlocks, threadsPerBlock,
							&execHistoryPerThread, &execHistoryPerBlock);
					block_reformMultiple(blockInfos, execHistoryPerBlock, numBlocks);
				);

				batch_list_convertToJobBatchDenseMultiple(
						execHistoryPerThread, numThreads, FALSE, &jobBatchArrPerThread, &jobBatchCountPerThread);


				vector_real_reset(yDwsTreeScatterSharedBlockWarmUp);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerBlock, numBlocks);
				block_resetMultiple(blockInfos, numBlocks);
				spmxv_dws_sharedBlock_rowWise_hybrid_JDS_CSR(
						x, yDwsTreeScatterSharedBlockWarmUp, blockInfos, numBlocks, threadsPerBlock,
						&execHistoryPerThread, &execHistoryPerBlock);
				input_reorderVector(&yDwsTreeScatterSharedBlockWarmUp, rowOrderLookupJDS);


				// Clean up dynamic scheduling data structures
				spmxv_dws_cleanUp(blockInfos, numBlocks);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerBlock, numBlocks);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_hybrid_JDS_CSR(x, yDwsTreeScatterSharedBlock,
							jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					dwsTreeScatterSharedBlockResults
				);

				vector_real_reset(yDwsTreeScatterSharedBlock);
				spmxv_static_rowWise_hybrid_JDS_CSR(
						x, yDwsTreeScatterSharedBlock, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				input_reorderVector(&yDwsTreeScatterSharedBlock, rowOrderLookupJDS);

				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}

		}

		// Correctness test
		// ------------------------------------------------------------------------------------------------

		vector_real_compare("SEQUENTIAL_CSR", ySequentialCSR, ySequential);
#ifdef __ICC
		vector_real_compare("MKL", ySequentialCSR, yMkl);
		vector_real_t* comparisonVector = yMkl;
#else
		vector_real_t* comparisonVector = ySequentialCSR;
#endif

		vector_real_compare("STATIC", comparisonVector, yStatic);
		vector_real_compare("STATIC_SCATTER", comparisonVector, yStaticScatter);
		vector_real_compare("OMP_TASK_WARMUP", comparisonVector, yOmpTaskWarmUp);
		vector_real_compare("OMP_TASK", comparisonVector, yOmpTask);
		vector_real_compare("OMP_TASK_SCATTER_WARMUP", comparisonVector, yOmpTaskScatterWarmUp);
		vector_real_compare("OMP_TASK_SCATTER", comparisonVector, yOmpTaskScatter);
		vector_real_compare("GWS_WARMUP", comparisonVector, yGwsWarmUp);
		vector_real_compare("GWS", comparisonVector, yGws);
		vector_real_compare("DWS_RING_WARMUP", comparisonVector, yDwsRingWarmUp);
		vector_real_compare("DWS_RING", comparisonVector, yDwsRing);
		vector_real_compare("DWS_RING_SCATTER_WARMUP", comparisonVector, yDwsRingScatterWarmUp);
		vector_real_compare("DWS_RING_SCATTER", comparisonVector, yDwsRingScatter);
		vector_real_compare("DWS_RING_SCATTER_SHARED_BLOCK_WARMUP", comparisonVector, yDwsRingScatterSharedBlockWarmUp);
		vector_real_compare("DWS_RING_SCATTER_SHARED_BLOCK", comparisonVector, yDwsRingScatterSharedBlock);
		vector_real_compare("DWS_TREE_WARMUP", comparisonVector, yDwsTreeWarmUp);
		vector_real_compare("DWS_TREE", comparisonVector, yDwsTree);
		vector_real_compare("DWS_TREE_SCATTER_WARMUP", comparisonVector, yDwsTreeScatterWarmUp);
		vector_real_compare("DWS_TREE_SCATTER", comparisonVector, yDwsTreeScatter);
		vector_real_compare("DWS_TREE_SCATTER_SHARED_BLOCK_WARMUP", comparisonVector, yDwsTreeScatterSharedBlockWarmUp);
		vector_real_compare("DWS_TREE_SCATTER_SHARED_BLOCK", comparisonVector, yDwsTreeScatterSharedBlock);


		// Clean up
		// ------------------------------------------------------------------------------------------------

		vector_real_delete(yMkl);
		vector_real_delete(yStatic);
		vector_real_delete(yStaticScatter);
		vector_real_delete(yOmpTask);
		vector_real_delete(yOmpTaskWarmUp);
		vector_real_delete(yOmpTaskScatter);
		vector_real_delete(yOmpTaskScatterWarmUp);
		vector_real_delete(yGws);
		vector_real_delete(yGwsWarmUp);
		vector_real_delete(yDwsRing);
		vector_real_delete(yDwsRingWarmUp);
		vector_real_delete(yDwsRingScatter);
		vector_real_delete(yDwsRingScatterWarmUp);
		vector_real_delete(yDwsRingScatterSharedBlock);
		vector_real_delete(yDwsRingScatterSharedBlockWarmUp);
		vector_real_delete(yDwsTree);
		vector_real_delete(yDwsTreeWarmUp);
		vector_real_delete(yDwsTreeScatter);
		vector_real_delete(yDwsTreeScatterWarmUp);
		vector_real_delete(yDwsTreeScatterSharedBlock);
		vector_real_delete(yDwsTreeScatterSharedBlockWarmUp);

		free(permutation);
		vector_int_delete(rowOrderLookupJDS);
		sub_mtx_tree_delete(&subMtxHead->node);
		spm_cmp_delete(spmCsr);
		spm_cmp_delete(spmCsrCounterpart);
	}

	strcat(results, sequentialResults);
#ifdef __ICC
	strcat(results, mklResults);
#endif
	strcat(results, staticResults);
	strcat(results, staticScatterResults);
	strcat(results, ompTaskResults);
	strcat(results, ompTaskScatterResults);
	strcat(results, gwsResults);
	strcat(results, dwsRingResults);
	strcat(results, dwsRingScatterResults);
	strcat(results, dwsRingScatterSharedBlockResults);
	strcat(results, dwsTreeResults);
	strcat(results, dwsTreeScatterResults);
	strcat(results, dwsTreeScatterSharedBlockResults);

	PRINTF(";%s\n", algorithmHeaders);
	PRINTF(";%s%s\n", statisticsHeaders, orderingHeaders);
	PRINTF("%s;%d;%d;%d;%0.0f;%0.1f;%0.0f;%0.0f;%0.1f;%0.1f;%0.1f;%d;%d;%s\n",
			options.mmfName, rowCount, colCount, nnz, rowMinNNZ, rowAvgNNZ, rowMaxNNZ,
			columnMinNNZ, columnAvgNNZ, columnMaxNNZ, options.targetedCacheSizeKB,
			unorderedPartitionCount, orderedPartitionCount, results);
	PRINTF("\n\n");


	vector_real_delete(ySequential);
	vector_real_delete(ySequentialCSR);

	free(triplets);
	cli_options_deleteNonPtr(&options);
	vector_real_delete(x);
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




