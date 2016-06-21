
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
#include "include/data_structure/spm_storage.h"
#include "include/scheduler/job_queue.h"
#include "include/scheduler/job_batch.h"
#include "include/scheduler/ring_scheduler.h"
#include "include/scheduler/tree_scheduler.h"
#include "include/scheduler/block_info.h"

#include "arch/generic/inc/spmxv_sequential.h"
#include "arch/generic/inc/spmxv_omp_loop.h"
#include "arch/generic/inc/spmxv_omp_task.h"
#include "arch/generic/inc/spmxv_static.h"
#include "arch/generic/inc/spmxv_gws.h"
#include "arch/generic/inc/spmxv_dws.h"

#ifdef __ICC
#include "arch/icc/inc/spmxv_mkl.h"
#endif

#define WARM_UP(runCount, func_calls)				\
	do{												\
		int i;										\
		for(i = 0; i < runCount; ++i)				\
		{											\
			func_calls								\
		}											\
	} while(0)

#define MEASURE(runCount, func_calls, timeElapsed)	\
	do{												\
		startTimer();								\
		for(i = 0; i < runCount; ++i)				\
		{											\
			func_calls								\
		}											\
		stopTimer(&timeElapsed);					\
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
	quintet_sortRowValue(triplets, nnz);

	// Don't print extra information
	// cli_print(&options);
	// PRINTF("\n");
	// mm_info_print(mmInfo);

	// Output buffers for CSV format
	// ---------------------------------------------------------------

	char statisticsHeaders[16384]; statisticsHeaders[0] = '\0';
	char algorithmHeaders[16384]; algorithmHeaders[0] = '\0';
	char orderingHeaders[16384]; orderingHeaders[0] = '\0';

	char results[32768]; results[0] = '\0';
	char sequentialResults[512]; sequentialResults[0] = '\0';
	char ompLoop32Results[512]; ompLoop32Results[0] = '\0';
	char ompLoop64Results[512]; ompLoop64Results[0] = '\0';
	char ompLoop128Results[512]; ompLoop128Results[0] = '\0';
	char mklResults[512]; mklResults[0] = '\0';
	char staticResults[512]; staticResults[0] = '\0';
	char staticScatterResults[512]; staticScatterResults[0] = '\0';
	char ompTaskResults[512]; ompTaskResults[0] = '\0';
	char ompTaskScatterResults[512]; ompTaskScatterResults[0] = '\0';
	char gwsResults[512]; gwsResults[0] = '\0';
	char dwsRingResults[512]; dwsRingResults[0] = '\0';
	char dwsRingScatterResults[512]; dwsRingScatterResults[0] = '\0';
	char dwsRingScatterSharedBlockResults[512]; dwsRingScatterSharedBlockResults[0] = '\0';
	char dwsTreeResults[512]; dwsTreeResults[0] = '\0';

	REAL rowMinNNZ = 0;
	REAL rowAvgNNZ = 0;
	REAL rowMaxNNZ = 0;
	REAL columnMinNNZ = 0;
	REAL columnAvgNNZ = 0;
	REAL columnMaxNNZ = 0;
	int unorderedPartitionCount;
	int orderedPartitionCount;

	// headers in CSV format in order
	// ---------------------------------------------------------------

	strcat(algorithmHeaders, ";;;;;;;;;;;;");
	strcat(algorithmHeaders, ALGORITHM_SEQUENTIAL_STR); strcat(algorithmHeaders, ";;");
	strcat(algorithmHeaders, "loop"); strcat(algorithmHeaders, "_32;;");
	strcat(algorithmHeaders, "loop"); strcat(algorithmHeaders, "_64;;");
	strcat(algorithmHeaders, "loop"); strcat(algorithmHeaders, "_128;;");
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

	strcat(statisticsHeaders, "rows;");
	strcat(statisticsHeaders, "columns;");
	strcat(statisticsHeaders, "nnz;");
	strcat(statisticsHeaders, "minr;");
	strcat(statisticsHeaders, "avgr;");
	strcat(statisticsHeaders, "maxr;");
	strcat(statisticsHeaders, "minc;");
	strcat(statisticsHeaders, "avgc;");
	strcat(statisticsHeaders, "maxc;");
	strcat(statisticsHeaders, "cache_size_KB;");
	strcat(statisticsHeaders, "k_un;");
	strcat(statisticsHeaders, "k_ord;");

	int i;
	for(i = 0; i < 13; ++i)
		strcat(orderingHeaders, "un;ord;");

	// ---------------------------------------------------------------

	// generate x vector randomly
	vector_real_t* x = vector_real_random(colCount);

	vector_real_t* ySequential = vector_real_new(rowCount);

	// Unordered runs
	// --------------------------------------------------------------------------------------------------

	{
		int i;
		vector_real_t* yOmpLoopDynamic32 = vector_real_new(rowCount);
		vector_real_t* yOmpLoopDynamic64 = vector_real_new(rowCount);
		vector_real_t* yOmpLoopDynamic128 = vector_real_new(rowCount);
		vector_real_t* yMkl = vector_real_new(rowCount);
		vector_real_t* yStatic = vector_real_new(rowCount);
		vector_real_t* yStaticScatter = vector_real_new(rowCount);
		vector_real_t* yOmpTask = vector_real_new(rowCount);
		vector_real_t* yOmpTaskWarmUp = vector_real_new(rowCount);
		vector_real_t* yOmpTaskScatter = vector_real_new(rowCount);
		vector_real_t* yOmpTaskScatterWarmUp = vector_real_new(rowCount);
		vector_real_t* yGws = vector_real_new(rowCount);
		vector_real_t* yGwsWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsTree = vector_real_new(rowCount);
		vector_real_t* yDwsTreeWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsRing = vector_real_new(rowCount);
		vector_real_t* yDwsRingWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsRingScatter = vector_real_new(rowCount);
		vector_real_t* yDwsRingScatterWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsRingScatterSharedBlock = vector_real_new(rowCount);
		vector_real_t* yDwsRingScatterSharedBlockWarmUp = vector_real_new(rowCount);

		spm_cmp_t* spm = NULL;

		REAL targetedCacheSizeInKB = options.targetedCacheSizeKB;

		converter_quintetToCSR(triplets, &spm, nnz, rowCount, colCount);

		// extract statistics
		spm_cmp_extractNonZeroStatistics_CSR(spm, &rowMinNNZ, &rowAvgNNZ, &rowMaxNNZ, &columnMinNNZ, &columnAvgNNZ, &columnMaxNNZ);

		// Sequential CSR Implementation
		// ------------------------------------------------------------------------------------------------
		{
			WARM_UP(
				options.runCountWarmUp,
				spmxv_sequential_CSR(spm, x, ySequential);
			);
			RUN(
				options.runCountMeasure,
				spmxv_sequential_CSR(spm, x, ySequential);,
				sequentialResults
			);

		}

		// Parallel Dynamic Implementation 32
		// ------------------------------------------------------------------------------------------------
		{
			int chunkSize = 32;

			WARM_UP(
				options.runCountWarmUp,
				spmxv_omp_loop_dynamic_CSR(spm, x, yOmpLoopDynamic32, chunkSize);
			);
			RUN(
				options.runCountMeasure,
				spmxv_omp_loop_dynamic_CSR(spm, x, yOmpLoopDynamic32, chunkSize);,
				ompLoop32Results
			);
		}

		// Parallel Dynamic Implementation 64
		// ------------------------------------------------------------------------------------------------
		{
			int chunkSize = 64;

			WARM_UP(
				options.runCountWarmUp,
				spmxv_omp_loop_dynamic_CSR(spm, x, yOmpLoopDynamic64, chunkSize);
			);
			RUN(
				options.runCountMeasure,
				spmxv_omp_loop_dynamic_CSR(spm, x, yOmpLoopDynamic64, chunkSize);,
				ompLoop64Results
			);
		}

		// Parallel Dynamic Implementation 128
		// ------------------------------------------------------------------------------------------------
		{
			int chunkSize = 128;

			WARM_UP(
				options.runCountWarmUp,
				spmxv_omp_loop_dynamic_CSR(spm, x, yOmpLoopDynamic128, chunkSize);
			);
			RUN(
				options.runCountMeasure,
				spmxv_omp_loop_dynamic_CSR(spm, x, yOmpLoopDynamic128, chunkSize);,
				ompLoop128Results
			);
		}

#ifdef __ICC
		// MKL Implementation
		// ------------------------------------------------------------------------------------------------
		{
			WARM_UP(
				options.runCountWarmUp,
				spmxv_mkl(spm, x, yMkl);
			);
			RUN(
				options.runCountMeasure,
				spmxv_mkl(spm, x, yMkl);,
				mklResults
			);
		}
		// ------------------------------------------------------------------------------------------------
#endif

		// Unordered algorithms with cache blocking
		// ------------------------------------------------------------------------------------------------

		{
			job_batch_t** jobBatchArrPerThread = NULL;
			int* jobBatchCountPerThread = NULL;
			job_batch_t** jobBatchArrPerBlock = NULL;
			int* jobBatchCountPerBlock = NULL;

			batch_list_t** execHistoryPerThread = NULL;
			batch_list_t** execHistoryPerBlock = NULL;

			sub_mtx_tree_t* subMtxHead = sub_mtx_tree_createUnorderedSubMtxTree_colNet_CSR(spm, targetedCacheSizeInKB);
			unorderedPartitionCount = tree_node_getLeafCount(&subMtxHead->node);

			// Static Implementation
			// ------------------------------------------------------------------------------------------------
			{
				spmxv_static_prepFromSubMtxTree_rowWise_CSR(spm, subMtxHead, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

				WARM_UP(
					options.runCountWarmUp,
					spmxv_static_rowWise_CSR(spm, x, yStatic, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				);
				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_CSR(spm, x, yStatic, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					staticResults
				);

				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}

			// Static Scatter Implementation
			// ------------------------------------------------------------------------------------------------
			{
				spmxv_static_prepFromSubMtxTreeScatter_rowWise_CSR(spm, subMtxHead, numBlocks, threadsPerBlock, &jobBatchArrPerThread, &jobBatchCountPerThread);

				WARM_UP(
					options.runCountWarmUp,
					spmxv_static_rowWise_CSR(spm, x, yStaticScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
				);
				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_CSR(spm, x, yStaticScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					staticScatterResults
				);

				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}

			// OMP victim any
			// ------------------------------------------------------------------------------------------------
			{
				job_batch_t** denseJobBatchArrPerThread = NULL;
				int* denseJobBatchCountPerThread = NULL;
				execHistoryPerThread = NULL;

				spmxv_omp_task_prepFromSubMtxTree_rowWise_CSR(spm, subMtxHead, numThreads, &denseJobBatchArrPerThread, &denseJobBatchCountPerThread);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					spmxv_omp_task_multiply_CSR(spm, x, yOmpTaskWarmUp, denseJobBatchArrPerThread, denseJobBatchCountPerThread, numThreads, &execHistoryPerThread);
				);

				batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

				// Clean up dynamic scheduling data structures
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_CSR(spm, x, yOmpTask, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					ompTaskResults
				);

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

				spmxv_omp_task_prepFromSubMtxTreeScatter_rowWise_CSR(
						spm, subMtxHead, numBlocks, threadsPerBlock, &denseJobBatchArrPerThread, &denseJobBatchCountPerThread);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					spmxv_omp_task_multiply_CSR(spm, x, yOmpTaskScatterWarmUp, denseJobBatchArrPerThread, denseJobBatchCountPerThread, numThreads, &execHistoryPerThread);
				);

				batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

				// Clean up dynamic scheduling data structures
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_CSR(spm, x, yOmpTaskScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					ompTaskScatterResults
				);

				// Clean up
				CLEAN_UP_JOB_BATCHES(denseJobBatchArrPerThread, denseJobBatchCountPerThread, numThreads);
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}

			// Global work stealing algorithm
			// ------------------------------------------------------------------------------------------------
			{
				execHistoryPerThread = NULL;
				job_queue_t* globalQueue = NULL;
				spmxv_gws_prepFromSubMtxTree_colNet_CSR(spm, subMtxHead, options.stealTreshold, &globalQueue);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					job_queue_reset(globalQueue);
					spmxv_gws_multiply_CSR(spm, x, yGwsWarmUp, globalQueue, numThreads, &execHistoryPerThread);
				);

				batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

				// Clean up dynamic scheduling data structures
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				job_queue_delete(globalQueue);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_CSR(spm, x, yGws, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					gwsResults
				);

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
				spmxv_dws_prepFromSubMtxTree_colNet_CSR(spm, subMtxHead, numThreads, options.stealTreshold, &blockInfos);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					spmxv_dws_rowWise_CSR(spm, x, yDwsRingWarmUp, blockInfos, numThreads, &execHistoryPerThread);
					block_reformMultiple(blockInfos, execHistoryPerThread, numThreads);
				);

				batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

				// Clean up dynamic scheduling data structures
				spmxv_dws_cleanUp(blockInfos, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_CSR(spm, x, yDwsRing, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					dwsRingResults
				);

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
				spmxv_dws_prepFromSubMtxTreeScatter_colNet_CSR(spm, subMtxHead, numBlocks, threadsPerBlock, options.stealTreshold, &blockInfos);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					spmxv_dws_rowWise_CSR(spm, x, yDwsRingScatterWarmUp, blockInfos, numThreads, &execHistoryPerThread);
					block_reformMultiple(blockInfos, execHistoryPerThread, numThreads);
				);

				batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

				// Clean up dynamic scheduling data structures
				spmxv_dws_cleanUp(blockInfos, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_CSR(spm, x, yDwsRingScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					dwsRingScatterResults
				);

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
				block_info_t** blockInfos = NULL;
				spmxv_dws_prepFromSubMtxTree_colNet_CSR(spm, subMtxHead, numBlocks, options.stealTreshold, &blockInfos);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					batch_list_deleteAllShallowMultiple(execHistoryPerBlock, numBlocks);
					spmxv_dws_sharedBlock_rowWise_CSR(spm, x, yDwsRingScatterSharedBlockWarmUp, blockInfos, numBlocks, threadsPerBlock, &execHistoryPerThread, &execHistoryPerBlock);
					block_reformMultiple(blockInfos, execHistoryPerBlock, numBlocks);
				);

				batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

				// Clean up dynamic scheduling data structures
				spmxv_dws_cleanUp(blockInfos, numBlocks);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerBlock, numBlocks);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_CSR(spm, x, yDwsRingScatterSharedBlock, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					dwsRingScatterSharedBlockResults
				);

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
				spmxv_dws_prepFromSubMtxTree_colNet_CSR(spm, subMtxHead, numThreads, options.stealTreshold, &blockInfos);

				WARM_UP(
					options.runCountWarmUp,
					batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
					spmxv_dws_rowWise_CSR(spm, x, yDwsTree, blockInfos, numThreads, &execHistoryPerThread);
					block_reformMultiple(blockInfos, execHistoryPerThread, numThreads);
				);

				batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

				// Clean up dynamic scheduling data structures
				spmxv_dws_cleanUp(blockInfos, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

				RUN(
					options.runCountMeasure,
					spmxv_static_rowWise_CSR(spm, x, yDwsTreeWarmUp, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
					dwsTreeResults
				);

				// Clean up
				CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			}

			sub_mtx_tree_delete(&subMtxHead->node);
		}

		// Correctness test
		// ------------------------------------------------------------------------------------------------

		vector_real_compare("OMP_LOOP_DYNAMIC_32", ySequential, yOmpLoopDynamic32);
		vector_real_compare("OMP_LOOP_DYNAMIC_64", ySequential, yOmpLoopDynamic64);
		vector_real_compare("OMP_LOOP_DYNAMIC_128", ySequential, yOmpLoopDynamic128);
#ifdef __ICC
		vector_real_compare("MKL", ySequential, yMkl);
#endif
		vector_real_compare("STATIC", ySequential, yStatic);
		vector_real_compare("STATIC_SCATTER", ySequential, yStaticScatter);
		vector_real_compare("OMP_TASK_WARMUP", ySequential, yOmpTaskWarmUp);
		vector_real_compare("OMP_TASK", ySequential, yOmpTask);
		vector_real_compare("OMP_TASK_SCATTER_WARMUP", ySequential, yOmpTaskScatterWarmUp);
		vector_real_compare("OMP_TASK_SCATTER", ySequential, yOmpTaskScatter);
		vector_real_compare("GWS_WARMUP", ySequential, yGwsWarmUp);
		vector_real_compare("GWS", ySequential, yGws);
		vector_real_compare("DWS_RING_WARMUP", ySequential, yDwsRingWarmUp);
		vector_real_compare("DWS_RING", ySequential, yDwsRing);
		vector_real_compare("DWS_RING_SCATTER_WARMUP", ySequential, yDwsRingScatterWarmUp);
		vector_real_compare("DWS_RING_SCATTER", ySequential, yDwsRingScatter);
		vector_real_compare("DWS_RING_SCATTER_SHARED_BLOCK_WARMUP", ySequential, yDwsRingScatterSharedBlockWarmUp);
		vector_real_compare("DWS_RING_SCATTER_SHARED_BLOCK", ySequential, yDwsRingScatterSharedBlock);
		vector_real_compare("DWS_TREE_WARMUP", ySequential, yDwsTreeWarmUp);
		vector_real_compare("DWS_TREE", ySequential, yDwsTree);

		// Clean up
		// ------------------------------------------------------------------------------------------------

		vector_real_delete(yOmpLoopDynamic32);
		vector_real_delete(yOmpLoopDynamic64);
		vector_real_delete(yOmpLoopDynamic128);
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

		spm_cmp_delete(spm);
	}

	// Ordered Runs
	// ---------------------------------------------------------------

	vector_real_t* ySequentialOrdered = vector_real_new(rowCount);

	{
		vector_real_t* yOmpLoopDynamic32 = vector_real_new(rowCount);
		vector_real_t* yOmpLoopDynamic64 = vector_real_new(rowCount);
		vector_real_t* yOmpLoopDynamic128 = vector_real_new(rowCount);
		vector_real_t* yOmpLoopGuided = vector_real_new(rowCount);
		vector_real_t* yMkl = vector_real_new(rowCount);
		vector_real_t* yStatic = vector_real_new(rowCount);
		vector_real_t* yStaticScatter = vector_real_new(rowCount);
		vector_real_t* yOmpTaskWarmUp = vector_real_new(rowCount);
		vector_real_t* yOmpTask = vector_real_new(rowCount);
		vector_real_t* yOmpTaskScatterWarmUp = vector_real_new(rowCount);
		vector_real_t* yOmpTaskScatter = vector_real_new(rowCount);
		vector_real_t* yGwsWarmUp = vector_real_new(rowCount);
		vector_real_t* yGws = vector_real_new(rowCount);
		vector_real_t* yDwsRingWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsRing = vector_real_new(rowCount);
		vector_real_t* yDwsRingScatterSharedBlockWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsRingScatterSharedBlock = vector_real_new(rowCount);
		vector_real_t* yDwsRingScatterWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsRingScatter = vector_real_new(rowCount);
		vector_real_t* yDwsTreeWarmUp = vector_real_new(rowCount);
		vector_real_t* yDwsTree = vector_real_new(rowCount);

		// order & partitioning
		vector_int_t* rowOrderLookup = NULL;
		vector_int_t* rowPartitioning = NULL;
		vector_int_t* colOrderLookup = NULL;
		vector_int_t* colPartitioning = NULL;

		input_kPatohOrder(
				options.mmfPath, options.mmfDir, options.mmfName, options.targetedCacheSizeKB, options.partitionTypeStr,
				triplets, nnz, rowCount, colCount, &rowOrderLookup, &rowPartitioning, &colOrderLookup, &colPartitioning);

		input_orderVector(&x, colOrderLookup);

		// In ordered SpMxV partition count is the length of one of two partitioning vectors (row or column)
		orderedPartitionCount = rowPartitioning->length;

		spm_cmp_t* spm = NULL;
		converter_quintetToCSR(triplets, &spm, nnz, rowCount, colCount);

		int i;

		// Sequential CSR Implementation
		// ------------------------------------------------------------------------------------------------
		{
			WARM_UP(
				options.runCountWarmUp,
				spmxv_sequential_CSR(spm, x, ySequentialOrdered);
			);

			RUN(
				options.runCountMeasure,
				spmxv_sequential_CSR(spm, x, ySequentialOrdered);,
				sequentialResults
			);
		}

		// Parallel Dynamic Implementation 32
		// ------------------------------------------------------------------------------------------------
		{
			int chunkSize = 32;

			WARM_UP(
				options.runCountWarmUp,
				spmxv_omp_loop_dynamic_CSR(spm, x, yOmpLoopDynamic32, chunkSize);
			);

			RUN(
				options.runCountMeasure,
				spmxv_omp_loop_dynamic_CSR(spm, x, yOmpLoopDynamic32, chunkSize);,
				ompLoop32Results
			);
		}

		// Parallel Dynamic Implementation 64
		// ------------------------------------------------------------------------------------------------
		{
			int chunkSize = 64;

			WARM_UP(
				options.runCountWarmUp,
				spmxv_omp_loop_dynamic_CSR(spm, x, yOmpLoopDynamic64, chunkSize);
			);

			RUN(
				options.runCountMeasure,
				spmxv_omp_loop_dynamic_CSR(spm, x, yOmpLoopDynamic64, chunkSize);,
				ompLoop64Results
			);
		}

		// Parallel Dynamic Implementation 128
		// ------------------------------------------------------------------------------------------------
		{
			int chunkSize = 128;

			WARM_UP(
				options.runCountWarmUp,
				spmxv_omp_loop_dynamic_CSR(spm, x, yOmpLoopDynamic128, chunkSize);
			);

			RUN(
				options.runCountMeasure,
				spmxv_omp_loop_dynamic_CSR(spm, x, yOmpLoopDynamic128, chunkSize);,
				ompLoop128Results
			);
		}

#ifdef __ICC
		// MKL Implementation
		// ------------------------------------------------------------------------------------------------
		{
			WARM_UP(
				options.runCountWarmUp,
				spmxv_mkl(spm, x, yMkl);
			);
			RUN(
				options.runCountMeasure,
				spmxv_mkl(spm, x, yMkl);,
				mklResults
			);
		}
		// ------------------------------------------------------------------------------------------------
#endif

		job_batch_t** jobBatchArrPerThread = NULL;
		int* jobBatchCountPerThread = NULL;
		job_batch_t** jobBatchArrPerBlock = NULL;
		int* jobBatchCountPerBlock = NULL;

		batch_list_t** execHistoryPerBlock = NULL;
		batch_list_t** execHistoryPerThread = NULL;

		// Static Implementation
		// ------------------------------------------------------------------------------------------------
		{
			spmxv_static_prepFromPartitionVector_rowWise_colNet_CSR(
					spm, numThreads, TRUE, rowPartitioning, colPartitioning, &jobBatchArrPerThread, &jobBatchCountPerThread);

			WARM_UP(
				options.runCountWarmUp,
				spmxv_static_rowWise_CSR(spm, x, yStatic, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			);

			RUN(
				options.runCountMeasure,
				spmxv_static_rowWise_CSR(spm, x, yStatic, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
				staticResults
			);

			// Clean up
			CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
		}

		// Static Scatter Implementation
		// ------------------------------------------------------------------------------------------------
		{
			spmxv_static_prepFromPartitionVectorScatter_rowWise_colNet_CSR(
					spm, numBlocks, threadsPerBlock, rowPartitioning, colPartitioning, &jobBatchArrPerThread, &jobBatchCountPerThread);

			WARM_UP(
				options.runCountWarmUp,
				spmxv_static_rowWise_CSR(spm, x, yStaticScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			);

			RUN(
				options.runCountMeasure,
				spmxv_static_rowWise_CSR(spm, x, yStaticScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
				staticScatterResults
			);

			// Clean up
			CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
		}

		// OMP victim any
		// ------------------------------------------------------------------------------------------------
		{
			execHistoryPerThread = NULL;
			job_batch_t** initialJobBatchArrPerThread = NULL;
			int* initialJobBatchCountPerThread = NULL;

			spmxv_omp_task_prepFromPartitionVector_rowWise_CSR(
					spm, numThreads, rowPartitioning, colPartitioning, &initialJobBatchArrPerThread, &initialJobBatchCountPerThread);

			WARM_UP(
				options.runCountWarmUp,
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				spmxv_omp_task_multiply_CSR(spm, x, yOmpTaskWarmUp, initialJobBatchArrPerThread, initialJobBatchCountPerThread, numThreads, &execHistoryPerThread);
			);

			batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

			// Clean up dynamic scheduling data structures
			batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

			RUN(
				options.runCountMeasure,
				spmxv_static_rowWise_CSR(spm, x, yOmpTask, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
				ompTaskResults
			);

			// Clean up
			CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
			CLEAN_UP_JOB_BATCHES(initialJobBatchArrPerThread, initialJobBatchCountPerThread, numThreads);
		}

		// OMP victim any scatter
		// ------------------------------------------------------------------------------------------------
		{
			job_batch_t** initialJobBatchArrPerThread = NULL;
			int* initialJobBatchCountPerThread = NULL;
			execHistoryPerThread = NULL;

			spmxv_omp_task_prepFromPartitionVectorScatter_rowWise_CSR(
					spm, numBlocks, threadsPerBlock, rowPartitioning, colPartitioning, &initialJobBatchArrPerThread, &initialJobBatchCountPerThread);

			WARM_UP(
				options.runCountWarmUp,
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				spmxv_omp_task_multiply_CSR(spm, x, yOmpTaskScatterWarmUp, initialJobBatchArrPerThread, initialJobBatchCountPerThread, numThreads, &execHistoryPerThread);
			);

			batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

			// Clean up dynamic scheduling data structures
			batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

			RUN(
				options.runCountMeasure,
				spmxv_static_rowWise_CSR(spm, x, yOmpTaskScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
				ompTaskScatterResults
			);

			// Clean up
			CLEAN_UP_JOB_BATCHES(initialJobBatchArrPerThread, initialJobBatchCountPerThread, numThreads);
			CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
		}

		// Global work stealing algorithm
		// ------------------------------------------------------------------------------------------------
		{
			execHistoryPerThread = NULL;
			job_queue_t* globalQueue = NULL;
			spmxv_gws_prepFromPartitionVector_colNet_CSR(spm, rowPartitioning, options.stealTreshold, &globalQueue);

			WARM_UP(
				options.runCountWarmUp,
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				job_queue_reset(globalQueue);
				spmxv_gws_multiply_CSR(spm, x, yGwsWarmUp, globalQueue, numThreads, &execHistoryPerThread);
			);

			batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

			// clean up scheduling history and intermediate data structures
			job_queue_delete(globalQueue);
			batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

			RUN(
				options.runCountMeasure,
				spmxv_static_rowWise_CSR(spm, x, yGws, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
				gwsResults
			);

			// clean up
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
			block_info_t** threadInfos = NULL;
			spmxv_dws_prepFromPartitionVector_rowWise_colNet_CSR(spm, numThreads, rowPartitioning, &threadInfos);

			WARM_UP(
				options.runCountWarmUp,
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				spmxv_dws_rowWise_CSR(spm, x, yDwsRingWarmUp, threadInfos, numThreads, &execHistoryPerThread);
				block_reformMultiple(threadInfos, execHistoryPerThread, numThreads);
			);

			batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

			// clean up scheduling history and intermediate data structures
			spmxv_dws_cleanUp(threadInfos, numThreads);
			batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

			RUN(
				options.runCountMeasure,
				spmxv_static_rowWise_CSR(spm, x, yDwsRing, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
				dwsRingResults
			);

			// clean up
			CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
		}

		// Distributed Work Stealing Algorithm (Ring Scheduler) Scatter
		// ----------------------------------------------------------------------------------------------------------------------------
		{
			block_schedType = SCHEDULER_RING;
			dws_scheduler_init = ring_scheduler_init;
			dws_scheduler_localityAwareStealHalf = ring_scheduler_localityAwareStealHalf;
			dws_scheduler_terminate = ring_scheduler_terminate;

			execHistoryPerThread = NULL;
			block_info_t** threadInfos = NULL;
			spmxv_dws_prepFromPartitionVectorScatter_rowWise_colNet_CSR(spm, numBlocks, threadsPerBlock, rowPartitioning, &threadInfos);

			WARM_UP(
				options.runCountWarmUp,
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				spmxv_dws_rowWise_CSR(spm, x, yDwsRingScatterWarmUp, threadInfos, numThreads, &execHistoryPerThread);
				block_reformMultiple(threadInfos, execHistoryPerThread, numThreads);
			);

			batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

			// clean up scheduling history and intermediate data structures
			spmxv_dws_cleanUp(threadInfos, numThreads);
			batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

			RUN(
				options.runCountMeasure,
				spmxv_static_rowWise_CSR(spm, x, yDwsRingScatter, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
				dwsRingScatterResults
			);

			// clean up
			CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
		}

		// Distributed Work Stealing Algorithm (Ring Scheduler) Shared Block
		// ----------------------------------------------------------------------------------------------------------------------------
		{
			block_schedType = SCHEDULER_RING;
			dws_scheduler_init = ring_scheduler_init;
			dws_scheduler_localityAwareStealHalf = ring_scheduler_localityAwareStealHalf;
			dws_scheduler_terminate = ring_scheduler_terminate;

			execHistoryPerThread = NULL;
			execHistoryPerBlock = NULL;
			block_info_t** blockInfos = NULL;
			spmxv_dws_prepFromPartitionVectorScatterSharedBlock_colNet_CSR(spm, numBlocks, rowPartitioning, &blockInfos);

			WARM_UP(
				options.runCountWarmUp,
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				batch_list_deleteAllShallowMultiple(execHistoryPerBlock, numBlocks);
				spmxv_dws_sharedBlock_rowWise_CSR(spm, x, yDwsRingScatterSharedBlockWarmUp, blockInfos, numBlocks, threadsPerBlock, &execHistoryPerThread, &execHistoryPerBlock);
				block_reformMultiple(blockInfos, execHistoryPerBlock, numBlocks);
			);

			batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

			// Clean up dynamic scheduling data structures
			spmxv_dws_cleanUp(blockInfos, numBlocks);
			batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
			batch_list_deleteAllShallowMultiple(execHistoryPerBlock, numBlocks);

			RUN(
				options.runCountMeasure,
				spmxv_static_rowWise_CSR(spm, x, yDwsRingScatterSharedBlock, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
				dwsRingScatterSharedBlockResults
			);

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
			block_info_t** threadInfos = NULL;
			spmxv_dws_prepFromPartitionVector_rowWise_colNet_CSR(spm, numThreads, rowPartitioning, &threadInfos);

			WARM_UP(
				options.runCountWarmUp,
				batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);
				spmxv_dws_rowWise_CSR(spm, x, yDwsTreeWarmUp, threadInfos, numThreads, &execHistoryPerThread);
				block_reformMultiple(threadInfos, execHistoryPerThread, numThreads);
			);

			batch_list_convertToJobBatchDenseMultiple(execHistoryPerThread, numThreads, TRUE, &jobBatchArrPerThread, &jobBatchCountPerThread);

			// clean up scheduling history and intermediate data structures
			spmxv_dws_cleanUp(threadInfos, numThreads);
			batch_list_deleteAllShallowMultiple(execHistoryPerThread, numThreads);

			RUN(
				options.runCountMeasure,
				spmxv_static_rowWise_CSR(spm, x, yDwsTree, jobBatchArrPerThread, jobBatchCountPerThread, numThreads);,
				dwsTreeResults
			);

			// clean up
			CLEAN_UP_JOB_BATCHES(jobBatchArrPerThread, jobBatchCountPerThread, numThreads);
		}

		// Correctness Test
		// ------------------------------------------------------------------------------------------------

		vector_real_compare("OMP_LOOP_DYNAMIC_32", ySequentialOrdered, yOmpLoopDynamic32);
		vector_real_compare("OMP_LOOP_DYNAMIC_64", ySequentialOrdered, yOmpLoopDynamic64);
		vector_real_compare("OMP_LOOP_DYNAMIC_128", ySequentialOrdered, yOmpLoopDynamic128);
#ifdef __ICC
		vector_real_compare("MKL", ySequentialOrdered, yMkl);
#endif
		vector_real_compare("STATIC", ySequentialOrdered, yStatic);
		vector_real_compare("STATIC_SCATTER", ySequentialOrdered, yStaticScatter);
		vector_real_compare("OMP_TASK_WARMUP", ySequentialOrdered, yOmpTaskWarmUp);
		vector_real_compare("OMP_TASK", ySequentialOrdered, yOmpTask);
		vector_real_compare("OMP_TASK_SCATTER_WARMUP", ySequentialOrdered, yOmpTaskScatterWarmUp);
		vector_real_compare("OMP_TASK_SCATTER", ySequentialOrdered, yOmpTaskScatter);
		vector_real_compare("GWS_WARMUP", ySequentialOrdered, yGwsWarmUp);
		vector_real_compare("GWS", ySequentialOrdered, yGws);
		vector_real_compare("DWS_RING_WARMUP", ySequentialOrdered, yDwsRingWarmUp);
		vector_real_compare("DWS_RING", ySequentialOrdered, yDwsRing);
		vector_real_compare("DWS_RING_SCATTER_WARMUP", ySequentialOrdered, yDwsRingScatterWarmUp);
		vector_real_compare("DWS_RING_SCATTER", ySequentialOrdered, yDwsRingScatter);
		vector_real_compare("DWS_RING_SCATTER_SHARED_BLOCK_WARMUP", ySequentialOrdered, yDwsRingScatterSharedBlockWarmUp);
		vector_real_compare("DWS_RING_SCATTER_SHARED_BLOCK", ySequentialOrdered, yDwsRingScatterSharedBlock);
		vector_real_compare("DWS_TREE_WARMUP", ySequentialOrdered, yDwsTreeWarmUp);
		vector_real_compare("DWS_TREE", ySequentialOrdered, yDwsTree);

		// reorder results vector
		input_reorderVector(&ySequentialOrdered, rowOrderLookup);

		// Clean up
		// ------------------------------------------------------------------------------------------------

		vector_real_delete(yOmpLoopGuided);
		vector_real_delete(yOmpLoopDynamic32);
		vector_real_delete(yOmpLoopDynamic64);
		vector_real_delete(yOmpLoopDynamic128);
		vector_real_delete(yMkl);
		vector_real_delete(yStatic);
		vector_real_delete(yStaticScatter);
		vector_real_delete(yOmpTaskWarmUp);
		vector_real_delete(yOmpTask);
		vector_real_delete(yOmpTaskScatterWarmUp);
		vector_real_delete(yOmpTaskScatter);
		vector_real_delete(yGwsWarmUp);
		vector_real_delete(yGws);
		vector_real_delete(yDwsRingWarmUp);
		vector_real_delete(yDwsRing);
		vector_real_delete(yDwsRingScatterWarmUp);
		vector_real_delete(yDwsRingScatter);
		vector_real_delete(yDwsRingScatterSharedBlockWarmUp);
		vector_real_delete(yDwsRingScatterSharedBlock);
		vector_real_delete(yDwsTreeWarmUp);
		vector_real_delete(yDwsTree);

		vector_int_delete(rowOrderLookup);
		vector_int_delete(colOrderLookup);
		vector_int_delete(rowPartitioning);
		vector_int_delete(colPartitioning);

		spm_cmp_delete(spm);
	}

	vector_real_compare("ORDERED_UNORDERED_CMP", ySequential, ySequentialOrdered);

	strcat(results, sequentialResults);
	strcat(results, ompLoop32Results);
	strcat(results, ompLoop64Results);
	strcat(results, ompLoop128Results);
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

	PRINTF(";%s\n", algorithmHeaders);
	PRINTF(";%s%s\n", statisticsHeaders, orderingHeaders);
	PRINTF("%s;%d;%d;%d;%0.0f;%0.1f;%0.0f;%0.0f;%0.1f;%0.0f;%0.1f;%d;%d;%s\n",
			options.mmfName, rowCount, colCount, nnz, rowMinNNZ, rowAvgNNZ, rowMaxNNZ, columnMinNNZ, columnAvgNNZ, columnMaxNNZ,
			options.targetedCacheSizeKB, unorderedPartitionCount, orderedPartitionCount, results);
	PRINTF("\n\n");

	vector_real_delete(ySequential);
	vector_real_delete(ySequentialOrdered);

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
