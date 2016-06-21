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
#include "include/control_unit/static.h"
#include "include/control_unit/omp_loop.h"
#include "include/control_unit/dws.h"
#include "include/control_unit/cu_options.h"
#include "include/control_unit/fast_run.h"

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

#define ALLOC alloc_if(1)
#define FREE free_if(1)
#define RETAIN free_if(0)
#define REUSE alloc_if(0)
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


int main(int argsc, char* argsv[])
{
	// read command line parameters
	// ---------------------------------------------------------------

	cli_options_t options;
	cli_options_init(&options);
	cli_parseInput(argsc, argsv, &options);

	if(!cli_isValid(&options))
	{
		cli_printUsage();
		exit(EXIT_FAILURE);
	}

	// Parameters that will be passed on to MICs
	// -------------------------------------------------------------------------------------------------------

	int numBlocks = options.numBlocks;
	int threadsPerBlock = options.threadsPerBlock;
	int numThreads = options.numBlocks * options.threadsPerBlock;

	int storageFormat = options.storageFormat;
	int orderingType = options.orderingType;
	int partitionType = options.partitionType;
	int partitionMethod = options.partitionMethod;

	int simdLength = options.simdLength;
	float targetedCacheSizeKB = options.targetedCacheSizeKB;
	int stealTreshold = options.stealTreshold;

	int runCountWarmUp = options.runCountWarmUp;
	int runCountMeasure = options.runCountMeasure;

	// read matrix
	// -------------------------------------------------------------------------------------------------------
	quintet_t* quintets = NULL;
	mm_info_t* mmInfo = NULL;
	DECIMAL nnz = 0;
	DECIMAL rowCount = 0;
	DECIMAL colCount = 0;

	// read mmf file
	input_readQuintets(options.mmfPath, &nnz, &rowCount, &colCount, &quintets, &mmInfo);

	// Using self-implemented hybrid sort instead of "input_sortQuintetRowValue(quintets, nnz);" for speed
	algorithm_parallelMergeQuickSort(quintets, nnz, numThreads, quintet_cmpRowValue);

	// generate x vector randomly
	vector_real_t* xProc = vector_real_random(colCount);
	REAL* xData = xProc->data;
	DECIMAL xLength = xProc->length;

	// -------------------------------------------------------------------------------------------------------

	PRINTF("%s;", options.mmfName);

	// Output buffers for CSV format
	// ---------------------------------------------------------------

	char sequentialResults[512]; sequentialResults[0] = '\0';
#ifdef __ICC
	char mklResults[512]; mklResults[0] = '\0';
#endif
	char staticResults[512]; staticResults[0] = '\0';
	char staticScatterResults[512]; staticScatterResults[0] = '\0';
//	char ompTaskResults[512]; ompTaskResults[0] = '\0';
//	char ompTaskScatterResults[512]; ompTaskScatterResults[0] = '\0';
	char gwsResults[512]; gwsResults[0] = '\0';
	char ompDynamicResults[512]; ompDynamicResults[0] = '\0';
	char ompGuidedResults[512]; ompGuidedResults[0] = '\0';
	char dwsRingResults[512]; dwsRingResults[0] = '\0';
	char dwsRingScatterResults[512]; dwsRingScatterResults[0] = '\0';
	char dwsRingSharedBlockResults[512]; dwsRingSharedBlockResults[0] = '\0';
	char dwsTreeResults[512]; dwsTreeResults[0] = '\0';
	char dwsTreeScatterResults[512]; dwsTreeScatterResults[0] = '\0';
	char dwsTreeSharedBlockResults[512]; dwsTreeSharedBlockResults[0] = '\0';


	// Unordered runs
	// --------------------------------------------------------------------------------------------------

#ifdef ARCH_MIC_OFFLOAD

	#pragma offload target(mic:3) \
		inout(sequentialResults, \
				mklResults, \
				staticResults, \
				staticScatterResults, \
				ompDynamicResults, \
				ompGuidedResults, \
				gwsResults, \
				dwsRingResults, \
				dwsRingScatterResults, \
				dwsRingSharedBlockResults, \
				dwsTreeResults, \
				dwsTreeScatterResults, \
				dwsTreeSharedBlockResults) \
		in(quintets[0:nnz], nnz, rowCount, colCount, numBlocks, numThreads, threadsPerBlock, \
				storageFormat, orderingType, partitionType, partitionMethod, \
				simdLength, targetedCacheSizeKB, stealTreshold, \
				runCountWarmUp, runCountMeasure, \
				xData[0:xLength], xLength)
#endif
	{

		// Sparse matrix statistics that will be outputed
		// ---------------------------------------------------------------

		REAL rowMinNNZ = 0;
		REAL rowAvgNNZ = 0;
		REAL rowMaxNNZ = 0;
		REAL columnMinNNZ = 0;
		REAL columnAvgNNZ = 0;
		REAL columnMaxNNZ = 0;

		// ----------------------------------------------------------------------------------------------------

		fast_run_t* fr = fast_run_new();
		fast_run_init(fr,
				quintets, rowCount, colCount, nnz, numBlocks, threadsPerBlock,
				targetedCacheSizeKB, simdLength,
				storageFormat, ORDERING_TYPE_NONE, partitionType, partitionMethod,
				NULL, NULL, NULL, NULL);

		int unorderedPartitionCount = fr->subMtxCount;

		// Hold x vector (created by processor itself)
		// ----------------------------------------------------------------------------------------------------
		vector_real_t* x = (vector_real_t*) malloc(sizeof(vector_real_t));
		x->data = xData;
		x->length = xLength;

		{
			spm_cmp_t* spmCsr = fr->spmCsr;
			vector_int_t* yVectorRowOrderLookup = fr->yVectorRowOrderLookup;

			// -------------------------------------------------------------------------------------------------

			// extract statistics
			spm_cmp_extractNonZeroStatistics_CSR(
					spmCsr, &rowMinNNZ, &rowAvgNNZ, &rowMaxNNZ, &columnMinNNZ, &columnAvgNNZ, &columnMaxNNZ);

			// Print statistics
			PRINTF("%d;%d;%d;%0.0f;%0.1f;%0.0f;%0.0f;%0.1f;%0.1f;%0.1f;%d;",
					rowCount, colCount, nnz, rowMinNNZ, rowAvgNNZ, rowMaxNNZ,
					columnMinNNZ, columnAvgNNZ, columnMaxNNZ, options.targetedCacheSizeKB,
					unorderedPartitionCount);

			// -------------------------------------------------------------------------------------------------

			vector_real_t* ySequentialCSR = vector_real_new(rowCount);
			spmxv_sequential_CSR(spmCsr, x, ySequentialCSR);

			int i;

			// Sequential CSR Implementation
			// ----------------------------------------------------------------------------------------------------
			{
				vector_real_t* ySequential = vector_real_new(rowCount);

				WARM_UP(
					runCountWarmUp,
					spmxv_sequential_CSR(spmCsr, x, ySequential);
				);
				RUN(
					runCountMeasure,
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
					runCountWarmUp,
					spmxv_mkl(spmCsr, x, yMkl);
				);
				RUN(
					runCountMeasure,
					spmxv_mkl(spmCsr, x, yMkl);,
					mklResults
				);

				vector_real_compare("MKL", ySequentialCSR, yMkl);
				vector_real_delete(yMkl);
			}
			// ----------------------------------------------------------------------------------------------------
	#endif

			// OpenMP default scheduling
			// ----------------------------------------------------------------------------------------------------
			{
				omp_loop_args_t* args = omp_loop_args_new();
				omp_loop_args_init(args, fr, 0);

				args->prep(args);

				// OpenMP dynamic scheduling
				// ------------------------------------------------------------------------------------------------
				{
					vector_real_t* yOmpDynamic = vector_real_new(rowCount);

					int chunkSize;
					for(chunkSize = 1; chunkSize < 5; ++chunkSize)
					{
						args->chunkSize = chunkSize;

						WARM_UP(
							runCountWarmUp,
							args->spmxv_omp_loop_dynamic(args, x, yOmpDynamic);
						);
						RUN(
							runCountMeasure,
							args->spmxv_omp_loop_dynamic(args, x, yOmpDynamic);,
							ompDynamicResults
						);

						vector_real_reset(yOmpDynamic);
						args->spmxv_omp_loop_dynamic(args, x, yOmpDynamic);
						input_reorderVector(&yOmpDynamic, yVectorRowOrderLookup);
					}

					vector_real_compare("OMP_Dynamic", ySequentialCSR, yOmpDynamic);
					vector_real_delete(yOmpDynamic);
				}

				// OpenMP guided scheduling
				// ------------------------------------------------------------------------------------------------
				{
					vector_real_t* yOmpGuided = vector_real_new(rowCount);

					WARM_UP(
						runCountWarmUp,
						args->spmxv_omp_loop_guided(args, x, yOmpGuided);
					);
					RUN(
						runCountMeasure,
						args->spmxv_omp_loop_guided(args, x, yOmpGuided);,
						ompGuidedResults
					);

					vector_real_reset(yOmpGuided);
					args->spmxv_omp_loop_guided(args, x, yOmpGuided);
					input_reorderVector(&yOmpGuided, yVectorRowOrderLookup);

					vector_real_compare("OMP_Guided", ySequentialCSR, yOmpGuided);
					vector_real_delete(yOmpGuided);
				}

				omp_loop_args_delete(args);
			}

			// Static Implementation
			// ------------------------------------------------------------------------------------------------
			{
				static_args_t* args = static_args_new();
				static_args_init(args, fr);

				args->prep(args);

				vector_real_t* yStatic = vector_real_new(rowCount);
				WARM_UP(
					runCountWarmUp,
					args->spmxv_static(args, x, yStatic);
				);

				RUN(
					runCountMeasure,
					args->spmxv_static(args, x, yStatic);,
					staticResults
				);

				vector_real_reset(yStatic);
				args->spmxv_static(args, x, yStatic);
				input_reorderVector(&yStatic, yVectorRowOrderLookup);
				vector_real_compare("STATIC", ySequentialCSR, yStatic);

				// Clean up
				static_args_delete(args);
				vector_real_delete(yStatic);
			}

			// Static Scatter Implementation
			// ------------------------------------------------------------------------------------------------
			{
				static_args_t* args = static_args_new();
				args->algorithmType = STATIC_ALGORITHM_TYPE_SCATTER;
				static_args_init(args, fr);

				args->prep(args);

				vector_real_t* yStatic = vector_real_new(rowCount);
				WARM_UP(
					runCountWarmUp,
					args->spmxv_static(args, x, yStatic);
				);

				RUN(
					runCountMeasure,
					args->spmxv_static(args, x, yStatic);,
					staticScatterResults
				);

				vector_real_reset(yStatic);
				args->spmxv_static(args, x, yStatic);
				input_reorderVector(&yStatic, yVectorRowOrderLookup);
				vector_real_compare("STATIC-SCATTER", ySequentialCSR, yStatic);

				// Clean up
				static_args_delete(args);
				vector_real_delete(yStatic);
			}

			// Global work stealing algorithm
			// ------------------------------------------------------------------------------------------------
			{
				gws_args_t* args = gws_args_new();
				gws_args_init(args, fr);

				args->prep(args);

				vector_real_t* yGwsWarmUp = vector_real_new(rowCount);
				WARM_UP(
					runCountWarmUp,
					gws_args_warmup(args, x, yGwsWarmUp);
				);

				vector_real_reset(yGwsWarmUp);
				gws_args_warmup(args, x, yGwsWarmUp);
				input_reorderVector(&yGwsWarmUp, yVectorRowOrderLookup);
				vector_real_compare("GWS_WARMUP", ySequentialCSR, yGwsWarmUp);

				gws_args_setup(args);

				vector_real_delete(yGwsWarmUp);

				vector_real_t* yGws = vector_real_new(rowCount);
				RUN(
					runCountMeasure,
					args->spmxv_gws(&args->staticArgs, x, yGws);,
					gwsResults
				);

				vector_real_reset(yGws);
				args->spmxv_gws(&args->staticArgs, x, yGws);
				input_reorderVector(&yGws, yVectorRowOrderLookup);

				vector_real_compare("GWS", ySequentialCSR, yGws);
				vector_real_delete(yGws);
				gws_args_delete(args);
			}

			// Distributed Work Stealing Algorithm (Ring Scheduler)
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				dws_args_t* args = dws_args_new();
				dws_args_init(args, fr, DWS_ALGORITHM_TYPE_DEFAULT, STEALING_SCHEME_RING, stealTreshold);

				args->prep(args);

				vector_real_t* yDwsRingWarmUp = vector_real_new(rowCount);
				WARM_UP(
					runCountWarmUp,
					dws_args_warmup(args, x, yDwsRingWarmUp);
				);

				vector_real_reset(yDwsRingWarmUp);
				dws_args_reset(args);
				dws_args_warmup(args, x, yDwsRingWarmUp);
				input_reorderVector(&yDwsRingWarmUp, yVectorRowOrderLookup);
				vector_real_compare("DWS_RING_WARMUP", ySequentialCSR, yDwsRingWarmUp);

				dws_args_setup(args);

				vector_real_delete(yDwsRingWarmUp);

				vector_real_t* yDwsRing = vector_real_new(rowCount);
				RUN(
					runCountMeasure,
					args->spmxv_dws(&args->staticArgs, x, yDwsRing);,
					dwsRingResults
				);


				vector_real_reset(yDwsRing);
				args->spmxv_dws(&args->staticArgs, x, yDwsRing);
				input_reorderVector(&yDwsRing, yVectorRowOrderLookup);
				vector_real_compare("DWS_RING", ySequentialCSR, yDwsRing);

				// Clean up
				dws_args_delete(args);
				vector_real_delete(yDwsRing);
			}

			// Distributed Work Stealing Scatter Algorithm (Ring Scheduler)
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				dws_args_t* args = dws_args_new();
				dws_args_init(args, fr, DWS_ALGORITHM_TYPE_SCATTER, STEALING_SCHEME_RING, stealTreshold);

				args->prep(args);

				vector_real_t* yDwsRingWarmUp = vector_real_new(rowCount);
				WARM_UP(
					runCountWarmUp,
					dws_args_warmup(args, x, yDwsRingWarmUp);
				);

				vector_real_reset(yDwsRingWarmUp);
				dws_args_reset(args);
				dws_args_warmup(args, x, yDwsRingWarmUp);
				input_reorderVector(&yDwsRingWarmUp, yVectorRowOrderLookup);
				vector_real_compare("DWS_RING_SCATTER_WARMUP", ySequentialCSR, yDwsRingWarmUp);

				dws_args_setup(args);

				vector_real_delete(yDwsRingWarmUp);

				vector_real_t* yDwsRing = vector_real_new(rowCount);
				RUN(
					runCountMeasure,
					args->spmxv_dws(&args->staticArgs, x, yDwsRing);,
					dwsRingResults
				);

				vector_real_reset(yDwsRing);
				args->spmxv_dws(&args->staticArgs, x, yDwsRing);
				input_reorderVector(&yDwsRing, yVectorRowOrderLookup);
				vector_real_compare("DWS_RING_SCATTER", ySequentialCSR, yDwsRing);

				// Clean up
				dws_args_delete(args);
				vector_real_delete(yDwsRing);
			}

			// Distributed Work Stealing Algorithm (Ring Scheduler) Shared Queue
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				dws_args_t* args = dws_args_new();
				dws_args_init(args, fr, DWS_ALGORITHM_TYPE_SHARED_BLOCK, STEALING_SCHEME_RING, stealTreshold);

				args->prep(args);

				vector_real_t* yDwsRingScatterSharedBlockWarmUp = vector_real_new(rowCount);
				WARM_UP(
					runCountWarmUp,
					dws_args_warmup(args, x, yDwsRingScatterSharedBlockWarmUp);
				);

				dws_args_setup(args);

				vector_real_reset(yDwsRingScatterSharedBlockWarmUp);
				dws_args_reset(args);
				dws_args_warmup(args, x, yDwsRingScatterSharedBlockWarmUp);
				input_reorderVector(&yDwsRingScatterSharedBlockWarmUp, yVectorRowOrderLookup);
				vector_real_compare("DWS_RING_SHARED_BLOCK_WARMUP", ySequentialCSR, yDwsRingScatterSharedBlockWarmUp);

				vector_real_delete(yDwsRingScatterSharedBlockWarmUp);

				vector_real_t* yDwsRingScatterSharedBlock = vector_real_new(rowCount);
				RUN(
					runCountMeasure,
					args->spmxv_dws(&args->staticArgs, x, yDwsRingScatterSharedBlock);,
					dwsRingSharedBlockResults
				);

				vector_real_reset(yDwsRingScatterSharedBlock);
				args->spmxv_dws(&args->staticArgs, x, yDwsRingScatterSharedBlock);
				input_reorderVector(&yDwsRingScatterSharedBlock, yVectorRowOrderLookup);
				vector_real_compare("DWS_RING_SHARED_BLOCK", ySequentialCSR, yDwsRingScatterSharedBlock);

				// Clean up
				dws_args_delete(args);
				vector_real_delete(yDwsRingScatterSharedBlock);
			}

			// Distributed Work Stealing Algorithm (Tree Scheduler)
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				dws_args_t* args = dws_args_new();
				dws_args_init(args, fr, DWS_ALGORITHM_TYPE_DEFAULT, STEALING_SCHEME_TREE, stealTreshold);

				args->prep(args);

				vector_real_t* yDwsRingWarmUp = vector_real_new(rowCount);
				WARM_UP(
					runCountWarmUp,
					dws_args_warmup(args, x, yDwsRingWarmUp);
				);

				vector_real_reset(yDwsRingWarmUp);
				dws_args_reset(args);
				dws_args_warmup(args, x, yDwsRingWarmUp);
				input_reorderVector(&yDwsRingWarmUp, yVectorRowOrderLookup);
				vector_real_compare("DWS_TREE_WARMUP", ySequentialCSR, yDwsRingWarmUp);

				dws_args_setup(args);

				vector_real_delete(yDwsRingWarmUp);

				vector_real_t* yDwsRing = vector_real_new(rowCount);
				RUN(
					runCountMeasure,
					args->spmxv_dws(&args->staticArgs, x, yDwsRing);,
					dwsRingResults
				);


				vector_real_reset(yDwsRing);
				args->spmxv_dws(&args->staticArgs, x, yDwsRing);
				input_reorderVector(&yDwsRing, yVectorRowOrderLookup);
				vector_real_compare("DWS_TREE", ySequentialCSR, yDwsRing);

				// Clean up
				dws_args_delete(args);
				vector_real_delete(yDwsRing);
			}

			// Distributed Work Stealing Algorithm (Tree Scheduler) Scatter
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				dws_args_t* args = dws_args_new();
				dws_args_init(args, fr, DWS_ALGORITHM_TYPE_SCATTER, STEALING_SCHEME_TREE, stealTreshold);

				args->prep(args);

				vector_real_t* yDwsRingWarmUp = vector_real_new(rowCount);
				WARM_UP(
					runCountWarmUp,
					dws_args_warmup(args, x, yDwsRingWarmUp);
				);

				vector_real_reset(yDwsRingWarmUp);
				dws_args_reset(args);
				dws_args_warmup(args, x, yDwsRingWarmUp);
				input_reorderVector(&yDwsRingWarmUp, yVectorRowOrderLookup);
				vector_real_compare("DWS_TREE_SCATTER_WARMUP", ySequentialCSR, yDwsRingWarmUp);

				dws_args_setup(args);

				vector_real_delete(yDwsRingWarmUp);

				vector_real_t* yDwsRing = vector_real_new(rowCount);
				RUN(
					runCountMeasure,
					args->spmxv_dws(&args->staticArgs, x, yDwsRing);,
					dwsRingResults
				);

				vector_real_reset(yDwsRing);
				args->spmxv_dws(&args->staticArgs, x, yDwsRing);
				input_reorderVector(&yDwsRing, yVectorRowOrderLookup);
				vector_real_compare("DWS_TREE_SCATTER", ySequentialCSR, yDwsRing);

				// Clean up
				dws_args_delete(args);
				vector_real_delete(yDwsRing);
			}

			// Distributed Work Stealing Algorithm (Tree Scheduler) Shared Queue
			// ----------------------------------------------------------------------------------------------------------------------------
			{
				dws_args_t* args = dws_args_new();
				dws_args_init(args, fr, DWS_ALGORITHM_TYPE_SHARED_BLOCK, STEALING_SCHEME_TREE, stealTreshold);

				args->prep(args);

				vector_real_t* yDwsRingScatterSharedBlockWarmUp = vector_real_new(rowCount);
				WARM_UP(
					runCountWarmUp,
					dws_args_warmup(args, x, yDwsRingScatterSharedBlockWarmUp);
				);

				dws_args_setup(args);

				vector_real_reset(yDwsRingScatterSharedBlockWarmUp);
				dws_args_reset(args);
				dws_args_warmup(args, x, yDwsRingScatterSharedBlockWarmUp);
				input_reorderVector(&yDwsRingScatterSharedBlockWarmUp, yVectorRowOrderLookup);
				vector_real_compare("DWS_TREE_SHARED_BLOCK_WARMUP", ySequentialCSR, yDwsRingScatterSharedBlockWarmUp);

				vector_real_delete(yDwsRingScatterSharedBlockWarmUp);

				vector_real_t* yDwsRingScatterSharedBlock = vector_real_new(rowCount);
				RUN(
					runCountMeasure,
					args->spmxv_dws(&args->staticArgs, x, yDwsRingScatterSharedBlock);,
					dwsRingSharedBlockResults
				);

				vector_real_reset(yDwsRingScatterSharedBlock);
				args->spmxv_dws(&args->staticArgs, x, yDwsRingScatterSharedBlock);
				input_reorderVector(&yDwsRingScatterSharedBlock, yVectorRowOrderLookup);
				vector_real_compare("DWS_TREE_SHARED_BLOCK", ySequentialCSR, yDwsRingScatterSharedBlock);

				// Clean up
				dws_args_delete(args);
				vector_real_delete(yDwsRingScatterSharedBlock);
			}

			vector_real_delete(ySequentialCSR);
			free(x);
		}

		fast_run_delete(fr);
	}

	// Ordered Runs
	// ----------------------------------------------------------------------------------------------------

	// Read ordering and partitioning info
	// ----------------------------------------------------------------------------------------------------

	vector_int_t* rowOrderLookupProc = NULL;
	vector_int_t* rowPartitioningProc = NULL;
	vector_int_t* colOrderLookupProc = NULL;
	vector_int_t* colPartitioningProc = NULL;

	input_powerOf2Order(
			options.mmfDir, options.mmfName, targetedCacheSizeKB,
			quintets, nnz, rowCount, colCount,
			&rowOrderLookupProc, &rowPartitioningProc, &colOrderLookupProc, &colPartitioningProc);

	input_orderVector(&xProc, colOrderLookupProc);
	xData = xProc->data;
	xLength = xProc->length;

	DECIMAL* rowOrderLookupData = rowOrderLookupProc->data;
	DECIMAL* rowPartitioningData = rowPartitioningProc->data;
//	TODO fix these
//	DECIMAL* colOrderLookupData = colOrderLookupProc->data;
//	DECIMAL* colPartitioningData = colPartitioningProc->data;

	int rowOrderLookupLength = rowOrderLookupProc->length;
	int rowPartitioningLength = rowPartitioningProc->length;
//	TODO fix these
//	int colOrderLookupLength = colOrderLookupProc->length;
//	int colPartitioningLength = colPartitioningProc->length;


	// TODO fix this (calculate orderedPartitionCount in an oblivious manner, probably in fast_run_t)
	int orderedPartitionCount = rowPartitioningProc->length - 1;
	PRINTF("%d;", orderedPartitionCount);

#ifdef ARCH_MIC_OFFLOAD

	#pragma offload target(mic:3) \
		inout(sequentialResults, \
				mklResults, \
				staticResults, \
				staticScatterResults, \
				ompDynamicResults, \
				ompGuidedResults, \
				gwsResults, \
				dwsRingResults, \
				dwsRingScatterResults, \
				dwsRingSharedBlockResults, \
				dwsTreeResults, \
				dwsTreeScatterResults, \
				dwsTreeSharedBlockResults) \
		in(quintets[0:nnz], nnz, rowCount, colCount, numBlocks, numThreads, threadsPerBlock, \
				storageFormat, orderingType, partitionType, partitionMethod, \
				simdLength, targetedCacheSizeKB, stealTreshold, \
				runCountWarmUp, runCountMeasure, \
				xData[0:xLength], xLength, \
				rowPartitioningData[0:rowPartitioningLength], rowPartitioningLength)
		// colPartitioningData[0:colPartitioningLength], colPartitioningLength)
#endif
	{
		// pack partitioning info sent from processor itself
		// -------------------------------------------------------------------------------------------------------

		vector_real_t* x = (vector_real_t*) malloc(sizeof(vector_real_t));
		x->data = xData;
		x->length = xLength;

		vector_int_t* rowPartitioning = (vector_int_t*) malloc(sizeof(vector_int_t));
		rowPartitioning->data = rowPartitioningData;
		rowPartitioning->length = rowPartitioningLength;

		// -------------------------------------------------------------------------------------------------------

		fast_run_t* fr = fast_run_new();
		fast_run_init(fr,
				quintets, rowCount, colCount, nnz, numBlocks, threadsPerBlock,
				targetedCacheSizeKB, simdLength,
				storageFormat, orderingType, partitionType, partitionMethod,
				NULL, rowPartitioning, NULL, NULL);

		spm_cmp_t* spmCsr = fr->spmCsr;
		vector_int_t* yVectorRowOrderLookup = fr->yVectorRowOrderLookup;

		vector_real_t* ySequentialCSR = vector_real_new(rowCount);
		spmxv_sequential_CSR(spmCsr, x, ySequentialCSR);

		int i;

		// Sequential CSR Implementation
		// ----------------------------------------------------------------------------------------------------
		{
			vector_real_t* ySequential = vector_real_new(rowCount);

			WARM_UP(
				runCountWarmUp,
				spmxv_sequential_CSR(spmCsr, x, ySequential);
			);
			RUN(
				runCountMeasure,
				spmxv_sequential_CSR(spmCsr, x, ySequential);,
				sequentialResults
			);


			vector_real_compare("SEQUENTIAL_CSR_ORD", ySequentialCSR, ySequential);
			vector_real_delete(ySequential);
		}

#ifdef __ICC
		// MKL Implementation
		// ----------------------------------------------------------------------------------------------------
		{
			vector_real_t* yMkl = vector_real_new(rowCount);

			WARM_UP(
				runCountWarmUp,
				spmxv_mkl(spmCsr, x, yMkl);
			);
			RUN(
				runCountMeasure,
				spmxv_mkl(spmCsr, x, yMkl);,
				mklResults
			);

			vector_real_compare("MKL_ORD", ySequentialCSR, yMkl);
			vector_real_delete(yMkl);
		}
		// ----------------------------------------------------------------------------------------------------
#endif

		// OpenMP default scheduling
		// ----------------------------------------------------------------------------------------------------
		{
			omp_loop_args_t* args = omp_loop_args_new();
			omp_loop_args_init(args, fr, 0);

			args->prep(args);

			// OpenMP dynamic scheduling
			// ------------------------------------------------------------------------------------------------
			{
				vector_real_t* yOmpDynamic = vector_real_new(rowCount);

				int chunkSize;
				for(chunkSize = 1; chunkSize < 5; ++chunkSize)
				{
					args->chunkSize = chunkSize;

					WARM_UP(
						runCountWarmUp,
						args->spmxv_omp_loop_dynamic(args, x, yOmpDynamic);
					);
					RUN(
						runCountMeasure,
						args->spmxv_omp_loop_dynamic(args, x, yOmpDynamic);,
						ompDynamicResults
					);

					vector_real_reset(yOmpDynamic);
					args->spmxv_omp_loop_dynamic(args, x, yOmpDynamic);
				}

				vector_real_compare("OMP_Dynamic_ORD", ySequentialCSR, yOmpDynamic);
				vector_real_delete(yOmpDynamic);
			}

			// OpenMP guided scheduling
			// ------------------------------------------------------------------------------------------------
			{
				vector_real_t* yOmpGuided = vector_real_new(rowCount);

				WARM_UP(
					runCountWarmUp,
					args->spmxv_omp_loop_guided(args, x, yOmpGuided);
				);
				RUN(
					runCountMeasure,
					args->spmxv_omp_loop_guided(args, x, yOmpGuided);,
					ompGuidedResults
				);

				vector_real_reset(yOmpGuided);
				args->spmxv_omp_loop_guided(args, x, yOmpGuided);

				vector_real_compare("OMP_Guided_ORD", ySequentialCSR, yOmpGuided);
				vector_real_delete(yOmpGuided);
			}

			omp_loop_args_delete(args);
		}

		// Static Implementation
		// ------------------------------------------------------------------------------------------------
		{
			static_args_t* args = static_args_new();
			static_args_init(args, fr);

			args->prep(args);

			vector_real_t* yStatic = vector_real_new(rowCount);
			WARM_UP(
				runCountWarmUp,
				args->spmxv_static(args, x, yStatic);
			);

			RUN(
				runCountMeasure,
				args->spmxv_static(args, x, yStatic);,
				staticResults
			);

			vector_real_reset(yStatic);
			args->spmxv_static(args, x, yStatic);
			vector_real_compare("STATIC_ORD", ySequentialCSR, yStatic);

			// Clean up
			static_args_delete(args);
			vector_real_delete(yStatic);
		}

		// Static Scatter Implementation
		// ------------------------------------------------------------------------------------------------
		{
			static_args_t* args = static_args_new();
			args->algorithmType = STATIC_ALGORITHM_TYPE_SCATTER;
			static_args_init(args, fr);

			args->prep(args);

			vector_real_t* yStatic = vector_real_new(rowCount);
			WARM_UP(
				runCountWarmUp,
				args->spmxv_static(args, x, yStatic);
			);

			RUN(
				runCountMeasure,
				args->spmxv_static(args, x, yStatic);,
				staticResults
			);

			vector_real_reset(yStatic);
			args->spmxv_static(args, x, yStatic);
			vector_real_compare("STATIC_SCATTER_ORD", ySequentialCSR, yStatic);

			// Clean up
			static_args_delete(args);
			vector_real_delete(yStatic);
		}

		// Global work stealing algorithm
		// ------------------------------------------------------------------------------------------------
		{
			gws_args_t* args = gws_args_new();
			gws_args_init(args, fr);

			args->prep(args);

			vector_real_t* yGwsWarmUp = vector_real_new(rowCount);
			WARM_UP(
				runCountWarmUp,
				gws_args_warmup(args, x, yGwsWarmUp);
			);

			vector_real_reset(yGwsWarmUp);
			gws_args_warmup(args, x, yGwsWarmUp);
			vector_real_compare("GWS_WARMUP_ORD", ySequentialCSR, yGwsWarmUp);

			gws_args_setup(args);

			vector_real_delete(yGwsWarmUp);

			vector_real_t* yGws = vector_real_new(rowCount);
			RUN(
				runCountMeasure,
				args->spmxv_gws(&args->staticArgs, x, yGws);,
				gwsResults
			);

			vector_real_reset(yGws);
			args->spmxv_gws(&args->staticArgs, x, yGws);

			vector_real_compare("GWS_ORD", ySequentialCSR, yGws);
			vector_real_delete(yGws);
			gws_args_delete(args);
		}

		// Distributed Work Stealing Algorithm (Ring Scheduler)
		// ----------------------------------------------------------------------------------------------------------------------------
		{
			dws_args_t* args = dws_args_new();
			dws_args_init(args, fr, DWS_ALGORITHM_TYPE_DEFAULT, STEALING_SCHEME_RING, stealTreshold);

			args->prep(args);

			vector_real_t* yDwsRingWarmUp = vector_real_new(rowCount);
			WARM_UP(
				runCountWarmUp,
				dws_args_warmup(args, x, yDwsRingWarmUp);
			);

			vector_real_reset(yDwsRingWarmUp);
			dws_args_reset(args);
			dws_args_warmup(args, x, yDwsRingWarmUp);
			vector_real_compare("DWS_RING_WARMUP_ORD", ySequentialCSR, yDwsRingWarmUp);

			dws_args_setup(args);

			vector_real_delete(yDwsRingWarmUp);

			vector_real_t* yDwsRing = vector_real_new(rowCount);
			RUN(
				runCountMeasure,
				args->spmxv_dws(&args->staticArgs, x, yDwsRing);,
				dwsRingResults
			);


			vector_real_reset(yDwsRing);
			args->spmxv_dws(&args->staticArgs, x, yDwsRing);
			vector_real_compare("DWS_RING_ORD", ySequentialCSR, yDwsRing);

			// Clean up
			dws_args_delete(args);
			vector_real_delete(yDwsRing);
		}

		// Distributed Work Stealing Algorithm (Ring Scheduler) Scatter
		// ----------------------------------------------------------------------------------------------------------------------------
		{
			dws_args_t* args = dws_args_new();
			dws_args_init(args, fr, DWS_ALGORITHM_TYPE_SCATTER, STEALING_SCHEME_RING, stealTreshold);

			args->prep(args);

			vector_real_t* yDwsRingWarmUp = vector_real_new(rowCount);
			WARM_UP(
				runCountWarmUp,
				dws_args_warmup(args, x, yDwsRingWarmUp);
			);

			vector_real_reset(yDwsRingWarmUp);
			dws_args_reset(args);
			dws_args_warmup(args, x, yDwsRingWarmUp);
			vector_real_compare("DWS_RING_SCATTER_WARMUP_ORD", ySequentialCSR, yDwsRingWarmUp);

			dws_args_setup(args);

			vector_real_delete(yDwsRingWarmUp);

			vector_real_t* yDwsRing = vector_real_new(rowCount);
			RUN(
				runCountMeasure,
				args->spmxv_dws(&args->staticArgs, x, yDwsRing);,
				dwsRingResults
			);

			vector_real_reset(yDwsRing);
			args->spmxv_dws(&args->staticArgs, x, yDwsRing);
			vector_real_compare("DWS_RING_SCATTER_ORD", ySequentialCSR, yDwsRing);

			// Clean up
			dws_args_delete(args);
			vector_real_delete(yDwsRing);
		}

		// Distributed Work Stealing Algorithm (Ring Scheduler) Shared Queue
		// ----------------------------------------------------------------------------------------------------------------------------
		{
			dws_args_t* args = dws_args_new();
			dws_args_init(args, fr, DWS_ALGORITHM_TYPE_SHARED_BLOCK, STEALING_SCHEME_RING, stealTreshold);

			args->prep(args);

			vector_real_t* yDwsRingScatterSharedBlockWarmUp = vector_real_new(rowCount);
			WARM_UP(
				runCountWarmUp,
				dws_args_warmup(args, x, yDwsRingScatterSharedBlockWarmUp);
			);

			dws_args_setup(args);

			vector_real_reset(yDwsRingScatterSharedBlockWarmUp);
			dws_args_reset(args);
			dws_args_warmup(args, x, yDwsRingScatterSharedBlockWarmUp);
			vector_real_compare("DWS_RING_SHARED_BLOCK_WARMUP_ORD", ySequentialCSR, yDwsRingScatterSharedBlockWarmUp);

			vector_real_delete(yDwsRingScatterSharedBlockWarmUp);

			vector_real_t* yDwsRingScatterSharedBlock = vector_real_new(rowCount);
			RUN(
				runCountMeasure,
				args->spmxv_dws(&args->staticArgs, x, yDwsRingScatterSharedBlock);,
				dwsRingSharedBlockResults
			);

			vector_real_reset(yDwsRingScatterSharedBlock);
			args->spmxv_dws(&args->staticArgs, x, yDwsRingScatterSharedBlock);
			vector_real_compare("DWS_RING_SHARED_BLOCK_ORD", ySequentialCSR, yDwsRingScatterSharedBlock);

			// Clean up
			dws_args_delete(args);
			vector_real_delete(yDwsRingScatterSharedBlock);
		}

		// Distributed Work Stealing Algorithm (Tree Scheduler)
		// ----------------------------------------------------------------------------------------------------------------------------
		{
			dws_args_t* args = dws_args_new();
			dws_args_init(args, fr, DWS_ALGORITHM_TYPE_DEFAULT, STEALING_SCHEME_TREE, stealTreshold);

			args->prep(args);

			vector_real_t* yDwsRingWarmUp = vector_real_new(rowCount);
			WARM_UP(
				runCountWarmUp,
				dws_args_warmup(args, x, yDwsRingWarmUp);
			);

			vector_real_reset(yDwsRingWarmUp);
			dws_args_reset(args);
			dws_args_warmup(args, x, yDwsRingWarmUp);
			vector_real_compare("DWS_TREE_WARMUP_ORD", ySequentialCSR, yDwsRingWarmUp);

			dws_args_setup(args);

			vector_real_delete(yDwsRingWarmUp);

			vector_real_t* yDwsRing = vector_real_new(rowCount);
			RUN(
				runCountMeasure,
				args->spmxv_dws(&args->staticArgs, x, yDwsRing);,
				dwsRingResults
			);


			vector_real_reset(yDwsRing);
			args->spmxv_dws(&args->staticArgs, x, yDwsRing);
			vector_real_compare("DWS_TREE_ORD", ySequentialCSR, yDwsRing);

			// Clean up
			dws_args_delete(args);
			vector_real_delete(yDwsRing);
		}

		// Distributed Work Stealing Algorithm (Tree Scheduler) Scatter
		// ----------------------------------------------------------------------------------------------------------------------------
		{
			dws_args_t* args = dws_args_new();
			dws_args_init(args, fr, DWS_ALGORITHM_TYPE_SCATTER, STEALING_SCHEME_TREE, stealTreshold);

			args->prep(args);

			vector_real_t* yDwsRingWarmUp = vector_real_new(rowCount);
			WARM_UP(
				runCountWarmUp,
				dws_args_warmup(args, x, yDwsRingWarmUp);
			);

			vector_real_reset(yDwsRingWarmUp);
			dws_args_reset(args);
			dws_args_warmup(args, x, yDwsRingWarmUp);
			vector_real_compare("DWS_TREE_SCATTER_WARMUP_ORD", ySequentialCSR, yDwsRingWarmUp);

			dws_args_setup(args);

			vector_real_delete(yDwsRingWarmUp);

			vector_real_t* yDwsRing = vector_real_new(rowCount);
			RUN(
				runCountMeasure,
				args->spmxv_dws(&args->staticArgs, x, yDwsRing);,
				dwsRingResults
			);

			vector_real_reset(yDwsRing);
			args->spmxv_dws(&args->staticArgs, x, yDwsRing);
			vector_real_compare("DWS_TREE_SCATTER_ORD", ySequentialCSR, yDwsRing);

			// Clean up
			dws_args_delete(args);
			vector_real_delete(yDwsRing);
		}


		// Distributed Work Stealing Algorithm (Tree Scheduler) Shared Queue
		// ----------------------------------------------------------------------------------------------------------------------------
		{
			dws_args_t* args = dws_args_new();
			dws_args_init(args, fr, DWS_ALGORITHM_TYPE_SHARED_BLOCK, STEALING_SCHEME_TREE, stealTreshold);

			args->prep(args);

			vector_real_t* yDwsRingScatterSharedBlockWarmUp = vector_real_new(rowCount);
			WARM_UP(
				runCountWarmUp,
				dws_args_warmup(args, x, yDwsRingScatterSharedBlockWarmUp);
			);

			dws_args_setup(args);

			vector_real_reset(yDwsRingScatterSharedBlockWarmUp);
			dws_args_reset(args);
			dws_args_warmup(args, x, yDwsRingScatterSharedBlockWarmUp);
			vector_real_compare("DWS_TREE_SHARED_BLOCK_WARMUP_ORD", ySequentialCSR, yDwsRingScatterSharedBlockWarmUp);

			vector_real_delete(yDwsRingScatterSharedBlockWarmUp);

			vector_real_t* yDwsRingScatterSharedBlock = vector_real_new(rowCount);
			RUN(
				runCountMeasure,
				args->spmxv_dws(&args->staticArgs, x, yDwsRingScatterSharedBlock);,
				dwsRingSharedBlockResults
			);

			vector_real_reset(yDwsRingScatterSharedBlock);
			args->spmxv_dws(&args->staticArgs, x, yDwsRingScatterSharedBlock);
			vector_real_compare("DWS_TREE_SHARED_BLOCK_ORD", ySequentialCSR, yDwsRingScatterSharedBlock);

			// Clean up
			dws_args_delete(args);
			vector_real_delete(yDwsRingScatterSharedBlock);
		}

		// reorder results vector
		// TODO check? how?

		// Print results
		// ------------------------------------------------------------------------------------------------

		PRINTF("%s", sequentialResults);
#ifdef __ICC
		PRINTF("%s", mklResults);
#endif
		PRINTF("%s", staticResults);
		PRINTF("%s", staticScatterResults);
		PRINTF("%s", gwsResults);
		PRINTF("%s", ompDynamicResults);
		PRINTF("%s", ompGuidedResults);
		PRINTF("%s", dwsRingResults);
		PRINTF("%s", dwsRingScatterResults);
		PRINTF("%s", dwsRingSharedBlockResults);
		PRINTF("%s", dwsTreeResults);
		PRINTF("%s", dwsTreeScatterResults);
		PRINTF("%s", dwsTreeSharedBlockResults);
		PRINTF("\n");

		// Clean up
		// ------------------------------------------------------------------------------------------------

		vector_real_delete(ySequentialCSR);
		free(x);
		fast_run_delete(fr);
	}

	vector_int_delete(rowOrderLookupProc);
	vector_int_delete(colOrderLookupProc);
	vector_int_delete(rowPartitioningProc);
	vector_int_delete(colPartitioningProc);


	free(quintets);
	cli_options_deleteNonPtr(&options);
	vector_real_delete(xProc);
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
