
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "include/config.h"
#include "include/algorithm/algorithm.h"
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

extern int isSameTriplet(quintet_t* t1, quintet_t* t2);
extern int compareTriplets(quintet_t* triplets1, quintet_t* triplets2, DECIMAL nnz);

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
	quintet_t* tripletsQSort = NULL;
	quintet_t* tripletsCustomSort = NULL;
	mm_info_t* mmInfo = NULL;
	mm_info_t* mmInfo2 = NULL;
	DECIMAL nnz = 0;
	DECIMAL rowCount = 0;
	DECIMAL colCount = 0;

	input_readQuintets(options.mmfPath, &nnz, &rowCount, &colCount, &tripletsQSort, &mmInfo);
	input_readQuintets(options.mmfPath, &nnz, &rowCount, &colCount, &tripletsCustomSort, &mmInfo2);

	int i;

	double timeElapsedQSort = 0;
	MEASURE(1,
			quintet_sortRowValue(tripletsQSort, nnz);,
			timeElapsedQSort
	);

	double timeElapsedCustomSort = 0;
	MEASURE(1,
			algorithm_parallelMergeQuickSort(tripletsCustomSort, nnz, numThreads, quintet_cmpRowValue);,
			timeElapsedCustomSort
	);
	printf("nnz: %d\n", nnz);

	printf("Are results same ?? : %d\n", compareTriplets(tripletsQSort, tripletsCustomSort, nnz));

	printf("Time elapsed for qsort: %lf\n", timeElapsedQSort);
	printf("Time elapsed for custom-sort: %lf\n", timeElapsedCustomSort);

	// clean up
	free(tripletsQSort);
	free(tripletsCustomSort);
	mm_info_delete(mmInfo);
	mm_info_delete(mmInfo2);
	cli_options_deleteNonPtr(&options);
}

int isSameTriplet(quintet_t* t1, quintet_t* t2)
{
	return (t1->i == t2->i) && (t1->j == t2->j) &&
			(t1->iVal == t2->iVal) && (t1->jVal == t2->jVal) &&
			(t1->val == t2->val);
}

void printTriplet(DECIMAL i, quintet_t* t)
{
	PRINTF("%d. i=%d j=%d ival=%d jval=%d val=%f\t", i, t->i, t->j, t->iVal, t->jVal, t->val);
}

int compareTriplets(quintet_t* triplets1, quintet_t* triplets2, DECIMAL nnz)
{
	DECIMAL i;
	for(i = 0; i < nnz; ++i)
	{
		if(isSameTriplet(&triplets1[i], &triplets2[i]) != TRUE)
		{
			printTriplet(i, &triplets1[i]);
			printTriplet(i, &triplets2[i]);
			PRINTF("\n");
		}
	}

	return TRUE;
}
