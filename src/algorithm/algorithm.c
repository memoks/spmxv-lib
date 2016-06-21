
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>
#include "omp.h"

#include "include/algorithm/algorithm.h"

void algorithm_parallelMergeQuickSort(
		quintet_t* quintets_inout, DECIMAL length, int threadCount,
		int (*quintet_cmpFunc) (const void* left, const void* right))
{
	DECIMAL* borders = (DECIMAL*) malloc(sizeof(DECIMAL) * (threadCount + 1));
	DECIMAL minCount = length / threadCount;
	DECIMAL extra = length % threadCount;

	DECIMAL i;
	borders[0] = 0;
	for(i = 1; i < threadCount; ++i)
	{
		borders[i] = borders[i - 1] + minCount;
		if(i < extra)
			++borders[i];
	}
	borders[i] = length;

	// Sort equal chunks simultaneously by quick sort
	#pragma omp parallel num_threads(threadCount)
	{
		int threadId = omp_get_thread_num();
		DECIMAL startInd = borders[threadId];
		DECIMAL endInd = borders[threadId + 1];

		quintet_sortCustom(&quintets_inout[startInd], endInd - startInd, quintet_cmpFunc);
	}

	// Merge results
	algorithm_parallelMerge(quintets_inout, borders, threadCount, quintet_cmpFunc);

	// clean up
	free(borders);
}

extern void algorithm_parallelMerge(
		quintet_t* quintets, DECIMAL* bordersInitial, DECIMAL borderCount,
		int (*quintet_cmpFunc) (const void* left, const void* right))
{
	if(borderCount < 2)
		return;

	quintet_t* left = NULL;
	DECIMAL leftLength = 0;
	quintet_t* right = NULL;
	DECIMAL rightLength = 0;

	// There are fewer elements (to sort) than there are threads. So, just quit.
	if(bordersInitial[borderCount] - bordersInitial[0] < 1)
		return;

	// recursive thread fork
	int borderCounts[2];
	borderCounts[0] = borderCount / 2 + borderCount % 2;
	borderCounts[1] = borderCount / 2;

	DECIMAL quintetStartIndex[3];
	quintetStartIndex[0] = bordersInitial[0];
	quintetStartIndex[1] = bordersInitial[borderCounts[0]];
	quintetStartIndex[2] = bordersInitial[borderCount];

	DECIMAL borderStartIndex[2];
	borderStartIndex[0] = 0;
	borderStartIndex[1] = borderCounts[0];

	leftLength = quintetStartIndex[1] - quintetStartIndex[0];
	left = &quintets[0];
	rightLength = quintetStartIndex[2] - quintetStartIndex[1];
	right = &quintets[leftLength];

	quintet_t* tQuintets[2] = {left, right};
	DECIMAL* tBorders[2] = {&bordersInitial[borderStartIndex[0]], &bordersInitial[borderStartIndex[1]]};

	// recursive fork
	if(borderCount > 2)
	{
		// Serial
		// algorithm_parallelMerge(tQuintets[0], tBorders[0], borderCounts[0], input_cmpQuintetFunc);
		// algorithm_parallelMerge(tQuintets[1], tBorders[1], borderCounts[1], input_cmpQuintetFunc);

		#pragma omp parallel sections
		{
			#pragma omp section
			{
				algorithm_parallelMerge(tQuintets[0], tBorders[0], borderCounts[0], quintet_cmpFunc);
			}
			#pragma omp section
			{
				algorithm_parallelMerge(tQuintets[1], tBorders[1], borderCounts[1], quintet_cmpFunc);
			}
		}

		/*
		#pragma omp parallel num_threads(2)
		{
			int tId = omp_get_thread_num();
			algorithm_parallelMerge(
					tQuintets[tId], tBorders[tId], borderCounts[tId], input_cmpQuintetFunc);
		}
		*/
	}

	// merge
	DECIMAL fullLength = leftLength + rightLength;
	quintet_t* merged = (quintet_t*) malloc(sizeof(quintet_t) * fullLength);
	algorithm_merge(
			left, leftLength, right, rightLength, merged,
			quintet_cmpFunc);
	quintet_copyQuintets(merged, quintets, fullLength);

	// clean up
	free(merged);
}

void algorithm_merge(
		quintet_t* quintets1, DECIMAL length1,
		quintet_t* quintets2, DECIMAL length2,
		quintet_t* quintetsMerged_inout,
		int (*quintet_cmpFunc) (const void*, const void*))
{
	DECIMAL index1 = 0;
	DECIMAL index2 = 0;
	DECIMAL indexMerged = 0;

	while((index1 < length1) && (index2 < length2))
	{
		if(quintet_cmpFunc((void*) &quintets1[index1], (void*) &quintets2[index2]) < 0)
			quintetsMerged_inout[indexMerged] = quintets1[index1++];
		else
			quintetsMerged_inout[indexMerged] = quintets2[index2++];

		++indexMerged;
	}

	while(index1 < length1)
		quintetsMerged_inout[indexMerged++] = quintets1[index1++];

	while(index2 < length2)
		quintetsMerged_inout[indexMerged++] = quintets2[index2++];
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
