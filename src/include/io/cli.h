/*
 * cli.h
 *
 *  Created on: Sep 2, 2013
 *      Author: matara
 */

#ifndef CLI_H_
#define CLI_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "include/config.h"

// Command line options
// --------------------------------------------------------------

#define CLO_MATRIX_PATH "matrix"

#define CLO_NUM_BLOCKS "num_blocks"
#define CLO_THREADS_PER_BLOCK "threads_per_block"

#define CLO_STORAGE_FORMAT "storage_format"
#define CLO_PARTITION_TYPE "partition_type"
#define CLO_PARTITION_METHOD "partition_method"
#define CLO_ORDERING_TYPE "ordering_type"
#define CLO_TARGETED_CACHE_SIZE_KB "targeted_cache_size_kb"
#define CLO_ROW_SLICE_LENGTH "row_slice_length"
#define CLO_STEAL_TRESHOLD "steal_treshold"
#define CLO_SIMD_LENGTH "simd_length"

#define CLO_RUN_COUNT_WARM_UP "run_count_warm_up"
#define CLO_RUN_COUNT_MEASURE "run_count_measure"


// Values of the command line options
// --------------------------------------------------------------

#define CLO_PARAM_PARTITION_METHOD_REGULAR "regular"
#define CLO_PARAM_PARTITION_METHOD_RECURSIVE_BIPARTITION "recursive_bipartition"

#define CLO_PARAM_STORAGE_FORMAT_CSR "CSR"
#define CLO_PARAM_STORAGE_FORMAT_JDS "JDS"
#define CLO_PARAM_STORAGE_FORMAT_HYBRID_JDS_CSR "HYBRID_JDS_CSR"

#define CLO_PARAM_PARTITION_TYPE_1D_ROW_STR "1D_row_wise"
#define CLO_PARAM_PARTITION_TYPE_1D_COLUMN_STR "1D_column_wise"
#define CLO_PARAM_PARTITION_TYPE_2D_CHECKERBOARD_STR "2D_checker_board"

#define CLO_PARAM_PARTITION_METHOD_REGULAR "regular"
#define CLO_PARAM_PARTITION_METHOD_RECURSIVE_BIPARTITION "recursive_bipartition"

#define CLO_PARAM_ORDERING_TYPE_NONE "none"
#define CLO_PARAM_ORDERING_TYPE_ROW_NET "row_net"
#define CLO_PARAM_ORDERING_TYPE_COLUMN_NET "column_net"

#define TRUE_STR "true"
#define FALSE_STR "false"



#define TRAVERSAL_IN_ORDER_STR "inorder"
#define TRAVERSAL_POST_ORDER_STR "postorder"

#define SCHEDULER_TREE_STR "tree"
#define SCHEDULER_RING_STR "ring"

// --------------------------------------------------------------

struct cli_options
{
	char* mmfPath;
	char* mmfDir;
	char* mmfName;

	int numBlocks;
	int threadsPerBlock;

	int storageFormat;
	int partitionType;
	int partitionMethod;
	int orderingType;

	int simdLength;
	float targetedCacheSizeKB;
	int stealTreshold;

	int runCountWarmUp;
	int runCountMeasure;
};

typedef struct cli_options cli_options_t;

extern void cli_options_init(cli_options_t* options);
extern void cli_options_deleteNonPtr(cli_options_t* options);

extern void cli_printUsage(void);
extern void cli_parseInput(
		int argsc, char* argsv[], cli_options_t* options);
extern int cli_isValid(cli_options_t* options);
extern void cli_print(cli_options_t* options);
extern char* cli_getTraversalTypeStr(cli_options_t* options);


#endif /* CLI_H_ */
