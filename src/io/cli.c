
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "include/config.h"
#include "include/control_unit/cu_options.h"

#include "include/io/cli.h"

// helper functions
// ----------------------------------------------------------------------------------------------------------------

static void cli_parseMMFInfo(char* mmfPath, char** mmfDir_out, char** mmfName_out);

// ----------------------------------------------------------------------------------------------------------------

void cli_options_init(cli_options_t* options)
{
	options->mmfPath = NULL;
	options->mmfDir = NULL;
	options->mmfName = NULL;

	options->numBlocks = -1;
	options->threadsPerBlock = -1;

	options->storageFormat = -1;
	options->partitionType = -1;
	options->partitionMethod = -1;

	options->targetedCacheSizeKB = -1;
	options->simdLength = DEFAULT_SIMD_LENGTH;
	options->stealTreshold = DEFAULT_STEAL_TRESHOLD;

	options->runCountWarmUp = DEFAULT_WARM_UP_RUN_COUNT;
	options->runCountMeasure = DEFAULT_MEASURE_RUN_COUNT;
}

void cli_options_deleteNonPtr(cli_options_t* options)
{
	if(options->mmfPath != NULL) free(options->mmfPath);
	if(options->mmfDir != NULL) free(options->mmfDir);
	if(options->mmfName != NULL) free(options->mmfName);
}

void cli_printUsage(void)
{
	PRINTF("NOTE: Set thread count affinity using OMP environment variables!!\n");
	PRINTF("\n");

	PRINTF("USAGE: ./bin\n");
	PRINTF("\t%s=<mmf_path>\n", CLO_MATRIX_PATH);
	PRINTF("\t%s=<integer>\n", CLO_NUM_BLOCKS);
	PRINTF("\t%s=<integer>\n", CLO_THREADS_PER_BLOCK);
	PRINTF("\t%s=<%s|%s|%s>\n", CLO_STORAGE_FORMAT,
			CLO_PARAM_STORAGE_FORMAT_CSR, CLO_PARAM_STORAGE_FORMAT_JDS, CLO_PARAM_STORAGE_FORMAT_HYBRID_JDS_CSR);
	PRINTF("\t%s=<%s|%s|%s>\n", CLO_PARTITION_TYPE,
			CLO_PARAM_PARTITION_TYPE_1D_ROW_STR, CLO_PARAM_PARTITION_TYPE_1D_COLUMN_STR, CLO_PARAM_PARTITION_TYPE_2D_CHECKERBOARD_STR);
	PRINTF("\t%s=<%s|%s>\n", CLO_PARTITION_METHOD,
			CLO_PARAM_PARTITION_METHOD_REGULAR, CLO_PARAM_PARTITION_METHOD_RECURSIVE_BIPARTITION);
	PRINTF("\t%s=<%s|%s|%s>\n", CLO_ORDERING_TYPE,
			CLO_PARAM_ORDERING_TYPE_NONE, CLO_PARAM_ORDERING_TYPE_COLUMN_NET, CLO_PARAM_ORDERING_TYPE_ROW_NET);
	PRINTF("\t%s=<integer>\n", CLO_TARGETED_CACHE_SIZE_KB);
	PRINTF("\t[%s=<integer:default=1>]\n", CLO_STEAL_TRESHOLD);
	PRINTF("\t[%s=<integer:default=(float_16 / double_8)>]\n", CLO_SIMD_LENGTH);
	PRINTF("\t[%s=<integer:default=10>]\n", CLO_RUN_COUNT_WARM_UP);
	PRINTF("\t[%s=<integer:default=100>]\n", CLO_RUN_COUNT_MEASURE);
	PRINTF("\n\n");


	PRINTF("EXAMPLE: ./bin\n");
	PRINTF("\t%s=input/wheel_601/wheel_601.mtx\n", CLO_MATRIX_PATH);
	PRINTF("\t%s=2 %s=2\n", CLO_NUM_BLOCKS, CLO_THREADS_PER_BLOCK);
	PRINTF("\t%s=%s\n", CLO_STORAGE_FORMAT, CLO_PARAM_STORAGE_FORMAT_CSR);
	PRINTF("\t%s=%s\n", CLO_PARTITION_TYPE, CLO_PARAM_PARTITION_TYPE_1D_ROW_STR);
	PRINTF("\t%s=%s\n", CLO_PARTITION_METHOD, CLO_PARAM_PARTITION_METHOD_REGULAR);
	PRINTF("\t%s=1024\n", CLO_TARGETED_CACHE_SIZE_KB);
	PRINTF("\n\n");
}

int cli_isValid(cli_options_t* options)
{
	return (options->mmfPath != NULL) &&
			(options->numBlocks > 0) &&
			(options->threadsPerBlock > 0) &&
			(options->storageFormat > 0) &&
			(options->partitionType > 0) &&
			(options->partitionMethod > 0) &&
			(options->orderingType > -1) &&
			(options->targetedCacheSizeKB > 0);
}

void cli_parseInput(int argsc, char* argsv[], cli_options_t* options)
{
	int inputSize = 12;
	char* inputNames[inputSize];
	inputNames[0] = CLO_MATRIX_PATH;
	inputNames[1] = CLO_NUM_BLOCKS;
	inputNames[2] = CLO_THREADS_PER_BLOCK;
	inputNames[3] = CLO_STORAGE_FORMAT;
	inputNames[4] = CLO_PARTITION_TYPE;
	inputNames[5] = CLO_PARTITION_METHOD;
	inputNames[6] = CLO_ORDERING_TYPE;
	inputNames[7] = CLO_TARGETED_CACHE_SIZE_KB;
	inputNames[8] = CLO_STEAL_TRESHOLD;
	inputNames[9] = CLO_RUN_COUNT_WARM_UP;
	inputNames[10] = CLO_RUN_COUNT_MEASURE;
	inputNames[11] = CLO_SIMD_LENGTH;

	int i;
	int j;
	for(i = 1; i < argsc; ++i)
	{
		for(j = 0; j < inputSize; ++j)
		{
			int inputNo = -1;
			if(strncmp(inputNames[j], argsv[i], strlen(inputNames[j])))
				continue;
			else
				inputNo = j;

			int offset = strlen(inputNames[j]) + 1;
			char* valuePtr = argsv[i] + offset;

			char* temp = (char *) malloc(sizeof(char) * (strlen(valuePtr) + 1));
			strcpy(temp, valuePtr);

			switch(inputNo)
			{
			case 0: // matrix file path
				options->mmfPath = temp;
				cli_parseMMFInfo(temp, &options->mmfDir, &options->mmfName);
				break;
			case 1: // numBlocks
				options->numBlocks = atoi(temp);
				free(temp);
				break;
			case 2: // thread count per block
				options->threadsPerBlock = atoi(temp);
				free(temp);
				break;
			case 3: // storage format
				if(!strcmp(temp, CLO_PARAM_STORAGE_FORMAT_CSR))
					options->storageFormat = SPM_STORAGE_CSR;
				else if(!strcmp(temp, CLO_PARAM_STORAGE_FORMAT_JDS))
					options->storageFormat = SPM_STORAGE_JDS;
				else if(!strcmp(temp, CLO_PARAM_STORAGE_FORMAT_HYBRID_JDS_CSR))
					options->storageFormat = SPM_STORAGE_HYBRID_JDS_CSR;
				else
				{
					PRINTF("Storage format \"%s\" not recognized!! ", temp);
					PRINTF("Available options: %s, %s, %s. ",
							CLO_PARAM_STORAGE_FORMAT_CSR,
							CLO_PARAM_STORAGE_FORMAT_JDS,
							CLO_PARAM_STORAGE_FORMAT_HYBRID_JDS_CSR);
					PRINTF("Exiting program.\n");
					exit(EXIT_FAILURE);
				}
				free(temp);
				break;
			case 4: // partition type
				if(!strcmp(temp, CLO_PARAM_PARTITION_TYPE_1D_COLUMN_STR))
					options->partitionType = PARTITION_TYPE_1D_COLUMN_SLICE;
				else if(!strcmp(temp, CLO_PARAM_PARTITION_TYPE_1D_ROW_STR))
					options->partitionType = PARTITION_TYPE_1D_ROW_SLICE;
				else if(!strcmp(temp, CLO_PARAM_PARTITION_TYPE_2D_CHECKERBOARD_STR))
					options->partitionType = PARTITION_TYPE_2D_CHECKER_BOARD;
				else
				{
					PRINTF("Partition type \"%s\" not recognized!! ", temp);
					PRINTF("Available options: %s, %s, %s. ",
							CLO_PARAM_PARTITION_TYPE_1D_COLUMN_STR,
							CLO_PARAM_PARTITION_TYPE_1D_ROW_STR,
							CLO_PARAM_PARTITION_TYPE_2D_CHECKERBOARD_STR);
					PRINTF("Exiting program.\n");
					exit(EXIT_FAILURE);
				}
				free(temp);
				break;
			case 5: // partition method
				if(!strcmp(temp, CLO_PARAM_PARTITION_METHOD_REGULAR))
					options->partitionMethod = PARTITION_METHOD_REGULAR;
				else if(!strcmp(temp, CLO_PARAM_PARTITION_METHOD_RECURSIVE_BIPARTITION))
					options->partitionMethod = PARTITION_METHOD_RECURSIVE_BIPARTITION;
				else
				{
					PRINTF("Partition method \"%s\" not recognized!! ", temp);
					PRINTF("Available options: %s, %s. ",
							CLO_PARAM_PARTITION_METHOD_REGULAR,
							CLO_PARAM_PARTITION_METHOD_RECURSIVE_BIPARTITION);
					PRINTF("Exiting program.\n");
					exit(EXIT_FAILURE);
				}
				free(temp);
				break;
			case 6: // ordering type
				if(!strcmp(temp, CLO_PARAM_ORDERING_TYPE_NONE))
					options->orderingType = ORDERING_TYPE_NONE;
				else if(!strcmp(temp, CLO_PARAM_ORDERING_TYPE_COLUMN_NET))
					options->orderingType = ORDERING_TYPE_COLUMN_NET;
				else if(!strcmp(temp, CLO_PARAM_ORDERING_TYPE_ROW_NET))
					options->orderingType = ORDERING_TYPE_ROW_NET;
				else
				{
					PRINTF("Ordering type \"%s\" not recognized!! ", temp);
					PRINTF("Available options: %s, %s, %s. ",
							CLO_PARAM_ORDERING_TYPE_NONE,
							CLO_PARAM_ORDERING_TYPE_COLUMN_NET,
							CLO_PARAM_ORDERING_TYPE_ROW_NET);
					PRINTF("Exiting program.\n");
					exit(EXIT_FAILURE);
				}
				free(temp);
				break;
			case 7: // targeted cache size in KB
				options->targetedCacheSizeKB = atof(temp);
				free(temp);
				break;
			case 8: // steal treshold for DWS algorithm
				options->stealTreshold = atoi(temp);
				free(temp);
				break;
			case 9: // run count warm up
				options->runCountWarmUp = atoi(temp);
				free(temp);
				break;
			case 10: // run count measure
				options->runCountMeasure = atoi(temp);
				free(temp);
				break;
			case 11: // SIMD length for vector processors
				options->simdLength = atoi(temp);
				free(temp);
				break;
			default:
				break;
			}

		}
	}
}

void cli_print(cli_options_t* options)
{
	PRINTF("COMMAND_LINE_INTERFACE OPTIONS\n");
	PRINTF("mmfPath: %s\n", options->mmfPath);
	PRINTF("mmfDir: %s\n", options->mmfDir);
	PRINTF("mmfName: %s\n", options->mmfName);

	PRINTF("numBlocks: %d\n", options->numBlocks);
	PRINTF("threadsPerBlock: %d\n", options->threadsPerBlock);

	PRINTF("storageFormat: %d\n", options->storageFormat);
	PRINTF("partitionType: %d\n", options->partitionType);
	PRINTF("partitionMethod: %d\n", options->partitionMethod);
	PRINTF("orderingType: %d\n", options->orderingType);

	PRINTF("targetedCacheSizeKB: %0.1f\n", options->targetedCacheSizeKB);
	PRINTF("simdLength: %d\n", options->simdLength);
	PRINTF("stealTreshold: %d\n", options->stealTreshold);

	PRINTF("runCountWarmUp: %d\n", options->runCountWarmUp);
	PRINTF("runCountMeasure: %d\n", options->runCountMeasure);
}

// Some helper functions
// ------------------------------------------------------------------------------------------------------

static void cli_parseMMFInfo(char* mmfPath, char** mmfDir_out, char** mmfName_out)
{
	int length = strlen(mmfPath);
	int delimeterIndex = length;
	int i;

	for(i = length - 1; i >= 0; --i)
	{
		if(mmfPath[i] == '/')
			break;

		--delimeterIndex;
	}

	char* mmfDir = (char*) malloc(sizeof(char) * delimeterIndex);
	char* mmfName = (char*) malloc(sizeof(char) * (length - delimeterIndex + 1));

	strncpy(mmfDir, mmfPath, delimeterIndex - 1);
	mmfDir[delimeterIndex - 1] = '\0';
	strcpy(mmfName, &mmfPath[delimeterIndex]);
	mmfName[length - delimeterIndex - 4] = '\0';

	*mmfDir_out = mmfDir;
	*mmfName_out = mmfName;
}
