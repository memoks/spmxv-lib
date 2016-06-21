
#include <stdlib.h>
#include <stdio.h>

#include "patoh.h"

#include "include/config.h"
#include "include/io/input_parser.h"
#include "include/io/cli.h"


int main(int argc, char** argv)
{
	if(argc < 10)
	{
		printf("Usage: %s mtx_path imbalance_ratio vertexOrder vertexDim netOrder netDim vertexPartVec netPartVec subMtxPathPrefix\n", argv[0]);
		printf("\nQuick start: %s t5-8-11.mtx 2  roworder rowdim colorder coldim vertexpartvector netpartvector\n", argv[0]);
		return EXIT_FAILURE;
	}

	char* mmfPath = argv[1];
	REAL imbalanceRatio = atof(argv[2]);
	char* vertexOrderingPath = argv[3];
	char* vertexDimensionsPath = argv[4];
	char* netOrderingPath = argv[5];
	char* netDimensionsPath[6];
	char* vertexPartVectorPath[7];
	char* netPartVectorPath[8];

	INTEGER NNZ = 0;
	INTEGER rowCount = 0;
	INTEGER colCount = 0;
	triplet_t* triplets = NULL;
	mm_info_t* mmInfo = NULL;
	input_readTriplets(mmfPath, &NNZ, &rowCount, &colCount, &triplets, &mmInfo);
	input_sortTripletRowIndex(triplets, NNZ);

	PaToH_Parameters args;
	PaToH_Initialize_Parameters(&args, PATOH_CONPART, PATOH_SUGPARAM_SPEED);
	args.final_imbal = imbalanceRatio;

	// column net model
	INTEGER cellCount = rowCount;
	INTEGER netCount = colCount;

	INTEGER* cellWeights = calloc(cellCount, sizeof(INTEGER));
	INTEGER* partVector = calloc(cellCount, sizeof(INTEGER));
	INTEGER* netWeights = calloc(cellCount, sizeof(INTEGER));

	INTEGER i;
	for(i = 0; i < NNZ; ++i)
	{
		triplet_t*
	}

	return EXIT_SUCCESS;
}
