
#include <stdlib.h>
#include <stdio.h>

#include "include/io/converter.h"
#include "include/io/input_parser.h"

#include "include/interface/common.h"


extern void common_routine_initialize(common_routine* self);
extern void common_routine_reorderVector(common_routine* self, vector_real_t* v);
extern void common_routine_terminate(common_routine* self);


common_routine* common_routine_new(cli_options_t* options)
{
	common_routine* cr = (common_routine*) malloc(sizeof(common_routine));
	cr->options = options;

	cr->initialize = common_routine_initialize;
	cr->reorderVector = common_routine_reorderVector;
	cr->terminate = common_routine_terminate;

	cr->nnz = 0;
	cr->rowCount = 0;
	cr->colCount = 0;
	cr->triplets = NULL;
	cr->rowOrderLookup = NULL;
	cr->colOrderLookup = NULL;
	cr->subMtxData = NULL;
	cr->spm = NULL;

	return cr;
}

void common_routine_cleanUp(common_routine* cr)
{
	cr->terminate(cr);
	free(cr);
}

// -------------------------------------------------------------------------------------------

void common_routine_initialize(common_routine* self)
{
	cli_options_t* opt = self->options;

	//self->triplets = input_readTriplets(self->options->mmfPath, &self->nnz, &self->rowCount, &self->colCount,
	//		NULL, self->options->isSymmetric, self->options->hasValues);
	quintet_sortRowIndex(self->triplets, self->nnz);

	if(self->options->useOrdering)
	{
//		input_readOrderLookup(opt->orderingPath,
//				self->triplets, self->nnz, self->rowCount, self->colCount,
//				&self->rowOrderLookup, &self->colOrderLookup);

		quintet_sortRowValue(self->triplets, self->nnz);

		// self->subMtxData = input_readBlockFile(opt->blockPath);
	}


	// self->spm = converter_tripletToCSR(self->triplets, self->nnz, self->rowCount, self->colCount);

	free(self->triplets);
}

void common_routine_reorderVector(common_routine* self, vector_real_t* v)
{
	input_reorderVector(&v, self->rowOrderLookup);
}

void common_routine_terminate(common_routine* self)
{
	vector_int_delete(self->rowOrderLookup);
	vector_int_delete(self->colOrderLookup);
	vector_int_delete(self->subMtxData);
	// spm_cmp_delete(self->spm);
}
