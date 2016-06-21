
#include <stdlib.h>
#include <stdio.h>

#include "arch/generic/inc/spmxv_sequential.h"

#ifdef __ICC
#include "arch/icc/inc/spmxv_mkl.h"
#endif

#include "include/interface/sequential.h"


// values to function pointers in sequential_routine struct
// -----------------------------------------------------------------------------------------------------------

extern void sequential_routine_initialize(sequential_routine* self);
extern void sequential_routine_warmUp(sequential_routine* self, vector_real_t* x);
extern void sequential_routine_copy(sequential_routine* self, vector_real_t* dest, vector_real_t* source);
extern void sequential_routine_add(sequential_routine* self, vector_real_t* r, vector_real_t* lhs, vector_real_t* rhs);
extern void sequential_routine_subtract(sequential_routine* self, vector_real_t* r, vector_real_t* lhs, vector_real_t* rhs);
extern void sequential_routine_accumulate(sequential_routine* self,
		vector_real_t* r, vector_real_t* v1, REAL scalar, vector_real_t* v2);
extern void sequential_routine_matrixMult(sequential_routine* self, vector_real_t* x, vector_real_t* y_out);
extern void sequential_routine_scalarMult(sequential_routine* self, vector_real_t* x, REAL scalar, vector_real_t* y_out);
extern REAL sequential_routine_innerProd(sequential_routine* self, vector_real_t* r1, vector_real_t* r2);

extern void sequential_routine_asolve(sequential_routine* self, vector_real_t* res, vector_real_t* b);

extern void sequential_routine_terminate(sequential_routine* self);

// -----------------------------------------------------------------------------------------------------------

sequential_routine* sequential_routine_new(cli_options_t* options)
{
	sequential_routine* sr = (sequential_routine*) malloc(sizeof(sequential_routine));

	sr->cr = common_routine_new(options);
	sr->cr->initialize(sr->cr);

	sr->spm = sr->cr->spm;

	sr->initialize = sequential_routine_initialize;
	sr->warmUp = sequential_routine_warmUp;
	sr->terminate = sequential_routine_terminate;

	sr->copy = sequential_routine_copy;
	sr->add = sequential_routine_add;
	sr->subtract = sequential_routine_subtract;
	sr->accumulate = sequential_routine_accumulate;
	sr->matrixMult = sequential_routine_matrixMult;
	sr->scalarMult = sequential_routine_scalarMult;
	sr->innerProd = sequential_routine_innerProd;

	sr->asolve = sequential_routine_asolve;

	return sr;
}

void sequential_routine_cleanUp(sequential_routine* sr)
{
	sr->terminate(sr);
	free(sr);
}


// sequential_routine struct function definitions
// -----------------------------------------------------------------------------------------------------------

void sequential_routine_initialize(sequential_routine* self)
{
}

void sequential_routine_warmUp(sequential_routine* self, vector_real_t* x)
{
	vector_real_t* yTemp = vector_real_new(self->spm->rowCount);
	vector_real_t* xTemp = vector_real_new(self->spm->colCount);
	vector_real_copy(x, xTemp);

	DECIMAL i;
	for(i = 0; i < self->cr->options->runCountWarmUp; ++i)
	{
		// spmxv_sequential_CSR(self->spm, xTemp, yTemp);
		vector_real_copy(yTemp, xTemp);
	}

	vector_real_delete(xTemp);
	vector_real_delete(yTemp);
}

void sequential_routine_copy(sequential_routine* self, vector_real_t* dest, vector_real_t* source)
{
	DECIMAL i;
	for(i = 0; i < dest->length; ++i)
		dest->data[i] = source->data[i];
}

void sequential_routine_add(sequential_routine* self, vector_real_t* r, vector_real_t* lhs, vector_real_t* rhs)
{
	DECIMAL i;
	for(i = 0; i < r->length; ++i)
		r->data[i] = lhs->data[i] + rhs->data[i];
}

void sequential_routine_subtract(sequential_routine* self, vector_real_t* r, vector_real_t* lhs, vector_real_t* rhs)
{
	DECIMAL i;
	for(i = 0; i < r->length; ++i)
		r->data[i] = lhs->data[i] - rhs->data[i];
}

void sequential_routine_accumulate(sequential_routine* self,
		vector_real_t* res, vector_real_t* lhs, REAL scalar, vector_real_t* rhs)
{
	DECIMAL i;
	for(i = 0; i < res->length; ++i)
		res->data[i] = lhs->data[i] + scalar * rhs->data[i];
}

void sequential_routine_matrixMult(sequential_routine* self, vector_real_t* x, vector_real_t* y_out)
{
#ifdef __ICC
	// printf("Using MKL\n");
	// spmxv_mkl(self->spm, x, y_out);
#else
	// spmxv_sequential_CSR(self->spm, x, y_out);
#endif
}

void sequential_routine_scalarMult(sequential_routine* self, vector_real_t* x, REAL scalar, vector_real_t* y_out)
{
	DECIMAL i;
	for(i = 0; i < x->length; ++i)
		y_out->data[i] = x->data[i] * scalar;
}

REAL sequential_routine_innerProd(sequential_routine* self, vector_real_t* r1, vector_real_t* r2)
{
	REAL result = 0.0;
	int i;
	for(i = 0; i < r1->length; ++i)
		result += r1->data[i] * r2->data[i];

	return result;
}

void sequential_routine_asolve(sequential_routine* self, vector_real_t* res, vector_real_t* b)
{
	DECIMAL i;
	for(i = 0; i < self->spm->rowCount; ++i)
	{
		res->data[i] = b->data[i];

		REAL diagonal = spm_cmp_get(self->spm, i, i);
		if(diagonal != 0)
			res->data[i] /= diagonal;
	}
}

void sequential_routine_terminate(sequential_routine* self)
{
	common_routine_cleanUp(self->cr);
}
