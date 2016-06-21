
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "include/config.h"
#include "include/io/cli.h"
#include "include/timer/custom_timer.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/spm.h"
#include "include/interface/sequential.h"

/*

void (*initialize)(struct sequential* self);
void (*warmUp)(struct sequential* self, vector_real_t* x);
void (*copy)(struct sequential* self, vector_real_t* dest, vector_real_t* source);
void (*add)(struct sequential* self, vector_real_t* r, vector_real_t* lhs, vector_real_t* rhs);
void (*subtract)(struct sequential* self, vector_real_t* r, vector_real_t* lhs, vector_real_t* rhs);
void (*accumulate)(struct sequential* self, vector_real_t* res, vector_real_t* lhs, REAL scalar, vector_real_t* rhs);
void (*matrixMult)(struct sequential* self,	vector_real_t* x, vector_real_t* y_out);
void (*scalarMult)(struct sequential* self,	vector_real_t* x, REAL scalar, vector_real_t* y_out);
REAL (*innerProd)(struct sequential* self, vector_real_t* r1, vector_real_t* r2);
void (*terminate)(struct sequential* self);

 */
int main(int argsc, char* argsv[])
{
	// read command line parameters
	// -----------------------------------------------------------------------------------------------------------------

	cli_options_t options;
	cli_options_init(&options);
	cli_parseInput(argsc, argsv, &options);
	cli_printUsage();
	cli_print(&options);

	sequential_routine* sr = sequential_routine_new(&options);
	sr->initialize(sr);


	// Basic Conjugate Gradient algorithm to solve Ax = b for x vector
	// -----------------------------------------------------------------------------------------------------------------

	// Choose an initial x vector
	vector_real_t* x = vector_real_random(sr->cr->colCount);
	vector_real_t* b = vector_real_new(sr->cr->rowCount);
	b->data[0] = 5; b->data[1] = 4; b->data[2] = 3; b->data[3] = 2; b->data[4] = 1;

	vector_real_t* r = vector_real_new(sr->cr->rowCount);
	vector_real_t* rr = vector_real_new(sr->cr->rowCount);
	vector_real_t* p = vector_real_new(sr->cr->rowCount);
	vector_real_t* pp = vector_real_new(sr->cr->rowCount);
	vector_real_t* z = vector_real_new(sr->cr->rowCount);
	vector_real_t* zz = vector_real_new(sr->cr->rowCount);

	REAL ak = 0.0;
	REAL akden = 0.0;
	REAL bk = 0.0;
	REAL bkden = 1.0;
	REAL bnrm = 0.0;
	REAL bknum = 0.0;

	REAL err;

	DECIMAL iter = 0;

	sr->matrixMult(sr, x, r);
	sr->subtract(sr, r, b, r);
	sr->copy(sr, rr, r);

	bnrm = sqrt(sr->innerProd(sr, b, b));
	sr->asolve(sr, z, r);

	while(iter < ITER_MAX)
	{
		++iter;
		sr->asolve(sr, zz, rr);

		bknum = sr->innerProd(sr, z, rr);

		if(iter == 1)
		{
			sr->copy(sr, p, z);
			sr->copy(sr, pp, zz);
		}
		else
		{
			bk = bknum / bkden;
			void (*accumulate)(struct sequential* self, vector_real_t* res, vector_real_t* lhs, REAL scalar, vector_real_t* rhs);
			sr->accumulate(sr, p, z, bk, p);
			sr->accumulate(sr, pp, zz, bk, pp);
		}

		bkden = bknum;
		sr->matrixMult(sr, z, p);

		akden = sr->innerProd(sr, z, pp);

		ak = bknum / akden;
		sr->matrixMult(sr, zz, pp);

		sr->accumulate(sr, x, x, ak, p);
		sr->accumulate(sr, r, r, -ak, z);
		sr->accumulate(sr, rr, rr, -ak, zz);

		sr->asolve(sr, z, r);

		err = sqrt(sr->innerProd(sr, r, r)) / bnrm;

		if(err <= RESULT_ERROR_TRESHOLD)
			break;
	}

	// Print result to file
	// -----------------------------------------------------------------------------------------------------------------
	// if(options.xPath != NULL)
		//input_printVectorToFile(options.xPath, x);


	vector_real_print("RESULT", x);

	// Clean-up
	// -----------------------------------------------------------------------------------------------------------------

	vector_real_delete(x);
	vector_real_delete(b);
	vector_real_delete(r);
	vector_real_delete(rr);
	vector_real_delete(p);
	vector_real_delete(pp);
	vector_real_delete(z);
	vector_real_delete(zz);

	sequential_routine_cleanUp(sr);
	cli_options_deleteNonPtr(&options);

	return EXIT_SUCCESS;
}
