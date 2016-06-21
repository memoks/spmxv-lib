
#include <stdlib.h>
#include <stdio.h>

#include "include/config.h"
#include "include/io/cli.h"
#include "include/timer/custom_timer.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/spm.h"
#include "include/interface/sequential.h"

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
	vector_real_t* p = vector_real_new(sr->cr->rowCount);
	vector_real_t* q = vector_real_new(sr->cr->colCount);

	REAL ro = 0.0;
	REAL roNew = 0.0;
	REAL pi = 0.0;
	REAL alpha = 0.0;
	REAL beta = 0.0;


	// r = p = b - Ax
	sr->matrixMult(sr, x, p);
	sr->subtract(sr, p, b, p);
	sr->copy(sr, r, p);

	// ro = <r, r>
	ro = sr->innerProd(sr, r, r);

	while(ro > RESULT_ERROR_TRESHOLD)
	{
		// q = Ap
		sr->matrixMult(sr, p, q);

		// pi = <p, q>
		pi = sr->innerProd(sr, p, q);

		alpha = ro / pi;

		// x = x + alpha * p
		sr->accumulate(sr, x, x, alpha, p);

		// r = r - alpha * q
		sr->accumulate(sr, r, r, -alpha, q);

		// roNew = <r, r>
		roNew = sr->innerProd(sr, r, r);

		beta = roNew / ro;
		ro = roNew;

		// p = r + beta * p
		sr->accumulate(sr, p, r, beta, p);
	}

	vector_real_print("x: ", x);

	// Print result to file
	// -----------------------------------------------------------------------------------------------------------------
	// if(options.xPath != NULL)
		//input_printVectorToFile(options.xPath, x);

	// Clean-up
	// -----------------------------------------------------------------------------------------------------------------

	vector_real_delete(x);
	vector_real_delete(b);
	vector_real_delete(r);
	vector_real_delete(p);
	vector_real_delete(q);

	sequential_routine_cleanUp(sr);
	cli_options_deleteNonPtr(&options);

	return EXIT_SUCCESS;
}
