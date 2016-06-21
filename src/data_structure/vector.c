
#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <alloca.h>
#include <omp.h>

#ifdef __ICC
#include <mkl.h>
#endif

#include "include/config.h"

#include "include/data_structure/vector.h"

vector_int_t* vector_int_new(int length)
{
	vector_int_t* vector = NULL;

#ifdef __ICC
	vector = (vector_int_t*) mkl_malloc(sizeof(vector_int_t), ALIGNMENT);
	vector->data = (DECIMAL*) mkl_malloc(length * sizeof(DECIMAL), ALIGNMENT);
#else
	posix_memalign(&vector, ALIGNMENT, sizeof(vector_int_t));
	posix_memalign(&vector->data, ALIGNMENT, sizeof(DECIMAL) * length);
#endif

	vector->length = length;

	int i;
	for(i = 0; i < length; ++i)
		vector->data[i] = 0;

	return vector;
}

vector_int_t* vector_int_newWithData(int *data, int length)
{
	vector_int_t* vector = (vector_int_t*) malloc(sizeof(vector_int_t));
	vector->data = data;
	vector->length = length;

	return vector;
}

void vector_int_delete(vector_int_t* vector)
{
	if(!vector)
		return;

#ifdef __ICC
	mkl_free(vector->data);
	mkl_free(vector);
#else
	free(vector->data);
	free(vector);
#endif

}

void vector_int_deleteMultiple(vector_int_t** vectors, int count)
{
	int i;
	for(i = 0; i < count; ++i)
		vector_int_delete(vectors[i]);

	free(vectors);
}

void vector_int_print(char* m, vector_int_t* vector)
{
	if(!vector)
		return;

	PRINTF("%s(%d):", m, vector->length);

	int i;
	for(i = 0; i < vector->length; ++i)
		PRINTF(" %d", vector->data[i]);

	PRINTF("\n");
}

void vector_int_copyLengthAmount(vector_int_t* dest, DECIMAL* contentsArr)
{
	int i;
	for(i = 0; i < dest->length; ++i)
		dest->data[i] = contentsArr[i];
}

vector_real_t* vector_real_new(int length)
{
	vector_real_t* vector = NULL;
#ifdef __ICC
	vector = (vector_real_t*) mkl_malloc(sizeof(vector_real_t), ALIGNMENT);
	vector->data = (REAL*) mkl_malloc(length * sizeof(REAL), ALIGNMENT);
#else
	posix_memalign(&vector, ALIGNMENT, sizeof(vector_real_t));
	posix_memalign(&vector->data, ALIGNMENT, sizeof(REAL) * length);
#endif

	vector->length = length;

	int i;
	for(i = 0; i < length; ++i)
		vector->data[i] = 0.0;

	return vector;
}

void vector_real_delete(vector_real_t* vector)
{
	if(!vector)
		return;

#ifdef __ICC
	mkl_free(vector->data);
	mkl_free(vector);
#else
	free(vector->data);
	free(vector);
#endif
}

void vector_real_deleteMultiple(vector_real_t** vectors, int count)
{
	int i;
	for(i = 0; i < count; ++i)
		vector_real_delete(vectors[i]);

	free(vectors);
}

vector_real_t* vector_real_random(int length)
{

	vector_real_t* randomVector = vector_real_new(length);
	int i;

	#pragma omp parallel for private(i)
	for(i = 0; i < length; ++i)
		randomVector->data[i] = rand() % 10;

	return randomVector;
}

void vector_real_print(char* m, vector_real_t* vector)
{
	if(!vector)
		return;

	PRINTF("%s(%d):", m, vector->length);

	int i;
	for(i = 0; i < vector->length; ++i)
		PRINTF(" %3.2f", vector->data[i]);

	PRINTF("\n");
}

int vector_real_compare(char* m, vector_real_t* v1, vector_real_t* v2)
{
	int isSame = 1;
	int i;
	int diffStart = 0;
	int diffEnd = 0;
	for(i = 0; i < v1->length; ++i)
	{
		if(((v1->data[i] - v2->data[i]) > 0.1) ||
			((v1->data[i] - v2->data[i]) < -0.1))
		{
			isSame = 0;
			break;
		}
	}

	if(!isSame)
	{
		if(strlen(m) > 0)
			PRINTF("%s...\n", m);
		int len = i + 10;
		for(; i < len && i < v1->length; ++i)
		{
			PRINTF("%d. %f\t%f\n", i, v1->data[i], v2->data[i]);
		}
	}

	return isSame;
}

int vector_real_copy(vector_real_t* source, vector_real_t* dest)
{
	if(source->length != dest->length)
		return FALSE;

	int i;
	for(i = 0; i < dest->length; ++i)
		dest->data[i] = source->data[i];

	return TRUE;
}

float vector_real_mult(vector_real_t* v1, vector_real_t* v2, int numThreads)
{
	REAL result = 0.0;
	int i;
	#pragma omp parallel for num_threads(numThreads)\
							schedule(guided)\
							reduction(+ : result)\
							private(i)
	for(i = 0; i < v1->length; ++i)
		result += v1->data[i] * v2->data[i];

	return result;
}

void vector_real_reset(vector_real_t* v)
{
	int i;
	for(i = 0; i < v->length; ++i)
		v->data[i] = 0.0;
}

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif
