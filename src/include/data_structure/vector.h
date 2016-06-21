
#ifndef VECTOR_H_
#define VECTOR_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"

struct vector_int
{
	int length;
	DECIMAL* data;
};

typedef struct vector_int vector_int_t;

extern vector_int_t* vector_int_new(int length);
extern vector_int_t* vector_int_newWithData(int *data, int length);
extern void vector_int_delete(vector_int_t* vector);
extern void vector_int_deleteMultiple(vector_int_t** vectors, int count);
extern void vector_int_print(char* m, vector_int_t* vector);
extern void vector_int_copyLengthAmount(vector_int_t* dest, DECIMAL* contentsArr);


struct vector_real
{
	int length;
	REAL* data;
};

typedef struct vector_real vector_real_t;


extern vector_real_t* vector_real_new(int length);
extern void vector_real_delete(vector_real_t* vector);
extern void vector_real_deleteMultiple(vector_real_t** vectors, int count);
extern vector_real_t* vector_real_random(int length);
extern void vector_real_print(char* m, vector_real_t* vector);
extern int vector_real_compare(char* m, vector_real_t* v1, vector_real_t* v2);
extern int vector_real_copy(vector_real_t* source, vector_real_t* dest);
extern float vector_real_mult(vector_real_t* v1, vector_real_t* v2, int numThreads);
extern void vector_real_reset(vector_real_t* v);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* VECTOR_H */
