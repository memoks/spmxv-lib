
#ifndef INPUT_PARSER_H_
#define INPUT_PARSER_H_

#include "include/io/mmio.h"
#include "include/data_structure/quintet.h"
#include "include/data_structure/vector.h"

/**
 * Prints given quintet array to file stream in Matrix Market Format
 */
extern void input_printQuintetsToFileInMMF(
		FILE *f, quintet_t* quintets,
		DECIMAL rowCount, DECIMAL colCount, DECIMAL nnz, int isBinary);

/**
 * Reads a matrix file written market matrix format and returns each
 * entry as quintet.
 * ASSUMPTION: By default, normal ordering of matrix entries are
 * copied to their manual order counterparts.
 */
extern void input_readQuintets(
		char* mmfPath,
		DECIMAL *nnz_out, DECIMAL* rowCount_out, DECIMAL* colCount_out,
		quintet_t** quintets_out, mm_info_t** mmInfo_out);

/**
 * Reads row & column ordering vector and injects this info into quintets.
 * @rowOrderLookup_out: a vector where (i-1)th element indicate ith row's
 * row number in the ordered format.
 * @colOrderLookup_out: a vector where (i-1)th element indicate ith column's
 * column number in the ordered format.
 */
extern void input_readOrderLookup(
		char* fileName, quintet_t* quintets_in_out,
		DECIMAL nnz, DECIMAL rowCount, DECIMAL columnCount,
		vector_int_t** rowOrderLookup_out,
		vector_int_t** colOrderLookup_out);

/**
 * Calculates and returns the line count in a given file.
 */
extern int input_getLineCount(char* fileName);

extern void input_readVector(
		char* fileName, vector_int_t* vector_inout);

/**
 * Just a wrapper function for reading & ordering matrix,
 * using PATOH's output.
 */
extern void input_kPatohOrder(
		char* mmfPath, char* mmfDirPath,
		char* mmfName, int cacheSizeKB,
		char* partitioningTypeStr,
		quintet_t* quintets, DECIMAL nnz,
		DECIMAL rowCount, DECIMAL colCount,
		vector_int_t** rowOrderLookup_out,
		vector_int_t** rowPartitioning_out,
		vector_int_t** colOrderLookup_out,
		vector_int_t** colPartitioning_out);

// TODO test
extern void input_powerOf2Order(
		char* mmfDirPath, char* mmfName, int cacheSizeKB,
		quintet_t* quintets_inout,
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount,
		vector_int_t** rowOrderLookup_out,
		vector_int_t** rowPartitioning_out,
		vector_int_t** colOrderLookup_out,
		vector_int_t** colPartitioning_out);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

/**
 * Orders a vector back by a given order lookup.
 * This is called on vector x (input) to order it
 * "if ordering is used".
 */
extern void input_orderVector(
		vector_real_t** inputAddress_inout, vector_int_t* orderLookup);

/**
 * Reorders a vector back to its original order (specified bu orderLookup).
 * This is called on vector y (results) to reorder it to its
 * original order "if ordering is used".
 */
extern void input_reorderVector(
		vector_real_t** inputAddress_inout, vector_int_t* orderLookup);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif


/**
 * Reads a vector of floating point numbers from a given file path.
 * It should be written like this (where x is the length of the vector)
 * x
 * 1.0
 * 2.0
 * 3.0
 * 4.0
 * 5.0
 * ...
 */
extern vector_real_t* input_readVectorRealFromFile(char* filePath);

/**
 * Reads a vector of integer numbers from a given file path.
 * It should be written like this (where x is the length of the vector)
 * x
 * 1
 * 2
 * 3
 * 4
 * 5
 * ...
 */
extern vector_int_t* input_readVectorIntFromFile(char* filePath);

/**
 * Reads partition vector from given file path.
 *
 * Suppose file has (first line is length);
 * 5
 * 4
 * 7
 * 5
 * 9
 * 5
 *
 * Then returned vector will be;
 * {0, 4, 11, 16, 25, 30}
 */
extern vector_int_t* input_readPartitionVectorFromFile(char* filePath);

/*
 * Dumps a vector of floating point numbers in a given file path.
 */
extern void input_printVectorRealToFile(char* filePath, vector_real_t* v);

/*
 * Dumps a vector of decimal numbers in a given file path.
 */
extern void input_printVectorIntToFile(char* filePath, vector_int_t* v);

extern void input_readPartitionStatFile(char* filePath, int* partCount_out);

extern vector_int_t* input_readPartVector(char* filePath);

#endif // INPUT_PARSER_H_
