
#ifndef QUINTET_H_
#define QUINTET_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"
#include "include/data_structure/vector.h"
#include "include/data_structure/sub_mtx.h"

// jds_row
// ------------------------------------------------------------------------------
/**
 * Structure that holds NNZ for each row of matrix. This is used
 * to convert a matrix from CSR to JDS format.
 */
struct jds_row
{
	DECIMAL nnz;
	DECIMAL rowId;
};

typedef struct jds_row jds_row_t;

/**
 * Compares two jds_rows. However, return value is the inverse of what
 * you would expect from a comparison function (in order to sort rows
 * in descending order).
 *
 * if left.nnz > right.nnz then return -1
 * if left.nnz < right.nnz then return 1
 * return 0 otherwise
 */
extern int jds_row_cmpInverse(
		const void* left, const void* right);

/**
 * Sorts jds_rows by their NNZ in descending order.
 */
extern void jds_row_sortDescending(
		jds_row_t* jdsRows, DECIMAL length);

extern void jds_row_accumulate(
		jds_row_t* jdsRows, DECIMAL length);

extern void jds_row_deaccumulate(
		jds_row_t* jdsRows, DECIMAL length);

extern void jds_row_printArr(
		jds_row_t* jdsRows, DECIMAL length);

// quintet
// ------------------------------------------------------------------------------

/**
 * Structure to represent each matrix entry
 */
struct quintet
{
	DECIMAL i; 		// row index when ordered
	DECIMAL iVal; 	// original row index
	DECIMAL j;		// column index when ordered
	DECIMAL jVal;	// original column index
	REAL val;		// entry value
};

typedef struct quintet quintet_t;

/**
 * Compare functions to be later used in sorting of quintets
 */
extern int quintet_cmpColIndex(const void* left, const void* right); // CSC
extern int quintet_cmpColValue(const void* left, const void* right); // CSC
extern int quintet_cmpRowIndex(const void* left, const void* right); // CSR
extern int quintet_cmpRowValue(const void* left, const void* right); // CSR

/**
 * Sorting functions for quintets, each has its own compare function.
 */
extern void quintet_sortColIndex(quintet_t* quintets, DECIMAL length); // CSC
extern void quintet_sortColValue(quintet_t* quintets, DECIMAL length); // CSC
extern void quintet_sortRowIndex(quintet_t* quintets, DECIMAL length); // CSR
extern void quintet_sortRowValue(quintet_t* quintets, DECIMAL length); // CSR
extern void quintet_sortCustom(
		quintet_t* quintets, DECIMAL length,
		int (*input_cmpTripletFunc) (const void*, const void*));

/**
 * Overwrites row & column ordering index of given quintets using
 * their ordering value.
 */
extern void quintet_overwriteIndex(
		quintet_t* quintets_inout, DECIMAL nnz);

/**
 * Generates only rowPtr array structure from quintets.
 *
 * ASSUMPTION: quintets must be ordered for CSR before
 * calling this function.
 */
extern void quintet_generateCSRRowPtr(
		quintet_t* quintets, DECIMAL length,
		DECIMAL** rowPtr_out, DECIMAL* rowPtrLength_out);

/**
 * Generates an array whose each element holds the # of
 * non-zero values in one of the rows of given quintets.
 *
 * @quintets: sparse-matrix elements.
 * @length: length of quintet array.
 * @rowCount: # of row in given quintet array.
 * @elemCountPerRowArr_out: (return value) array of numbers
 *  where each number represents nnz of a particular row.
 */
extern void quintet_countElementsPerRow(
		quintet_t* quintets, DECIMAL length, DECIMAL rowCount,
		DECIMAL** elemCountPerRowArr_out);

/**
 * For binary matrices, instead of making value of every nnz 1,
 * this function calculates values of each element in the same
 * row according to this formula:
 * (nnz_count_this_row / col_count)
 *
 * Input;
 * 1.0 0.0 0.0 1.0 0.0
 * 0.0 1.0 1.0 0.0 1.0
 * 1.0 1.0 0.0 1.0 1.0
 * 0.0 0.0 0.0 0.0 0.0
 * 1.0 0.0 0.0 0.0 0.0
 *
 * Output:
 * 0.4 0.0 0.0 0.4 0.0
 * 0.0 0.6 0.6 0.0 0.6
 * 0.8 0.8 0.0 0.8 0.8
 * 0.0 0.0 0.0 0.0 0.0
 * 0.2 0.0 0.0 0.0 0.0
 */
extern void quintet_adjustBinaryMatrixCoefficients(
		quintet_t* quintets_inout, int numTriplets);

/**
 * Transposes quintets, row count, column count.
 */
extern void quintet_transpose(quintet_t* quintets_inout, DECIMAL nnz,
		DECIMAL* rowCount_inout, DECIMAL* colCount_inout);


/**
 * Injects row ordering info to each of the quintets.
 * It is assumed that quintets are ordered by their
 * row index beforehand.
 */
extern void quintet_addOrderingInfoByIndex(
		quintet_t* quintets_inout, DECIMAL nnz,
		vector_int_t* rowOrderLookup,
		vector_int_t* columnOrderLookup);

/**
 * Injects row ordering info to each of the quintets.
 * It is assumed that quintets are ordered by their
 * row value beforehand.
 */
extern void quintet_addRowOrderingInfoByValue(
		quintet_t* quintets_inout, DECIMAL nnz,
		vector_int_t* rowOrderLookup);

/**
 * Copies row & column indices of all non-zeros on to
 * row & column values respectively.
 */
extern void quintet_copyDefaultOrderingInfoToValue(
		quintet_t* quintets_inout, DECIMAL nnz);

/*
 * Just a console print method for debugging purposes.
 */
extern void quintet_printQuintets(
		char* m, quintet_t* quintets, int length);

extern void quintet_copyQuintet(
		quintet_t* source, quintet_t* destination);

extern void quintet_copyQuintets(
		quintet_t* source, quintet_t* destination, DECIMAL length);

// TODO test & comment
extern void quintet_orderForJDSCounterpart(
		quintet_t** quintets_inout,
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount,
		sub_mtx_dim_t* subMtxArr, int subMtxArrLength,
		vector_int_t** rowOrderLookup_out,
		DECIMAL** elemCountPerRow_out);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* QUINTET_H_ */
