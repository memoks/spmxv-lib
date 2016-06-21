
#ifndef CONVERSION_H_
#define CONVERSION_H_

#include "include/io/input_parser.h"

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include "include/config.h"
#include "include/data_structure/spm_storage.h"
#include "include/scheduler/job_batch.h"


extern void converter_CSR_to_ICSR(
		spm_cmp_t* spmCsr, spm_inc_t** spmIcsr_out);
extern void converter_CSR_to_ICSR_multiple(
		spm_cmp_t** spmCsrs, spm_inc_t*** spmIcsrs_out, int numSpms);
extern void converter_ICSR_to_CSR(
		spm_inc_t* spmIcsr, spm_cmp_t** spmCsr_out);

// quintet conversion functions
// ----------------------------------------------------------------------------------

/**
 * Converts quintets to CSR structure.
 *
 * ASSUMPTION: Quintets must be ordered accordingly (to
 * match with CSR structure) before this functions is called.
 *
 * ASSUMPTION: when reading quintets for the first time
 * their original index must be copied to index values
 * even though there is no ordering scheme used.
 * As a result iVal and jVal must never be empty.
 */
extern void converter_quintetToCSR_alloc(
		quintet_t* quintets, spm_cmp_t** spmCsr_out,
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount);

/**
 * Converts quintets to JDS structure.
 *
 * ASSUMPTION: Quintets must be ordered accordingly (to
 * match with CSR structure) before this functions is called.
 *
 * ASSUMPTION: when reading quintets for the first time
 * their original index must be copied to index values
 * even though there is no ordering scheme used.
 * As a result iVal and jVal must never be empty.
 */
extern void converter_quintetToCSR(
		quintet_t* quintets, spm_cmp_t* spmCsr_inout);

/*
 * Extracts NNZ from CSR ordered quintets of given length.
 * Can execute in parallel fashion.
 *
 * ASSUMPTION: Quintets must be ordered accordingly (to
 * match with CSR structure) before this functions is called.
 *
 * ASSUMPTION: when reading quintets for the first time
 * their original index must be copied to index values
 * even though there is no ordering scheme used.
 * As a result iVal and jVal must never be empty.
 *
 * ASSUMPTION: nnzPerRow_inout array must be initialized
 * to zero before calling this function.
 */
extern void converter_extractCsrNnzPerRow(
		quintet_t* quintets, DECIMAL length,
		DECIMAL* nnzPerRow_inout);

/**
 * Converts quintets to JDS structure.
 *
 * ASSUMPTION: Quintets must be ordered accordingly (to
 * match with CSR structure) before this functions is called.
 *
 * ASSUMPTION: when reading quintets for the first time
 * their original index must be copied to index values
 * even though there is no ordering scheme used.
 * As a result iVal and jVal must never be empty.
 */
// TODO test
extern void converter_quintetToJDS_alloc(
		quintet_t* quintets, spm_jds_t** spmJds_out,
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount,
		DECIMAL idiagLength);

/**
 * Converts quintets to JDS structure.
 *
 * ASSUMPTION: Quintets must be ordered accordingly (to
 * match with CSR structure) before this functions is called.
 * And CSR has to have its rows ordered by NNZs in
 * descending order.
 *
 * ASSUMPTION: when reading quintets for the first time
 * their original index must be copied to index values
 * even though there is no ordering scheme used.
 * As a result iVal and jVal must never be empty.
 */
// TODO test
extern void converter_quintetToJDS(
		quintet_t* quintets, spm_jds_t* spmJds_inout);

/*
 * Extracts NNZ per "jds column" (not normal column).
 * Can execute in parallel fashion.
 *
 * ASSUMPTION: Quintets must be ordered accordingly (to
 * match with CSR structure) before this functions is called.
 *
 * ASSUMPTION: when reading quintets for the first time
 * their original index must be copied to index values
 * even though there is no ordering scheme used.
 * As a result iVal and jVal must never be empty.
 *
 * ASSUMPTION: nnzPerColumn_inout must be initialized
 * to zero before calling this function.
 */
// TODO test
extern void converter_extractJdsNnzPerColumn(
		quintet_t* quintets, DECIMAL length,
		DECIMAL* nnzPerColumn_inout);

/*
 * Used to convert nnz per row to csr row_ptr array.
 */
extern void converter_accumulate(
		DECIMAL* arr_inout, DECIMAL length);

/**
 * Converts quintets to CSR (ordered for JDS) structure.
 *
 * ASSUMPTION: Quintets must be sorted partially for
 * partial JDS structure.
 * Call "converter_quintetToJDSOrderedQuintet" before
 * calling this one.
 *
 * ASSUMPTION: when reading quintets for the first time
 * their original index must be copied to index values
 * even though there is no ordering scheme used.
 * As a result iVal and jVal must never be empty.
 */
extern void converter_quintetToJDSCounterpartCSR(
		quintet_t* quintets,
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount,
		sub_mtx_dim_t* subMtxArr, int subMtxCount,
		spm_cmp_t** spmCsrCounterpart_out,
		DECIMAL** permutation_out,
		vector_int_t** rowOrderLookup_out);

extern void converter_quintetToCSC(
		quintet_t* quintets, spm_cmp_t** spmCsc_out,
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount);

/*
 * ASSUMPTION: when reading quintets for the first time
 * their original index must be copied to index values
 * even though there is no ordering scheme used.
 * As a result iVal and jVal must never be empty.
 */
// TODO needs to be gone or changed at least
extern void converter_quintetToJDSOrderedQuintet(
		quintet_t* quintets_inout,
		DECIMAL nnz, DECIMAL rowCount, DECIMAL colCount,
		sub_mtx_dim_t* subMtxArr, int subMtxCount,
		vector_int_t** rowOrderLookup_out);

/**
 * Converts given CSR structure to JDS structure.
 */
extern void converter_CSRToJDS(
		spm_cmp_t* spmCsr, DECIMAL* permutation,
		spm_jds_t** spmJds_out);

// TODO needs to be gone / changed
extern void converter_CSRPartToJDS(
		spm_cmp_t* spmCsr, sub_mtx_dim_t* subMtx,
		DECIMAL* permutation, spm_jds_t** spmJds_out);

// TODO needs to be gone / changed
extern void converter_extractPartialJDSArr(
		spm_cmp_t* spmCsrCounterpart,
		job_batch_t* subMtxArr, int subMtxCount,
		DECIMAL* permutation,
		job_batch_t** partialJdsBatchArr_out);

// TODO needs to be gone / changed
extern void converter_extractPartialJDSArrMultiple(
		spm_cmp_t* spmCsrCounterpart, int numBlocks,
		job_batch_t** subMtxArrPerBlock, int* subMtxCountPerBlock,
		DECIMAL* permutation,
		job_batch_t*** partialJdsBatchArrPerBlock_out);

// TODO needs to be gone / changed
extern void converter_extractHybridJDSCSR(
		job_batch_t* partialJdsBatch,
		job_batch_t* hybridJdsCsr_inout);

/**
 * ASSUMPTION: Given batch in partial_JDS format,
 * must have its CSR counterpart in it for this function
 * to work.
 */
// TODO needs to be gone / changed
extern void converter_extractHybridJDSCSR_JDS(
		job_batch_t* partialJdsBatch,
		job_batch_t* hybridJdsCsr_inout);

/**
 * ASSUMPTION: Given batch in partial_JDS format,
 * must have its CSR counterpart in it for this function
 * to work.
 */
// TODO needs to be gone / changed
extern void converter_extractHybridJDSCSR_CSR(
		job_batch_t* partialJdsBatch,
		job_batch_t* hybridJdsCsr_inout);

// TODO needs to be gone / changed
extern void converter_extractHybridJDSCSRArr(
		job_batch_t* partialJdsBatchArr, int length,
		job_batch_t** hybridJdsCsrArr_out);

// TODO needs to be gone / changed
extern void converter_extractHybridJDSCSRArrMultiple(
		int numBlocks, job_batch_t** partialJdsBatchArrPerBlock,
		int* partialJdsBatchCountPerBlock,
		job_batch_t*** hybridJdsCsrArrPerBlock_out);

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* CONVERSION_H_ */
