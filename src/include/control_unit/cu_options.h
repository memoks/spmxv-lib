
#ifndef SRC_INCLUDE_CONTROL_UNIT_CU_OPTIONS_H_
#define SRC_INCLUDE_CONTROL_UNIT_CU_OPTIONS_H_


// Sparse Matrix Storage Formats
// ----------------------------------------------------------------------
// There are 3 storage formats for which data preparation and SpMxV vary.
//  - CSR
//  - JDS
//  - Hybrid JDS-CSR

#define SPM_STORAGE_NONE 0
#define SPM_STORAGE_CSR 1
#define SPM_STORAGE_JDS 2
#define SPM_STORAGE_HYBRID_JDS_CSR 3
#define SPM_STORAGE_ICSC 4
#define SPM_STORAGE_CSC 5
#define SPM_STORAGE_ICSR 6

// Work stealing schemes
// -----------------------------------------------------------------------------
// There are 2 work stealing schemes for DWS
//  - tree: where each block uses the bipartitioning tree to decide the
//    next victim
//  - ring: where each block uses ring structure to choose the next victim
#define STEALING_SCHEME_RING 1
#define STEALING_SCHEME_TREE 2

#define DEFAULT_STEALING_SCHEME STEALING_SCHEME_RING

// Stealing Schemes and parameters
// -----------------------------------------------------------------------------
#define DEFAULT_STEAL_TRESHOLD 1


// Ordering type
// -----------------------------------------------------------------------------
// There are 3 ordering models
//  - unordered
//  - colNet
//  - rowNet (not implemented)
#define ORDERING_TYPE_NONE 0
#define ORDERING_TYPE_COLUMN_NET 1
#define ORDERING_TYPE_ROW_NET 2


// Partitioning type
// -----------------------------------------------------------------------------
// There are 3 partitioning types
//  - 1D row slice
//  - 1D column slice
//  - 2D checker-board
#define PARTITION_TYPE_1D_ROW_SLICE 1
#define PARTITION_TYPE_1D_COLUMN_SLICE 2
#define PARTITION_TYPE_2D_CHECKER_BOARD 3

#define DEFAULT_PARTITION_TYPE PARTITION_TYPE_1D_ROW_SLICE


// Partitioning method
// -----------------------------------------------------------------------------
// There are 2 ways to partition matrices into smaller sub-matrices of
// targeted size
//  - REGULAR: Starting from the first row / column, each consecutive row
//    is added up while being traversed. When their size exceeds targeted
//    size then those rows are grouped in a sub-matrix.
//  - RECURSIVE-BIPARTITION: Matrix is divided into 2 (row / column wise)
//    at every step until it drops below targeted size. Then sub-matrices
//    are created from each of those matrix parts.
#define PARTITION_METHOD_REGULAR 1
#define PARTITION_METHOD_RECURSIVE_BIPARTITION 2

#define DEFAULT_PARTITION_METHOD PARTITION_METHOD_REGULAR

// Variants of static algorithm
// -----------------------------------------------------------------------------
// There are 2 variants of static algorithm
//  - STATIC (original): where each execution-context has its own block
//  - STATIC-Scatter: the difference from the original is initial partition of
//    job-batches
// Depending on the chosen type, data-preparation routines vary.
#define STATIC_ALGORITHM_TYPE_DEFAULT 1
#define STATIC_ALGORITHM_TYPE_SCATTER 2


// Variants of DWS algorithm
// -----------------------------------------------------------------------------
// There are 3 variants of DWS algorithm
//  - DWS (original): where each execution-context has its own block
//  - DWS-Scatter: the difference from the original is initial partition of
//    job-batches
//  - DWS-Shared-Block n: where n execution context shares a block together
// Depending on the chosen type, data-preparation and SpMxV routines vary.
#define DWS_ALGORITHM_TYPE_DEFAULT 1
#define DWS_ALGORITHM_TYPE_SCATTER 2
#define DWS_ALGORITHM_TYPE_SHARED_BLOCK 3



#endif /* SRC_INCLUDE_CONTROL_UNIT_CU_OPTIONS_H_ */
