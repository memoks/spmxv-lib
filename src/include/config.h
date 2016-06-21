
#ifndef CONFIG_H_
#define CONFIG_H_

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>


// Some useful macros to help readability
// ----------------------------------------------------------------------

#define FALSE 0
#define TRUE 1
#define DEFAULT_STR_BUFF_SIZE 2048

// Some headers for Intel math kernel library
// ----------------------------------------------------------------------

// Some data types and alignment rules for flexibility (float, double, etc...)
// ----------------------------------------------------------------------

// Single Precision or not (double precision otherwise)
// #define SINGLE_PRECISION

// #define SINGLE_PRECISION

#ifdef SINGLE_PRECISION
#define REAL float
#define POW_FUNC_PREC(base, power) powf(base, power);
#define DEFAULT_SIMD_LENGTH 16
#define DEFAULT_ROW_SLICE_LENGTH SIMD_LENGTH
#else
#define REAL double
#define POW_FUNC_PREC(base, power) powl(base, power);
#define DEFAULT_SIMD_LENGTH 8
#define DEFAULT_ROW_SLICE_LENGTH SIMD_LENGTH
#endif

#define DECIMAL int

// Cache alignment (32 or 64 byte)
#define ALIGN_32 32
#define ALIGN_64 64

#define ALIGNMENT ALIGN_64

// Compilation flag that is passed by makefile (for code not to break)
// ----------------------------------------------------------------------

#ifdef STD_C99
#define RESTRICT restrict
#define REGISTER register
#else
#define RESTRICT
#define REGISTER
#endif

// Some dynamic pragmas to implement routines
// using macros instead of functions.
// Although not used, let's keep them just in case
// ----------------------------------------------------------------------

#define OMP_PARAM_(...) #__VA_ARGS__
#define OMP_PARAM(...) OMP_PARAM_(__VA_ARGS__)
#define __OMP__ omp
#define PRAGMA_OMP_PARALLEL(...) _Pragma(OMP_PARAM(__OMP__ parallel __VA_ARGS__))
#define PRAGMA_OMP_PARALLEL_FOR(...) _Pragma(OMP_PARAM(__OMP__ parallel for __VA_ARGS__))
#define PRAGMA_OMP_ATOMIC _Pragma("omp atomic")
#define PRAGMA_IVDEP _Pragma("ivdep")
#define PRAGMA_SIMD _Pragma("simd")

// Conjugate Gradient API
// ----------------------------------------------------------------------

#define RESULT_ERROR_TRESHOLD 0.01
#define ITER_MAX 10000

// Default values to parameters of omp-dynanmic scheduling algorithms
// ----------------------------------------------------------------------

#define DEFAULT_CHUNK_SIZE 64

// Job Queue states
// ----------------------------------------------------------------------

#define JOB_QUEUE_STATE_DONE 0
#define JOB_QUEUE_STATE_EMPTY 1
#define JOB_QUEUE_STATE_STEALING 2
#define JOB_QUEUE_STATE_INITIAL 3

// How many times will kernel code run
// ----------------------------------------------------------------------

#define DEFAULT_WARM_UP_RUN_COUNT 10
#define DEFAULT_MEASURE_RUN_COUNT 100


// TODO about hybrid data extraction methods, design a better way,
// use parameter or something instead of macro and delete all these macros

// Hybrid Data Extraction Method
#define JDS_BASED_HYBRID_DATA_EXTRACTION
// If not defined that it is CSR_BASED_HYBRID_DATA_EXTRACTION
#define TRUE_HYBRID_FORM

// Program mode
// ----------------------------------------------------------------------

#define NO_OUTPUT_MODE 0 // Does not really print anything
#define MEASURE_MODE 1 // Prints timing info
#define DEBUG_MODE 2 // Prints enough messages to watch overall process
#define TRACE_MODE 3 // Prints literally everything, useful for debug

// Program mode is always passed through makefile
// This is just a backup otherwise
#ifndef PROGRAM_MODE
#define PROGRAM_MODE TRACE_MODE // Actual flag that holds program mode
#endif

// A Macro for printing toString returnee
// ----------------------------------------------------------------------

#define PRINTF(format, ...)					\
	do{										\
		printf(format, ##__VA_ARGS__);		\
		fflush(0);							\
	} while(0)

#define PRINT_BUFF(format, ...)				\
	do{										\
		printf(format, ##__VA_ARGS__);		\
		fflush(0);							\
		free(buff);							\
	} while(0)

#define DEBUG(format, ...)					\
	do{										\
		if(PROGRAM_MODE >= DEBUG_MODE)		\
		{									\
			printf(format, ##__VA_ARGS__);	\
			fflush(0);						\
		}									\
	}										\
	while(0)

#define DEBUG_BUFF(format, ...)				\
	do{										\
		if(PROGRAM_MODE >= DEBUG_MODE)		\
		{									\
			printf(format, ##__VA_ARGS__);	\
			fflush(0);						\
			free(buff);						\
		}									\
	}										\
	while(0)

// ----------------------------------------------------------------------

#ifdef ARCH_MIC_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif /* CONFIG_H_ */
