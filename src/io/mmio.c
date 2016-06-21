/* 
*   Matrix Market I/O library for ANSI C
*
*   See http://math.nist.gov/MatrixMarket for details.
*
*
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "include/config.h"
#include "include/io/mmio.h"

int mm_read(FILE* f,
		DECIMAL* rowCount_out, DECIMAL* colCount_out, DECIMAL *nnz_out,
		REAL** values_out, DECIMAL** rowInds_out, DECIMAL** colInds_out,
		mm_info_t** mmInfo_out)
{
    mm_info_t* mmInfo = (mm_info_t*) malloc(sizeof(mm_info_t));
    mm_info_initDefault(mmInfo);
    MM_typecode matcode;
    DECIMAL rowCount = 0;
    DECIMAL colCount = 0;
    DECIMAL nnz = 0;
    DECIMAL* rowIndices = NULL;
    DECIMAL* colIndices = NULL;
    REAL* values = NULL;

    if(f == NULL)
    	return -1;

    if(mm_read_banner(f, &matcode) != 0)
    {
        printf("mm_read_unsymetric: Could not process Matrix Market banner.\n");
        return -1;
    }

    mm_info_init(mmInfo, &matcode);

    if(!(mm_is_matrix(matcode) && mm_is_sparse(matcode)))
    {
        fprintf(stderr, "Sorry, this application does not support ");
        fprintf(stderr, "Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        return -1;
    }

    if (mm_read_mtx_crd_size(f, &rowCount, &colCount, &nnz) != 0)
    {
        fprintf(stderr, "mm_read(): could not parse matrix size.\n");
        return -1;
    }

    int tempNnz = nnz;
    if(mmInfo->isSymmetric)
    	nnz = nnz * 2 - rowCount;

    rowIndices = (DECIMAL*) malloc(sizeof(DECIMAL) * nnz);
    colIndices = (DECIMAL*) malloc(sizeof(int) * nnz);
    values = (REAL*) malloc(sizeof(REAL) * nnz);

    DECIMAL index = 0;
    if(mmInfo->isPattern)
    {
    	for(index = 0; index < tempNnz; index++)
		{
			fscanf(f, "%d %d\n", &rowIndices[index], &colIndices[index]);
			rowIndices[index]--;  /* adjust from 1-based to 0-based */
			colIndices[index]--;
			values[index] = 1.0;
		}
    }
    else
    {
    	for(index = 0; index < tempNnz; index++)
		{

#ifdef SINGLE_PRECISION
    		fscanf(f, "%d %d %f\n", &rowIndices[index], &colIndices[index], &values[index]);
#else
    		fscanf(f, "%d %d %lf\n", &rowIndices[index], &colIndices[index], &values[index]);
#endif

			rowIndices[index]--;  /* adjust from 1-based to 0-based */
			colIndices[index]--;
		}
    }

    if(mmInfo->isSymmetric)
    {
    	DECIMAL k = 0;
		while(index < nnz)
		{
			if(rowIndices[k] != colIndices[k])
			{
				rowIndices[index] = colIndices[k];
				colIndices[index] = rowIndices[k];
				values[index] = values[k];
				++index;
			}

			++k;
		}
    }

    // return values
    *mmInfo_out = mmInfo;
    *rowCount_out = rowCount;
	*colCount_out = colCount;
	*nnz_out = nnz;
    *rowInds_out = rowIndices;
    *colInds_out = colIndices;
    *values_out = values;

    mmInfo->nnz = nnz;
    mmInfo->rowCount = rowCount;
    mmInfo->columnCount = colCount;

    return 0;
}

int mm_is_valid(MM_typecode matcode)
{
    if (!mm_is_matrix(matcode)) return 0;
    if (mm_is_dense(matcode) && mm_is_pattern(matcode)) return 0;
    if (mm_is_real(matcode) && mm_is_hermitian(matcode)) return 0;
    if (mm_is_pattern(matcode) && (mm_is_hermitian(matcode) || 
                mm_is_skew(matcode))) return 0;
    return 1;
}

int mm_read_banner(FILE *f, MM_typecode *matcode)
{
    char line[MM_MAX_LINE_LENGTH];
    char banner[MM_MAX_TOKEN_LENGTH];
    char mtx[MM_MAX_TOKEN_LENGTH]; 
    char crd[MM_MAX_TOKEN_LENGTH];
    char data_type[MM_MAX_TOKEN_LENGTH];
    char storage_scheme[MM_MAX_TOKEN_LENGTH];
    char *p;


    mm_clear_typecode(matcode);  

    if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL) 
        return MM_PREMATURE_EOF;

    if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type, 
        storage_scheme) != 5)
        return MM_PREMATURE_EOF;

    for (p=mtx; *p!='\0'; *p=tolower(*p),p++);  /* convert to lower case */
    for (p=crd; *p!='\0'; *p=tolower(*p),p++);  
    for (p=data_type; *p!='\0'; *p=tolower(*p),p++);
    for (p=storage_scheme; *p!='\0'; *p=tolower(*p),p++);

    /* check for banner */
    if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
        return MM_NO_HEADER;

    /* first field should be "mtx" */
    if (strcmp(mtx, MM_MTX_STR) != 0)
        return  MM_UNSUPPORTED_TYPE;
    mm_set_matrix(matcode);


    /* second field describes whether this is a sparse matrix (in coordinate
            storgae) or a dense array */


    if (strcmp(crd, MM_SPARSE_STR) == 0)
        mm_set_sparse(matcode);
    else
    if (strcmp(crd, MM_DENSE_STR) == 0)
            mm_set_dense(matcode);
    else
        return MM_UNSUPPORTED_TYPE;
    

    /* third field */

    if (strcmp(data_type, MM_REAL_STR) == 0)
        mm_set_real(matcode);
    else
    if (strcmp(data_type, MM_COMPLEX_STR) == 0)
        mm_set_complex(matcode);
    else
    if (strcmp(data_type, MM_PATTERN_STR) == 0)
        mm_set_pattern(matcode);
    else
    if (strcmp(data_type, MM_INT_STR) == 0)
        mm_set_integer(matcode);
    else
        return MM_UNSUPPORTED_TYPE;
    

    /* fourth field */

    if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
        mm_set_general(matcode);
    else
    if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
        mm_set_symmetric(matcode);
    else
    if (strcmp(storage_scheme, MM_HERM_STR) == 0)
        mm_set_hermitian(matcode);
    else
    if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
        mm_set_skew(matcode);
    else
        return MM_UNSUPPORTED_TYPE;
        

    return 0;
}

int mm_write_mtx_crd_size(FILE *f, int M, int N, int nz)
{
    if (fprintf(f, "%d %d %d\n", M, N, nz) != 3)
        return MM_COULD_NOT_WRITE_FILE;
    else 
        return 0;
}

int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz )
{
    char line[MM_MAX_LINE_LENGTH];
    int num_items_read;

    /* set return null parameter values, in case we exit with errors */
    *M = *N = *nz = 0;

    /* now continue scanning until you reach the end-of-comments */
    do 
    {
        if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL) 
            return MM_PREMATURE_EOF;
    }while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
    if (sscanf(line, "%d %d %d", M, N, nz) == 3)
        return 0;
        
    else
    do
    { 
        num_items_read = fscanf(f, "%d %d %d", M, N, nz); 
        if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 3);

    return 0;
}


int mm_read_mtx_array_size(FILE *f, int *M, int *N)
{
    char line[MM_MAX_LINE_LENGTH];
    int num_items_read;
    /* set return null parameter values, in case we exit with errors */
    *M = *N = 0;
	
    /* now continue scanning until you reach the end-of-comments */
    do 
    {
        if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL) 
            return MM_PREMATURE_EOF;
    }while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
    if (sscanf(line, "%d %d", M, N) == 2)
        return 0;
        
    else /* we have a blank line */
    do
    { 
        num_items_read = fscanf(f, "%d %d", M, N); 
        if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 2);

    return 0;
}

int mm_write_mtx_array_size(FILE *f, int M, int N)
{
    if (fprintf(f, "%d %d\n", M, N) != 2)
        return MM_COULD_NOT_WRITE_FILE;
    else 
        return 0;
}



/*-------------------------------------------------------------------------*/

/******************************************************************/
/* use when I[], J[], and val[]J, and val[] are already allocated */
/******************************************************************/

int mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, int I[], int J[],
        double val[], MM_typecode matcode)
{
    int i;
    if (mm_is_complex(matcode))
    {
        for (i=0; i<nz; i++)
            if (fscanf(f, "%d %d %lg %lg", &I[i], &J[i], &val[2*i], &val[2*i+1])
                != 4) return MM_PREMATURE_EOF;
    }
    else if (mm_is_real(matcode))
    {
        for (i=0; i<nz; i++)
        {
            if (fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i])
                != 3) return MM_PREMATURE_EOF;

        }
    }

    else if (mm_is_pattern(matcode))
    {
        for (i=0; i<nz; i++)
            if (fscanf(f, "%d %d", &I[i], &J[i])
                != 2) return MM_PREMATURE_EOF;
    }
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;
        
}

int mm_read_mtx_crd_entry(FILE *f, int *I, int *J,
        double *real, double *imag, MM_typecode matcode)
{
    if (mm_is_complex(matcode))
    {
            if (fscanf(f, "%d %d %lg %lg", I, J, real, imag)
                != 4) return MM_PREMATURE_EOF;
    }
    else if (mm_is_real(matcode))
    {
            if (fscanf(f, "%d %d %lg\n", I, J, real)
                != 3) return MM_PREMATURE_EOF;

    }

    else if (mm_is_pattern(matcode))
    {
            if (fscanf(f, "%d %d", I, J) != 2) return MM_PREMATURE_EOF;
    }
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;
        
}


/************************************************************************
    mm_read_mtx_crd()  fills M, N, nz, array of values, and return
                        type code, e.g. 'MCRS'

                        if matrix is complex, values[] is of size 2*nz,
                            (nz pairs of real/imaginary values)
************************************************************************/

int mm_read_mtx_crd(char *fname, int *M, int *N, int *nz, int **I, int **J, 
        double **val, MM_typecode *matcode)
{
    int ret_code;
    FILE *f;

    if (strcmp(fname, "stdin") == 0) f=stdin;
    else
    if ((f = fopen(fname, "r")) == NULL)
        return MM_COULD_NOT_READ_FILE;


    if ((ret_code = mm_read_banner(f, matcode)) != 0)
        return ret_code;

    if (!(mm_is_valid(*matcode) && mm_is_sparse(*matcode) && 
            mm_is_matrix(*matcode)))
        return MM_UNSUPPORTED_TYPE;

    if ((ret_code = mm_read_mtx_crd_size(f, M, N, nz)) != 0)
        return ret_code;


    *I = (int *)  malloc(*nz * sizeof(int));
    *J = (int *)  malloc(*nz * sizeof(int));
    *val = NULL;

    if (mm_is_complex(*matcode))
    {
        *val = (double *) malloc(*nz * 2 * sizeof(double));
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val, 
                *matcode);
        if (ret_code != 0) return ret_code;
    }
    else if (mm_is_real(*matcode))
    {
        *val = (double *) malloc(*nz * sizeof(double));
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val, 
                *matcode);
        if (ret_code != 0) return ret_code;
    }

    else if (mm_is_pattern(*matcode))
    {
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val, 
                *matcode);
        if (ret_code != 0) return ret_code;
    }

    if (f != stdin) fclose(f);
    return 0;
}

int mm_write_banner(FILE *f, MM_typecode matcode)
{
    char *str = mm_typecode_to_str(matcode);
    int ret_code;

    ret_code = fprintf(f, "%s %s\n", MatrixMarketBanner, str);
    free(str);
    if (ret_code !=2 )
        return MM_COULD_NOT_WRITE_FILE;
    else
        return 0;
}

int mm_write_mtx_crd(char fname[], int M, int N, int nz, int I[], int J[],
        double val[], MM_typecode matcode)
{
    FILE *f;
    int i;

    if (strcmp(fname, "stdout") == 0) 
        f = stdout;
    else
    if ((f = fopen(fname, "w")) == NULL)
        return MM_COULD_NOT_WRITE_FILE;
    
    /* print banner followed by typecode */
    fprintf(f, "%s ", MatrixMarketBanner);
    fprintf(f, "%s\n", mm_typecode_to_str(matcode));

    /* print matrix sizes and nonzeros */
    fprintf(f, "%d %d %d\n", M, N, nz);

    /* print values */
    if (mm_is_pattern(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d\n", I[i], J[i]);
    else
    if (mm_is_real(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d %20.16g\n", I[i], J[i], val[i]);
    else
    if (mm_is_complex(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d %20.16g %20.16g\n", I[i], J[i], val[2*i], 
                        val[2*i+1]);
    else
    {
        if (f != stdout) fclose(f);
        return MM_UNSUPPORTED_TYPE;
    }

    if (f !=stdout) fclose(f);

    return 0;
}
  

/**
*  Create a new copy of a string s.  mm_strdup() is a common routine, but
*  not part of ANSI C, so it is included here.  Used by mm_typecode_to_str().
*
*/
char *mm_strdup(const char *s)
{
	int len = strlen(s);
	char *s2 = (char *) malloc((len+1)*sizeof(char));
	return strcpy(s2, s);
}

char  *mm_typecode_to_str(MM_typecode matcode)
{
    char buffer[MM_MAX_LINE_LENGTH];
    char *types[4];
	char *mm_strdup(const char *);
    int error =0;

    /* check for MTX type */
    if (mm_is_matrix(matcode)) 
        types[0] = MM_MTX_STR;
    else
        error=1;

    /* check for CRD or ARR matrix */
    if (mm_is_sparse(matcode))
        types[1] = MM_SPARSE_STR;
    else
    if (mm_is_dense(matcode))
        types[1] = MM_DENSE_STR;
    else
        return NULL;

    /* check for element data type */
    if (mm_is_real(matcode))
        types[2] = MM_REAL_STR;
    else
    if (mm_is_complex(matcode))
        types[2] = MM_COMPLEX_STR;
    else
    if (mm_is_pattern(matcode))
        types[2] = MM_PATTERN_STR;
    else
    if (mm_is_integer(matcode))
        types[2] = MM_INT_STR;
    else
        return NULL;


    /* check for symmetry type */
    if (mm_is_general(matcode))
        types[3] = MM_GENERAL_STR;
    else
    if (mm_is_symmetric(matcode))
        types[3] = MM_SYMM_STR;
    else 
    if (mm_is_hermitian(matcode))
        types[3] = MM_HERM_STR;
    else 
    if (mm_is_skew(matcode))
        types[3] = MM_SKEW_STR;
    else
        return NULL;

    sprintf(buffer,"%s %s %s %s", types[0], types[1], types[2], types[3]);
    return mm_strdup(buffer);

}

// structure to be included in SPM
// ----------------------------------------------------------------------------------------------------------------------

void mm_info_initDefault(mm_info_t* mmInfo)
{
	mmInfo->isMatrix = -1;
	mmInfo->isSparse = -1;
	mmInfo->isCoordinate = -1;
	mmInfo->isDense = -1;
	mmInfo->isArray = -1;
	mmInfo->isComplex = -1;
	mmInfo->isReal = -1;
	mmInfo->isPattern = -1;
	mmInfo->isInteger = -1;
	mmInfo->isSymmetric = -1;
	mmInfo->isGeneral = -1;
	mmInfo->isSkew = -1;
	mmInfo->isHermitian = -1;

	mmInfo->nnz = 0;
	mmInfo->rowCount = 0;
	mmInfo->columnCount = 0;
}


void mm_info_init(mm_info_t* mmInfo, MM_typecode* typecode)
{
	mmInfo->isMatrix = mm_is_matrix(*typecode);
	mmInfo->isSparse = mm_is_sparse(*typecode);
	mmInfo->isCoordinate = mm_is_coordinate(*typecode);
	mmInfo->isDense = mm_is_dense(*typecode);
	mmInfo->isArray = mm_is_array(*typecode);
	mmInfo->isComplex = mm_is_complex(*typecode);
	mmInfo->isReal = mm_is_real(*typecode);
	mmInfo->isPattern = mm_is_pattern(*typecode);
	mmInfo->isInteger = mm_is_integer(*typecode);
	mmInfo->isSymmetric = mm_is_symmetric(*typecode);
	mmInfo->isGeneral = mm_is_general(*typecode);
	mmInfo->isSkew = mm_is_skew(*typecode);
	mmInfo->isHermitian = mm_is_hermitian(*typecode);
}

void mm_info_print(mm_info_t* mmInfo)
{
	printf("MM_INFO\n");
	printf("isMatrix: %d\n", mmInfo->isMatrix);
	printf("isSparse: %d\n", mmInfo->isSparse);
	printf("isCoordinate: %d\n", mmInfo->isCoordinate);
	printf("isDense: %d\n", mmInfo->isDense);
	printf("isArray: %d\n", mmInfo->isArray);
	printf("isComplex: %d\n", mmInfo->isComplex);
	printf("isReal: %d\n", mmInfo->isReal);
	printf("isPattern: %d\n", mmInfo->isPattern);
	printf("isInteger: %d\n", mmInfo->isInteger);
	printf("isSymmetric: %d\n", mmInfo->isSymmetric);
	printf("isGeneral: %d\n", mmInfo->isGeneral);
	printf("isSkew: %d\n", mmInfo->isSkew);
	printf("isHermitian: %d\n", mmInfo->isHermitian);
}

void mm_info_deleteNonPtr(mm_info_t* mmInfo)
{
}

void mm_info_delete(mm_info_t* mmInfo)
{
	if(mmInfo == NULL)
		return;

	mm_info_deleteNonPtr(mmInfo);
	free(mmInfo);
}

// ----------------------------------------------------------------------------------------------------------------------
