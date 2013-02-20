/*
 *  CLPKComplexOpsCM.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 1/13/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef ___CLPKComplexOpsCM_Header_
#define ___CLPKComplexOpsCM_Header_
 
#include "SharedFEData.h"
 
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <g2cEdited.h>

#include "SharedFEData.h"

#include <vecLib/vectorOps.h>
#include <vecLib/clapack.h>
#include <vecLib/vBLAS.h>

#define kHasAltiVecMask    ( 1 << gestaltPowerPCHasVectorInstructions )  // used in looking for a g4

#pragma mark ____SingleToComplex____

/*!
 convert a square complex matrix of size into a real matrix
 */
inline void realComplexMakeSingleVector(__CLPK_real* vecOut,
                                        __CLPK_complex* vecIn, 
                                        int size		)
{	int i;						
	for ( i=0; i<size; i++){		
		vecOut[i] = vecIn[i].r;
		//printf("wr[%d]: %lf\n", i, vecOut[i]);
	}
}


#pragma mark ___SingleOps___

/*!
 return the complex conjugate of a complex number
 */
inline __CLPK_complex	complexSingleGetConjugate(__CLPK_complex number) 
{
	__CLPK_complex result;
	result.r = number.r; 
	result.i = -number.i;
	return result;
}

/*!
 return the phase of a complex number
 */
inline double	complexSingleGetPhase(__CLPK_complex number) 
{
	return atan2(number.i, number.r);
}

/*!
 return the magnitude of a complex number
 */
inline double	complexSingleGetMagnitude(__CLPK_complex number)
{
	return sqrtf(pow(number.i,2) + pow(number.r,2));
}

/*!
 return the a complex number raised to a power
 */
inline __CLPK_complex complexSinglePow(__CLPK_complex number, float power)
{
#ifdef DEBUG_PRINT
	printf("number: %lf + %lfi power: %lf\n", number.r, number.i, power);
#endif

	double mag = complexSingleGetMagnitude(number); //printf("mag: %lf\n", mag);
	double phase = complexSingleGetPhase(number); //printf("phase: %lf\n", phase);
	
	__CLPK_complex result;
	result.r = pow(mag, power)*cos(phase*power);
	result.i = pow(mag, power)*sin(phase*power);
	return result;	
}

/*!
 add three square complex matricies and return the result in d
 */
inline void complexSingleAdd3Matricies(	__CLPK_complex* d, 
										__CLPK_complex* a, __CLPK_complex* b, __CLPK_complex* c, 
										int size)
{	int i,j;
	for ( i=0; i<size; i++){
		for ( j=0; j< size; j++){
						
			d[j*size+i].r = a[j*size+i].r + b[j*size+i].r + c[j*size+i].r;
			d[j*size+i].i = a[j*size+i].i + b[j*size+i].i + c[j*size+i].i;
			
#ifdef DEBUG_PRINT
	printf(" d[%ld][%ld]: %lf + %lfi", i,j, d[j*size+i].r,  d[j*size+i].i );
#endif

		}
	}

	return;
}

/*!
 multiply a complex matrix by a complex number and return it in another matrix
 */
inline void complexSingleMultMatrix(__CLPK_complex* matOut, __CLPK_complex scale, __CLPK_complex* matIn, int size)
{
	int i,j;
	for ( i=0; i<size; i++){
		for ( j=0; j<size; j++){
		
			matOut[j*size + i].r = matIn[j*size+i].r*scale.r - matIn[j*size+i].i*scale.i;
			matOut[j*size + i].i = matIn[j*size+i].r*scale.i + matIn[j*size+i].i*scale.i;
			
#ifdef DEBUG_PRINT
	printf("matIn[%ld][%ld]: %lf + %lfi scale: %lf + %lfi matOut[%ld][%ld]: %lf + %lfi", 
			i,j, matIn[i*size+j].r, matIn[i*size+j].i,
			scale.r, scale.i, i,j, 
			matOut[i*size+j].r, matOut[i*size+j].i);	
#endif
		}
	}
/*  
	//enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
	//enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
	
	void cblas_cgemm(CblasColMajor, CblasNoTrans,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc); 
	
*/	
}

/*!
 normalize the columns of a square complex matrix in place
 */
inline void complexSingleNormalizeMatrixColumns(__CLPK_complex* mat, int size)
{	int i,j;
	for ( i=0; i<size; i++){
		float norm = 0;
		for ( j=0; j<size; j++){
			norm += mat[j*size + i].r*mat[j*size + i].r + mat[j*size + i].i*mat[j*size + i].i;
		}
		norm = sqrtf(norm);
#ifdef DEBUG_PRINT
	printf("norm: %lf\n", norm);
#endif		
		for ( j=0; j<size; j++){
			mat[j*size + i].r /= norm;
			mat[j*size + i].i /= norm;
#ifdef DEBUG_PRINT
	printf("mat[%ld][%ld]: %lf + %lfi\n",j,i, mat[j*size + i].r, mat[j*size + i].i);
#endif

		}
	}
}

/*!
 out of place reorder the contents of a matrix from column major to row major format
 */
inline void complexSingleMakeRowMajor(__CLPK_complex* matOut, __CLPK_complex* matIn, int size)
{	
	int i,j;
	for ( i=0; i<size; i++){
		for ( j=0; j<size; j++){
			matOut[j*size + i] = matIn[i*size +j];		
		}		
	}
}

#pragma mark ___Eigen___
/*!
 perform an eigendecomposition of a complex matrix
 */
inline void complexSingleGeneralEigenDecompose(__CLPK_complex* wr, 
												__CLPK_complex* vr, __CLPK_complex*vl, 
												__CLPK_complex* mat, int size)
{	// in column mmajor
  integer matSize = (integer) size;

  /* Set up temporary variables. */ 
  __CLPK_complex	*work;                /* Working scratch area */
  __CLPK_real		*work2;
  __CLPK_integer	info;                 /* Return value of algorithm */
 
  /* Allocate space for LAPACK to do its computation. */
   __CLPK_integer lwork = 32*matSize;
  work = (__CLPK_complex*)malloc( lwork*sizeof(__CLPK_complex) );
  work2 = (__CLPK_real*)malloc( lwork*sizeof(__CLPK_real) );

#ifdef DEBUG_PRINT
	int row, col;
	printf("\nMatrix whose eigendecomposition is desired:\n");
	for (row = 0; row < matSize; row++)  {
	for (col = 0; col < matSize; col++)
	  printf(" mat[%ld][%ld]: %lf + %lfi\n", row, col, mat[col*matSize+row].r,  mat[col*matSize+row].i);
	printf("\n");
	}
#endif
 
#ifdef DEBUG_PRINT
	double time = currentTime();
#endif	
	// void cgeev(char jobvl, char jobvr, int n,  complex  *a,  int
	//          lda,  complex  *w,  complex *vl, int ldvl, complex
	//         *vr, int ldvr, int *info);
	//http://docs.sun.com/source/819-3691/cgeev.html

	cgeev_ (		&N, &V,  
					&matSize,  mat, 
					&matSize,  wr,  
					vl,  &matSize, vr, &matSize, 
					work, &lwork, work2, &info);

#ifdef DEBUG_PRINT
	time = currentTime() - time;
	printf("time to full decompose: %lf ms\n", time*1E3);


	 ///Check the output flag. 
	if (info==0)
		printf("\nComputation successful. Optimal lwork: %f\n",work[0]);
	else{
		printf("\nFailed to compute: %d\n", info);
		exit(1);
	}

	//Show the output. 
	int indx;
	printf("\nEigenvalues:\n");
	for ( indx = 0; indx < matSize; indx++)
		printf("  % 5.4e \n", wr[indx].r );

	printf("\nRight eigenvectors:\n");
	for ( col = 0; col < matSize; col++){
		printf( "\n  Eigenvector %d\n", col);
		for ( row = 0; row < matSize; row++)
			printf("    % 5.4e + % 5.4e\n", vr[col*matSize+row].r, vr[col*matSize+row].i );
	}
#endif
	// Free up the memory we allocated. 
	free(work), free(work2);

	/* Return success. */
}

/*!
 computes all eigenvalues and, optionally, eigenvectors of a
    complex Hermitian matrix mat
 */
inline void complexSingleHermetianEigenDecompose(	__CLPK_real* wr, 
													__CLPK_complex* vr, 
													__CLPK_complex* mat, int size)
{	// in column mmajor
  integer matSize = (integer) size;

  /* Set up temporary variables. */ 
  __CLPK_complex	*work;                /* Working scratch area */
  __CLPK_real		*work2;
  __CLPK_integer	info;                 /* Return value of algorithm */
 
  /* Allocate space for LAPACK to do its computation. */
   __CLPK_integer lwork = 32*matSize;
  work = (__CLPK_complex*)malloc( lwork*sizeof(__CLPK_complex) );
  work2 = (__CLPK_real*)malloc( lwork*sizeof(__CLPK_real) );

#ifdef DEBUG_PRINT
	int row, col;
	printf("\nMatrix whose eigendecomposition is desired:\n");
	for (row = 0; row < matSize; row++)  {
	for (col = 0; col < matSize; col++)
	  printf(" mat[%ld][%ld]: %lf + %lfi\n", row, col, mat[col*matSize+row].r,  mat[col*matSize+row].i);
	printf("\n");
	}
#endif

#ifdef DEBUG_PRINT 
	double time = currentTime();
#endif	

// Subroutine  int cheev_(char *jobz, char *uplo, __CLPK_integer *n, __CLPK_complex *a, 
//	__CLPK_integer *lda, __CLPK_real *w, __CLPK_complex *work, __CLPK_integer *lwork, __CLPK_real *rwork, 
//	__CLPK_integer *info);
//	http://docs.sun.com/source/819-3691/cheev.html
	
	cheev_ (		&V, &U, 
					&matSize,  mat, &matSize,
					wr, work, &lwork, work2, &info);

	 memcpy(vr, mat, matSize*matSize*sizeof(float));
					
#ifdef DEBUG_PRINT
	time = currentTime() - time;
	printf("time to full decompose: %lf ms\n", time*1E3);


	 ///Check the output flag. 
	if (info==0)
		printf("\nComputation successful. Optimal lwork: %f\n",work[0]);
	else{
		printf("\nFailed to compute: %d\n", info);
		exit(1);
	}

	//Show the output. 
	int indx;
	printf("\nEigenvalues:\n");
	for ( indx = 0; indx < matSize; indx++)
		printf("  % 5.4e \n", wr[indx] );

	printf("\nRight eigenvectors:\n");
	for ( col = 0; col < matSize; col++){
		printf( "\n  Eigenvector %d\n", col);
		for ( row = 0; row < matSize; row++)
			printf("    % 5.4e + % 5.4e\n", vr[col*matSize+row].r, vr[col*matSize+row].i );
	}
#endif
	// Free up the memory we allocated. 
	free(work), free(work2);

	/* Return success. */
}



#pragma mark ____Double____

inline __CLPK_doublecomplex	complexDoubleGetConjugate(__CLPK_doublecomplex number) 
{
	__CLPK_doublecomplex result;
	result.r = number.r; 
	result.i = -number.i;
	return result;
}

inline double	complexDoubleGetPhase(__CLPK_doublecomplex number) 
{
	return atan2(number.i, number.r);
}

inline double	complexDoubleGetMagnitude(__CLPK_doublecomplex number) 
{
	return sqrtf(pow(number.i,2) + pow(number.r,2));
}

inline __CLPK_doublecomplex complexDoublePow(__CLPK_doublecomplex number, float power)
{
	double phase = complexDoubleGetPhase(number);
	double mag = complexDoubleGetMagnitude(number);
	
	__CLPK_doublecomplex result;
	result.r = pow(mag, power)*cos(phase*power);
	result.i = pow(mag, power)*sin(phase*power);
	return result;	
}

inline void complexDoubleAdd3Matricies(	__CLPK_doublecomplex* d, 
										__CLPK_doublecomplex* a, __CLPK_doublecomplex* b, __CLPK_doublecomplex* c, 
										int size)
{	int i,j;
	for ( i=0; i<size; i++){
		for ( j=0; j< size; j++){
						
			d[j*size+i].r = a[j*size+i].r + b[j*size+i].r + c[j*size+i].r;
			d[j*size+i].i = a[j*size+i].i + b[j*size+i].i + c[j*size+i].i;
			
#ifdef DEBUG_PRINT
	if (i==0 || i== 65){
			printf(" [%ld][%ld] - a: %lf + %lfi b %lf + %lfi c: %lf + %lfi d: %lf + %lfi\n", i,j, 
			a[j*size+i].r,  a[j*size+i].i, b[j*size+i].r,  b[j*size+i].i, c[j*size+i].r,  c[j*size+i].i, d[j*size+i].r,  d[j*size+i].i );
	}
#endif

		}
	}

	return;
}

inline void complexDoubleMultMatrix(__CLPK_doublecomplex* matOut, __CLPK_doublecomplex scale, __CLPK_doublecomplex* matIn, int size)
{	int i,j;
	for ( i=0; i<size; i++){
		for ( j=0; j<size; j++){
		
			matOut[j*size + i].r = matIn[j*size+i].r*scale.r - matIn[j*size+i].i*scale.i;
			matOut[j*size + i].i = matIn[j*size+i].r*scale.i + matIn[j*size+i].i*scale.i;
			
#ifdef DEBUG_PRINT
			printf("matIn[%ld][%ld]: %lf + %lfi scale: %lf + %lfi matOut[%ld][%ld]: %lf + %lfi", 
			i,j, matIn[i*size+j].r, matIn[i*size+j].i,
			scale.r, scale.i, i,j, 
			matOut[i*size+j].r, matOut[i*size+j].i);			
#endif
		}
	}
}

inline void complexDoubleNormalizeMatrixColumns(__CLPK_doublecomplex* mat, int size)
{	int i,j;
	for ( i=0; i<size; i++){
		float norm = 0;
		for ( j=0; j<size; j++){
			norm += mat[j*size + i].r*mat[j*size + i].r + mat[j*size + i].i*mat[j*size + i].i;
		}
		norm = sqrtf(norm);
#ifdef DEBUG_PRINT
		printf("norm: %lf\n", norm);
#endif		
		for ( j=0; j<size; j++){
			mat[j*size + i].r /= norm;
			mat[j*size + i].i /= norm;
			
#ifdef DEBUG_PRINT
	printf("mat[%ld][%ld]: %lf + %lfi\n",j,i, mat[j*size + i].r, mat[j*size + i].i);
#endif

		}
	}
}

#pragma mark ___Eigen___

inline void complexDoubleGeneralEigenDecompose(	__CLPK_doublecomplex* wr, __CLPK_doublecomplex* vr, __CLPK_doublecomplex*vl, 
											__CLPK_doublecomplex* mat, 
											int size)
{
  integer matSize = (integer) size;

  /* Set up temporary variables. */ 
  __CLPK_doublecomplex *work;                   /* Working scratch area */
  __CLPK_doublereal *work2;
  __CLPK_integer info;                 /* Return value of algorithm */
 
  /* Allocate space for LAPACK to do its computation. */
   __CLPK_integer lwork = 32*matSize;
  work = (__CLPK_doublecomplex*)malloc( lwork*sizeof(__CLPK_doublecomplex) );
  work2 = (__CLPK_doublereal*)malloc( lwork*sizeof(__CLPK_doublereal) );
	
#ifdef DEBUG_PRINT
	int row, col;
	printf("\nMatrix whose eigendecomposition is desired:\n");
	for ( row = 0; row < matSize; row++)  {
	for ( col = 0; col < matSize; col++)
	  printf("mat[%ld][%ld]:  %lf + %lf\n", row, col, mat[col*matSize+row].r,  mat[col*matSize+row].i);
	printf("\n");
	}
#endif

	zgeev_ (	&N, &V,  
				&matSize,  mat, 
				&matSize,  wr,  
				vl,  &matSize, vr, &matSize, 
				work, &lwork, work2, &info);

	
#ifdef DEBUG_PRINT
	// Check the output flag. 
	if (info==0)
		printf("\nComputation successful. Optimal lwork: %f\n",work[0]);
	else{
		printf("\nFailed to compute: %d\n", info);
		exit(1);
	}
	
	// Show the output. 
	printf("\nEigenvalues:\n");
	int indx;
	for ( indx = 0; indx < matSize; indx++)
		printf("  % 5.4e \n", wr[indx] );
	
	printf("\nRight eigenvectors:\n");
	for ( col = 0; col < matSize; col++){
		printf( "\n  Eigenvector %d\n", col);
		for ( row = 0; row < matSize; row++)
			printf("    % 5.4e + % 5.4e\n", vr[col*matSize+row].r, vr[col*matSize+row].i );
	}
#endif
	// Free up the memory we allocated. 
	free(work), free(work2);
}


#endif