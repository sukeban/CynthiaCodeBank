/*
 *  CLPKRealMatrixOpsCM.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 1/11/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef __CLPKRealMatrixOpsCM_Header_
#define __CLPKRealMatrixOpsCM_Header_

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

//#include <g2c.h>
#include "g2cEdited.h"

#include "SharedFEData.h"

#include <vecLib/vecLibTypes.h>
#include <vecLib/clapack.h>
#include <vecLib/vDSP.h>
#include <vecLib/vectorOps.h>
#include <vecLib/vBLAS.h>

#define kHasAltiVecMask    ( 1 << gestaltPowerPCHasVectorInstructions )  // used in looking for a g4


#pragma mark -
#pragma mark change

/*!
 convert a rsquare eal matrix of size into a complex matrix
 */
void realSingleMakeSingleComplex(__CLPK_complex*   matOut,
                                  __CLPK_real*      matIn,
                                  int				size	)
{	int i,j;						
	for ( i=0; i<size; i++){
		for ( j=0; j<size; j++){
			matOut[j*size + i].r = matIn[j*size+i];
		}
	}
}

/*!
 convert a real matrix to a complex matrix times a scale
*/
void realSingleToComplexSingleMultMatrix(__CLPK_complex*   matOut, 
                                          __CLPK_complex    scale, 
                                          __CLPK_real*      matIn, 
                                          int               size	)				// used for conveting the real stiffness matrix into complex fourier representation
{	
	int i,j;
	for ( i=0; i<size; i++){
		for ( j=0; j<size; j++){
		
			matOut[j*size + i].r = matIn[j*size+i]*scale.r;
			matOut[j*size + i].i = matIn[j*size+i]*scale.i;
			
#ifdef DEBUG_PRINT
	printf("matIn[%ld][%ld]: %lf + %lfi scale: %lf + %lfi matOut[%ld][%ld]: %lf + %lfi\n", 
			i,j, matIn[j*size+i], 0.0,
			scale.r, scale.i, i,j, 
			matOut[j*size+i].r, matOut[j*size+i].i);	
#endif		
		}
	}
}

#pragma mark -

#pragma mark vector
/*!
 scale a single precision vector
 */
void realSingleVectorScale(float* B, float* A, float beta, int length)
{
	float b = 1.0/beta;
	vDSP_vsdiv (A, 1, &b, B, 1, length);
} 

/*!
 scale adouble precision vector
 */
void realDoubleVectorScale(double* B, double* A, double beta, int length)
{
	double b = 1.0/beta;
	vDSP_vsdivD (A, 1, &b, B, 1, length);
} 

/*!
 add a single precision vector
 */
void realSingleVectorAdd(float* C, float* A, float* B, int length)
{
	vDSP_vadd (A, 1, B, 1, C, 1, length);
}

/*!
 add a double precision vector
 */
void realDoubleVectorAdd(double* C, double* A, double* B, int length)
{
	vDSP_vaddD (A, 1, B, 1, C, 1, length);
}

/*!
 subtract a single precision vector
 */
void realSingleVectorSubtract(float* C, float* A, float* B, int length)
{
	vDSP_vsub (B, 1, A, 1, C, 1, length);
}

/*!
 subtract a double precision vector
 */
void realDoubleVectorSubtract(double* C, double* A, double* B, int length)
{
	vDSP_vsubD (B, 1, A, 1, C, 1, length);
}

/*!
 solve A*x = B for single precision
 */
void realSingleSolve(float* kcopy, float* fcopy, int cdofs)
{
	integer w = 1;
	integer ipiv[cdofs];
	integer info;
		
	integer n = (integer) cdofs;
		
	//  star timer
//	double time = currentTime();
	
	//http://www.netlib.org/clapack/single/sgesv.c
	sgesv_(&n, &w, kcopy, &n, ipiv, fcopy, &n, &info);
	
	// end timer
//	time = currentTime() - time;
//	printf("time to solve %d equations =  %lf ms\n", cdofs, time*1E3);
}

/*!
 solve A*x = B for double precision
 */
void realDoubleSolve(double* kcopy, double* fcopy, int cdofs)
{
	integer w = 1;
	integer ipiv[cdofs];
	integer info;
		
	integer n = (integer) cdofs;
		
//	double time = currentTime();
	
	//http://www.netlib.org/clapack/double/dgesv.c
	dgesv_(&n, &w, kcopy, &n, ipiv, fcopy, &n, &info);
	
	// end timer
//	time = currentTime() - time;
//	printf("time to solve %d equations =  %lf ms\n", cdofs, time*1E3);
}

#pragma mark -
#pragma mark matrix
/*!
 normalize the rows of a matrix
 */
void realSingleNormalizeMatrixRows(__CLPK_real* mat, int size	)
{	
	int i,j;
	for ( i=0; i<size; i++){
		__CLPK_real norm = 0;
		for ( j=0; j<size; j++){
			norm += pow(mat[i*size + j],2);
		}
		norm = sqrtf(norm);
		
#ifdef DEBUG_PRINT
		printf("norm: %lf\n", norm);
#endif		
		for ( j=0; j<size; j++){
			mat[i*size + j] /= norm;
			
#ifdef DEBUG_PRINT
			printf("mat[%ld][%ld]: %lf\n",j,i, mat[i*size + j]);
#endif
		}
	}
}

/*!
 normalize the columns of a matrix
 */
void realSingleNormalizeMatrixColumns(	__CLPK_real* mat, int size	)
{	
	int i,j;
	for ( i=0; i<size; i++){
		__CLPK_real norm = 0;
		for ( j=0; j<size; j++){
			norm += pow(mat[j*size + i],2);
		}
		norm = sqrtf(norm);
		
#ifdef DEBUG_PRINT
		printf("norm: %lf\n", norm);
#endif		
		for ( j=0; j<size; j++){
			mat[j*size + i] /= norm;
			
#ifdef DEBUG_PRINT
			printf("mat[%ld][%ld]: %lf\n",j,i, mat[j*size + i]);
#endif
		}
	}
}

/*!
 multiply matrix a an b and put the result in c in single precision
 */
void realSingleMatrixVectorMult(__CLPK_real* c, __CLPK_real* a, __CLPK_real* b, int size)
{
	for (int i=0; i< size; i++){
		c[i] = 0.0;
		for (int j=0; j< size; j++){
			c[i] += a[i + j*size]*b[j];	
		}
	}
}

/*!
 multiply matrix a an b and put the result in c in double precision
 */
void realDoubleMatrixVectorMult(__CLPK_doublereal* c, __CLPK_doublereal* a, __CLPK_doublereal* b, int size)
{
	for (int i=0; i< size; i++){
		c[i] = 0.0;
		for (int j=0; j< size; j++){
 			c[i] += a[i + j*size]*b[j];	
		}
	}
}

/*!
 scale matrix a by scalar b and put the result in c in single precision
 */
void realSingleMatrixScale(__CLPK_real* c, __CLPK_real* a, __CLPK_real beta, int size)
{
//	printf("\nm_scale:\n");

	for (int i=0; i< size; i++){
		for (int j=0; j<size; j++){
			c[i+j*size] = beta*a[i+j*size];	
			//printf("%lf\t", c[i+j*size]);
		}
	//	printf("\n");
	}	
}

/*!
 scale matrix a by scalar b and put the result in c in double precision
 */
void realDoubleMatrixScale(__CLPK_doublereal* c, __CLPK_doublereal* a, __CLPK_doublereal beta, int size)
{
//	printf("\nm_scale:\n");

	for (int i=0; i< size; i++){
		for (int j=0; j<size; j++){
			c[i+j*size] = beta*a[i+j*size];	
			//printf("%lf\t", c[i+j*size]);
		}
	//	printf("\n");
	}	
}

/*!
 add square matrix a to square matrix b and put the result in squre matrix c in single precision
 */
void realSingleAdd2Matricies(float* c, float* a, float* b, int size) // square matrix
{	

//	const vFloat* aa = (const vFloat*) a;
//	const vFloat* bb = (const vFloat*) b;
//	vFloat* cc = (vFloat*) c;
//	vSgeadd(size, size, aa,'n', bb,'n', cc);   //non pow4 execute does not

	//printf("\nadded \n");

	int i,j;
	for ( i=0; i<size; i++){
		for ( j=0; j< size; j++){
			c[j*size+i] = a[j*size+i]+b[j*size+i]; 
			//printf(" %lf\t", c[j*size+i]);
		}
		//printf("\n");
	}
}

/*!
 add square matrix a to square matrix b and put the result in squre matrix c in double precision
 */
void realDoubleAdd2Matricies(__CLPK_doublereal* c, __CLPK_doublereal* a, __CLPK_doublereal* b, int size) // square matrix
{	

//	const vFloat* aa = (const vFloat*) a;
//	const vFloat* bb = (const vFloat*) b;
//	vFloat* cc = (vFloat*) c;
//	vSgeadd(size, size, aa,'n', bb,'n', cc);   //non pow4 execute does not

	//printf("\nadded \n");

	int i,j;
	for ( i=0; i<size; i++){
		for ( j=0; j< size; j++){
			c[j*size+i] = a[j*size+i]+b[j*size+i]; 
			//printf(" %lf\t", c[j*size+i]);
		}
		//printf("\n");
	}

}

/*!
 add three square matricies in single precision, put the result in matrix d
 */
void realSingleAdd3Matricies(float* d, float* a, float* b, float*c, int size	) // square matrix
{	
	float* temp = (float*) calloc(size*size, sizeof(float));
	realSingleAdd2Matricies(temp,a,b,size);   
	realSingleAdd2Matricies(d,temp,c,size);
/*	int i,j;
	for ( i=0; i<size; i++)
		for ( j=0; j< size; j++)
			d[j*size+i] = a[j*size+i]+b[j*size+i]+c[j*size+i];
*/			
	free(temp);
}

/*!
 add three square matricies in double precision, put the result in matrix d
 */
void realDoubleAdd3Matricies(double* d, double* a, double* b, double*c, int size	) // square matrix
{	
	double* temp = (double*) calloc(size*size, sizeof(double));
	realDoubleAdd2Matricies(temp,a,b,size);   
	realDoubleAdd2Matricies(d,temp,c,size);
/*	int i,j;
	for ( i=0; i<size; i++)
		for ( j=0; j< size; j++)
			d[j*size+i] = a[j*size+i]+b[j*size+i]+c[j*size+i];
*/			
	free(temp);
}

/*!
 subtract square matrix a from square matrix b and put the result in squre matrix c in single precision
 */
void realSingleSubtract2Matricies(float* c, float* a, float* b, int size) // square matrix
{	
//	const vFloat* aa = (const vFloat*) a;
//	const vFloat* bb = (const vFloat*) b;
//	vFloat* cc = (vFloat*) c;
	
//	vSgesub(size,size, aa,'n', bb,'n', cc); // non pow4 does not execute
  
	int i,j;
	for ( i=0; i<size; i++){
		for ( j=0; j< size; j++){
			c[j*size+i] = a[j*size+i]-b[j*size+i]; 
			//printf("c[%d][%d]: %lf\n", i,j, c[j*size+i]);
		}
	}
}

/*!
 orthonormalize matrix m of size n x m stored in array sizeM in single precision
 */
void realSingleOrthonormalizeMatrix(__CLPK_real* m, __CLPK_integer sizeM[])
{
	__CLPK_integer info;
	__CLPK_integer lwork = 32*sizeM[1];
	__CLPK_real *tau = (__CLPK_real*) calloc(sizeM[1], sizeof(__CLPK_real));
	__CLPK_real* work = (__CLPK_real*) calloc(lwork, sizeof(__CLPK_real));

	//http://www.csee.umbc.edu/help/fortran/LAPACK/SRC/sgeqrf.f
	int ret = sgeqrf_(
						&sizeM[0], &sizeM[1], 
						m, 
						&sizeM[0], 
						tau, 
						work, 
						&lwork, &info);
	
	
	//http://www.csee.umbc.edu/help/fortran/LAPACK/SRC/sorg2r.f
	ret = sorg2r_(	
					&sizeM[0], &sizeM[1], 
					&sizeM[1], 
					m,
					&sizeM[0], 
					tau,
					work, 
					&info);
						 
	free(tau);
	free(work);
}

/*!
 orthonormalize matrix m of size n x m stored in array sizeM in double precision
 */
void realDoubleOrthonormalizeMatrix(__CLPK_doublereal* m, __CLPK_integer sizeM[])
{
	__CLPK_integer info;
	__CLPK_integer lwork = 32*sizeM[1];
	__CLPK_doublereal *tau = (__CLPK_doublereal*) calloc(sizeM[1], sizeof(__CLPK_doublereal));
	__CLPK_doublereal* work = (__CLPK_doublereal*) calloc(lwork, sizeof(__CLPK_doublereal));
	
		
	//http://www.csee.umbc.edu/help/fortran/LAPACK/SRC/dgeqrf.f
	int ret = dgeqrf_(
						&sizeM[0], &sizeM[1], 
						m, 
						&sizeM[0], 
						tau, 
						work, 
						&lwork, &info);
	
	
	//http://www.csee.umbc.edu/help/fortran/LAPACK/SRC/dorg2r.f
	//http://www.netlib.org/lapack/explore-html/dorg2r.f.html
	ret = dorg2r_(	
					&sizeM[0], &sizeM[1], 
					&sizeM[1], 
					m,
					&sizeM[0], 
					tau,
					work, 
					&info);
						 
	free(tau);
	free(work);
}

/*!
 transpose matrix a of size n x m stored in array size and return it in matrix b, in single precision
 */
void realSingleMatrixOuOfPlaceTranspose(float* a, int size[], float* b)
{
	//vSgetmi inplace
	//vSgetmo outof place

	for (int i=0; i<size[0]; i++){
		for (int j=0; j<size[1]; j++){
			b[j+i*size[1]] = a[i+j*size[0]];
			//printf("a^T[%d][%d]: %lf\n", j, i, b[j+i*size[1]]);
		}
	}

	//const vFloat* aa = (const vFloat*) a;
	//vFloat* bb = (vFloat*) b;
	
	//vSgetmo( size[0], size[1], aa, bb); // if not a multiple of 4 it doesn't execute the scalar code
}

/*!
 transpose matrix a of size n x m stored in array size and return it in matrix b, in double precision
 */
void realDoubleMatrixOutOfPlaceTranspose(double* a, int size[], double* b)
{
	//vSgetmi inplace
	//vSgetmo outof place

	for (int i=0; i<size[0]; i++){
		for (int j=0; j<size[1]; j++){
			b[j+i*size[1]] = a[i+j*size[0]];
#ifdef DEBUG_PRINT
			printf("a^T[%d][%d]: %lf\n", j, i, b[j+i*size[1]]);
#endif
		}
	}

	//const vFloat* aa = (const vFloat*) a;
	//vFloat* bb = (vFloat*) b;
	
	//vSgetmo( size[0], size[1], aa, bb); // if not a multiple of 4 it doesn't exect the scalar code
}

/*!
 compute the  Cholesky decomposition of matrix a of size n x m stored in array size and return it in place, in single precision
 */
void realSingleMatrixCholesky(float* a, int size)
{
	//http://www.netlib.org/clapack/single/spotrf.c
	//int spotrf_(char *uplo, integer *n, real *a, integer *lda, integer *info)
	integer info;
	integer s = size;
	spotrf_(&L, &s, a, &s,&info);
}

/*!
 compute the  inverse of matrix a of size n x m stored in array size and return it in place, in single precision
 */
void realSingleMatrixInverse(float* a, int size)
{
	integer info;
	integer s = size;
	integer ipiv[size];
	float work[size];
	//http://www.netlib.org/clapack/single/sgetrf.c
	sgetrf_(&s,&s,a,&s,ipiv,&info);     
	//http://www.netlib.org/clapack/single/sgetri.c
	sgetri_(&s,a,&s,ipiv,work,&s,&info);	
	
	//printf("\ninv_a:\n");
	//for (int i=0; i< size; i++){
	//	for (int j=0; j<size; j++){
	//		printf(" %lf\t", a[i+j*size]);
	//	}
	//	printf("\n");
	//}
}

/*!
 compute the  inverse of matrix a of size n x m stored in array size and return it in place, in double precision
 */
void realDoubleMatrixInverse(double* a, int size)
{
	integer info;
	integer s = size;
	integer ipiv[size];
	double work[size];
	//http://www.netlib.org/clapack/single/sgetrf.c
	dgetrf_(&s,&s,a,&s,ipiv,&info);     
	//http://www.netlib.org/clapack/single/sgetri.c
	dgetri_(&s,a,&s,ipiv,work,&s,&info);	
	
	//printf("\ninv_a:\n");
	//for (int i=0; i< size; i++){
	//	for (int j=0; j<size; j++){
	//		printf(" %lf\t", a[i+j*size]);
	//	}
	//	printf("\n");
	//}
}



#pragma mark -
#pragma mark matrix eigen

/*!
 compute the eigen decomposition for a non-symmatric matrix in single precision
 wr and wi are the eigenvalues, vr and vl are the eigenvectors for the square matrix mat of size sizek
 */
void realSingleNonSymmetricEigenDecomp(float* wr, float* wi, float* vr, float*vl, float* mat, int sizeK)
{
	integer matSize = (integer) sizeK;

	// Set up temporary variables. 
	integer lwork;                // Size of working scratch area 
	float *work;                  // Working scratch area 
	integer info;                 // Return value of algorithm 

	// Allocate space for LAPACK to do its computation. 
	lwork = 32*matSize;
	work = (float*)malloc( lwork*sizeof(float) );
	if ( (work==NULL) )	{
		char *memString = "\nBuy more memory!\n";
		printf(memString);
		exit(1);
	}
	
#ifdef DEBUG_PRINT
	printf("\nMatrix whose eigendecomposition is desired:\n");
	int row, col;
	for (row = 0; row < matSize; row++)  {
		for (col = 0; col < matSize; col++)
			printf(" %lf", mat[col*matSize+row] );
		printf("\n");
	}
	printf("\nComputing eigenvectors and eigenvalues.\n");
#endif

	/* 	
	 LAPACK libraries require that all arguments be passed as pointers. 
	 sgeev - compute for an N-by-N real  nonsymmetric  matrix  A,
     the  eigenvalues  and,  optionally,  the  left  and/or right
     eigenvectors
	 http://docs.sun.com/source/819-3691/sgeev.html
	*/ 
	 sgeev_(&N,&V,&matSize,mat,&matSize,wr,wi,vl,&matSize,vr,&matSize,work,&lwork,&info);

	int indx;
	for ( indx = 0; indx < matSize; indx++)
		if (wr[indx] < 0.0) wr[indx] = 0.0;

	//  but now you have to sort this
	//float* wr_sorted = (float*) calloc(matSize, sizeof(float));
	//reclaimRealEigenValuesRealEigenVectors((__CLPK_real*) wr_sorted, (__CLPK_real*) wr, (__CLPK_real*) vl, matSize, 1);
	//memcpy(wr, wr_sorted, matSize*sizeof(float));
	//free(wr_sorted);
	
#ifdef DEBUG_PRINT	
	// Check the output flag. 
	if (info==0)
		printf("\nComputation successful. Optimal lwork: %f\n",work[0]);
	else{
		printf("\nFailed to compute: %d\n", info);
		exit(1);
	}
	
	printf("\nEigenvalues:\n");
	for ( indx = 0; indx < matSize; indx++)
		printf("  % 5.4e+% 5.4ei \n", wr[indx], wi[indx] );
	
//#ifdef DEBUG_PRINT	// Show the output. 	
	printf("\nRight eigenvectors:\n");
	
	for ( col = 0; col < matSize; col++){
		printf( "\n  Eigenvector %d\n", col+1 );
		for ( row = 0; row < matSize; row++)
			printf("    % 5.4e\n", vr[col*matSize+row] );
	}
#endif

	//Free up the memory we allocated. 
	free(work);

	// Return success. 
}

/*!
 compute the eigen decomposition for a generalized symmetric-definite eigenproblem in single precision
 wr are eigenvalues, on entry K is the stiffness matrix and on exit K contains the eigenvectors, the size of the matrix is K
 on entry M is the mass matrix, on exit, it contains the cholesky decomposition of the mass matrix
 */
void realSingleSymmetricGeneralEigenDecomp(float* wr, float* K, float* M, int sizeK)
{
	integer matSize = (integer) sizeK;

	// Set up temporary variables. 
	integer lwork;                // Size of working scratch area 
	float *work;                  // Working scratch area 
	integer info;                 // Return value of algorithm 

	// Allocate space for LAPACK to do its computation. 
	lwork = 32*matSize;
	work = (float*)malloc( lwork*sizeof(float) );
	if ( (work==NULL) )	{
		char *memString = "\nBuy more memory!\n";
		printf(memString);
		exit(1);
	}
	
#ifdef DEBUG_PRINT
	int row, col;
	printf("\nMatrix whose eigendecomposition is desired:\n");

	for (row = 0; row < matSize; row++)  {
		for (col = 0; col < matSize; col++)
		printf(" %lf\n", K[col*matSize+row] );
		printf("\n");
	}
	printf("\nComputing eigenvectors and eigenvalues.\n");
#endif

	/* 	
	LAPACK libraries require that all arguments be passed as pointers. 
	ssygv - compute all the  eigenvalues,  and  optionally,  the
	eigenvectors of a real generalized symmetric-definite eigen-
	problem, of the form A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or
	B*A*x=(lambda)*x
	http://docs.sun.com/source/819-3691/ssygv.html
	*/ 
	__CLPK_integer type  = 1;
	ssygv_(&type, &V, &U, &matSize, K, &matSize, M, &matSize, wr, work, &lwork, &info);	 
	
	
#ifdef DEBUG_PRINT	
	printf("n: %d, info: %d\n", matSize,info);
	// Check the output flag. 
	if (info==0)
		printf("\nComputation successful. Optimal lwork: %f\n",work[0]);
	else{
		printf("\nFailed to compute: %d\n", info);
		exit(1);
	}
	
	printf("\nEigenvalues:\n");
	for (int  indx = 0; indx < matSize; indx++)
		printf("  % 5.4e \n", wr[indx] );
	
//#ifdef DEBUG_PRINT	// Show the output. 	
	printf("\nRight eigenvectors:\n");
	for ( col = 0; col < matSize; col++){
		printf( "\n  Eigenvector %d\n", col+1 );
		for ( row = 0; row < matSize; row++)
			printf("    % 5.4e\n", K[col*matSize+row] );
	}
#endif

	for (int indx = 0; indx < matSize; indx++)
		if (wr[indx] < 0.0) wr[indx] = 0.0;

	//Free up the memory we allocated. 
	free(work);
}

/*!
 compute the eigen decomposition for a symmetric matrix in single precision
 wr are eigenvalues, on entry lkl is the stiffness matrix and on exit lkl contains the eigenvectors, the size of the matrix is sizeK
 */
void realSingleSymmetricEigenDecomp(float* wr, float* lkl, int sizeK)
{
	integer matSize = (integer) sizeK;

	// Set up temporary variables. 
	integer lwork;                // Size of working scratch area 
	float *work;                  // Working scratch area 
	integer info;                 // Return value of algorithm 

	// Allocate space for LAPACK to do its computation. 
	lwork = 32*matSize;
	work = (float*)malloc( lwork*sizeof(float) );
	if ( (work==NULL) )	{
		char *memString = "\nBuy more memory!\n";
		printf(memString);
		exit(1);
	}
	
	
#ifdef DEBUG_PRINT
	int row, col;
	printf("\nMatrix whose eigendecomposition is desired:\n");

	for (row = 0; row < matSize; row++)  {
		for (col = 0; col < matSize; col++)
		printf(" %lf\n", lkl[col*matSize+row] );
		printf("\n");
	}
	printf("\nComputing eigenvectors and eigenvalues.\n");
#endif


	//if (num == sizeK){
		/* 	
		 LAPACK libraries require that all arguments be passed as pointers. 
		 sgeev - compute for an N-by-N real  nonsymmetric  matrix  A,
		 the  eigenvalues  and,  optionally,  the  left  and/or right
		 eigenvectors
		//	http://docs.sun.com/source/819-3691/ssyev.html
		*/ 
		 ssyev_(&V,&U,&matSize, lkl, &matSize, wr, work,&lwork,&info);	 
	//}
	//else {

		/*  sgeevx - compute for an N-by-N real nonsymmetric  matrix  A,
		the  eigenvalues  and,  optionally,  the  left  and/or right
		eigenvectors
		http://docs.sun.com/source/819-3691/sgeevx.html
		*/

	//}

	// come back sorted
	
#ifdef DEBUG_PRINT	
	// Check the output flag. 
	if (info==0)
		printf("\nComputation successful. Optimal lwork: %f\n",work[0]);
	else{
		printf("\nFailed to compute: %d\n", info);
		exit(1);
	}
	
	
	printf("\nEigenvalues:\n");
	for (int indx = 0; indx < matSize; indx++)
		printf("  % 5.4e \n", wr[indx] );
	
//#ifdef DEBUG_PRINT	// Show the output. 	
	printf("\nRight eigenvectors:\n");
	
	for ( col = 0; col < matSize; col++){
		printf( "\n  Eigenvector %d\n", col+1 );
		for ( row = 0; row < matSize; row++)
			printf("    % 5.4e\n", lkl[col*matSize+row] );
	}
#endif

	for (int indx = 0; indx < matSize; indx++)
		if (wr[indx] < 0.0) wr[indx] = 0.0;

	//Free up the memory we allocated. 
	free(work);

}

/*!
 compute the eigen decomposition for a non-symmatric matrix in double precision
 wr and wi are the eigenvalues, vr and vl are the eigenvectors for the square matrix mat of size sizek
 */
void realDoubleNonSymmetricEigenDecomp(double* wr, double* wi, double* vr, double*vl, double* mat, int sizeK)
{
	integer matSize = (integer) sizeK;

	// Set up temporary variables. 
	integer lwork;                // Size of working scratch area 
	double *work;                  // Working scratch area 
	integer info;                 // Return value of algorithm 

	// Allocate space for LAPACK to do its computation. 
	lwork = 32*matSize;
	work = (double*)malloc( lwork*sizeof(double) );
	if ( (work==NULL) )	{
		char *memString = "\nBuy more memory!\n";
		printf(memString);
		exit(1);
	}
	
#ifdef DEBUG_PRINT
	printf("\nMatrix whose eigendecomposition is desired:\n");
	int row, col;
	for (row = 0; row < matSize; row++)  {
		for (col = 0; col < matSize; col++)
		printf(" %lf", mat[col*matSize+row] );
		printf("\n");
	}
	printf("\nComputing eigenvectors and eigenvalues.\n");
#endif


	/* 	
	 LAPACK libraries require that all arguments be passed as pointers. 
	 sgeev - compute for an N-by-N real  nonsymmetric  matrix  A,
     the  eigenvalues  and,  optionally,  the  left  and/or right
     eigenvectors
	 http://docs.sun.com/source/819-3691/dgeev.html
	*/ 
	 dgeev_(&N,&V,&matSize,mat,&matSize,wr,wi,vl,&matSize,vr,&matSize,work,&lwork,&info);

	int indx;
	for ( indx = 0; indx < matSize; indx++)
		if (wr[indx] < 0.0) wr[indx] = 0.0;

	//  but now you have to sort this
	//double* wr_sorted = (double*) calloc(matSize, sizeof(double));
	//reclaimRealEigenValuesRealEigenVectors((__CLPK_real*) wr_sorted, (__CLPK_real*) wr, (__CLPK_real*) vl, matSize, 1);
	//memcpy(wr, wr_sorted, matSize*sizeof(double));
	//free(wr_sorted);
	
#ifdef DEBUG_PRINT	
	// Check the output flag. 
	if (info==0)
		printf("\nComputation successful. Optimal lwork: %f\n",work[0]);
	else{
		printf("\nFailed to compute: %d\n", info);
		exit(1);
	}
	
	
	printf("\nEigenvalues:\n");
	for ( indx = 0; indx < matSize; indx++)
		printf("  % 5.4e+% 5.4ei \n", wr[indx], wi[indx] );
	
//#ifdef DEBUG_PRINT	// Show the output. 	
	printf("\nRight eigenvectors:\n");
	
	for ( col = 0; col < matSize; col++){
		printf( "\n  Eigenvector %d\n", col+1 );
		for ( row = 0; row < matSize; row++)
			printf("    % 5.4e\n", vr[col*matSize+row] );
	}
#endif

	//Free up the memory we allocated. 
	free(work);

	// Return success. 
}


/*!
 compute the eigen decomposition for a generalized symmetric-definite eigenproblem in double precision
 wr are eigenvalues, on entry K is the stiffness matrix and on exit K contains the eigenvectors, the size of the matrix is K
 on entry M is the mass matrix, on exit, it contains the cholesky decomposition of the mass matrix
 */
void realDoubleSymmetricGeneralEigenDecomp(double* wr, double* K, double* M, int sizeK)
{
	integer matSize = (integer) sizeK;

	// Set up temporary variables. 
	integer lwork;
	double *work;

	integer info;					// Return value of algorithm 

	// Allocate space for LAPACK to do its computation. 
	lwork =  32*matSize;
	work = (double*)malloc( lwork*sizeof(double) );
	if ( (work==NULL) )	{
		char *memString = "\nBuy more memory!\n";
		printf(memString);
		exit(1);
	}

	
#ifdef DEBUG_PRINT	
	int row, col;
	printf("K=[\n");
	for (row = 0; row < matSize; row++)  {
		for (col = 0; col < matSize; col++){
			printf(" %10.10e\t", K[col*matSize+row] );
		}
		printf(";\n");
	}
	printf("];\n");
	

	printf("M=[\n");
	for (row = 0; row < matSize; row++)  {
		for (col = 0; col < matSize; col++){
			printf(" %10.10e\t", M[col*matSize+row] );
		}
		printf(";\n");
	}
	printf("];\n");
#endif

	/* 	
	LAPACK libraries require that all arguments be passed as pointers. 
	ssygv - compute all the  eigenvalues,  and  optionally,  the
	eigenvectors of a real generalized symmetric-definite eigen-
	problem, of the form A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or
	B*A*x=(lambda)*x
	http://docs.sun.com/source/819-3691/dsygv.html
	*/ 
	__CLPK_integer type  = 1;
	dsygv_(&type, &V, &U, &matSize, K, &matSize, M, &matSize, wr, work, &lwork, &info);	 
	
	
// pp 179 and 186 in demmel
// or divide and conquer

// int dsygvd_(__CLPK_integer *itype, char *jobz, char *uplo, __CLPK_integer *
//	n, __CLPK_doublereal *a, __CLPK_integer *lda, __CLPK_doublereal *b, __CLPK_integer *ldb, 
//	__CLPK_doublereal *w, __CLPK_doublereal *work, __CLPK_integer *lwork, __CLPK_integer *iwork, 
//	__CLPK_integer *liwork, __CLPK_integer *info); 
//	dsygvd_(&type, &V, &U, &matSize, K, &matSize, M, &matSize, wr, work, &lwork, iwork, &liwork, &info);	 


#ifdef DEBUG_PRINT	
	printf("n: %d, info: %d\n", matSize,info);

	// Check the output flag. 
	if (info==0)
		printf("\nComputation successful. Optimal lwork: %f\n",work[0]);
	else{
		printf("\nFailed to compute: %d\n", info);
		exit(1);
	}	
#endif

#ifdef DEBUG_PRINT	
	printf("D=[\n");
	for (int  indx = 0; indx < matSize; indx++)
		printf("  % 10.10e \n", wr[indx] );
	printf("];\n");

//#ifdef DEBUG_PRINT
	//int row, col;	
	printf("\nV=[\n");
	for ( row = 0; row < matSize; row++){
		for ( col = 0; col < matSize; col++){		
			printf("%10.10e\t", K[col*matSize+row] );
		}
		printf(";\n");
	}
	printf("];\n");
#endif

	for (int indx = 0; indx < matSize; indx++)
		if (wr[indx] < 0.0) wr[indx] = 0.0;

	//Free up the memory we allocated. 
	free(work);

}

/*!
 same as above but instead computes a specified number of values and vectors
 */
void realDoubleSymmetricGeneralEigenDecompComputeN(double* wr, double* K, double* M, int sizeK, int* range)
{
	integer matSize = (integer) sizeK;

	// Set up temporary variables. 
	integer lwork;
	double *work;
	__CLPK_integer *iwork;

	__CLPK_integer info;					// Return value of algorithm 
	__CLPK_integer *ifail = (__CLPK_integer*) calloc(matSize, sizeof(__CLPK_integer));
	
	// Allocate space for LAPACK to do its computation. 
	lwork =  32*matSize;
	work = (double*)malloc( lwork*sizeof(double) );
	if ( (work==NULL) )	{
		char *memString = "\nBuy more memory!\n";
		printf(memString);
		exit(1);
	}
	
	iwork = (__CLPK_integer*) malloc(5*matSize*sizeof(__CLPK_integer));
	if ( (iwork==NULL) )	{
		char *memString = "\nBuy more memory!\n";
		printf(memString);
		exit(1);
	}	
	
//#ifdef DEBUG_PRINT
	int row, col;
	printf("K=[\n");
	for (row = 0; row < matSize; row++)  {
		for (col = 0; col < matSize; col++){
			printf(" %10.10e\t", K[col*matSize+row] );
		}
		printf(";\n");
	}
	printf("];\n");
	
	printf("M=[\n");
	for (row = 0; row < matSize; row++)  {
		for (col = 0; col < matSize; col++){
			printf(" %10.10e\t", M[col*matSize+row] );
		}
		printf(";\n");
	}
	printf("];\n");
//#endif

	/* 	
	dsygvx computes selected eigenvalues, and optionally, eigen-
     vectors  of a real generalized symmetric-definite eigenprob-
     lem, of the  form  A*x=(lambda)*B*x,   A*Bx=(lambda)*x,   or
     B*A*x=(lambda)*x.   Here A and B are assumed to be symmetric
     and B is also positive definite.  Eigenvalues and  eigenvec-
     tors  can be selected by specifying either a range of values
     or a range of indices for the desired eigenvalues.	 
	//http://docs.sun.com/source/819-3691/dsygvx.html
	*/ 
	__CLPK_integer type  = 1;
	
	__CLPK_integer IL = range[0];
	__CLPK_integer IU = range[1];
	__CLPK_integer m = 0;
	
	double VL = 0.0;
	double VU = 0.0;
	
	double ABSTOL = 2.0*dlamch_(&S);
	
	dsygvx_(	&type, &V, &I, &U, &matSize, 
				K, &matSize, M, &matSize, 
				&VL, &VU, &IL, &IU, &ABSTOL, &m,
				wr, 
				K, &matSize,
				work, &lwork, iwork,
				ifail, &info);	 
		
#ifdef DEBUG_PRINT	
	printf("m: %d, info: %d\n", m,info);

	// Check the output flag. 
	if (info==0)
		printf("\nComputation successful. Optimal lwork: %f\n",work[0]);
	else{
		printf("\nFailed to compute: %d\n", info);
		exit(1);
	}	
#endif

//#ifdef DEBUG_PRINT	
	printf("D=[\n");
	for (int  indx = 0; indx < m; indx++)
		printf("  %11.11e; \n", wr[indx] );
	printf("];\n");

//#ifdef DEBUG_PRINT
	// Show the output. 	
	printf("V = [\n");
	for ( row = 0; row < matSize; row++){
		for ( col = 0; col < m; col++){
			printf("% 10.10e\t", K[col*matSize+row] );
		}
		printf(";\n");
	}
	printf("];\n");
//#endif

	for (int indx = 0; indx < m; indx++)
		if (wr[indx] < 0.0) wr[indx] = 0.0;

	//Free up the memory we allocated. 
	free(work);
	free(iwork);
	free(ifail);
}

/*!
 compute the eigen decomposition for a symmetric matrix in double precision
 wr are eigenvalues, on entry lkl is the stiffness matrix and on exit lkl contains the eigenvectors, the size of the matrix is sizeK
 */
void realDoubleSymmetricEigenDecomp(double* wr, double* lkl, int sizeK)
{
	integer matSize = (integer) sizeK;

	// Set up temporary variables. 
	integer lwork;                // Size of working scratch area 
	double *work;                  // Working scratch area 
	integer info;                 // Return value of algorithm 

	// Allocate space for LAPACK to do its computation. 
	lwork = 32*matSize;
	work = (double*)malloc( lwork*sizeof(double) );
	if ( (work==NULL) )	{
		char *memString = "\nBuy more memory!\n";
		printf(memString);
		exit(1);
	}
	
	
#ifdef DEBUG_PRINT
	int row, col;
	printf("\nMatrix whose eigendecomposition is desired:\n");

	for (row = 0; row < matSize; row++)  {
		for (col = 0; col < matSize; col++)
		printf(" %lf\n", lkl[col*matSize+row] );
		printf("\n");
	}
	printf("\nComputing eigenvectors and eigenvalues.\n");
#endif

	//if (num == sizeK){
		/* 	
		 LAPACK libraries require that all arguments be passed as pointers. 
		 sgeev - compute for an N-by-N real  nonsymmetric  matrix  A,
		 the  eigenvalues  and,  optionally,  the  left  and/or right
		 eigenvectors
		//	http://docs.sun.com/source/819-3691/dsyev.html
		*/ 
		 dsyev_(&V,&U,&matSize, lkl, &matSize, wr, work,&lwork,&info);	 
	//}
	//else {

		/*  sgeevx - compute for an N-by-N real nonsymmetric  matrix  A,
		the  eigenvalues  and,  optionally,  the  left  and/or right
		eigenvectors
		http://docs.sun.com/source/819-3691/sgeevx.html
		*/

	//}

	// come back sorted
	
#ifdef DEBUG_PRINT	
	// Check the output flag. 
	if (info==0)
		printf("\nComputation successful. Optimal lwork: %f\n",work[0]);
	else{
		printf("\nFailed to compute: %d\n", info);
		exit(1);
	}
	
	printf("\nEigenvalues:\n");
	for (int indx = 0; indx < matSize; indx++)
		printf("  % 5.4e \n", wr[indx] );
	
//#ifdef DEBUG_PRINT	// Show the output. 	
	printf("\nRight eigenvectors:\n");
	
	for ( col = 0; col < matSize; col++){
		printf( "\n  Eigenvector %d\n", col+1 );
		for ( row = 0; row < matSize; row++)
			printf("    % 5.4e\n", lkl[col*matSize+row] );
	}
#endif

	for (int indx = 0; indx < matSize; indx++)
		if (wr[indx] < 0.0) wr[indx] = 0.0;

	//Free up the memory we allocated. 
	free(work);

}

#pragma mark ___Accummulate___

/*!
 compute the accumulation operation c =  a^T b a
//kb = kb + kinmtsb'*matmtsb*kinmtsb*wtx*wty*detjacob
*/
inline void	accumulate(			double* c, int size,	// cbnote only for square matrix??
								double* a, int sizea[], 
								double* b, int sizeb,
								double scale)
{
	//	c =  a^T b a	
	// (m x k)(k x k) (k x m)

	double* temp = (double*) calloc(sizea[0]*sizea[1], sizeof(double));	

	int m = sizea[1];
	int k = sizea[0];

/*
	http://www.netlib.org/blas/dgemm.f
*/										
	// (m x k)(k x k)																																												
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,	// temp = a^T b
				m,											// m rows op(A)
				k,											// n cols op(B)
				k,											// k cols op(A), rows op(B)
				scale,										// alpha
				a,											// a
				k,											// k
				b,											// b
				k,											// k
                0.0,										// beta
				temp,										// alpha*op( A )*op( B ) + beta*C
				m);											// m
	
	
	// (m x k) (k x m) 
															
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,	// K = temp*B		
				 m,											// m rows op(A)
				 m,											// n cols op(B)
				 k,											// k cols op(A), rows op(B)
                 1.0,										// alpha
				 temp,										// a 
				 m,											// m
				 a,											// b
				 k,											// k
                 1.0,										// beta
				 c,											// alpha*op( A )*op( B ) + beta*C
				 m);										// m


#ifdef DEBUG_PRINT
	int i,j;
	printf("accumulate=[\n");
	for ( i=0; i<size; i++){
		for ( j=0; j<size; j++){
			printf("%lf\t", c[j+size*i]);
		}
		printf("];\n");
	}
	printf(";\n");		
#endif	
	free(temp);
}

/*!
 test that V^T K V produces a diagonal matrix
 */
void testDiagonalizes(double* V, int sizeV[], double* K, int dofs)
{
	// check V^T K V is a diagonal matrix
	//  
	// (cols x dofs)(dofs x dofs)(dofs x cols)

	double* temp = (double*) calloc(sizeV[1]*sizeV[1], sizeof(double));

	accumulate(	temp, sizeV[1], 
				V, sizeV, 
				K, dofs,
				1.0);	
	
//#ifdef DEBUG_PRINT
	int i,j;
	printf("DD=[");
	for ( i=0; i<sizeV[1]; i++){
		for ( j=0; j<sizeV[1]; j++){
			printf("%10.10e\t",  temp[j*sizeV[1] + i]);
		}
		printf("\n");
	}
	printf("];\n");
//#endif
																								
	free(temp);
}


#endif