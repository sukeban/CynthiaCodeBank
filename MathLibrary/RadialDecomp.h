/*
 *  RadialDecomp.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 1/9/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 
 
 *	Post discretization eigendecomposition
 
 
 */


#include <math.h>

#include <vecLib/clapack.h>
#include <vecLib/vDSP.h>


#include "CLPKComplexOpsCM.h"
#include "CLPKRealMatrixOpsCM.h"

/*!
 build the permutation matrix
 */
inline void buildPermutationMatrix(__CLPK_real* perm, int rankPerm, int numLong, int numLat)
{
	// make the permutaion matrix
	int i;
	for ( i=0; i<rankPerm-1; i++){
		//printf("i: %ld\n", i);
		perm[i*rankPerm+i+rankPerm] = 1.0;
	}
	
	i= rankPerm-1;
	perm[i] = 1.0;
	
	/*for ( i=0; i<rankPerm; i++)
		for ( j=0; j<rankPerm; j++)
			printf("perm[%ld][%ld]: %lf\n", i,j , perm[j*rankPerm+i]);
	*/
}

/*!
 build the eigenvector permutation matrix
 */
inline void buildEvecPermutationMatrix(__CLPK_complex* evecPerm, int rankPerm, __CLPK_complex omega)
{
	int i,j,k;
	float normRow;
	__CLPK_complex* row = (__CLPK_complex*) calloc(rankPerm, sizeof(__CLPK_complex));	
	for ( i=0; i<rankPerm ;i++){
	
		__CLPK_complex omega_k = complexSinglePow(omega,i);
		
		//printf("\n\n\n");
		//printf("omega_k[%ld]: %lf + %lfi\n", i, omega_k.r, omega_k.i);
			
		normRow = 0.0;
		for ( j=0; j<rankPerm; j++){
			row[j] = complexSinglePow(omega_k,j);
			normRow += pow(row[j].r,2)+pow(row[j].i,2);
			//printf("row[%ld]: %lf + %lfi\n", j, row[j].r, row[j].i);
		}		
		normRow = sqrtf(normRow);
		//printf("normRow: %lf\n", normRow);
		
		for ( k=0; k<rankPerm; k++){
			evecPerm[k*rankPerm + i].r = row[k].r/normRow;
			evecPerm[k*rankPerm + i].i = row[k].i/normRow;
		}
	}
	
	/*for ( i=0; i<rankPerm; i++)
		for ( j=0; j<rankPerm; j++)
			printf("evecPerm[%ld][%ld]: %lf + %lfi\n", i,j,evecPerm[j*rankPerm+i].r, evecPerm[j*rankPerm+i].i);
	*/	
	free(row);
	
}

/*!
 build a rotation matrix
 */
inline void buildRotationMatrix(__CLPK_real* Q, int partitions, int rankBlock)
{
	// build rotation matrix

	float an = (M_PI/2.0)/(2.0*partitions);
	float c  = cos(an);
	float s  = sin(an);
	
	Q[0] = c;
	Q[3] = -s;
	Q[6] = 0;
	
	Q[1] = s;
	Q[4] = c;
	Q[7] = 0;
	
	Q[2] = 0;
	Q[5] = 0;
	Q[8] = 1;
	
}

/*!
 rotate a stiffness matrix in a given direction
 */
inline void realSingleRotateMatrixDiag(	__CLPK_real* Knew, 
										__CLPK_real* Kold, int sizeK, 
										__CLPK_real* rot, int sizeR,
										int direction	)
{	int i,j, k;

#ifdef DEBUG_PRINT
	printf("\n\n\n");
	for ( i=0; i<sizeK; i++)
		for ( j=0; j< sizeK; j++)
			printf("Kold[%ld][%ld]: %lf\n", i,j, Kold[j*sizeK+i]);
#endif
	
	int kindex, rindex;
	if (direction > 0){ //right

		for ( i=0; i<sizeK; i++){		
			int groupr = -1;
			
#ifdef DEBUG_PRINT
			printf("i: %ld\n", i);
#endif
			for ( j=0; j<sizeK; j++){
			
#ifdef DEBUG_PRINT
				printf("j: %ld\n", j);
#endif
				if (j%3 == 0){					
					groupr++; 
					
#ifdef DEBUG_PRINT
					printf("groupr: %ld\n", groupr);	//groupr to jump down in K
#endif
				}
						
				for ( k=0; k<sizeR; k++){	
				
					kindex = i+ k*sizeK  + groupr*sizeK*sizeR; 
					rindex = k + (j%3)*sizeR;
							
#ifdef DEBUG_PRINT
					printf("filling in K[%ld][%ld] with  K[%ld]:%lf with R[%ld]:%lf\n", i, j, kindex, Kold[kindex], rindex, rot[rindex]);
#endif
					Knew[j*sizeK + i] += Kold[kindex] *rot[rindex];
					
				}					
			}			
		}
	
	}
	else {			//left
		int groupr = -1;
		for ( i=0; i<sizeK; i++){
			
#ifdef DEBUG_PRINT
			printf("i: %ld\n", i);
#endif
			if (i%3 == 0){					
				groupr++;
#ifdef DEBUG_PRINT
			printf("groupr: %ld\n", groupr);	//groupr to jump down in K
#endif
			}
				
			for ( j=0; j<sizeK; j++){
#ifdef DEBUG_PRINT
				printf("j: %ld\n", j);
#endif
				for ( k=0; k<sizeR; k++){	
				
					kindex = groupr*sizeR + k  + j*sizeK;
					rindex = k + (i%3)*sizeR;
							
#ifdef DEBUG_PRINT
				printf("filling in K[%ld][%ld]:%lf with R[%ld]:%lf with Kold[%ld]:%lf\n", i, j ,Knew[j*sizeK + i], rindex, rot[rindex], kindex, Kold[kindex]);
#endif
				Knew[j*sizeK + i] += rot[rindex]*Kold[kindex];
					
				}							
			}
		}
	}
	
#ifdef DEBUG_PRINT
	printf("\n\n\n");
	for ( i=0; i<sizeK; i++)
		for ( j=0; j< sizeK; j++)
			printf("Knew[%ld][%ld]: %lf\n", i,j, Knew[i*sizeK+j]);
#endif
}

/*!
 recreate the larger eigenvalue vector and the real eigenvectos from the radial decomposition
 */
inline void reclaimRealEigenValuesRealEigenVectors( __CLPK_real* evals, __CLPK_real* evalBLOCKS, __CLPK_real *evecBLOCKS, int rankBlock, int rankPerm)
{	
	int chunkV = rankBlock*rankBlock; // size of each clump		
	// create space for evevblockcopy
	__CLPK_real* evecBLOCKSCopy = (__CLPK_real*) calloc(rankPerm*chunkV,sizeof(__CLPK_real));
#ifdef DEBUG_PRINT
	if (evecBLOCKSCopy == NULL){
		printf("could not allocate that much memory!\n");	
	}
	else {
		printf("allocated %ld blocks\n", rankPerm*chunkV);
	}	
#endif
		
	int i;
	vDSP_Length length = rankBlock*rankPerm;
	vDSP_Length *inds;        
    inds = (vDSP_Length*)malloc(length*sizeof(vDSP_Length) );
	for (i=0; i<length; i++)
		inds[i] = i;
	
	vDSP_vsorti(evalBLOCKS, inds, NULL, length, 1); 	// this orders the eigenvalues smallest to largest
	
	__CLPK_real* evecBLOCKSCopyPtr; 
	__CLPK_real* evecBLOCKSPtr;

	int bytesChunkV = rankBlock*sizeof(__CLPK_real);
	for (i=0; i<length; i++){
	
		evals[i] = evalBLOCKS[inds[i]];
	
		evecBLOCKSCopyPtr = evecBLOCKSCopy + (i*rankBlock);
		evecBLOCKSPtr = evecBLOCKS + (inds[i]*rankBlock);

		//printf("i: %d index[i]: %d evecBLOCKSCopyPtr: 0x%X evecBLOCKSPtr: 0x%X\n", i, inds[i], evecBLOCKSCopyPtr, evecBLOCKSPtr);
	
		memcpy(evecBLOCKSCopyPtr, evecBLOCKSPtr, bytesChunkV);// get from evecBLOCK with index, write to new location
	}	

	memcpy(evecBLOCKS, evecBLOCKSCopy, rankPerm*bytesChunkV);
	
	free(evecBLOCKSCopy);
	free(inds);	
}

/*!
 recreate the larger eigenvalue vector and the complex eigenvectos from the radial decomposition
 */
inline void reclaimRealEigenValuesComplexEigenVectors( __CLPK_real* evals, __CLPK_real* evalBLOCKS, __CLPK_complex *evecBLOCKS, int rankBlock, int rankPerm)
{	
	int chunkV = rankBlock*rankBlock; // size of each clump		
	// create space for evevblockcopy
	__CLPK_complex* evecBLOCKSCopy = (__CLPK_complex*) calloc(rankPerm*chunkV,sizeof(__CLPK_complex));
#ifdef DEBUG_PRINT
	if (evecBLOCKSCopy == NULL){
		printf("could not allocate that much memory!\n");	
	}
	else {
		printf("allocated %ld blocks\n", rankPerm*chunkV);
	}	
#endif
		
	int i;
	vDSP_Length length = rankBlock*rankPerm;
	vDSP_Length *inds;        
    inds = (vDSP_Length*)malloc(length*sizeof(vDSP_Length) );
	for (i=0; i<length; i++)
		inds[i] = i;
	
	__CLPK_real* evalBLOCKSCopy = (__CLPK_real*) calloc(length, sizeof(__CLPK_real));
	for (i=0; i<length; i++){
		evalBLOCKSCopy[i] = evalBLOCKS[i];
	}

	vDSP_vsorti(evalBLOCKSCopy, inds, NULL, length, 1); 	// this orders the eigenvalues smallest to largest
	
	free(evalBLOCKSCopy);
	
	__CLPK_complex* evecBLOCKSCopyPtr; 
	__CLPK_complex* evecBLOCKSPtr;

	int bytesChunkV = rankBlock*sizeof(__CLPK_complex);
	for (i=0; i<length; i++){
	
		evals[i] = evalBLOCKS[inds[i]]; 
		printf("evals[%d]: %lf\n", i, evals[i]);
	
		evecBLOCKSCopyPtr = evecBLOCKSCopy + (i*rankBlock);
		evecBLOCKSPtr = evecBLOCKS + (inds[i]*rankBlock);

		//printf("i: %d index[i]: %d evecBLOCKSCopyPtr: 0x%X evecBLOCKSPtr: 0x%X\n", i, inds[i], evecBLOCKSCopyPtr, evecBLOCKSPtr);
	
		memcpy(evecBLOCKSCopyPtr, evecBLOCKSPtr, bytesChunkV);// get from evecBLOCK with index, write to new location
	}	

	memcpy(evecBLOCKS, evecBLOCKSCopy, rankPerm*bytesChunkV);
	
	free(evecBLOCKSCopy);
	free(inds);	
}

/*!
 recreate the larger eigenvector matrix from the axis-symmetric represenation
 */
inline void reclaimEigenVectors(__CLPK_real* evecs, __CLPK_complex *evecBLOCKS, int rankBlock, __CLPK_complex* evecPerm, int rankPerm)
{
	// get the eigenvectors
	// evecPERM = kron(evecPerm, eye(rankBlock));
	// V_Kpolar	= evecPERM * evecBLOCKED;
	
	int chunkV = rankBlock*rankBlock;
	int rankV = rankPerm*rankBlock;

	int groupp = -1;
	int indexPerm;
	int indexBlock;
	
	int length = rankV;
	int i,j;

	for (i=0; i<length; i++){
	
		int groupv=-1;
		if (i%rankBlock ==0){
			groupp++;
			//printf("groupp: %ld\n", groupp);
		}
		
		for (j=0; j<length; j++){
		
			if (j%rankBlock ==0){
				groupv++;
				//printf("groupv: %ld\n", groupv);
			}		
			indexPerm = groupp + groupv*rankPerm; 
			indexBlock = (j%rankBlock)*rankBlock + groupv*chunkV + (i%rankBlock);
		
			//printf("filling in evec[%ld][%ld] (or %ld)  from evecPerm[%ld]: %lf and evecBLOCK[%ld]%lf\n", 
			//i, j, j*rankV +i, indexPerm, evecPerm[indexPerm].r, indexBlock, evecBLOCKS[indexBlock].r);
			
			evecs[j*rankV + i] = evecPerm[indexPerm].r*evecBLOCKS[indexBlock].r - evecPerm[indexPerm].i*evecBLOCKS[indexBlock].i;
						
			//printf("evec[%ld][%ld] : %lf\n", i, j, evecs[j*rankV + i]);
		
		}	
	}// no need to noramlize them
}

/*!
 decompose an axis-symmetric object from the stiffness matricies
 */
inline void radialDecompose(float* K11, float* K12, float* K21, int sizeK, 
							float* evecs, float* evals, 
							int numLong, int numLat) //num along length
{	
	int rankPerm    = numLat;		// how many across radius
	int rankBlock   = numLong*6;	// number down
	int chunkV		= rankBlock*rankBlock;

	float angle     = 2.0*M_PI/rankPerm;
	__CLPK_complex omega;
	omega.r = cos(angle);
	omega.i = sin(angle);
	//printf("omega: %lf + %lfi\n", omega.r, omega.i);
	
	__CLPK_complex* evecPerm = (__CLPK_complex*) calloc(rankPerm*rankPerm,sizeof(__CLPK_complex));
	buildEvecPermutationMatrix(evecPerm, rankPerm, omega);
			
	__CLPK_real* Q = (__CLPK_real*) calloc(9, sizeof(__CLPK_real));
	int rankQ = 3;
	buildRotationMatrix(Q, rankPerm, rankBlock);
	
	__CLPK_real* K12Q = (__CLPK_real*) calloc(chunkV, sizeof(__CLPK_real));
	realSingleRotateMatrixDiag(K12Q, K12, rankBlock, Q, rankQ, 1);

	__CLPK_real* QK21 = (__CLPK_real*) calloc(chunkV, sizeof(__CLPK_real));
	realSingleRotateMatrixDiag(QK21, K21, rankBlock, Q, rankQ, -1); 
			
	__CLPK_complex* pencilBlock = (__CLPK_complex*) calloc(chunkV, sizeof(__CLPK_complex));	
	__CLPK_complex* OK11 = (__CLPK_complex*) calloc(chunkV, sizeof(__CLPK_complex));
	__CLPK_complex* OK12 = (__CLPK_complex*) calloc(chunkV, sizeof(__CLPK_complex));
	__CLPK_complex* OK21 = (__CLPK_complex*) calloc(chunkV, sizeof(__CLPK_complex));
	
/*	__CLPK_complex *wr;					// float & imaginary parts of eigenvalues 
    wr = (__CLPK_complex*)malloc( rankBlock*sizeof(__CLPK_complex) );
	__CLPK_complex* evalBLOCKS = (__CLPK_complex*) calloc(rankPerm*rankBlock,sizeof(__CLPK_complex));
	if (evalBLOCKS == NULL){
		printf("could not allocate that much memory!\n");	
	}
	else {
		printf("allocated %ld blocks\n", rankPerm*rankBlock);
	}
*/	

	__CLPK_real *wr;					
    wr = (__CLPK_real*)malloc( rankBlock*sizeof(__CLPK_real) );
	__CLPK_real* evalBLOCKS = (__CLPK_real*) calloc(rankPerm*rankBlock,sizeof(__CLPK_real));
#ifdef DEBUG_PRINT
	if (evalBLOCKS == NULL){
		printf("could not allocate that much memory!\n");	
	}
	else {
		printf("allocated %ld blocks\n", rankPerm*rankBlock);
	}
#endif
	__CLPK_complex *vl, *vr;              /* Left and right eigenvector matrices */
    vr = (__CLPK_complex*)malloc(chunkV*sizeof(__CLPK_complex) );
	vl = NULL; 
	
	__CLPK_complex* evecBLOCKS = (__CLPK_complex*) calloc(rankPerm*chunkV,sizeof(__CLPK_complex));
#ifdef DEBUG_PRINT
	if (evecBLOCKS == NULL){
		printf("could not allocate that much memory!\n");	
	}
	else {
		printf("allocated %ld blocks\n", rankPerm*chunkV);
	}
#endif	
	realSingleMakeSingleComplex(OK11, K11, rankBlock);

	int evalsSize = rankBlock*sizeof(__CLPK_real); 
	int evecsSize = chunkV*sizeof(__CLPK_complex);
	
	double time = currentTime();
	
	int i;
	for ( i=0; i<rankPerm; i++){
		__CLPK_complex omega_k = complexSinglePow(omega,i);
		
#ifdef DEBUG_PRINT
		printf("\n\n\n\nomega_k[%ld]: %lf + %lfi\n", i, omega_k.r, omega_k.i);
#endif 
#ifdef DEBUG_PRINT	
		printf("\n\n\n\n\nOK12: ");
#endif
		realSingleToComplexSingleMultMatrix(OK12, omega_k, K12Q, rankBlock);	
		
#ifdef DEBUG_PRINT
		printf("\n\n\n\n\nOK21: ");
#endif
		__CLPK_complex omega_k_conj = complexSingleGetConjugate(omega_k);
		realSingleToComplexSingleMultMatrix(OK21, omega_k_conj, QK21, rankBlock);
			
		complexSingleAdd3Matricies(pencilBlock, OK11, OK12, OK21, rankBlock);	
		
		complexSingleHermetianEigenDecompose(wr, vr, pencilBlock, rankBlock);

		memcpy(evalBLOCKS+(i*rankBlock), wr, evalsSize);		
		memcpy(evecBLOCKS+(i*chunkV), vr, evecsSize); 
	}

	reclaimRealEigenValuesComplexEigenVectors(evals, evalBLOCKS, evecBLOCKS, rankBlock, rankPerm);
	printf("done getting eigenvalues\n");

	reclaimEigenVectors(evecs, evecBLOCKS, rankBlock, evecPerm, rankPerm);
	printf("done getting eigenvectors\n");
	time = currentTime() - time;
	printf(" time to radial decompose: %lf ms\n", time*1E3); 

	free(pencilBlock);
	
	free(evecPerm);
	free(Q);
	
	free(K12Q), free(QK21);
	free(OK11), free(OK12), free(OK21);
	
	free(evecBLOCKS), free(evalBLOCKS);
	
	free(vl), free(vr),	free(wr);
}
