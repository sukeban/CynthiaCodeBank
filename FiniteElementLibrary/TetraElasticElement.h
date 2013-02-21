/*
 *  TetraElasticElement.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 9/5/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _TetraElastic_Header_
#define _TetraElastic_Header_

#include <stdio.h>
#include <stdlib.h> 

#include "ShapeFunctions.h"
#include "GeneralFE.h"
#include "3dElements.h"

/*!
 compute the shape functions for a tetrehedral element given a set of coordinates
 return them in N
 
 x,y,z are unused at the moment
 */
void makeN(double* N, double* coords, double x, double y, double z)
{
	double x1=coords[0+0]; double y1=coords[0+4]; double z1=coords[0+8];
	double x2=coords[1+0]; double y2=coords[1+4]; double z2=coords[1+8];
	double x3=coords[2+0]; double y3=coords[2+4]; double z3=coords[2+8];
	double x4=coords[3+0]; double y4=coords[3+4]; double z4=coords[3+8];
	
	double a1 = -(x4*y3*z2) + x3*y4*z2 + x4*y2*z3 - x2*y4*z3 - x3*y2*z4 +  x2*y3*z4;
	double a2 = x4*y3*z1 - x3*y4*z1 - x4*y1*z3 + x1*y4*z3 + x3*y1*z4 - x1*y3*z4;
	double a3 = -(x4*y2*z1) + x2*y4*z1 + x4*y1*z2 - x1*y4*z2 - x2*y1*z4 + x1*y2*z4;
	double a4 = x3*y2*z1 - x2*y3*z1 - x3*y1*z2 + x1*y3*z2 + x2*y1*z3 - x1*y2*z3;

	double b1 = y4*(-z2 + z3) + y3*(z2 - z4) + y2*(-z3 + z4);
	double b2 = y4*(z1 - z3) + y1*(z3 - z4) + y3*(-z1 + z4);
	double b3 = y4*(-z1 + z2) + y2*(z1 - z4) + y1*(-z2 + z4);
	double b4 = y3*(z1 - z2) + y1*(z2 - z3) + y2*(-z1 + z3);
	
	double c1 = x4*(z2 - z3) + x2*(z3 - z4) + x3*(-z2 + z4);
	double c2 = x4*(-z1 + z3) + x3*(z1 - z4) + x1*(-z3 + z4);
	double c3 = x4*(z1 - z2) + x1*(z2 - z4) + x2*(-z1 + z4);
	double c4 = x3*(-z1 + z2) + x2*(z1 - z3) + x1*(-z2 + z3);
	
	double d1 = x4*(-y2 + y3) + x3*(y2 - y4) + x2*(-y3 + y4);
	double d2 = x4*(y1 - y3) + x1*(y3 - y4) + x3*(-y1 + y4);
	double d3 = x4*(-y1 + y2) + x2*(y1 - y4) + x1*(-y2 + y4);
	double d4 = x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3);

	double V = tetraVolume(coords);
	double scale =(1.0/(6.0*V));

/*	
	double c[3];
	getTetraCenter(coords, c);
	double x = c[0]; 
	double y = c[1]; 
	double z = c[2];
*/	
	double n[4];
	
	n[0] = scale*(a1+x*b1+y*c1+z*d1);
	n[1] = scale*(a2+x*b2+y*c2+z*d2);
	n[2] = scale*(a3+x*b3+y*c3+z*d3);
	n[3] = scale*(a4+x*b4+y*c4+z*d4);
	
#ifdef DEBUG_PRINT
	printf("n[0]: %lf n[1]: %lf n[2]: %lf n[3]: %lf\n", n[0], n[1], n[2], n[3]);
#endif
	
	int i;
	for (i=0; i<4; i++){
		N[i*9  ] = n[i];
		N[i*9+4] = n[i];
		N[i*9+8] = n[i];
	}
	
#ifdef DEBUG_PRINT
	int j;
	for ( i=0; i<3; i++){
		for (j=0; j< 12; j++){
			printf("%lf\t",  N[i+j*3]);
		}
		printf(";\n");
	}
#endif
}

/*!
 compute the matrix B which satifies the equation K = B^T D B
 */
void makeB(double* B, double* coords)
{
	double x1=coords[0+0]; double y1=coords[0+4]; double z1=coords[0+8];
	double x2=coords[1+0]; double y2=coords[1+4]; double z2=coords[1+8];
	double x3=coords[2+0]; double y3=coords[2+4]; double z3=coords[2+8];
	double x4=coords[3+0]; double y4=coords[3+4]; double z4=coords[3+8];

	double b1 = y4*(-z2 + z3) + y3*(z2 - z4) + y2*(-z3 + z4);
	double b2 = y4*(z1 - z3) + y1*(z3 - z4) + y3*(-z1 + z4);
	double b3 = y4*(-z1 + z2) + y2*(z1 - z4) + y1*(-z2 + z4);
	double b4 = y3*(z1 - z2) + y1*(z2 - z3) + y2*(-z1 + z3);
	
	double c1 = x4*(z2 - z3) + x2*(z3 - z4) + x3*(-z2 + z4);
	double c2 = x4*(-z1 + z3) + x3*(z1 - z4) + x1*(-z3 + z4);
	double c3 = x4*(z1 - z2) + x1*(z2 - z4) + x2*(-z1 + z4);
	double c4 = x3*(-z1 + z2) + x2*(z1 - z3) + x1*(-z2 + z3);
	
	double d1 = x4*(-y2 + y3) + x3*(y2 - y4) + x2*(-y3 + y4);
	double d2 = x4*(y1 - y3) + x1*(y3 - y4) + x3*(-y1 + y4);
	double d3 = x4*(-y1 + y2) + x2*(y1 - y4) + x1*(-y2 + y4);
	double d4 = x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3);
	
	int size = 12*6;
	B[0+0] = b1;	B[0+12] = 0.0;	B[0+12*2] = 0.0;	B[0+12*3] = c1;		B[0+12*4] = 0.0;	B[0+12*5] = d1;
	B[1+0] = 0.0;	B[1+12] = c1;	B[1+12*2] = 0.0;	B[1+12*3] = b1;		B[1+12*4] = d1;		B[1+12*5] = 0.0;
	B[2+0] = 0.0;	B[2+12] = 0.0;	B[2+12*2] = d1;		B[2+12*3] = 0.0;	B[2+12*4] = c1;		B[2+12*5] = b1;
	B[3+0] = b2;	B[3+12] = 0.0;	B[3+12*2] = 0.0;	B[3+12*3] = c2;		B[3+12*4] = 0.0;	B[3+12*5] = d2;
	B[4+0] = 0.0;	B[4+12] = c2;	B[4+12*2] = 0.0;	B[4+12*3] = b2;		B[4+12*4] = d2;		B[4+12*5] = 0.0;
	B[5+0] = 0.0;	B[5+12] = 0.0;	B[5+12*2] = d2;		B[5+12*3] = 0.0;	B[5+12*4] = c2;		B[5+12*5] = b2;
	B[6+0] = b3;	B[6+12] = 0.0;	B[6+12*2] = 0.0;	B[6+12*3] = c3;		B[6+12*4] = 0.0;	B[6+12*5] = d3;
	B[7+0] = 0.0;	B[7+12] = c3;	B[7+12*2] = 0.0;	B[7+12*3] = b3;		B[7+12*4] = d3;		B[7+12*5] = 0.0;
	B[8+0] = 0.0;	B[8+12] = 0.0;	B[8+12*2] = d3;		B[8+12*3] = 0.0;	B[8+12*4] = c3;		B[8+12*5] = b3;
	B[9+0] = b4;	B[9+12] = 0.0;	B[9+12*2] = 0.0;	B[9+12*3] = c4;		B[9+12*4] = 0.0;	B[9+12*5] = d4;
	B[10+0] = 0.0;	B[10+12] = c4;	B[10+12*2] = 0.0;	B[10+12*3] = b4;	B[10+12*4] = d4;	B[10+12*5] = 0.0;
	B[11+0] = 0.0;	B[11+12] = 0.0;	B[11+12*2] = d4;	B[11+12*3] = 0.0;	B[11+12*4] = c4;	B[11+12*5] = b4;
	
	
	double V = tetraVolume(coords);
	double scale = (1.0/(6.0*V));
	
	int i;
	for (i=0; i<size; i++){
		B[i] *= scale;
	}
	
	
#ifdef DEBUG_PRINT
	int j;
	for ( i=0; i<12; i++){
		for ( j=0; j<6; j++){
			printf("B[%d][%d]: %lf", i,j, B[i+j*12]);
		}
		printf("\n");
	}
#endif

}

/*!
 compute the stiffness matrix of a tetrahedral element
 */
void tetraElasticStiffnessMatrix(	double* K,
									double* coords,
									AnalysisInfo info)
{		
	MaterialInfo	minfo = info.materialinfo;	
	
	int sizeD = 6;
	double* D = (double*) calloc(sizeD*sizeD,sizeof(double));
	materialMatrix(D, minfo);	
		
	int sizeB[2];	
	sizeB[0] = 12; //before
	sizeB[1] = 6;
	
	double* B = (double*) calloc(sizeB[0]*sizeB[1],sizeof(double));
	makeB(B, coords);
	
	double*BT = (double*) calloc(sizeB[1]*sizeB[0], sizeof(double));
	realDoubleMatrixOutOfPlaceTranspose(B, sizeB, BT); // this changes the size, so this to be like the other elements and reuse accumulate code
	sizeB[0] = 6;
	sizeB[1] = 12;
	
	double V = tetraVolume(coords);
	
	int sizeK = 12;
	accumulate(			K, sizeK, 
						BT, sizeB, 
						D, sizeD,
						V);
	
#ifdef DEBUG_PRINT
	int i,j;
	printf("\n\n\n");
	for ( i=0; i<sizeK; i++){
		for ( j=0; j<sizeK; j++){
			printf("K[%ld][%ld]: %lf\n", i,j,K[j*sizeK+i]);
		}
		printf("\n");
	}		
#endif	

	free(D);
	free(B);
	free(BT);
}

/*!
 
 compute the consistent stiffness matrix of a tetrahedral ekement
 
//http://icl.cs.utk.edu/lapack-forum/viewtopic.php?p=769&sid=51047c7c05a5be0be278c20cc118eecf
// page 362 rao
 */
void	tetraElasticConsistentMassMatrix(double* M, double* coords, MaterialInfo minfo)
{
	double V = tetraVolume(coords);
	
	double diag = minfo.density*V/10.0;
	double offdiag = minfo.density*V/20.0;
	
	int i,j,k;
	int rowblock,colblock;
	
	
	int sizeM = 12;
	
	for ( i=0; i<4; i++){
		//printf("i: %d\n", i);		
		rowblock = i*3;
		
		for ( j=0; j<4; j++){
			//printf("j: %d\n", j);
			
			colblock = j*3;
			
			for ( k=0; k<3; k++){
				//printf("k: %d\n", k);
				
					if (i==j)
						M[(rowblock + k) + (colblock + k)*sizeM] = diag;
					else
						M[(rowblock + k) + (colblock + k)*sizeM] = offdiag;
				
			}		
		}	
	}

#ifdef DEBUG_PRINT
	printf("\n\n\n");
	for ( i=0; i<sizeM; i++){
		for ( j=0; j<sizeM; j++){
			printf("M[%ld][%ld]: %lf\t\t", i,j,M[j*sizeM+i]);
		}
		printf("\n");
	}
#endif

}

/*!
 compute the lumped mass matrix for a tetrahedral element
*/
void	tetraElasticLumpedMassMatrix(double* M, double* coords, MaterialInfo minfo)
{
	double vol = tetraVolume(coords);
	
	double mass = vol*minfo.density/12.;  
	
	int i;	
	int length = 4*3; // cbnote numverts * numdofs
	for ( i=0; i<length; i++) 
		M[i + i*length] = mass; 
		
#ifdef DEBUG_PRINT
	int j;
	printf("\n\n\n");
	for ( i=0; i<length; i++){
		for ( j=0; j<length; j++){
			printf("M[%ld][%ld]: %lf\n", i,j,M[j*length+i]);
		}
		printf("\n");
	}
#endif
}

/*!
 assemble the stiffness and mass matricies for an object made of tetrahedral elements
 */
void assembleTetraElasticMatricies(	double* K, double* M, 
									double* verts, int* inds, 
									AnalysisInfo info)
{
	MaterialInfo	minfo = info.materialinfo;
	GeomInfo		ginfo = info.geominfo;
	IntegrationInfo iinfo = info.intinfo;
	
	int edofs = ginfo.numelnodes*ginfo.nodedofs;
	int chunk = edofs*edofs;
	double* Ke = (double*) calloc(chunk, sizeof(double));
	double* Me = (double*) calloc(chunk, sizeof(double));	
	
	double* coords = (double*) calloc(ginfo.numelnodes*ginfo.dimension,sizeof(double));
	
	int el;
	for ( el=0; el<ginfo.numelements; el++){ 
	
		memset(Ke,0,chunk*sizeof(double));
		memset(Me,0,chunk*sizeof(double));
			
		get3dCoords(coords,el,verts,inds, ginfo); 
		//if (checkTetraElement(coords)){
		//	//printf("tetra: %d is facing inwards\n", el);
		//}
			
#ifdef DEBUG_PRINT
		printf("\n\n\nelement: %ld\n", el);
#endif	
		tetraElasticStiffnessMatrix(  Ke, coords, info);	
								
		if (iinfo.lumpedmass)
			tetraElasticLumpedMassMatrix( Me, coords, minfo); 
		else
			tetraElasticConsistentMassMatrix( Me, coords, minfo);
					
#ifdef DEBUG_PRINT
		int i,j;
		printf("\n\n\n\n\n");
		for ( i=0; i<edofs; i++)
			//for ( j=0; j<edofs; j++)
				if (Ke[i + i*edofs] == 0.0){
					printf("hey tetra: %d is bad\n", el);
					printf("Ke[%ld][%ld]: %lf\n", i,j, Ke[i + j*edofs]);		
				}
																					
		for ( i=0; i<edofs; i++)
			//for ( j=0; j<edofs; j++)
				if (Me[i + i*edofs] == 0.0){
					printf("hey tetra: %d is bad\n", el);
					printf("Me[%ld][%ld]: %lf\n", i,j, Me[i + j*edofs]);		
				}					
#endif	
		assemble(K, M, el, Ke, Me, ginfo, inds);
	}
	
	//printf("testing original matrix\n");
	//int dofs = ginfo.numnodes*ginfo.nodedofs;
	//symmetricTest(K, dofs);

	
	free(Ke), free(Me);
	free(coords);
}


#endif