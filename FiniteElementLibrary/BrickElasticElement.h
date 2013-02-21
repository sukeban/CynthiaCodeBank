/*
 *  BrickElasticElement.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 9/5/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _BrickElastic_Header_
#define _BrickElastic_Header_

#include <stdio.h>
#include <stdlib.h> 

#include "ShapeFunctions.h"
#include "GeneralFE.h"
#include "3dElements.h"

/*!
 compute the stiffness matrix for a 3D brick element
 */
void	brickElasticStiffnessMatrix(	float* K,
										float* coords,
										AnalysisInfo info) 
{	
	MaterialInfo minfo = info.materialinfo;
	GeomInfo ginfo = info.geominfo;
	IntegrationInfo iinfo = info.intinfo;
		
	int sizeK = ginfo.numelnodes*ginfo.nodedofs; 
	
	float* gp  = (float*) calloc(iinfo.numgauss, sizeof(float));
	float* w = (float*) calloc(iinfo.numgauss, sizeof(float));
	getGaussPoints( gp, w, iinfo.numgauss);
	
	int sizeD = 6;
	float* D = (float*) calloc(sizeD*sizeD,sizeof(float));
	materialMatrix(D, minfo);	
	
	float* J = (float*) calloc(3*3, sizeof(float));
	float* invJ = (float*) calloc(3*3, sizeof(float));
	float detJ = 0.;

	float* N   = (float*) calloc(8, sizeof(float));
	float* dNdSdN = (float*) calloc(8*3, sizeof(float));
		
	int sizeB[2];
	sizeB[0] = 6;	sizeB[1] = 3*8;
	float* B = (float*) calloc(sizeB[0]*sizeB[1],sizeof(float));
	
	float* dNdXdY = (float*) calloc(sizeB[1], sizeof(float)); 
	
	int i,j,k;
	for ( i=0; i<iinfo.numgauss; i++){
		for ( j=0; j<iinfo.numgauss; j++){
			for ( k=0; k<iinfo.numgauss; k++){
		
				brickShapeFunctions(N, dNdSdN, gp[i], gp[j], gp[k],iinfo.order);
				
				brickJacobian(J, invJ, detJ, dNdSdN, coords);			
				brickPhysicalDerivs(dNdXdY, invJ, dNdSdN);		
					
				brickConstMatrix(B, dNdXdY);	
				
				accumulateStiffness(K, sizeK, B, sizeB, D, sizeD, w[i], w[j], w[j], detJ);
			}
		}
	}
	
	free(gp), free(w);
	free(D), free(B);
	free(J), free(invJ);
	free(N), free(dNdSdN);
	free(dNdXdY);
}

/*!
 compute the lumped mass matrix for a 3D elastic element
 */
void	brickElasticLumpedMassMatrix(float* M, float* coords, MaterialInfo minfo)
{
	float vol = brickVolume(coords);
	float mass = vol*minfo.density/8.;
	int i;
	for ( i=0; i<8*3; i++)
		M[i] += mass;
}

/*!
 compute the mass and stiffness matricies for a 3d object
 */
void assembleBrickElasticMatricies(	float* K, float* M,
									float* verts, int* inds, 
									AnalysisInfo info)
{
	MaterialInfo	minfo = info.materialinfo;
	GeomInfo		ginfo = info.geominfo;
	IntegrationInfo    iinfo = info.intinfo;
	
	int edofs = ginfo.numelnodes*ginfo.nodedofs;
	float* Ke = (float*) calloc(edofs*edofs, sizeof(float));
	float* Me = (float*) calloc(edofs, sizeof(float));	
	
	float* coords = (float*) calloc(ginfo.numelnodes*ginfo.dimension,sizeof(float));
	
	int el;
	for ( el=0; el<ginfo.numelements; el++){ 
	
		memset(Ke,0,edofs*edofs*sizeof(float));
		memset(Me,0,edofs*sizeof(float));
			
		get3dCoords(coords,el,verts,inds, ginfo); 
		
#ifdef DEBUG_PRINT
		printf("\n\n\nelement: %ld\n", el);
#endif	
		brickElasticStiffnessMatrix(  Ke, coords, info);							
		brickElasticLumpedMassMatrix( Me, coords, minfo);
					
#ifdef DEBUG_PRINT
		int i,j;
		printf("\n\n\n\n\n");
		for ( i=0; i<edofs; i++)
			for ( j=0; j<edofs; j++)
				printf("Ke[%ld][%ld]: %lf\n", i,j, Ke[i*edofs + j]);					
																					
		for ( i=0; i<edofs; i++)
			printf("Me[%ld]: %lf\n", i, Me[i]);	
					
#endif	
		assemble(K, M, el, Ke, Me, ginfo, inds);
	}

	//int dofs = ginfo.numnodes*ginfo.nodedofs;	// cbnote not needed for symmetric eigendecomp
	//symmetricTest(K, dofs);

	free(Ke), free(Me);
	free(coords);
}


#endif