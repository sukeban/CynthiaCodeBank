/*
 *  QuadAcousticElement.h
 *  AcousticView
 *
 *  Created by Cynthia Bruyns on 4/23/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef _QuadAcousticElement_Header_
#define _QuadAcousticElement_Header_


#include <stdio.h>
#include <stdlib.h> 


#include "ShapeFunctions.h"
#include "GeneralFE.h"

#include "2dElements.h"

#include "CLPKRealMatrixOpsCM.h"


#pragma mark ___Mass___
void	quadAcousticLumpedMassMatrix(double* M, double* coords, MaterialInfo minfo)
{
	double area = quadArea(coords);
#ifdef DEBUG_PRINT	
	printf("area: %lf\n", area);
#endif
	double mass = area/4.0;
	int i;
	for ( i=0; i<4; i++)
		M[i+i*4] += mass;
}


#pragma mark ___Stiffness___

void	quadAcousticStiffnessMatrix(	double* K, double* M,
										double* coords,
										AnalysisInfo info)
{	
	MaterialInfo	minfo = info.materialinfo;
	GeomInfo		ginfo = info.geominfo;
	IntegrationInfo iinfo = info.intinfo;
	
	int sizeK = ginfo.numelnodes*ginfo.nodedofs;
	
	double* gp  = (double*) calloc(iinfo.numgauss, sizeof(double));
	double* w = (double*) calloc(iinfo.numgauss, sizeof(double));
	getGaussPoints( gp, w, iinfo.numgauss);
	
	int sizeD = 2;
	double* D = (double*) calloc(sizeD*sizeD,sizeof(double));
	materialMatrix(D, minfo);	// acoustic
	
	double* J = (double*) calloc(2*2, sizeof(double));
	double* invJ = (double*) calloc(2*2, sizeof(double));
	double detJ = 0.;

	double* N   = (double*) calloc(4, sizeof(double));
	int sizeN[2];
	sizeN[0] = 1;
	sizeN[1] = 4;
	
	double* dNdSdN = (double*) calloc(4*2, sizeof(double));
		
	int sizeB[2];
	sizeB[0] = 2;	sizeB[1] = 4;
	double* B = (double*) calloc(2*4,sizeof(double));
	
	int sizedNdXdY[2];
	sizedNdXdY[0] = 4;	sizedNdXdY[1] = 2;
	double* dNdXdY = (double*) calloc(4*2, sizeof(double)); 
	

	int i,j;
	for ( i=0; i<iinfo.numgauss; i++){
		for ( j=0; j<iinfo.numgauss; j++){
		
			quadrilateralShapeFunctions(N, dNdSdN, gp[i], gp[j], iinfo.order);
			
			twoDJacobian(4, J, invJ, detJ, dNdSdN, coords);			
			twoDPhysicalDerivs(4, dNdXdY, invJ, dNdSdN);
			
			realDoubleMatrixOutOfPlaceTranspose(dNdXdY, sizedNdXdY, B); 
			// this changes the size, so this to be like the other elements and reuse accumulate code
							
			accumulateStiffness(K, sizeK, B, sizeB, D, sizeD, w[i], w[j], 1.0, pow(minfo.soundspeed,2)*detJ);

			if (!iinfo.lumpedmass){
				accumulateMass(M, 4, N, sizeN, w[i], w[j], 1.0, detJ);
				
			}
		}		
	}
	
	free(gp), free(w);
	free(D);
	free(J), free(invJ);
	free(N), free(dNdSdN);
	free(B), free(dNdXdY);
	
}

#pragma mark ___Assemble___

void assembleQuadAcousticMatricies(	double* K, double* M, 
									double* verts, int* inds, 
									AnalysisInfo info)
{
	MaterialInfo	minfo = info.materialinfo;
	GeomInfo		ginfo = info.geominfo;
	IntegrationInfo iinfo = info.intinfo;
	
	int edofs = ginfo.numelnodes*ginfo.nodedofs;
	double* Ke = (double*) calloc(edofs*edofs, sizeof(double));
	double* Me = (double*) calloc(edofs*edofs, sizeof(double));	
	
	double* coords = (double*) calloc(ginfo.numelnodes*ginfo.dimension,sizeof(double));
	
	int el;
	
	for ( el=0; el<ginfo.numelements; el++){ 
	
		memset(Ke,0,edofs*edofs*sizeof(double));
		memset(Me,0,edofs*edofs*sizeof(double));
			
		get2dCoords(coords,el,verts,inds, ginfo); 
		
#ifdef DEBUG_PRINT
		printf("\n\n\nelement: %ld\n", el);
#endif	
		quadAcousticStiffnessMatrix(  Ke, Me, coords, info);
		if (iinfo.lumpedmass)
			quadAcousticLumpedMassMatrix( M, coords, minfo);
			
						
#ifdef DEBUG_PRINT
		int i,j;
		printf("\n\n\n\n\n");
		for ( i=0; i<edofs; i++)
			for ( j=0; j<edofs; j++)
				printf("Ke[%ld][%ld]: %lf\n", i,j, Ke[i*edofs + j]);					
																					
		for ( i=0; i<edofs; i++)
			for ( j=0; j<edofs; j++)
				printf("Me[%ld][%ld]: %lf\n", i, j, Me[i+i*edofs]);	
					
#endif	

		assemble(K, M, el, Ke, Me, ginfo, inds);
	}

	free(Ke), free(Me);
	free(coords);
}



#endif
