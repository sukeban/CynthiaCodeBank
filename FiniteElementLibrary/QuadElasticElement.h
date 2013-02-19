/*
 *  QuadElasticElement.h
 *  Finite Element
 *
 *  Created by Cynthia Bruyns on 1/8/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _QuadElastic_Header_
#define _QuadElastic_Header_

#include <stdio.h>
#include <stdlib.h> 

#include "ShapeFunctions.h"
#include "GeneralFE.h"

#include "2dElements.h"


void	quadElasticStiffnessMatrix(		double* K,
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
	
	int sizeD = 3;
	double* D = (double*) calloc(sizeD*sizeD,sizeof(double));
	materialMatrix(D, minfo);	
	
	double* J = (double*) calloc(2*2, sizeof(double));
	double* invJ = (double*) calloc(2*2, sizeof(double));
	double detJ = 0.;

	double* N   = (double*) calloc(4, sizeof(double));
	double* dNdSdN = (double*) calloc(4*2, sizeof(double));
		
	int sizeB[2];
	sizeB[0] = 3;	sizeB[1] = 4*2;
	double* B = (double*) calloc(sizeB[0]*sizeB[1],sizeof(double));
	
	double* dNdXdY = (double*) calloc(2*4, sizeof(double)); 
	
	int i,j;
	for ( i=0; i<iinfo.numgauss; i++){
		for ( j=0; j<iinfo.numgauss; j++){
		
			quadrilateralShapeFunctions(N, dNdSdN, gp[i], gp[j], iinfo.order);
			
			twoDJacobian(4, J, invJ, detJ, dNdSdN, coords);			
			twoDPhysicalDerivs(4, dNdXdY, invJ, dNdSdN);		
				
			twoDConstMatrix(4, B, dNdXdY);	
			
			accumulateStiffness(K, sizeK, B, sizeB, D, sizeD, w[i], w[j], 1.0, detJ);
		}		
	}
	
	free(gp), free(w);
	free(D);
	free(J), free(invJ);
	free(N), free(dNdSdN);
	free(B), free(dNdXdY);
}



void assembleQuadElasticMatricies(	double* K, double* M, 
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
	int chunk = edofs*edofs*sizeof(double);
	for ( el=0; el<ginfo.numelements; el++){ 
	
		memset(Ke,0,chunk);
		memset(Me,0,chunk);
			
		get2dCoords(coords,el,verts,inds, ginfo); 
		
#ifdef DEBUG_PRINT
		printf("\n\n\nelement: %ld\n", el);
#endif	
		quadElasticStiffnessMatrix(  Ke, coords, info);							
		if (iinfo.lumpedmass)
			quadElasticLumpedMassMatrix( Me, coords, minfo);
		else
			quadElasticConsistentMassMatrix( Me, coords, minfo);
			
					
#ifdef DEBUG_PRINT
		int i,j;
		printf("\n\n\n\n\n");
		for ( i=0; i<edofs; i++){
			for ( j=0; j<edofs; j++){
				printf("Ke[%ld][%ld]: %lf\n", i,j, Ke[i*edofs + j]);		
			}
			printf("\n");
		}
																					
		for ( i=0; i<edofs; i++){
			for ( j=0; j<edofs; j++){
				printf("Me[%ld][%ld]: %lf\n", i,j, Me[i*edofs + j]);
			}
			printf("\n");
		}
					
#endif	
		assemble(K, M, el, Ke, Me, ginfo, inds);
	}
	
	//int dofs = ginfo.numnodes*ginfo.nodedofs;

	free(Ke), free(Me);
	free(coords);
}



#endif