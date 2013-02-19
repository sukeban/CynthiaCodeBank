/*
 *  TriAxisSymBaryElasticElement.h
 *  AxisSymmetricSolidApp
 *
 *  Created by Cynthia Bruyns on 11/19/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
//http://he-cda.wiley.com/WileyCDA/HigherEdTitle/productCd-0471648078.html
//http://bcs.wiley.com/he-bcs/Books?action=index&itemId=0471648078&bcsId=3105

#include "GeneralFE.h"
#include "2dElements.h"


#pragma mark ___Stiffness___
void	triAxisSymBaryElasticStiffnessMatrix(	double* K, double* M,
												double* coords,
												AnalysisInfo info,
												double theta,
												int n)
{	
	MaterialInfo minfo = info.materialinfo;
	GeomInfo ginfo = info.geominfo;
	IntegrationInfo iinfo = info.intinfo;
	
	int sizeK = ginfo.numelnodes*ginfo.nodedofs; 
	
	int sizeD = 6;
	double* D = (double*) calloc(sizeD*sizeD,sizeof(double));
	materialMatrix(D, minfo);	 // solid
	
	double* N   = (double*) calloc(3, sizeof(double));
	double* dNdRdZ = (double*) calloc(3*2, sizeof(double)); 
	
	int sizeB[2];
	sizeB[0] = 6;	sizeB[1] = 9;
	double* B = (double*) calloc(sizeB[0]*sizeB[1],sizeof(double));
	
	double center[2];
	getTriCenter(coords, center[0], center[1]);	
	triangularBaryShapeFunctions(N, dNdRdZ, center[0], center[1], coords, iinfo.order);	
	double r = N[0]* coords[0] + N[1]*coords[1] + N[2]*coords[2];	
		
	if (n==0) 
		axisSymConstMatrix(3, B, dNdRdZ, N, r, theta, n);	
	else
		antiAxisSymConstMatrix(3, B, dNdRdZ, N, r, theta, n);	
		
		
	if (iinfo.lumpedmass)
		triAxisSymElasticLumpedMassMatrix(M, coords, minfo);
	else
		triAxisSymElasticConsistentMassMatrix(M, coords, minfo);
	
		
	double area = triArea(coords);
	accumulateStiffness(K, sizeK, B, sizeB, D, sizeD, 1.0, 1.0, 1.0,r*2.0*M_PI*area);

	free(D);
	free(N), free(dNdRdZ);
	free(B);

}

#pragma mark ___Assemble___
void assembleTriAxisSymBaryElasticMatricies(	double* K, double* M,		
												double* verts, int* inds, 
												AnalysisInfo info,
												double theta,
												int n)
{
	MaterialInfo		minfo = info.materialinfo;
	GeomInfo			ginfo = info.geominfo;
	IntegrationInfo		iinfo = info.intinfo;
	
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
		triAxisSymBaryElasticStiffnessMatrix(Ke, Me, coords, info, theta, n);							
				
					
#ifdef DEBUG_PRINT
		
		printf("\n\n\n\n\n");
		for ( i=0; i<edofs; i++)
			for ( j=0; j<edofs; j++)
				printf("Ke[%ld][%ld]: %lf\n", i,j, Ke[i*edofs + j]);					

	int i,j;
																					
		for ( i=0; i<edofs; i++)
			for ( j=0; j<edofs; j++)
				printf("Me[%ld][%d]: %lf\n", i, j, Me[i + j*edofs]);						
#endif	
		assemble(K, M, el, Ke, Me, ginfo, inds);
	}

	free(Ke), free(Me);
	free(coords);
	
}