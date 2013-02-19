/*
 *  QuadAxisSymElasticElement.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 11/29/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */


//http://he-cda.wiley.com/WileyCDA/HigherEdTitle/productCd-0471648078.html
//http://bcs.wiley.com/he-bcs/Books?action=index&itemId=0471648078&bcsId=3105

#include "GeneralFE.h"
#include "ShapeFunctions.h"



#pragma mark __Mass___
void	quadAxisSymElasticLumpedMassMatrix(double* M, double* coords, MaterialInfo minfo)
{
	double area = quadArea(coords);		
	double rbar = (coords[0] + coords[1] + coords[2] + coords[3])/4.0;
	double mass = 2.0*M_PI*rbar*area*minfo.density/4.0; 

	int i;	
	int sizeM = 12;
	for ( i=0; i<sizeM; i++)
		M[i+i*sizeM] += mass;
		
	return;
}

void	quadAxisSymElasticConsistentMassMatrix(double* M, double* N, double* coords, MaterialInfo minfo)
{
	double area = quadArea(coords);
	double rbar = (coords[0] + coords[1] + coords[2] + coords[3])/4.0;
	double prefix = 2.0*M_PI*rbar*minfo.density*area;
//	printf("prefix: %lf\n", prefix);
	
	double* NL = (double*) calloc(3*12, sizeof(double));
	int sizeNL[2];
	sizeNL[0] = 3;
	sizeNL[1] = 12;
	
/*	for (int i=0; i< 4; i++)
		printf("N[%d]: %lf\n", i, N[i]);
*/	
	for (int i=0; i< 4; i++){
		
		NL[i*9  ] = N[i];  // 0x 3  6    9x 12  15    18x 21  24    27x 30  33
		NL[i*9+4] = N[i];  // 1  4x 7    10 13x 16    19  22x 25    28  31x 34
		NL[i*9+8] = N[i];  // 2  5  8x   11 14  17x   20  23  26x   29  32  35x
		
	}
	
/*	printf("NL\n");
	for (int i=0; i<3; i++){
		for (int j=0; j< 12; j++){
			printf(" %lf\t", NL[i+j*3]);
		}
		printf("\n");
	}
	
	printf("Befoer Mass\n");
	for (int i=0; i<12; i++){
		for (int j=0; j< 12; j++){
			printf(" %lf\t", M[i+j*12]);
		}
		printf("\n");
	}
*/		
	accumulateMass(M, 12, NL, sizeNL, 1, 1, 1, prefix);

/*	printf("Mass\n");
	for (int i=0; i<12; i++){
		for (int j=0; j< 12; j++){
			printf(" %lf\t", M[i+j*12]);
		}
		printf("\n");
	}
*/
	
	free(NL);
	return;
}

#pragma mark ___Stiffness___
void	quadAxisSymElasticStiffnessMatrix(	double* K, double* M,
											double* coords,
											AnalysisInfo info,
											double theta,
											int n)
{	
	MaterialInfo	minfo = info.materialinfo;
	GeomInfo		ginfo = info.geominfo;
	IntegrationInfo iinfo = info.intinfo;
	
	int sizeK = ginfo.numelnodes*ginfo.nodedofs;
	
	double* gp  = (double*) calloc(iinfo.numgauss, sizeof(double));
	double* w = (double*) calloc(iinfo.numgauss, sizeof(double));
	getGaussPoints( gp, w, iinfo.numgauss);
	
	int sizeD = 6;
	double* D = (double*) calloc(sizeD*sizeD,sizeof(double));
	materialMatrix(D, minfo);	//solid
	
	double* J = (double*) calloc(2*2, sizeof(double));
	double* invJ = (double*) calloc(2*2, sizeof(double));
	double detJ = 0.0;

	double* N   = (double*) calloc(4, sizeof(double));
	double* dNdSdN = (double*) calloc(4*2, sizeof(double));
		
	int sizeB[2];
	sizeB[0] = 6;	
	sizeB[1] = 12;

	double* B = (double*) calloc(sizeB[0]*sizeB[1],sizeof(double));
	
	double* dNdXdY = (double*) calloc(4*2, sizeof(double)); 
	
	double r;

	int i,j;
	for ( i=0; i<iinfo.numgauss; i++){
		for ( j=0; j<iinfo.numgauss; j++){

			quadrilateralShapeFunctions(N, dNdSdN, gp[i], gp[j], iinfo.order);	
			r = N[0]* coords[0] + N[1]*coords[1] + N[2]*coords[2]+ N[3]*coords[3];	
			
			twoDJacobian(4, J, invJ, detJ, dNdSdN, coords);			
			twoDPhysicalDerivs(4, dNdXdY, invJ, dNdSdN);		

			if (n==0) 
				axisSymConstMatrix(4, B, dNdXdY, N, r, theta, n);	
			else
				antiAxisSymConstMatrix(4, B, dNdXdY, N, r, theta, n);
			

			if (iinfo.lumpedmass)
				quadAxisSymElasticLumpedMassMatrix(M, coords, minfo);
			else
				quadAxisSymElasticConsistentMassMatrix(M, N, coords, minfo);
				
			accumulateStiffness(K, sizeK, B, sizeB, D, sizeD, w[i], w[j], 1.0, detJ*r*M_PI);
			
	
			}
	}

#ifdef DEBUG_PRINT
		
		printf("\n\n\n\n\n");
		for ( i=0; i<sizeK; i++)
			for ( j=0; j<sizeK; j++)
				printf("Ke[%ld][%ld]: %lf\n", i,j, K[i*sizeK + j]);					
#endif	

#ifdef DEBUG_PRINT
																						
		for ( i=0; i<sizeK; i++)
			for ( j=0; j<sizeK; j++)
				printf("Me[%ld][%ld]: %lf\n", i, j, M[i + j*sizeK]);						
#endif


	free(gp), free(w);
	free(D), free(B);
	free(J), free(invJ);
	free(N), free(dNdSdN), free(dNdXdY);

}

#pragma mark ___Assemble___
void assembleQuadAxisSymElasticMatricies(	double* K, double* M,		
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
		
		quadAxisSymElasticStiffnessMatrix(Ke, Me, coords, info, theta, n);							
							
#ifdef DEBUG_PRINT
		int i,j;
		printf("\n\n\n\n\n");
		for ( i=0; i<edofs; i++)
			for ( j=0; j<edofs; j++)
				printf("Ke[%ld][%ld]: %lf\n", i,j, Ke[i*edofs + j]);					
#endif	
#ifdef DEBUG_PRINT
																						
		for ( i=0; i<edofs; i++)
			for ( j=0; j<edofs; j++)
				printf("Me[%ld][%ld]: %lf\n", i, j, Me[i + j*edofs]);						
#endif
	
		assemble(K, M, el, Ke, Me, ginfo, inds);
	}

	free(Ke), free(Me);
	free(coords);
	
}