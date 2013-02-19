/*
 *  TriAxisSymMapElasticElement.h
 *  AxisSymmetricSolidApp
 *
 *  Created by Cynthia Bruyns on 11/20/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */


//http://he-cda.wiley.com/WileyCDA/HigherEdTitle/productCd-0471648078.html
//http://bcs.wiley.com/he-bcs/Books?action=index&itemId=0471648078&bcsId=3105

#include "GeneralFE.h"
#include "ShapeFunctions.h"

// cbnote go and check chapter 5 of the fundamentals book for shape functions

#pragma mark __Mass___
void	triAxisSymElasticLumpedMassMatrix(double* M, double* coords, MaterialInfo minfo)
{
	double area = triArea(coords);		
	double rbar = (coords[0] + coords[1] + coords[2])/3.0;
	double mass = 2.0*M_PI*rbar*area*minfo.density/3.0; 

	int i;	
	int sizeM = 9;
	for ( i=0; i<sizeM; i++)
		M[i+i*sizeM] += mass;
		
	return;
}

void	triAxisSymElasticConsistentMassMatrix(double* M, double* coords, MaterialInfo minfo)
{
	double area = triArea(coords);
	double prefix = minfo.density*M_PI*area/10.0;

	double rbar = (coords[0] + coords[1] + coords[2])/3.0;
	
	double r413 = prefix*(2.0*rbar + 4.0*coords[0]/3.0);
	double r423 = prefix*(2.0*rbar + 4.0*coords[1]/3.0);
	double r433 = prefix*(2.0*rbar + 4.0*coords[2]/3.0);

	double r113 = prefix*(2.0*rbar - coords[0]/3.0);
	double r123 = prefix*(2.0*rbar - coords[1]/3.0);
	double r133 = prefix*(2.0*rbar - coords[2]/3.0); 	

	M[0] = r413;	M[9] = 0.0;		M[18] = 0.0;		M[27] = r133;	M[36] = 0.0;	M[45] = 0.0;		M[54] = r123;	M[63] = 0.0;	M[72] = 0.0;
	M[1] = 0.0;		M[10] = r413;	M[19] = 0.0;		M[28] = 0.0;	M[37] = r133;	M[46] = 0.0;		M[55] = 0.0;	M[64] = r123;	M[73] = 0.0;
	M[2] = 0.0;		M[11] = 0.0;	M[20] = r413;		M[29] = 0.0;	M[38] = 0.0;	M[47] = r133;		M[56] = 0.0;	M[65] = 0.0;	M[74] = r123;
	
	M[3] = r133;	M[12] = 0.0;	M[21] = 0.0;		M[30] = r423;	M[39] = 0.0;	M[48] = 0.0;		M[57] = r113;	M[66] = 0.0;	M[75] = 0.0;
	M[4] = 0.0;		M[13] = r133;	M[22] = 0.0;		M[31] = 0.0;	M[40] = r423;	M[49] = 0.0;		M[58] = 0.0;	M[67] = r113;	M[76] = 0.0;
	M[5] = 0.0;		M[14] = 0.0;	M[23] = r133;		M[32] = 0.0;	M[41] = 0.0;	M[50] = r423;		M[59] = 0.0;	M[68] = 0.0;	M[77] = r113;
	
	M[6] = r123;	M[15] = 0.0;	M[24] = 0.0;		M[33] = r113;	M[42] = 0.0;	M[51] = 0.0;		M[60] = r433;	M[69] = 0.0;	M[78] = 0.0;
	M[7] = 0.0;		M[16] = r123;	M[25] = 0.0;		M[34] = 0.0;	M[43] = r113;	M[52] = 0.0;		M[61] = 0.0;	M[70] = r433;	M[79] = 0.0;
	M[8] = 0.0;		M[17] = 0.0;	M[26] = r123;		M[35] = 0.0;	M[44] = 0.0;	M[53] = r113;		M[62] = 0.0;	M[71] = 0.0;	M[80] = r433;
		
	return;
}

#pragma mark ___Stiffness___
void	triAxisSymMapElasticStiffnessMatrix(	double* K, double* M,
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

	double* N   = (double*) calloc(3, sizeof(double));
	double* dNdSdN = (double*) calloc(3*2, sizeof(double));
		
	int sizeB[2];
	sizeB[0] = 6;	sizeB[1] = 9;

	double* B = (double*) calloc(sizeB[0]*sizeB[1],sizeof(double));
	
	double* dNdXdY = (double*) calloc(2*3, sizeof(double)); 
	
	double r;
	
	//double area = triArea(coords);
	//printf("twice area: %lf\n", 2*area);

	int i,j;
	for ( i=0; i<iinfo.numgauss; i++){
		for ( j=0; j<iinfo.numgauss; j++){

			triangularShapeFunctions(N, dNdSdN, gp[i], gp[j], iinfo.order);	
			r = N[0]* coords[0] + N[1]*coords[1] + N[2]*coords[2];	
			
			twoDJacobian(3, J, invJ, detJ, dNdSdN, coords);			
			twoDPhysicalDerivs(3, dNdXdY, invJ, dNdSdN);		

			if (n==0) 
				axisSymConstMatrix(3, B, dNdXdY, N, r, theta, n);	
			else
				antiAxisSymConstMatrix(3, B, dNdXdY, N, r, theta, n);
			
			if (iinfo.lumpedmass)
				triAxisSymElasticLumpedMassMatrix(M, coords, minfo);
			else
				triAxisSymElasticConsistentMassMatrix(M, coords, minfo);

			accumulateStiffness(K, sizeK, B, sizeB, D, sizeD, w[i], w[j], 1.0, detJ*r*M_PI);
		}
	}

	free(gp), free(w);
	free(D), free(B);
	free(J), free(invJ);
	free(N), free(dNdSdN), free(dNdXdY);

}

#pragma mark ___Assemble___
void assembleTriAxisSymMapElasticMatricies(	double* K, double* M,		
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
		
		triAxisSymMapElasticStiffnessMatrix(Ke, Me, coords, info, theta, n);							
		
					
#ifdef DEBUG_PRINT
		int i,j;
		printf("\n\n\n\n\n");
		for ( i=0; i<edofs; i++)
			for ( j=0; j<edofs; j++)
				printf("Ke[%ld][%ld]: %lf\n", i,j, Ke[i*edofs + j]);					
																					
		for ( i=0; i<edofs; i++)
			for ( j=0; j<edofs; j++)
				printf("Me[%ld][%d]: %lf\n", i, j, Me[i + j*edofs]);						
#endif	
		assemble(K, M, el, Ke, Me, ginfo, inds);
	}

	free(Ke), free(Me);
	free(coords);
	
}