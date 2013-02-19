/*
 *  TriPlateElement.h
 *  ModelChangingApp
 *
 *  Created by Cynthia Maxwell on 8/18/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "GeneralFE.h"
#include "2dElements.h"


void	triPlateLumpedMassMatrix(double* M, double* coords, MaterialInfo minfo)
{
	double area = triArea(coords);
	double mass = area*minfo.density/9.0;
	int i;
	for ( i=0; i<9; i++)
		M[i+i*9] += mass; 
}



void getPlateGaussPoints(double* gp, double* w, int num)
{
	gp[0] = 1./3.; gp[4] = 1./3.;
	gp[1] = 1./5.; gp[5] = 1./5.;
	gp[2] = 3./5.; gp[6] = 1./5.;
	gp[3] = 1./5.; gp[7] = 3./5.;
	
	 
	w[0] = -9./32.;
	w[1] = 25./96.;
	w[2] = 25./96.;
	w[3] = 25./96.;
}

void plateMaterialMatrix(double* D,  PlateInfo constants)
{
	double frac = pow(constants.thickness,3.0)/12.0;
	for (int i=0; i< 9; i++){
		D[i] *= frac;
	}
//	for (int i=0; i<3; i++)
//		for (int j=0; j<3; j++)
//			printf("D[%d][%d]: %lf\n", i, j, D[i+j*3]);
}


void	triPlateStiffnessMatrix(		double*			K,
										double*			coords,
										AnalysisInfo	info,
										PlateInfo		pinfo)
{	
	MaterialInfo	minfo = info.materialinfo;
	GeomInfo		ginfo = info.geominfo;
	IntegrationInfo iinfo = info.intinfo;
	

	
	int sizeK = ginfo.numelnodes*ginfo.nodedofs;

	iinfo.numgauss = 4;
	
	double* gp  = (double*) calloc(iinfo.numgauss*2, sizeof(double));
	double* w = (double*) calloc(iinfo.numgauss, sizeof(double));
	getPlateGaussPoints( gp, w, iinfo.numgauss);
	
	int sizeD = 3;
	double* D = (double*) calloc(sizeD*sizeD,sizeof(double));
	materialMatrix(D, minfo);
	plateMaterialMatrix(D, pinfo);	
	
	double detJ = 0.;
		
	int sizeB[2];
	sizeB[0] = 3;	sizeB[1] = 9;
	double* B = (double*) calloc(sizeB[0]*sizeB[1],sizeof(double));
	
	double x1 = coords[0]; 
	double x2 = coords[1];
	double x3 = coords[2];
	
	double y1 = coords[3]; 
	double y2 = coords[4];
	double y3 = coords[5];
	
	double f1 = x2*y3 - x3*y2; 
	double f2 = x3*y1 - x1*y3; 
	double f3 = x1*y2 - x2*y1;
	
	double b1 = y2 - y3; 
	double b2 = y3 - y1; 
	double b3 = y1 - y2;

	double c1 = x3 - x2; 
	double c2 = x1 - x3; 
	double c3 = x2 - x1; 
	
	double A = (f1 + f2 + f3)/2.0;
	
	double x,y;
	double A1, A2, A3;

	int i;
	for ( i=0; i<iinfo.numgauss; i++){		
		
		x = (1. - gp[i] - gp[i+4])*x1 + gp[i]*x2 + gp[i+4]*x3;
		y = (1. - gp[i] - gp[i+4])*y1 + gp[i]*y2 + gp[i+4]*y3;
		
		
		A1 = 1./(2.*A)* (f1 + x*b1 + y*c1);
		A2 = 1./(2.*A)* (f2 + x*b2 + y*c2);
		A3 = 1./(2.*A)* (f3 + x*b3 + y*c3);
		
		B[0] = -(-(A*b1*(3*b1 + 2*b2)) + 3*pow(b1,3)*x + pow(b2,2)*(f1 + c1*y) + pow(b1,2)*(3*f1 + f2 + 3*b2*x + 3*c1*y + c2*y) + b1*b2*(2*f1 + 2*f2 + 3*b2*x + 2*c1*y + 2*c2*y))/(2*pow(A,3));
		B[3] = -(b3*(-2*A*b1*b2 + b2*b3*(f1 + c1*y) + pow(b1,2)*(f2 + 3*b2*x + c2*y) + b1*b2*(2*f1 + 2*f3 + 3*b3*x + 2*c1*y + 2*c3*y)))/(4*pow(A,3));
		B[6] = -(-2*A*b1*b3*c2 + pow(b3,2)*c2*(f1 + 3*b1*x + c1*y) + b1*c3*(2*b2*f1 + b1*f2 + 3*b1*b2*x + 2*b2*c1*y + b1*c2*y) + 2*b1*b3*c2*(f3 + c3*y))/(4*pow(A,3));
		B[9] = -(-(A*b2*(2*b1 + 3*b2)) + pow(b1,2)*(f2 + 3*b2*x + c2*y) + b1*b2*(2*f1 + 2*f2 + 3*b2*x + 2*c1*y + 2*c2*y) + pow(b2,2)*(f1 + 3*f2 + 3*b2*x + c1*y + 3*c2*y))/(2*pow(A,3));
		B[12] = -(b1*(-2*A*b2*b3 + b1*b3*(f2 + c2*y) + b2*b3*(2*f1 + 2*f2 + 3*b1*x + 2*c1*y + 2*c2*y) + pow(b2,2)*(f3 + 3*b3*x + c3*y)))/(4*pow(A,3));
		B[15] = -(-2*A*b1*b2*c3 + 2*b1*b2*c3*(f1 + c1*y) + pow(b1,2)*c3*(f2 + 3*b2*x + c2*y) + b2*c1*(2*b3*f2 + b2*f3 + 3*b2*b3*x + 2*b3*c2*y + b2*c3*y))/(4*pow(A,3));
		B[18] = -(-(A*b3*(2*b1 + 3*b3)) + pow(b1,2)*(f3 + 3*b3*x + c3*y) + b1*b3*(2*f1 + 2*f3 + 3*b3*x + 2*c1*y + 2*c3*y) + pow(b3,2)*(f1 + 3*f3 + 3*b3*x + c1*y + 3*c3*y))/(2*pow(A,3));
		B[21] = -(b2*(-2*A*b1*b3 + pow(b3,2)*(f1 + 3*b1*x + c1*y) + b1*b2*(f3 + c3*y) + b1*b3*(2*f2 + 2*f3 + 3*b2*x + 2*c2*y + 2*c3*y)))/(4*pow(A,3));
		B[24] = -(-2*A*b2*b3*c1 + pow(b3,2)*c2*(f1 + 3*b1*x + c1*y) + pow(b2,2)*c1*(f3 + c3*y) + b3*(3*pow(b2,2)*c1*x + 2*b2*c1*(f2 + c2*y) + 2*b1*c2*(f3 + c3*y)))/(4*pow(A,3));
 
 
		B[1] = -(-(A*c1*(3*c1 + 2*c2)) + pow(c2,2)*(f1 + b1*x) + 3*pow(c1,3)*y + pow(c1,2)*(3*f1 + f2 + 3*b1*x + b2*x + 3*c2*y) + c1*c2*(2*f1 + 2*f2 + 2*b1*x + 2*b2*x + 3*c2*y))/(2*pow(A,3));
		B[4] = -(-2*A*b2*c1*c3 + b2*c3*(2*c1*f3 + c3*(f1 + b1*x + 3*c1*y)) + b3*c1*(c1*f2 + b2*c1*x + 2*b2*c3*x + c2*(2*f1 + 2*b1*x + 3*c1*y)))/(4*pow(A,3));
		B[7] = -(c3*(-2*A*c1*c2 + c2*c3*(f1 + b1*x) + pow(c1,2)*(f2 + b2*x + 3*c2*y) + c1*c2*(2*f1 + 2*f3 + 2*b1*x + 2*b3*x + 3*c3*y)))/(4*pow(A,3));
		B[10] = -(-(A*c2*(2*c1 + 3*c2)) + pow(c1,2)*(f2 + b2*x + 3*c2*y) + c1*c2*(2*f1 + 2*f2 + 2*b1*x + 2*b2*x + 3*c2*y) + pow(c2,2)*(f1 + 3*f2 + b1*x + 3*b2*x + 3*c2*y))/(2*pow(A,3));
		B[13] = -(-2*A*b3*c1*c2 + b3*(b1*pow(c2,2)*x + 2*c1*c2*(f1 + b1*x) + pow(c1,2)*(f2 + b2*x + 3*c2*y)) + b1*c2*(c2*f3 + c3*(2*f2 + 2*b2*x + 3*c2*y)))/(4*pow(A,3));
		B[16] = -(c1*(-2*A*c2*c3 + c1*c3*(f2 + b2*x) + c2*c3*(2*f1 + 2*f2 + 2*b1*x + 2*b2*x + 3*c1*y) + pow(c2,2)*(f3 + b3*x + 3*c3*y)))/(4*pow(A,3));
		B[19] = -(-(A*c3*(2*c1 + 3*c3)) + pow(c1,2)*(f3 + b3*x + 3*c3*y) + c1*c3*(2*f1 + 2*f3 + 2*b1*x + 2*b3*x + 3*c3*y) + pow(c3,2)*(f1 + 3*f3 + b1*x + 3*b3*x + 3*c3*y))/(2*pow(A,3));
		B[22] = -(-2*A*b1*c2*c3 + b1*c2*(2*c3*f2 + c2*f3 + b3*c2*x + 3*c2*c3*y) + b2*c3*(2*(c1*f3 + b3*c1*x + b1*c2*x) + c3*(f1 + b1*x + 3*c1*y)))/(4*pow(A,3));
		B[25] = -(c2*(-2*A*c1*c3 + c1*c2*(f3 + b3*x) + pow(c3,2)*(f1 + b1*x + 3*c1*y) + c1*c3*(2*f2 + 2*f3 + 2*b2*x + 2*b3*x + 3*c2*y)))/(4*pow(A,3));
		
		B[2] = -((-(A*(b2*c1 + b1*(3*c1 + c2))) + pow(b1,2)*(3*c1 + c2)*x + b2*(c2*f1 + pow(c1,2)*y + c1*(f1 + f2 + b2*x + 2*c2*y)) + b1*(3*pow(c1,2)*y + c2*(f1 + f2 + 2*b2*x + c2*y) + c1*(3*f1 + f2 + 2*b2*x + 2*c2*y)))/pow(A,3));
		B[5] = -(-(A*b2*(b3*c1 + b1*c3)) + b1*b3*(c1*f2 + c2*(f1 + b1*x + 2*c1*y)) + b2*(pow(b3,2)*c1*x + b1*c3*(f3 + c3*y) + b3*(c3*(f1 + 2*b1*x) + pow(c1,2)*y + c1*(f1 + f3 + 2*b1*x + 2*c3*y))))/(2*pow(A,3));
		B[8] = -(-(A*c2*(b3*c1 + b1*c3)) + b1*c2*c3*f1 + b3*c2*c3*f1 + b1*c1*c3*f2 + b3*c1*c2*f3 + b1*c2*c3*f3 + pow(b3,2)*c1*c2*x + pow(b1,2)*c2*c3*x + 2*b1*b3*c2*c3*x + 2*b1*c1*c2*c3*y + 2*b3*c1*c2*c3*y + b1*c2*pow(c3,2)*y + b2*c1*c3*(f1 + 2*b1*x + c1*y))/(2*pow(A,3));
		B[11] = -((-(A*(b1*c2 + b2*(c1 + 3*c2))) + pow(b2,2)*(c1 + 3*c2)*x + b1*(c1*f2 + pow(c2,2)*y + c2*(f1 + f2 + b1*x + 2*c1*y)) + b2*(pow(c1,2)*y + c1*(f1 + f2 + 2*b1*x + 2*c2*y) + c2*(f1 + 3*f2 + 2*b1*x + 3*c2*y)))/pow(A,3));
		B[14] = -(-(A*b3*(b2*c1 + b1*c2)) + b1*pow(b2,2)*c3*x + b1*b3*(c1*f2 + pow(c2,2)*y + c2*(f1 + f2 + b1*x + 2*c1*y)) + b2*(b3*(c1*f1 + 2*b1*c1*x + 2*b1*c2*x + pow(c1,2)*y) + b1*(c3*f2 + c2*f3 + 2*c2*c3*y)))/(2*pow(A,3));
		B[17] = -(-(A*(b2*c1 + b1*c2)*c3) + b1*c2*c3*f1 + b3*c1*c2*f2 + b1*c1*c3*f2 + pow(b2,2)*c1*c3*x + pow(b1,2)*c2*c3*x + b3*c1*pow(c2,2)*y + 2*b1*c1*c2*c3*y + b2*c1*(c2*(f3 + 2*b3*x) + c3*(f1 + f2 + 2*b1*x + c1*y + 2*c2*y)))/(2*pow(A,3));
		B[20] = -((-(A*(b1*c3 + b3*(c1 + 3*c3))) + pow(b3,2)*(c1 + 3*c3)*x + b1*(c1*f3 + pow(c3,2)*y + c3*(f1 + f3 + b1*x + 2*c1*y)) + b3*(pow(c1,2)*y + c1*(f1 + f3 + 2*b1*x + 2*c3*y) + c3*(f1 + 3*f3 + 2*b1*x + 3*c3*y)))/pow(A,3));
		B[23] = -(-(A*b1*(b3*c2 + b2*c3)) + b1*pow(b2,2)*c3*x + b1*b3*c2*(f2 + c2*y) + b2*(pow(b3,2)*c1*x + b3*(c1*f3 + 2*b1*c2*x + c3*(f1 + 2*b1*x + 2*c1*y)) + b1*(c2*f3 + pow(c3,2)*y + c3*(f2 + f3 + 2*c2*y))))/(2*pow(A,3));
		B[26] = -(-(A*c1*(b3*c2 + b2*c3)) + b2*c1*c3*f2 + b2*c1*c2*f3 + b1*c2*c3*f3 + pow(b3,2)*c1*c2*x + pow(b2,2)*c1*c3*x + 2*b2*c1*c2*c3*y + b1*c2*pow(c3,2)*y + b3*c2*(c3*(f1 + 2*b1*x + 2*c1*y) + c1*(f2 + f3 + 2*b2*x + c2*y)))/(2*pow(A,3));
 
	    detJ = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
					
		accumulate(K, sizeK, B, sizeB, D, sizeD, w[i]*detJ);
				
	}
	
	free(gp), free(w);
	free(D);
	free(B);
}



void assembleTriPlateMatricies(		double* K, double* M, 
									double* verts, int* inds, 
									AnalysisInfo info,
									PlateInfo pinfo)
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
		triPlateStiffnessMatrix(  Ke, coords, info, pinfo);							
		triPlateLumpedMassMatrix( Me, coords, minfo);
			
#ifdef DEBUG_PRINT					
		int i,j;
		printf("\n\n\n\n\nKe=[");
		for ( i=0; i<edofs; i++){
			for ( j=0; j<edofs; j++){
				printf(" %lf\t",  Ke[i*edofs + j]);		
			}
			printf("\n");
		}
		printf("];");
																					

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