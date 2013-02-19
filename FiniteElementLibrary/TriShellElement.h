/*
 *  TriShellElement.h
 *  ModelViewer
 *
 *  Created by Cynthia Bruyns on 3/19/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _TriShellElement_H_
#define _TriShellElement_H_
 
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <string.h>

#include "SharedFEData.h"
#include "GeneralFE.h"
#include "ShapeFunctions.h"

#include "2dElements.h"

#include "MatrixOps.h"

#include "CLPKRealMatrixOpsCM.h"

#pragma mark ___ConstMatrix___
#pragma mark ___ConstMatrix___
inline void shellConstBending(double* Bbend, double* dNdXdY)
{
	int cr = 3;
	int i,i1;
	for ( i=0; i<3; i++){
		  i1 = i*18;  
		 
		 Bbend[i1+cr*4  ] =  dNdXdY[i  ];		
		 Bbend[i1+cr*3+1] = -dNdXdY[i+3];
		 Bbend[i1+cr*4+2] =  dNdXdY[i+3];
		 Bbend[i1+cr*3+2] = -dNdXdY[i  ];
		
	}
#ifdef DEBUG_PRINT
	printf("\n\n\n\n");
	for ( i=0; i<18; i++){
		for ( j=0; j<3; j++){
			printf("Bbend[%ld][%ld]: %lf\n",j,i, Bbend[i*3+j]);
		}
		printf("\n");
	}
#endif	
}
inline void shellConstMembrane(double* Bmemb, double* dNdXdY)
{
	int cr = 3;
	int i,i1;
	for ( i=0; i<3; i++){
		  i1 = i*18;  
		
		 Bmemb[i1     ] = dNdXdY[i  ];		
		 Bmemb[i1+cr+1] = dNdXdY[i+3];
		 Bmemb[i1+2   ] = dNdXdY[i+3];
		 Bmemb[i1+cr+2] = dNdXdY[i  ];
	}

#ifdef DEBUG_PRINT
	int j;
	for ( i=0; i<18; i++){
		for ( j=0; j<3; j++){
			printf("Bmemb[%ld][%ld]: %lf\n",j,i,Bmemb[i*3+j]);
		}
		printf("\n");
	}
#endif
}

inline void shellConstShear(double* Bshear, double* N, double* dNdXdY)
{
	int cr=2;
	int i, i1;
	 for ( i=0; i<3; i++){
		 i1 = i*12;  
		 		 		
		 Bshear[i1+cr*2] = dNdXdY[i];
		 Bshear[i1+cr*4] = N[i];
		 
		 Bshear[i1+cr*2+1] = dNdXdY[i+3];
		 Bshear[i1+cr*3+1] = -1.0*N[i];
	}
#ifdef DEBUG_PRINT
	int j;
	for ( i=0; i<18; i++)
		for ( j=0; j<2; j++)
			printf("Bshear[%ld][%ld]: %lf\n", j,i, Bshear[i*2+j]);
#endif
}



#pragma mark ___Transform___

inline void computeDirectionCosines(double* transformMatrix, double* prime, double* coords, GeomInfo ginfo)
{
	//  compute direction cosines
	double v12[3];
	v12[0] = coords[1]-coords[0];  //2 - 1
	v12[1] = coords[4]-coords[3];
	v12[2] = coords[7]-coords[6];

	double l12 = v3_length(v12);
	
	double v23[3];
	v23[0] = coords[2]-coords[1]; // 3 - 2
	v23[1] = coords[5]-coords[4];
	v23[2] = coords[8]-coords[7];

	double l23 = v3_length(v23);

	double v13[3];
	v13[0] = coords[0]-coords[2]; // 1 - 3
	v13[1] = coords[3]-coords[5];
	v13[2] = coords[6]-coords[8];

	double l13 = v3_length(v13);
	
	double z[3];
	v3_cross(v12,v13,z);
	double lz = v3_length(z);
	
	double y[3];
	v3_cross(z,v12,y);
	double ly = v3_length(y);

	double vxx = v12[0]/l12;
	double vxy = v12[1]/l12;
	double vxz = v12[2]/l12;

	double vyx = y[0]/ly;
	double vyy = y[1]/ly;
	double vyz = y[2]/ly;
	
	double vzx = z[0]/lz;
	double vzy = z[1]/lz;
	double vzz = z[2]/lz;

	//   transformation matrix	
	int edof = ginfo.numelnodes*ginfo.nodedofs;
	int i1,i;										
	for ( i=0; i< ginfo.numelnodes*(ginfo.dimension-1); i++){ 
	
		i1 = i*ginfo.dimension*edof + i*ginfo.dimension;
		
		transformMatrix[i1       ] = vxx; 
		transformMatrix[i1+edof  ] = vxy; 
		transformMatrix[i1+2*edof] = vxz; 
		
		transformMatrix[i1+1	   ] = vyx; 
		transformMatrix[i1+edof+1  ] = vyy; 
		transformMatrix[i1+2*edof+1] = vyz; 
		
		transformMatrix[i1+2       ] = vzx; 
		transformMatrix[i1+edof+2  ] = vzy; 
		transformMatrix[i1+2*edof+2] = vzz; 
	}
	
#ifdef DEBUG_PRINT	
	for( i=0; i<edof; i++)
		for ( j=0; j<edof; j++)
			printf("transformMatrix[%ld][%ld]: %lf\n", i,j, transformMatrix[j*edof+i]);
#endif	
	 
	// compute nodal values in terms of local axes
	double alpha = acos(v3_dot(v13, v12)/(v3_length(v13)*v3_length(v12)));

	prime[0] = 0;					prime[3] = 0;
	prime[1] = l12;					prime[4] = 0;
	prime[2] = l13*cos(alpha);		prime[5] = l13*sin(alpha);
	
#ifdef DEBUG_PRINT
	for ( i=0; i<6; i++)
		printf("prime[%ld]: %lf\n", i, prime[i]);
#endif	
}

inline void  transformToGlobal(double* Kout, double* K, double* T, int size) 
{
	double* temp = (double*) calloc(size*size, sizeof(double));		
	
/*#ifdef DEBUG_PRINT
	int i,j;
	printf("\n\ntransform to global\n");
	for ( i=0; i<size; i++)
		for ( j=0; j<size; j++)
			printf("Kin[%ld][%ld]: %lf\n", i,j,K[j*size+i]);		
#endif
*/
/*
	http://www.netlib.org/blas/sgemm.f
*/															//				m k k k	k m
	cblas_dgemm(CblasColMajor,CblasTrans, CblasNoTrans,		//temp = T^T K  9x9 9x9 9x9
				size, 
				size,										
				size,									 
				 1.0,
				 T,											
				 size,										 
				 K, 
				 size,										
                 1.0, 
				 temp, 
				 size);										
				 
															//				m n k m 
	cblas_dgemm(CblasColMajor,CblasNoTrans, CblasNoTrans,	// K = temp*T		9x9 9x9 9x9
				 size, 
				 size, 
				 size,
                 1.0,  
				 temp, 
				 size,
				 T, 
				 size,
                 1.0, 
				 Kout, 
				 size);


#ifdef DEBUG_PRINT
	int i,j;
	printf("\n\ntransform to global\n");
	for ( i=0; i<size; i++){
		for ( j=0; j<size; j++){
			printf("Kout[%ld][%ld]: %lf\n", i,j,Kout[j*size+i]);	
		}
	}
#endif	

	free(temp);
}

#pragma mark ___OneElement___

inline void	shellStiffnessMatrix(	double* K,
									 double* coords, 
									 AnalysisInfo info, 
									 ShellInfo sinfo)
{
	int i,j;
	
	MaterialInfo minfo = info.materialinfo;
	GeomInfo ginfo = info.geominfo;
	IntegrationInfo iinfo = info.intinfo;
	
	int edof = ginfo.numelnodes*ginfo.nodedofs;
	//printf("edof: %ld\n", edof);

	double* transformMatrix = (double*) calloc(edof*edof,sizeof(double));
	double  prime[6];
	computeDirectionCosines(transformMatrix, prime, coords, ginfo);
		
	double* J = (double*) calloc(2*2, sizeof(double));
	double* invJ = (double*) calloc(2*2, sizeof(double));
	double  detJ = 0.;
	
	// planar element

	double* N      = (double*) calloc(ginfo.numelnodes  , sizeof(double));
	double* dNdSdN = (double*) calloc(ginfo.numelnodes*2, sizeof(double));
	double* dNdXdY = (double*) calloc(ginfo.numelnodes*2, sizeof(double));
	
	//bending	
	int sizeD = 3;	
	double* Dmem = (double*) calloc(sizeD*sizeD,sizeof(double));
	materialMatrix(Dmem, minfo);
	for ( i=0; i<sizeD*sizeD; i++){
		Dmem[i] *= sinfo.thickness;
	}
#ifdef DEBUG_PRINT
	for ( i=0; i<3; i++)
		for ( j=0; j<3; j++)
			printf("Dmem[%ld][%ld]: %lf\n", i, j,Dmem[j*3 + i]);
#endif
	
	double* Dbend = (double*) calloc(sizeD*sizeD,sizeof(double));
	materialMatrix(Dbend, minfo);
	for ( i=0; i<sizeD*sizeD; i++){
		Dbend[i] *= pow(sinfo.thickness,3)/12.0;
	}
#ifdef DEBUG_PRINT
	for ( i=0; i<3; i++)
		for ( j=0; j<3; j++)
			printf("Dbend[%ld][%ld]: %lf\n", i, j,Dbend[j*3 + i]);
#endif			
	
		
	double* Bbend = (double*) calloc(3*edof,sizeof(double));
	double* Bmem = (double*) calloc(3*edof,sizeof(double));
	int sizeB[2];
	sizeB[0] = ginfo.dimension;	
	sizeB[1] = edof;
	
	double* Kb = (double*) calloc(edof*edof, sizeof(double));
	double* Km =(double*) calloc(edof*edof, sizeof(double));
	int sizeK = edof;

	double* gpb  = (double*) calloc(iinfo.numgauss, sizeof(double));
	double* wb = (double*) calloc(iinfo.numgauss, sizeof(double));
	getGaussPoints(gpb, wb, iinfo.numgauss);

	for ( i=0; i<iinfo.numgauss; i++){
		for ( j=0; j<iinfo.numgauss; j++){
		
			triangularShapeFunctions(N, dNdSdN, gpb[i], gpb[j], iinfo.order);				
			twoDJacobian(3, J, invJ, detJ, dNdSdN, prime);
			twoDPhysicalDerivs(3, dNdXdY, invJ, dNdSdN);
			
			shellConstBending(Bbend, dNdXdY); 
			shellConstMembrane(Bmem, dNdXdY);
			
#ifdef DEBUG_PRINT	
			printf("Kb\n");
#endif
			accumulateStiffness(Kb, sizeK, Bbend, sizeB, Dbend, sizeD, wb[i], wb[j], 1.0, detJ);	
			
#ifdef DEBUG_PRINT	
			printf("Km\n");
#endif			
			accumulateStiffness(Km, sizeK, Bmem , sizeB, Dmem , sizeD, wb[i], wb[j], 1.0, detJ);
			
		}		
	}
	
#ifdef DEBUG_PRINT	
	printf("\n\n\n");
	for ( i=0; i<sizeK; i++)
		for ( j=0; j<sizeK; j++)
			printf("Kb[%ld][%ld]: %lf\n", i,j,Kb[j*sizeK+i]);

	printf("\n\n\n");
	for ( i=0; i<sizeK; i++)
		for ( j=0; j<sizeK; j++)
			printf("Km[%ld][%ld]: %lf\n", i,j,Km[j*sizeK+i]);
#endif
	
	//shear	
	double shearModulus = 0.5*minfo.elasticModulus/(1.0+minfo.poisson);	// shear modulus
	double shearCorrection = 5.0/6.0;									// shear correction factor					
	double* Dshear = (double*) calloc(2*2,sizeof(double));					// shear material property
	Dshear[0] = shearModulus*shearCorrection*sinfo.thickness;
	Dshear[2] = 0.;
	Dshear[1] = 0.;
	Dshear[3] = shearModulus*shearCorrection*sinfo.thickness;
	sizeD= 2;
#ifdef DEBUG_PRINT
	for ( i=0; i<2; i++)
		for( j=0; j<2; j++)
			printf("Dshear[%ld][%ld]: %lf\n", i,j, Dshear[j*2+i]);
#endif	
	
	double* Bshear = (double*) calloc(2*edof,sizeof(double));
	sizeB[0] = 2;	sizeB[1] = edof;

	double* gps  = (double*) calloc(sinfo.numgauss_s, sizeof(double));
	double* ws = (double*) calloc(sinfo.numgauss_s, sizeof(double));
	getGaussPoints(gps, ws, sinfo.numgauss_s);
	
	double* Ks =(double*) calloc(edof*edof, sizeof(double));

	for (i=0; i<sinfo.numgauss_s; i++){
		for ( j=0; j<sinfo.numgauss_s; j++){
		
			triangularShapeFunctions(N, dNdSdN, gps[i], gps[j], iinfo.order);	

			twoDJacobian(3, J, invJ, detJ, dNdSdN, prime);			
			twoDPhysicalDerivs(3, dNdXdY, invJ, dNdSdN);

			shellConstShear(Bshear, N, dNdXdY);
			
#ifdef DEBUG_PRINT	
			printf("Ks\n");
#endif
			accumulateStiffness(Ks, sizeK, Bshear, sizeB, Dshear, sizeD, ws[i], ws[j], 1.0, detJ);		
		}		
	}
#ifdef DEBUG_PRINT
	printf("\n\n\n");
	for ( i=0; i<edof; i++)
		for ( j=0; j<edof; j++)
			printf("Ks[%ld][%ld]: %lf\n", i,j,Ks[j*sizeK+i]); 
#endif
		
	double* Kt = (double*) calloc(edof*edof, sizeof(double));
	realDoubleAdd3Matricies(Kt, Kb, Km, Ks, edof);			
	
#ifdef DEBUG_PRINT
	printf("\n\n\n");
	for ( i=0; i<edof; i++)
		for ( j=0; j<edof; j++)
			printf("Kt[%ld][%ld]: %lf\n", i,j,Kt[j*sizeK+i]);
#endif

	transformToGlobal(K, Kt, transformMatrix, edof);
	
#ifdef DEBUG_PRINT
	printf("\n\n\n");
	for ( i=0; i<edof; i++)
		for ( j=0; j<edof; j++)
			printf("Kfin[%ld][%ld]: %lf\n", i,j,K[j*sizeK+i]);
#endif

	free(gpb), free(wb), free(gps), free(ws);
	
	free(transformMatrix);
	free(J), free(invJ);
	free(N), free(dNdSdN), free(dNdXdY);
	
	free(Dbend), free(Dmem), free(Dshear);
	free(Bbend), free(Bmem), free(Bshear);
	free(Kb), free(Km), free(Ks), free(Kt);
}

inline void	shellLumpedMassMatrix(double* M, double* coords, AnalysisInfo info, ShellInfo sinfo)
{
	GeomInfo ginfo = info.geominfo;
	MaterialInfo minfo = info.materialinfo;
	
	double area = quadArea(coords);
	double volume = sinfo.thickness*area;
	//printf("quad area: %lf\n", area);
	double mass = volume*minfo.density/3.0;
		
	int edof = 	ginfo.numelnodes*ginfo.nodedofs;
	int i;
	for ( i=0; i<edof; i++)
		M[i+i*edof] += mass;
}

#pragma mark ___ManyElement___

inline void assembleShellMatricies(	double* K, double* M, 
								double* verts, int* inds, 
								AnalysisInfo info, ShellInfo sinfo)
{
	MaterialInfo	minfo = info.materialinfo;
	GeomInfo		ginfo = info.geominfo;
	
	int edofs = ginfo.numelnodes*ginfo.nodedofs;
	double* Ke = (double*) calloc(edofs*edofs, sizeof(double));
	double* Me = (double*) calloc(edofs*edofs, sizeof(double));	
	
	double* coords = (double*) calloc(ginfo.numelnodes*ginfo.dimension,sizeof(double));
	

	int el;
	for ( el=0; el<ginfo.numelements; el++){ 
	
		memset(Ke,0,edofs*edofs*sizeof(double));
		memset(Me,0,edofs*edofs*sizeof(double));
			
		get3dCoords(coords,el,verts,inds, ginfo); 
		
#ifdef DEBUG_PRINT
		printf("\n\n\nelement: %ld\n", el);
#endif	
		shellStiffnessMatrix(  Ke, coords, info, sinfo);							
		shellLumpedMassMatrix( Me, coords, info, sinfo);
					
#ifdef DEBUG_PRINT
		printf("\n\n\n\n\n");
		for ( i=0; i<edofs; i++)
			for ( j=0; j<edofs; j++)
				printf("Ke[%ld][%ld]: %lf\n", i,j, Ke[i*edofs + j]);					
																					
		for ( i=0; i<edofs; i++)
			printf("Me[%ld]: %lf\n", i, Me[i + i*edofs]);	
					
#endif	
		assemble(K, M, el, Ke, Me, ginfo, inds);
	}

	//int dofs = ginfo.numnodes*ginfo.nodedofs;

	free(Ke), free(Me);
	free(coords);
}

inline void	checkDrilling(double* K, int size)
{	
	int i,j;
	for ( i=0; i<size; i++){
	
		if(dabs(K[i*size + i]) < 1e-5){

			double sum=0.0;
			for ( j=0; j< size; j++){
				sum += dabs(K[j*size + i]);
			}
			
#ifdef DEBUG_PRINT
			printf("sum[%d]: %lf\n", i, sum);
#endif
			if (sum < 1e-5){
				K[i*size + i] = 1.0;
			}

		}
	}	
}

#endif
