/*
 *  LineAcousticElasticCouplingElement.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 1/8/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef _LineAcousticElasticCouplingElement_Header_
#define _LineAcousticElasticCouplingElement_Header_


#include <stdio.h>
#include <stdlib.h> 


#include "ShapeFunctions.h"
#include "GeneralFE.h"


double getLineLength(double* coords)
{
	double dx = -coords[0] + coords[1];
	double dy = -coords[2] + coords[3];
#ifdef DEBUG_PRINT
	printf("length: %lf\n", abs(dx)+abs(dy));
#endif	
	return (abs(dx)+abs(dy));

}

void	couplingMatrix(double* B, double* N, double* normals)
{
	// B = [N1] *[ N1 {n1}		 N2{n1} ]
	//	   [N2]	 [    {n2}		   {n2} ]
	//	   
	
	double n[2];
	n[0] = normals[0];
	n[1] = normals[1];

    B[0] = N[0]*n[0]*N[0]; B[4] = N[0]*n[0]*N[1];
    B[1] = N[0]*n[1]*N[0]; B[5] = N[0]*n[1]*N[1];
    B[2] = N[1]*n[0]*N[0]; B[6] = N[1]*n[0]*N[1];
    B[3] = N[1]*n[1]*N[0]; B[7] = N[1]*n[1]*N[1];

#ifdef DEBUG_PRINT
	printf("make H\n");
	for ( int ii=0; ii<4; ii++)
			for ( int jj=0; jj<2; jj++)
				printf("B[%ld][%ld]: %lf\n", ii, jj, B[jj*4 + ii]);					
#endif
}


void 		accumulateCouplingMatrix(double* H, double* B, double w, double detJ)
{
	// H += H*w*detJ;
	
	for (int i=0; i<4; i++)
		for (int j=0; j<2; j++)
			H[i+j*4] += B[i+j*4]*w*detJ;	

#ifdef DEBUG_PRINT
	printf("accum\n");
	for ( int ii=0; ii<4; ii++)
			for ( int jj=0; jj<2; jj++)
				printf("H[%ld][%ld]: %lf\n", ii, jj, H[jj*4 + ii]);					
#endif
}

void	lineAcousticElasticCouplingMatrix(		double* H,
												double* coords,
												double* normals,
												AnalysisInfo info)
{	

	GeomInfo		ginfo = info.geominfo;
	IntegrationInfo iinfo = info.intinfo;	

	double* gp  = (double*) calloc(iinfo.numgauss, sizeof(double));
	double* w = (double*) calloc(iinfo.numgauss, sizeof(double));
	getGaussPoints( gp, w, iinfo.numgauss);
	
	double* N   = (double*) calloc(2, sizeof(double));	
	double* dNdSdN = (double*) calloc(2, sizeof(double));
	
	double* B = (double*) calloc(4*2,sizeof(double));

	double detJ;
	
	int i;
	for ( i=0; i<iinfo.numgauss; i++){
			
		lineShapeFunctions(N, dNdSdN, gp[i], iinfo.order);
		detJ = 0.5*getLineLength(coords);
#ifdef DEBUG_PRINT
		printf("detJ: %lf\n", detJ);
#endif		
		couplingMatrix(B, N, normals);
		
		accumulateCouplingMatrix(H, B, w[i], detJ);
			
	}
	
	free(gp), free(w);
	free(N), free(dNdSdN);
	free(B);
}


inline void assembleH(double* H, int id, double* He, GeomInfo ginfo, CoupleInfo cinfo, int* inds, int* map)
{	
	int i,j,ii,jj;

	int bdofs = cinfo.numbeamverts*2;
//	int fdofs = cinfo.numfluidverts;
	
#ifdef DEBUD_PRINT
	printf("\n\n\nbdofs: %d fdofs: %d\n", bdofs, fdofs);
	
	for (int i=0; i<4; i++)
		for (int j=0; j<2; j++)
			printf("He[%ld][%ld]: %lf\n", i,j, He[i + j*4]);	
#endif	
	// how to tell this, a lookuptable of
	//[ couple_node_index beam_node_index fluid__node_index]

	for ( i=0; i<2; i++){
	
		int iindex = inds[id + i*ginfo.numelements];
#ifdef DEBUG_PRINT
		printf("iindex: %d\n", iindex);
#endif		
		ii = map[iindex+ 1*ginfo.numnodes]*2;
#ifdef DEBUG_PRINT
		printf("ii: %d\n", ii);
#endif		
		for ( j=0; j<2; j++){
		
			int jindex = inds[id + j*ginfo.numelements];
#ifdef DEBUG_PRINT
			printf("jindex: %d\n", jindex);
#endif		
			jj = map[jindex + 2*ginfo.numnodes];
#ifdef DEBUG_PRINT
			printf("jj: %d\n", jj);
#endif			
			H[jj*bdofs + ii]   += He[j*4 + i*2    ];
			H[jj*bdofs + ii+1] += He[j*4 + (i*2)+1];
						
		}
	}

#ifdef DEBUG_PRINT
	printf("\n");
	
	for (int i=0; i<bdofs; i++)
		for (int j=0; j<fdofs; j++)
			printf("H[%ld][%ld]: %lf\n", i,j, H[i + j*bdofs]);					
#endif	

}


void assembleLineAcousticElasticCouplingMatricies(	double* H, 
													double* verts, int* inds, double* norms, 
													int* map,
													AnalysisInfo info)
{
	GeomInfo		ginfo = info.geominfo;
	CoupleInfo		cinfo = info.coupleinfo;	
	
	double* He = (double*) calloc(4*2, sizeof(double));
	
	double* coords = (double*) calloc(ginfo.numelnodes*ginfo.dimension,sizeof(double));
	double* normals = (double*) calloc(ginfo.dimension,sizeof(double));	

	int el;
	for ( el=0; el<ginfo.numelements; el++){ 
	
		memset(He,0,4*2*sizeof(double));
			
		get2dCoords( coords  ,el, verts,inds, ginfo); 
		
		normals[0] = norms[el];
		normals[1] = norms[el+	ginfo.numelements];
		
#ifdef DEBUG_PRINT
		printf("\n\n\nelement: %ld\n", el);
#endif	
		lineAcousticElasticCouplingMatrix(He, coords, normals, info);							
					
#ifdef DEBUG_PRINT
		printf("\n\n\n\n\n");
		for ( i=0; i<4; i++)
			for ( j=0; j<2; j++)
				printf("He[%ld][%ld]: %lf\n", i,j, He[i + j*4]);					
#endif	
		assembleH(H, el, He, ginfo, cinfo, inds, map);
	}
	

	free(He);
	free(coords), free(normals);
}




#endif