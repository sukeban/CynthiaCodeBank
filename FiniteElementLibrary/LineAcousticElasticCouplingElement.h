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



void	lineAcousticElasticCouplingMatrix(		double* H,
												double* coords,
												AnalysisInfo info)
{	
	double L = getLineLength(coords);
	
	double thickness = 1.0;
	
	H[1] = 7.0/20.0;	H[7] = 3.0/20.0;
	H[2] = L/20.0;		H[8] = L/30.0;
	
	H[4] = 3.0/20.0;	H[10] = 7.0/20.0;
	H[5] = -L/30.0;		H[11] = -L/20.0;
	
	for (int i=0; i<12; i++){
		H[i] *= L;
	}
	
#ifdef DEBUG_PRINT
	printf("h=[\n");
	for (int i=0; i<6; i++){
		for (int j=0; j<2; j++){
			printf("%lf\t", H[i+j*6]);
		}
		printf(";\n");
	}
	printf("];\n");
#endif	
}


inline void assembleH(double* H, int id, double* He, GeomInfo ginfo, CoupleInfo cinfo, int* inds, int* map)
{	
	int i,j,ii,jj;

	int bdofs = cinfo.numbeamverts*3;	
	int fdofs = cinfo.numfluidverts;
	
#ifdef DEBUG_PRINT
	printf("\n\n\nbdofs: %d fdofs: %d\n", bdofs, fdofs);
	
	for (int i=0; i<6; i++)
		for (int j=0; j<2; j++)
			printf("He[%ld][%ld]: %lf\n", i,j, He[i + j*6]);	
#endif	
	// how to tell this, a lookuptable of
	//[ couple_node_index beam_node_index fluid__node_index]

	for ( i=0; i<2; i++){				// for each node in the couple element
	
		int iindex = inds[id + i*ginfo.numelements];
#ifdef DEBUG_PRINT
		printf("iindex: %d\n", iindex);
#endif		
		ii = map[iindex+ ginfo.numnodes]*3;	// look up the associated beam node and get its dof
#ifdef DEBUG_PRINT
		printf("ii: %d\n", ii);
#endif		
		for ( j=0; j<2; j++){				// for each node in the couple element
		
			int jindex = inds[id + j*ginfo.numelements];
#ifdef DEBUG_PRINT
			printf("jindex: %d\n", jindex);
#endif		
			jj = map[jindex + 2*ginfo.numnodes];	// get the fluid node
			
			
#ifdef DEBUG_PRINT
			printf("jj: %d\n", jj);
#endif		


			//printf("wriring to H[%d][%d]=%lf\n", ii, jj, He[j*6 + (i*3)  ]);	
			H[jj*bdofs + ii  ] += He[j*6 + (i*3)  ];
			
			
			//printf("wriring to H[%d][%d]=%lf\n", ii+1, jj, He[j*6 + (i*3)+1  ]);	
			H[jj*bdofs + ii+1] += He[j*6 + (i*3)+1];		


			//printf("wriring to H[%d][%d]=%lf\n", ii+2, jj, He[j*6 + (i*3)+2  ]);	
			H[jj*bdofs + ii+2] += He[j*6 + (i*3)+2];		
						
		}
	}

#ifdef DEBUG_PRINT
	printf("H=[\n");	
	for (int i=0; i<bdofs; i++){
		for (int j=0; j<fdofs; j++){
			printf(" %lf\t",  H[i + j*bdofs]);	
		}
		printf(";\n");
	}
	printf("];");
#endif	


}


void assembleLineAcousticElasticCouplingMatricies(	double* H, 
													double* verts, int* inds,
													int* map,
													AnalysisInfo info,
													CoupleInfo cinfo)
{
	GeomInfo		ginfo = info.geominfo;
		
	
	
	int bdofs = cinfo.numbeamverts*3;
	int fdofs = cinfo.numfluidverts;
	
	double* He = (double*) calloc(6*2, sizeof(double));
	
	double* coords = (double*) calloc(ginfo.numelnodes*ginfo.dimension,sizeof(double));

	int el;
	for ( el=0; el<ginfo.numelements; el++){ 
	
		memset(He,0,6*2*sizeof(double));
			
		get2dCoords(coords, el, verts, inds, ginfo); 
				
#ifdef DEBUG_PRINT
		printf("\n\n\nelement: %ld\n", el);
#endif	
		lineAcousticElasticCouplingMatrix(He, coords, info);							
		assembleH(H, el, He, ginfo, cinfo, inds, map);
	}
	

#ifdef DEBUG_PRINT
	printf("H=[\n");	
	for (int i=0; i<bdofs; i++){
		for (int j=0; j<fdofs; j++){
			printf(" %lf\t",  H[i + j*bdofs]);	
		}
		printf(";\n");
	}
	printf("];");
#endif	

	free(He);
	free(coords);
}




#endif