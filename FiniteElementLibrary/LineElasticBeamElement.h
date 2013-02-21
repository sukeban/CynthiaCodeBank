/*
 *  LineElasticBeamElement.h
 *  FiniteElement
 *
 *  Created by Cynthia Maxwell on 9/24/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef _LineElasticBeamElement_Header_
#define _LineElasticBeamElement_Header_


#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

#include "SharedFEData.h"
#include "GeneralFE.h"

/*!
 compute the stiffness matrix for a linear elastic beam element
*/
void	lineBeamElasticStiffnessMatrix(	double* K,
										double* coords, 
										AnalysisInfo info, 
										BeamInfo binfo) 
{	
	MaterialInfo minfo = info.materialinfo;
	double E = minfo.elasticModulus;
	
	double I = binfo.momentOfInertia;
	double A = binfo.crossSectionalArea;

	double L = getLineLength(coords);	
	double L2  = pow(L,2.0);
	double L3 = pow(L,3.0);
	
	K[0] = E*A/L;
	K[3] = -E*A/L;
	
	K[7] = 12*E*I/L3;
	K[8] = 6*E*I/L2;	
	K[10] = -12*E*I/L3;
	K[11] = 6*E*I/L2;
	
	K[13] = 6*E*I/L2;
	K[14] = 4*E*I/L;	
	K[16] = -6*E*I/L2;
	K[17] = 2*E*I/L;
	
	K[18] = -E*A/L;
	K[21] = E*A/L;
		
	K[25] = -12*E*I/L3;
	K[26] = -6*E*I/L2;	
	K[28] = 12*E*I/L3;
	K[29] = -6*E*I/L2;
	
	K[31] = 6*E*I/L2;
	K[32] = 2*E*I/L;	
	K[34] = -6*E*I/L2;
	K[35] = 4*E*I/L;

}

/*!
 compute the consistent mass matrix for a linear elastic beam element
 */
void	lineBeamElasticConsistentMassMatrix(double* M, double* coords, MaterialInfo minfo, BeamInfo binfo)
{
	double mass = minfo.density*binfo.crossSectionalArea;
	
	double L = getLineLength(coords);	
	double L2  = pow(L,2.0);
	
	M[0] = 140.0;
	M[3] = 70.0;	
	
	M[7] = 156.0;
	M[8] = 22.0*L;	
	M[10] = 54.0;
	M[11] = -13*L;
	
	M[13] = 22.0*L;
	M[14] = 4.0*L2;	
	M[16] = 13*L;
	M[17] = -3.0*L2;
	
	M[18] = 70.0;
	M[21] = 140.0;
	
	M[25] = 54.0;
	M[26] = 13.0*L;
	M[28] = 156.0;
	M[29] = -22.0*L;
	
	M[31] = -13.0*L;
	M[32] = -3.0*L2;
	M[34] = -22.0*L;
	M[35] = 4.0*L2;
	
	for (int i=0; i<36; i++)
		M[i] *= (mass*L/420.0);
	
}

/*!
 compute the mass and stiffness matricies for a linear elastic beam elmement
 */
void assembleLineBeamElasticMatricies(	double* K, double* M, 
										double* verts, int* inds, 
										AnalysisInfo info,
										BeamInfo binfo)
{
	MaterialInfo		minfo = info.materialinfo;
	GeomInfo			ginfo = info.geominfo;
	
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
		lineBeamElasticStiffnessMatrix(  Ke, coords, info, binfo);							
		lineBeamElasticConsistentMassMatrix( Me, coords, minfo, binfo);
					
#ifdef DEBUG_PRINT
		int i,j;
		printf("\n\n\n\n\n");
		for ( i=0; i<edofs; i++)
			for ( j=0; j<edofs; j++)
				printf("Ke[%ld][%ld]: %lf\n", i,j, Ke[i*edofs + j]);					
																					
		for ( i=0; i<edofs; i++)
			for ( j=0; j<edofs; j++)
				printf("Me[%ld][%ld]: %lf\n", i,j, Me[i*edofs + j]);					
					
#endif	
		assemble(K, M, el, Ke, Me, ginfo, inds);
	}

	//int dofs = ginfo.numnodes*ginfo.nodedofs;	// cbnote not needed for symmetric eigendecomp
	//symmetricTest(K, dofs);

	free(Ke), free(Me);
	free(coords);
}

#endif