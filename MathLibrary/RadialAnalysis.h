/*
 *  RadialAnalysis.h
 *  AxisSymmetricSolidApp


 * Pre-discretization eigendecomposition

 *  Created by Cynthia Bruyns on 12/20/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>

#include <vecLib/clapack.h>
#include <vecLib/vDSP.h>

#include "ModelSharedData.h"

#include "CLPKComplexOpsCM.h"
#include "CLPKRealMatrixOpsCM.h"

#import "TriAxisSymMapElasticElement.h"
#import "TriAxisSymBaryElasticElement.h"
#import "QuadAxisSymElasticElement.h"

/*!
 compute the eigeninformation for a radially axis-symmetric object
 */
inline void radialAnalyze(GeomData* geomdata, EigenData* eigendata, AnalysisInfo info, int rankPerm, int slice)
							
{	
	double theta = 2.0*M_PI/ (double) rankPerm;

	int rankBlock = geomdata->numverts*3;	
	int evalsSize = rankBlock*sizeof(__CLPK_doublereal);

	int chunkV	  = rankBlock*rankBlock;
	int evecSize = chunkV*sizeof(__CLPK_doublereal);
	
	double* K = (double*) malloc(evecSize);
	double* M = (double*) malloc(evecSize);
	//double* LKL = (double*) malloc(evecSize);
	//double* invL = (double*) malloc(evecSize);
	
	//double* VTL = (double*) malloc(evecSize);	
	//double* LTV = (double*) malloc(evecSize);

	__CLPK_doublereal *wr;					
    wr = (__CLPK_doublereal*)malloc(evalsSize);

	__CLPK_doublereal *vr;              /* Left and right eigenvector matrices */
    vr = (__CLPK_doublereal*)malloc(evecSize);
	
	//double time = currentTime();
	
	int i,j;
	for ( i=slice; i<(slice+1); i++){ 
	
		memset(K, 0, evecSize);
		memset(M, 0, evecSize);
		//memset(LKL,0,evecSize);
		//memset(invL, 0, evecSize);
	
		if (geomdata->vertexperface ==3)
			//assembleTriAxisSymBaryElasticMatricies
			assembleTriAxisSymMapElasticMatricies(	K, M, 		
													geomdata->vertex_pos, geomdata->faces,
													info,
													theta,
													i);
		else
			assembleQuadAxisSymElasticMatricies(	K, M, 
													geomdata->vertex_pos, geomdata->faces,
													info,
													theta,
													i);	

																								
		applyBoundaryConditions(K, M, info.constinfo, rankBlock); 		
		
		//double time = currentTime();
		//printf("made stiffness matrix: %lf ms\n", time*1E3);
		
		//realSingleReduceToStandardEigenProblem( LKL, K, M, invL, rankBlock);
		
		//printf("testing cholesky\n");
		//symmetricTest(LKL, rankBlock);	
		
		//memcpy(vr, LKL, evecSize);
		//realSingleSymmetricEigenDecomp(wr, vr, rankBlock); 
		//projectOntoCholesky(vr, VTL, LTV, invL, rankBlock);
				
		symmetricTest(K, rankBlock);
		symmetricTest(M, rankBlock);
		memcpy(vr, K, evecSize);
		realDoubleSymmetricGeneralEigenDecomp(wr, vr, M, rankBlock);

		//printf("m: %d\n", i);
		for (j=0; j<rankBlock; j++)
			printf("%lf\n", sqrt(wr[j])/2.0/M_PI);

		//time = currentTime() - time;
		//printf("time to full decompose: %lf ms\n", time*1E3);	
		
		//memcpy(eigendata->EigenValues+(i*rankBlock), wr, evalsSize);		
		//memcpy(eigendata->EigenVectors+(i*chunkV), vr, evecSize); 
		//memcpy(eigendata->VTL+(i*chunkV), VTL, evecSize);
		//memcpy(eigendata->LTV+(i*chunkV), LTV, evecSize);
		
	}
		
	//time = currentTime() - time;
	//printf(" time to radial analyze: %lf ms\n", time*1E3); 
	
	
	if (vr) free(vr); if (wr) free(wr);									
	//free(VTL), free(LTV);
	//free(invL), free(LKL);
	if (K) free(K); if (M) free(M); 
}
