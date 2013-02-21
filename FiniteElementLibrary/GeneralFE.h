/*
 *  GeneralFE.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 1/9/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
 #ifndef _GeneralFE_Header_
 #define _GeneralFE_Header_
 
#include <vecLib/vectorOps.h>
#include <vecLib/clapack.h>
#include <vecLib/vBLAS.h>

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>

#include <math.h>
#include <g2cEdited.h>

#include "SharedFEData.h"
#include "CLPKRealMatrixOpsCM.h"

/*!
 compute the length of a line given an array of coordinates
 */
double getLineLength(double* coords)
{
	double dx = -coords[0] + coords[1];
	double dy = -coords[2] + coords[3];
#ifdef DEBUG_PRINT
	printf("length: %lf\n", dabs(dx)+dabs(dy));
#endif	
	return (dabs(dx)+dabs(dy));

}

#pragma mark ___Constants____
/*!
 compute the material matrix given the material constants for an object
 */
inline void	materialMatrix(double* D, MaterialInfo constants)
{	
	int i;
	int sizeD;
	int sizeDD;
	
	double c;
	
	switch (constants.state){
		case planestress:{
			// plane stress, going left to right then down
			c = constants.elasticModulus/(1.0 - powf(constants.poisson,2.0));
			
			D[0] = c*1.0;
			D[1] = c*constants.poisson;
			D[2] = 0.0;
			
			D[3] = c*constants.poisson;
			D[4] = c*1.0;
			D[5] = 0.0;
			
			D[6] = 0.0;
			D[7] = 0.0;
			D[8] = c*(1.0-constants.poisson)/2.0;
			
			sizeD = 3;
			sizeDD = sizeD*sizeD;

		break;
		}   
		case planestrain:{
			// plane strain
			
			c = constants.elasticModulus/((1.0 + constants.poisson)*(1.0 - 2.0*constants.poisson));
			
			D[0] = c*(1.0-constants.poisson);
			D[1] = c*constants.poisson;
			D[2] = 0;
			
			D[3] = c*constants.poisson;
			D[4] = c*(1.0-constants.poisson);
			D[5] = 0.0;
			
			D[6] = 0.0;
			D[7] = 0.0;
			D[8] = c*(1.0-2.0*constants.poisson)/2.0;
			
			sizeD = 3;
			sizeDD = sizeD*sizeD;

			
		 break;
		 }
		 case pure_axissymmetric:{
			D[0] = (1.0-constants.poisson);
			D[1] = constants.poisson;
			D[2] = constants.poisson;
			D[3] = 0.0;
			
			D[4] = constants.poisson;
			D[5] = (1.0-constants.poisson);
			D[6] = constants.poisson;
			D[7] = 0.0;
			
			D[8] = constants.poisson;
			D[9] = constants.poisson;
			D[10] = (1.0-constants.poisson);
			D[11] = 0.0;
			
			D[12] = 0.0;
			D[13] = 0.0;
			D[14] = 0.0;
			D[15] = (1.0-2.0*constants.poisson)/2.0;

			sizeD = 4;
			sizeDD = sizeD*sizeD;

			for ( i=0; i<sizeDD; i++){
				D[i] *= constants.elasticModulus/((1.0 + constants.poisson)*(1.0 - 2.0*constants.poisson));
			}
		 
		 break;
		}
		 case solid:{
			D[0] = (1.0-constants.poisson);
			D[1] = constants.poisson;
			D[2] = constants.poisson;
			D[3] = 0.0;			
			D[4] = 0.0;
			D[5] = 0.0;
			
			D[6] = constants.poisson;	
			D[7] = (1.0-constants.poisson);	
			D[8] = constants.poisson;
			D[9] = 0.0;
			D[10] = 0.0;
			D[11] = 0.0;
			
			D[12] = constants.poisson;
			D[13] = constants.poisson;
			D[14] = (1.0-constants.poisson);	
			D[15] = 0.0;
			D[16] = 0.0;
			D[17] = 0.0;
			
			D[18] = 0.0;
			D[19] = 0.0;
			D[20] = 0.0;
			D[21] = (1.0-2.0*constants.poisson)/2.0;
			D[22] = 0.0;
			D[23] = 0.0;

			D[24] = 0.0;
			D[25] = 0.0;
			D[26] = 0.0;
			D[27] = 0.0;
			D[28] = (1.0-2.0*constants.poisson)/2.0;
			D[29] = 0.0;

			D[30] = 0.0;
			D[31] = 0.0;
			D[32] = 0.0;
			D[33] = 0.0;
			D[34] = 0.0;
			D[35] = (1.0-2.0*constants.poisson)/2.0;
			
			sizeD = 6;
			sizeDD = sizeD*sizeD;

			for ( i=0; i<sizeDD; i++){
				D[i] *= constants.elasticModulus/((1.0 + constants.poisson)*(1.0 - 2.0*constants.poisson));
			}
					 
		 break;
		 }
		 case acoustic:{
			D[0] = 1.0;
			D[1] = 0.0;
			
			D[2] = 0.0;
			D[3] = 1.0;	
			
			sizeD = 2;		
		 
		 break;
		 }
	}
	
#ifdef DEBUG_PRINT
	int j;
	for ( i=0; i<sizeD; i++)
		for ( j=0; j<sizeD; j++)
			printf("D[%ld][%ld]: %lf\n", i, j, D[j*sizeD + i]);
#endif

}


#pragma mark ___Test___
/*!
 test whether a stiffness matrix came back symmetric or not
 */
inline void symmetricTest(double* K, int size)
{
	//double eps = 1e-2;
	int i,j;
	for ( i=0; i<size; i++){
		for ( j=0; j<size; j++){
		
			//if (dabs(K[j*size+i] - K[i*size+j]) > eps){
#ifdef DEBUG_PRINT
				printf("not symmetric K[%d][%d] : %lf K[%d][%d] :%lf\n",i,j, K[j*size+i], j,i, K[i*size+j]);
#endif				
				K[j*size+i] = (K[j*size+i]+K[i*size+j])/2.0;
				K[i*size+j] = K[j*size+i];
			//}
				
		}		
	}
}



#pragma mark ___GetCoords___

/*!
 grab 3d element coordinates from large coordinates list
 */
inline void get3dCoords(double* coords, int id, double* verts, int* inds, GeomInfo ginfo)
{	
	int i;
	int idx;


	for ( i=0; i<ginfo.numelnodes; i++){
		
		idx = id + i*ginfo.numelements;

#ifdef DEBUG_PRINT
		printf("id: %d i: %ld idx: %ld vertindex: %ld\n", id, i, idx, inds[idx]);
#endif
		coords[i			         ] = ginfo.scale*verts[inds[idx]                ]; 		
		coords[i + ginfo.numelnodes	 ] = ginfo.scale*verts[inds[idx]+ginfo.numnodes ];		
		coords[i + ginfo.numelnodes*2] = ginfo.scale*verts[inds[idx]+ginfo.numnodes*2];
	}
	
#ifdef DEBUG_PRINT
	for ( i=0; i<ginfo.numelnodes; i++)
		printf("getting [%ld]: %lf %lf %lf\n",i, coords[i], coords[ginfo.numelnodes+i], coords[2*ginfo.numelnodes+i]);
#endif
}

/*!
 grab 2d element coordinates from large coordinates list
 */
inline void get2dCoords(double* coords, int id, double* verts, int* inds, GeomInfo ginfo)
{	

	int i;
	int idx;
	for ( i=0; i<ginfo.numelnodes; i++){
		
		idx = id + i*ginfo.numelements;

#ifdef DEBUG_PRINT
		printf("i: %ld idx: %ld vertindex: %ld\n", i, idx, inds[idx]);
#endif
		coords[i			         ] = ginfo.scale*verts[inds[idx]                ]; 		
		coords[i + ginfo.numelnodes	 ] = ginfo.scale*verts[inds[idx]+ ginfo.numnodes ];		
	}
	
#ifdef DEBUG_PRINT
	for ( i=0; i<ginfo.numelnodes; i++)
		printf("getting coords[%ld]: %lf %lf\n",i, coords[i], coords[ginfo.numelnodes+i]);
#endif
}

#pragma mark ___GetNormals___

/*!
 grab 2d element normals from large normals list
 */
inline void get2dNormals(double* normals, int id, double* norms, int* inds, GeomInfo ginfo)
{	
	printf("\n");
	int i;
	int idx;
	for ( i=0; i<ginfo.numelnodes; i++){
		
		idx = id + i*ginfo.numelements;

#ifdef DEBUG_PRINT
		printf("i: %ld idx: %ld vertindex: %ld\n", i, idx, inds[idx]);
#endif
		normals[i			         ] = norms[inds[idx]                ]; 		
		normals[i + ginfo.numelnodes ] = norms[inds[idx]+ ginfo.numnodes ];		
	}
	
#ifdef DEBUG_PRINT
	for ( i=0; i<ginfo.numelnodes; i++)
		printf("getting normals[%ld]: %lf %lf\n",i, normals[i], normals[ginfo.numelnodes+i]);
#endif
}


#pragma mark ___Assembly___

/*!
 determine where to place the infomration in the larger matrix based on vertex index and degrees of freedom
 */
inline void findDofs(int* index, int id, int* inds, GeomInfo ginfo)
{
	
   int k =0 ;
   int start;
   int i,j;
   for ( i=0; i<ginfo.numelnodes; i++){
      start = inds[id + i*ginfo.numelements]*ginfo.nodedofs;
	  
#ifdef DEBUG_PRINT
	printf("start: %ld\n", start);
#endif	  
	   for ( j=0; j<ginfo.nodedofs; j++){        
         index[k] = start+j; k++;			
		}
   }
   
#ifdef DEBUG_PRINT
	for ( i=0; i<ginfo.numelnodes*ginfo.nodedofs; i++)
		printf("index[%ld]: %ld\n", i, index[i]);
#endif
	
}

/*!
 assemble the smaller matricies stiffness and mass matricies into the larger ones
 */
inline void assemble(double* K, double* M, int id, double* Ke, double* Me, GeomInfo ginfo, int* inds)
{
	int dofs = ginfo.numnodes*ginfo.nodedofs;
	int edofs = ginfo.numelnodes*ginfo.nodedofs;
	
	int* index = (int*) calloc(edofs, sizeof(int));
	findDofs(index, id, inds, ginfo);
	
	int i,j,ii,jj;
	for ( i=0; i<edofs; i++){
		ii = index[i];
	
		for ( j=0; j<edofs; j++){
			jj = index[j];
			
			//printf("adding %lf to %lf to make: %lf\n",  
			//Ke[j*edofs + i], K[jj*dofs + ii],K[jj*dofs + ii] + Ke[j*edofs + i]);
			
			
			K[jj*dofs + ii] += Ke[j*edofs + i];
			
			
			M[jj*dofs + ii] += Me[j*edofs + i];

		}
	}


#ifdef DEBUG_PRINT
		printf("\n\n\n\nassemble:\n");
		for ( i=0; i<dofs; i++){
			for ( j=0; j<dofs; j++){
				printf("K[%ld][%ld]: %lf\n", j,i, K[j*dofs + i]);
			}
			printf("\n");
		}
																					

		for ( i=0; i<dofs; i++){
			for ( j=0; j<dofs; j++){
				printf("M[%ld][%ld]: %lf\n", j,i, M[j*dofs + i]);
			}
			printf("\n");
		}					
#endif	

	free(index);
}

/*!
 compute the mass matrix for a quadrilateral element using consistent mass matrix approach
 */
inline void	accumulateMass(	double* M, int sizeM, 
							double* N, int sizeN[],
							double wx, double wy, double wz,
							double rho)
{
	

/*
	http://www.netlib.org/blas/sgemm.f
*/				 
															//				      m n k m 
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,		//M = rho*N^T N  12x3 3x12
				sizeN[1],   //m
				sizeN[1],	//n						
				sizeN[0],	//k							
				 wx*wy*rho, //alpha
				 N,											
				 sizeN[0],	//k								
				 N, 
				 sizeN[0],	//k							
                 1.0, 
				 M,			//
				 sizeN[1]);	//m								


#ifdef DEBUG_PRINT
	int i,j;
	printf("\n\naccumulate mass \n");
	for ( i=0; i<sizeM; i++){
		for ( j=0; j<sizeM; j++){
			printf("%lf\t", M[j*sizeM+i]);
		}
		printf(";\n");
	}
#endif	
	
}

#pragma mark ___Constraints___

/*!
 knock out the rows and columns of the stiffness matrix to apply the appropriate boundary conditions
 */
void applyBoundaryConditions(double* K, double* M, ConstraintInfo cinfo, int size)
{
	int n = cinfo.numconstraints;
	
	int i,j,c;
	
	for( i=0; i< n; i++){
	
		c = cinfo.constraints[i];
		
		for (j=0; j<size; j++){
			K[c + j*size] = 0.0;						
			K[j + c*size] = 0.0;				
		}
		M[c+c*size] = 1.0;
	}
	
/*	for ( i=0; i<size; i++)
		for ( j=0; j< size; j++)
			printf("Kc[%ld][%ld]: %lf\n", i,j, K[j*size+i]);
*/
}

/*!
 compress the matricies to remove zero elements from stiffness matrix
*/
void removeRowCols(double* K, double* M, double* Knew, double* Mnew, ConstraintInfo cinfo, int size)
{
	int n = cinfo.numconstraints;
	
	int i,j,k;
	int c;
	
	int newsize = size-n;
	
	int ii, jj;
	
	bool found ;
	
	ii = 0;
	for( i=0; i< size; i++){
		jj = 0;
		
		found = false;
		for (k=0; k<n; k++){
			c = cinfo.constraints[k];
			if (i == c)
				found = true;
		
		}
		if (found) continue;
		
		
		for (j=0; j<size; j++){
		
			found = false;	
			for (k=0; k<n; k++){
				c = cinfo.constraints[k];
				if (j == c)
					found = true;		
			}
			if (found) continue;
			
			Knew[ii + jj*newsize] = K[i+j*size];
			Mnew[ii + jj*newsize] = M[i+j*size];
			jj++;					
			
		}
		
		ii++;
	}
	
#ifdef DEBUG_PRINT
	for (i=0; i<n; i++)
		printf("c[%d]: %d\n", i, cinfo.constraints[i]);
	
	for ( i=0; i<size; i++)
		for ( j=0; j< size; j++)
			printf("K[%ld][%ld]: %lf\n", i,j, K[j*size+i]);

	
	for ( i=0; i<newsize; i++)
		for ( j=0; j< newsize; j++)
			printf("Knew[%ld][%ld]: %lf\n", i,j, Knew[j*newsize+i]);
#endif

}

/*!
 knock out the rows and columns of the hessian matrix to apply the appropriate boundary conditions
 */
void applyBoundaryConditionsH(double* H, ConstraintInfo cinfo, int size[])
{
	int n = cinfo.numconstraints;
	
	int i,j,c;
	for( i=0; i< n; i++){
	
		c = cinfo.constraints[i];
		
		for (j=0; j<size[1]; j++){
			H[c + j*size[0]] = 0.0;						
		}
		
	}
	
/*	for ( i=0; i<size[0]; i++)
		for ( j=0; j< size[1]; j++)
			printf("Hc[%ld][%ld]: %lf\n", i,j, H[j*size[0]+i]);
*/
}

/*!
 compress the hessian matrix to remove zero elements from stiffness matrix
 */
void removeRowColsH(double* H, double* Hnew,  ConstraintInfo cinfo, int size[])
{
	int n = cinfo.numconstraints;
	
	int i,j,k;
	int c;
	
	int newsize = size[0]-n;
	
	
	int ii;
	
	bool found ;
	
	ii = 0;
	for( i=0; i< size[0]; i++){
		
		found = false;
		for (k=0; k<n; k++){
			c = cinfo.constraints[k];
			if (i == c)
				found = true;
		
		}
		if (found) continue;
		
		
		for (j=0; j<size[1]; j++){
		
			Hnew[ii + j*newsize] = H[i+j*size[0]];
								
			
		}
		
		ii++;
	}
	
#ifdef DEBUG_PRINT
	for (i=0; i<n; i++)
		printf("c[%d]: %d\n", i, cinfo.constraints[i]);
	
	for ( i=0; i<size[0]; i++)
		for ( j=0; j< size[1]; j++)
			printf("H[%ld][%ld]: %lf\n", i,j, H[j*size[0]+i]);

	
	for ( i=0; i<newsize; i++)
		for ( j=0; j< size[1]; j++)
			printf("Hnew[%ld][%ld]: %lf\n", i,j, Hnew[j*newsize+i]);
#endif

}


#endif
