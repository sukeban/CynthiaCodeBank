/*
 *  3dElements.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 9/5/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _3dElements_H_
#define _3dElements_H_

#include "MatrixOps.h"

#pragma mark ____Solid___

double brickVolume(double* coords) 
{
	//cbnote compute this
	return 0;
}
 

#pragma mark ___Tetra___
void getTetraCenter(double* coords, double* center)
{
	center[0] = center[1] = center[2] = 0.0;
	int i;
	for (i=0; i<4; i++){
		center[0] += coords[i];
		center[1] += coords[i+4];
		center[2] += coords[i+8];
	}
	center[0] /= 4.0;
	center[1] /= 4.0;
	center[2] /= 4.0;
}


double tetraVolume(double* coords) 
{
	double x[16];
	int i;
	x[0] = x[1] = x[2] = x[3] = 1.0;
	for ( i=0; i<12;i++){
		x[i+4] = coords[i];
		//printf("x: %lf\n", x[i+4]);
	}
	
	
	double vol = (1.0/6.0)*dabs(m4_det(x));
	
	if (vol == 0.0){
		vol = 1e-6;
	}
	
#ifdef DEBUG_PRINT
	printf("tetra vol: %lf\n", vol);
#endif


	return vol;
}


void tetraBaseNormal(double* coords, double* normal)
{
   double u[3];
   u[0] = coords[1] - coords[0]; // 1-0
   u[1] = coords[5] - coords[4];
   u[2] = coords[9] - coords[8];

   double v[3];
   v[0] = coords[2 ] - coords[0]; // 2-0
   v[1] = coords[6 ] - coords[4];
   v[2] = coords[10] - coords[8];

   double dir[3];
   v3_cross(u,v, dir);
   v3_normalize(dir);
}
 
void tetraBaseCenter(double* coords, double* center)
{
	int i;
	center[0] = center[1] = center[2] = 0.0;
	for ( i=0; i<3; i++){
		   center[0]+= coords[i];
		   center[1]+= coords[i+4];
		   center[2]+= coords[i+8];
	}
	center[0]/=3.0;
	center[1]/=3.0;
	center[2]/=3.0;
}
 
bool tetraBaseIsFacingInwards(double* coords) 
{
   double normal[3];
   tetraBaseNormal(coords, normal);

   double center[3];
   tetraBaseCenter(coords, center);

   double w[3];
   w[0] = coords[3 ] - center[0]; //3-center
   w[1] = coords[7 ] - center[1];
   w[2] = coords[11] - center[2];
   v3_normalize(w);

   double dist = v3_dot(normal,w)/v3_dot(w,w);// angle between the top node and the bottom face normal

   if (dist < 0.0) return 1;
   return 0;
}

bool checkTetraElement(double* coords)
{
	bool inwards = false;
   double tempcoords[4*3];
   memcpy(tempcoords, coords, 4*3*sizeof(double));

   // advaned p. 47
   // make sure 1,2,3 is cw wrt 4
   if ( tetraBaseIsFacingInwards(coords) ){

		   tempcoords[1] = coords[2];
		   tempcoords[5] = coords[6];
		   tempcoords[9] = coords[10];

		   tempcoords[2] = coords[1];
		   tempcoords[6] = coords[5];
		   tempcoords[10] = coords[9];

		  // printf("is facing inwards swapping node 1 and 2\n");
		   inwards = true;
   }

   memcpy(coords, tempcoords, 4*3*sizeof(double));
   return inwards;
}
 
 

#pragma mark ___General___

void threeDJacobian(int numnodes, double* J, double* invJ, double &detJ, double* dNdSdN, double* x)
 {
	int i;
	int dim = 3;
	
	memset(J,0,(dim*dim)*sizeof(double));	
	
	for ( i=0; i<numnodes; i++){
		J[0]+=dNdSdN[i]*x[i           ];
		J[3]+=dNdSdN[i]*x[i+numnodes  ];
		J[6]+=dNdSdN[i]*x[i+2*numnodes];
	 
		J[1]+=dNdSdN[i+numnodes]*x[i           ];
		J[4]+=dNdSdN[i+numnodes]*x[i+numnodes  ];
		J[7]+=dNdSdN[i+numnodes]*x[i+2*numnodes];
 
		J[2]+=dNdSdN[i+2*numnodes]*x[i           ];
		J[5]+=dNdSdN[i+2*numnodes]*x[i+numnodes  ];
		J[8]+=dNdSdN[i+2*numnodes]*x[i+2*numnodes];		
	}

	detJ =  m3_det(J);
	
	m3_inverse(invJ,J);
		
#ifdef DEBUG_PRINT
	int j;
	for ( i=0; i<dim; i++)
		for ( j=0; j<dim; j++)
			printf("J[%ld][%ld]: %llf\n",i,j, J[j*dim+i]);
		
	for ( i=0; i<dim; i++)
		for ( j=0; j<dim; j++)
			printf("invJ[%ld][%ld]: %llf\n",i,j, invJ[j*dim+i]);
		
	printf("detJ: %lf\n", detJ);
#endif
 }



 void threeDPhysicalDerivs(int numnodes, double* dNdXdY, double* invJ, double* dNdSdN)
{	
	int i;
	
	for ( i=0; i<numnodes; i++){	
		dNdXdY[i           ] = dNdSdN[i   ]*invJ[0] + dNdSdN[numnodes+i]*invJ[3] + dNdSdN[2*numnodes+i]*invJ[6]; //dNi / dx		
		dNdXdY[numnodes+i  ] = dNdSdN[i   ]*invJ[1] + dNdSdN[numnodes+i]*invJ[4] + dNdSdN[2*numnodes+i]*invJ[7]; //dNi / dy
		dNdXdY[2*numnodes+i] = dNdSdN[i   ]*invJ[2] + dNdSdN[numnodes+i]*invJ[5] + dNdSdN[2*numnodes+i]*invJ[8]; //dNi / dz
	}
	
#ifdef DEBUG_PRINT
	int j;
	int dim = 3;
	for ( i=0; i< numnodes; i++)
		for ( j=0; j<dim; j++)
			printf("dNdXdY[%ld][%ld]: %llf\n",i,j, dNdXdY[j*numnodes+i]);	
#endif

}

void	threeDConstMatrix(int numnodes, double* B, double* dNdXdY)
{	
	int i,j;
	int sizeB1 = 6;
	int dim = 3;
	int sizeBB = dim*sizeB1;
	
	for ( i=0; i< numnodes; i++){	
	
		j = i*sizeBB; 
		 
		B[j   ]	= dNdXdY[i           ];	
		B[j+7 ]	= dNdXdY[numnodes+i  ];	
		B[j+14]	= dNdXdY[2*numnodes+i];	
		
		B[j+3 ]	= dNdXdY[numnodes+i];
		B[j+9 ]	= dNdXdY[i         ];
		
		B[j+10]	= dNdXdY[2*numnodes+i];
		B[j+16]	= dNdXdY[numnodes+i  ];
		
		B[j+5 ]	= dNdXdY[2*numnodes+i];
		B[j+17]	= dNdXdY[i           ];
			
	}
	

#ifdef DEBUG_PRINT

	for ( i=0; i<sizeB1; i++)
		for ( j=0; j<numnodes*dim; j++)
			printf("B[%ld][%ld]: %lf\n",i,j, B[j*sizeB1+i]);
#endif			
}


#endif
