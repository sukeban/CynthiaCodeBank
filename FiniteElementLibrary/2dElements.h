/*
 *  2dElements.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 1/9/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _2dElements_Header_
#define _2dElements_Header_
 
#include <math.h>
#include <stdlib.h>

#include "SharedFEData.h"
#include "g2cEdited.h"

#include "MatrixOps.h"


#pragma mark _Quads___

/*!
 compute the area of a quadrilateral from the list of coordinates
 */
inline double quadArea(double* coords) //assuming all same z ccords
{
	double p[3];
	double q[3];
	double w[3];
	
	p[0] = coords[2] - coords[0];
	p[1] = coords[6] - coords[4];
	p[2] = 0;
	
	q[0] = coords[3] - coords[1];
	q[1] = coords[7] - coords[5];
	q[2] = 0;
	
	v3_cross(p,q,w);
	double area = 0.5*v3_length(w);
	
#ifdef DEBUG_PRINT
	printf("area: %lf\n", area);
#endif	

	return area;
}

/*!
 compute the center of a quadrilateral from the list of coordinates and return them in x and y
 */
inline void getQuadCenter(double* coords, double &x, double &y)
{
	double x1,x2,x3,x4, y1,y2,y3,y4;
	x1 = coords[0]; y1 = coords[4];
	x2 = coords[1]; y2 = coords[5];
	x3 = coords[2]; y3 = coords[6];
	x4 = coords[3]; y4 = coords[7];
  
    x=(x1+x2+x3+x4)/4;                      // x-centroid of triangle
    y=(y1+y2+y3+y4)/4;                      // y-centroid of triangle
}

/*!
 compute the lumped mass matrix (http://en.wikipedia.org/wiki/Mass_matrix) from a quadrilateral given the coordinates and the material info 
 */
void	quadElasticLumpedMassMatrix(double* M, double* coords, MaterialInfo minfo)
{
	double area = quadArea(coords);
	double mass = area*minfo.density/4;
	int i;
	for ( i=0; i<8; i++)
		M[i+i*8] += mass;
}

/*!
 compute the consistent mass matrix (http://www.softeng.rl.ac.uk/st/projects/felib4/Docs/html/Intro/intro-node55.html) from a quadrilateral given the coordinates and the material info
 */
void	quadElasticConsistentMassMatrix(double* M, double* coords, MaterialInfo minfo)
{
	double area = quadArea(coords);
	double mass = area*minfo.density/36.0;
	int i;
	for ( i=0; i<8; i++){
		M[((i+0)%8+i*8)] = mass*4.0;
		M[((i+2)%8+i*8)] = mass*2.0;
		M[((i+4)%8+i*8)] = mass*1.0;
		M[((i+6)%8+i*8)] = mass*2.0;
	}
	
#ifdef DEBUG_PRINT
	printf("consistent\n");
	
	for ( i=0; i<8; i++){
		for (int j=0; j<8; j++){
			printf("M[%d][%d]: %lf\n", i,j, M[i+j*8]);		
		}
	}	
#endif
	
}

#pragma mark ___Tri___
/*!
 compute the area of a triangle from the list of coordinates
 */
inline double triArea(double* coords) //assuming all same z ccords
{
	double x1,x2,x3,y1,y2,y3,area;
	
	x1 = coords[0]; y1 = coords[3];
	x2 = coords[1]; y2 = coords[4];
	x3 = coords[2]; y3 = coords[5];
	
	area = 0.5*((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1));
	
#ifdef DEBUG_PRINT
	printf("area: %lf\n", area);
#endif	

	return dabs(area);
}

/*!
 compute the center of a triangle from the list of coordinates and return them in x and y
 */
inline void getTriCenter(double* coords, double &x, double &y)
{
	double x1,x2,x3,y1,y2,y3;
	x1 = coords[0]; y1 = coords[3];
	x2 = coords[1]; y2 = coords[4];
	x3 = coords[2]; y3 = coords[5];
  
    x=(x1+x2+x3)/3.0;                      // x-centroid of triangle
    y=(y1+y2+y3)/3.0;                      // y-centroid of triangle
}

#pragma mark ___General___
/*!
 compute the Jacobian matrix for a 2d isoparametric element
 */
inline void twoDJacobian(int numnodes, double* J, double* invJ, double &detJ, double* dNdSdN, double* x)
 {
	int i;
	int dim = 2;
	
	memset(J,0,(dim*dim)*sizeof(double));	

	for ( i=0; i<numnodes; i++){
		J[0] += x[i         ]*dNdSdN[i];			
		J[1] += x[i         ]*dNdSdN[numnodes+i];	
		J[2] += x[numnodes+i]*dNdSdN[i];	
		J[3] += x[numnodes+i]*dNdSdN[numnodes+i];			
	}

	detJ = m2_det(J);

	m2_inverse(invJ,J);
		
#ifdef DEBUG_PRINT
	int j;
	for ( i=0; i<dim; i++)
		for ( j=0; j<dim; j++)
			printf("J[%ld][%ld]: %llf\n",i,j, J[j*dim+i]);
		
	for ( i=0; i< dim; i++)
		for ( j=0; j<dim; j++)
			printf("invJ[%ld][%ld]: %llf\n",i,j, invJ[j*dim+i]);
		
	printf("detJ: %lf\n", detJ);
#endif
 }
 
/*!
 compute the global derivitive matrix for a 2d isoparametric element
 */
inline void twoDPhysicalDerivs(int numnodes, double* dNdXdY, double* invJ, double* dNdSdN)
{	
	int i;
	for ( i=0; i<numnodes; i++){	
		dNdXdY[i         ] = dNdSdN[i   ]*invJ[0] + dNdSdN[numnodes+i]*invJ[2]; //dNi / dx		
		dNdXdY[numnodes+i] = dNdSdN[i   ]*invJ[1] + dNdSdN[numnodes+i]*invJ[3]; //dNi / dy
	}
	
#ifdef DEBUG_PRINT
	//printf("\n");
	int dim = 2;
	int j;
	for ( i=0; i< dim; i++){
		for ( j=0; j<numnodes; j++){
			printf("dNdXdY[%ld][%ld]: %llf\n",j,i, dNdXdY[i*numnodes+j]);
		}
		//printf("\n");
	}
#endif

}

/*!
  compute the kinematic equation for a 2d isoparametric element
 */
/inline void	twoDConstMatrix(int numnodes, double* B, double* dNdXdY)
{	
	int sizeB1 = 3;
	int dim = 2;
	int sizeBB = sizeB1*dim;
	
	int i,j;
	
	for ( i=0; i< numnodes; i++){	
	
		j = i*sizeBB; 
		 
		B[j    ]	= dNdXdY[i         ];		
		B[j + 4]	= dNdXdY[numnodes+i];
		
		B[j + 2]	= dNdXdY[numnodes+i];
		B[j + 5]	= dNdXdY[i         ];			
	}
	
#ifdef DEBUG_PRINT	
	for ( i=0; i<sizeB1; i++)
		for ( j=0; j<numnodes*dim; j++)
			printf("B[%ld][%ld]: %lf\n",i,j, B[j*sizeB1+i]);
#endif			
}

#pragma mark ___AxisSym___

/*!
 compute the kinematic equation for an even axis-symmetric radial element
 */
void	axisSymConstMatrix(int numnodes, double* B, double* dNdRdZ, double* N, double r, double theta, int n)
{
	int sizeB1 = 6;
	int dim =3;
	int sizeBB = sizeB1*dim;
	
	int i,j;
	
	double cosnt = cos(n*theta);
	double sint = sin(n*theta);
	
	for ( i=0; i< numnodes; i++){	
	
		j = i*sizeBB; 
		 
		B[j    ]	= dNdRdZ[i]*cosnt;	
		//	
		B[j + 2]	= N[i]*cosnt/r;
		B[j + 3]	= dNdRdZ[numnodes+i]*cosnt;
		//
		B[j + 5]	= -n*N[i]*sint/r;
		//
		//
		B[j + 8]	= n*N[i]*cosnt/r;
		//
		B[j + 10]	= dNdRdZ[numnodes+i]*sint;	
		B[j + 11]	= (dNdRdZ[i] - N[i]/r)*sint;
		//
		B[j + 13]	= dNdRdZ[numnodes+i]*cosnt;		
		//		
		B[j + 15]	= dNdRdZ[i]*cosnt;		
		B[j + 16]	= -n*N[i]*sint/r;		
		//
	}
	
#ifdef DEBUG_PRINT	
	for ( i=0; i<sizeB1; i++){
		for ( j=0; j<numnodes*dim; j++){
			printf("B[%ld][%ld]: %lf\n",i,j, B[j*sizeB1+i]);
		}
		printf("\n");
	}
#endif			
}

/*!
 compute the kinematic equation for an odd axis-symmetric radial element
 */
void	antiAxisSymConstMatrix(int numnodes, double* B, double* dNdRdZ, double* N, double r, double theta, int n)
{
	int sizeB1 = 6;
	int dim =3;
	int sizeBB = sizeB1*dim;
	
	int i,j;
	
	double cosnt =cos(n*theta);
	double sint = sin(n*theta);
	
	for ( i=0; i< numnodes; i++){	
	
		j = i*sizeBB; 
		 
		B[j    ]	= dNdRdZ[i]*sint;	
		//	
		B[j + 2]	= N[i]*sint/r;
		B[j + 3]	= dNdRdZ[numnodes+i]*sint;
		//
		B[j + 5]	= n*N[i]*cosnt/r;
		//
		//
		B[j + 8]	= n*N[i]*sint/r;
		//
		B[j + 10]	= -dNdRdZ[numnodes+i]*cosnt;	
		B[j + 11]	= (-dNdRdZ[i] + N[i]/r)*cosnt;
		//
		B[j + 13]	= dNdRdZ[numnodes+i]*sint;		
		//		
		B[j + 15]	= dNdRdZ[i]*sint;		
		B[j + 16]	= n*N[i]*cosnt/r;		
		//
	}
	
#ifdef DEBUG_PRINT	
	for ( i=0; i<sizeB1; i++)
		for ( j=0; j<numnodes*dim; j++)
			printf("B[%ld][%ld]: %lf\n",i,j, B[j*sizeB1+i]);
#endif			
}



#endif
