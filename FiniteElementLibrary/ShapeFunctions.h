/*
 *  ShapeFunctions.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 1/7/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _ShapeFunction_Header_
#define _ShapeFunction_Header_
 
#include <math.h>
#include <string.h>

#include "2dElements.h"
#include "SharedFEData.h" 

#pragma mark ___GaussPoints___

/*!
 compute the guass points for a quadrilateral element given an n point quadrature rule
 */
inline void getGaussPoints(double* points, double* weights, int number)
 { //page 351 reddy
	switch(number){
	
	default:
	case 1:
		// 1.0  - point quadrature rule
		points[0] = 0.0;
		weights[0]= 2.0;
		break;
		
	case 2:
		// 2-point quadrature rule
		points[0] =-0.577350269189626;
		points[1] = 0.577350269189626;
		
		weights[0]= 1.0;
		weights[1]= 1.0;
		break;

	 case 3:
		//3-point quadrature rule
		points[0]=-0.774596669241483;
		points[1]= 0.0;
		points[2]= 0.774596669241483;
		
		weights[0]=0.555555555555556;
		weights[1]=0.888888888888889;
		weights[2]=0.555555555555556;
		break;

	case 4:
		//4-point quadrature rule
		points[0]=-0.861136311594053;
		points[1]=-0.339981043584856;
		points[2]= 0.339981043584856;
		points[3]= 0.861136311594053;
		
		weights[0]=0.347854845137454;
		weights[1]=0.652145154862546;
		weights[2]=0.652145154862546;
		weights[3]=0.347854845137454;
		break;
	 
	case 5:
		// 5-point quadrature rule
		points[0]=-0.906179845938664;
		points[1]=-0.538469310105683;
		points[2]= 0.0;
		points[3]= 0.538469310105683;
		points[4]= 0.906179845938664;
		
		weights[0]=0.236926885056189;
		weights[1]=0.478628670499366;
		weights[2]=0.568888888888889;
		weights[3]=0.478628670499366;
		weights[4]=0.236926885056189;
		break;
		
	case 6:
		points[0]=-0.932469514203512;
		points[1]=-0.661209386466265;
		points[2]=-0.238619186083197;
		points[3]= 0.238619186083197;
		points[4]= 0.661209386466265;
		points[5]= 0.932469514203512;

		weights[0]=0.17132449237917;
		weights[1]=0.360761573048139;
		weights[2]=0.467913934572691;
		weights[3]=0.467913934572691;
		weights[4]=0.360761573048139;
		weights[5]=0.17132449237917;	
		break;
	
	}
#ifdef DEBUG_PRINT
	int i;
	for ( i=0; i<number; i++)
		printf("point[%ld]: %lf weight[%ld]: %lf\n", i, points[i], i, weights[i]);
#endif	
}
 
 #pragma mark ___ShapeFunctions___

/*!
 return the shape functions for a line given a polynomial approximation order
 */
 inline void lineShapeFunctions(double* N, double* dNdSdN, double xi, PolynomialOrderEnum order)
 { // page 346 reddy
  
	switch (order){
		
		case cubic:
			N[0] = -0.5625*(1.0-xi  )*(0.1111  - pow(xi,2.0));
			N[1] =  1.6875*(1.0-pow(xi,2.0))*(0.3333-xi  );
			N[2] =  1.6875*(1.0-pow(xi,2.0))*(0.3333+xi  );
			N[3] = -0.5625*(1.0+xi  )*(0.1111  - pow(xi,2.0));
			
			dNdSdN[0]=   0.0625 - 0.5625*pow(xi,2.0) - 2.0*(-0.5625 + 0.5625*xi)*xi;
			dNdSdN[1]=           -3.3750*xi*(0.3333-xi) - 1.6875 + 1.6875*pow(xi,2.0);
			dNdSdN[2]=		     -3.3750*xi*(0.3333+xi) + 1.6875 - 1.6875*pow(xi,2.0);
			dNdSdN[3] = -0.0625 + 0.5625*pow(xi,2.0) - 2.0*(-0.5625 - 0.5625*xi)*xi;

		break;
		case quadratic:
			N[0] = -0.5*xi*(1.0-xi  );
			N[1] =   (1.0-pow(xi,2.0)); 
			N[2] =  0.5*xi*(1.0+xi  );
			
			dNdSdN[0]= -0.5+xi;
			dNdSdN[1]= -2.0*xi;
			dNdSdN[2]=  0.5+xi;
			
		break;
		case linear:
		default:
			N[0] = 0.5*(1.0-xi);
			N[1] = 0.5*(1.0+xi);
			
			dNdSdN[0]= -0.5;
			dNdSdN[1]=  0.5;
		break;
	}
#ifdef DEBUG_PRINT
			printf("xi: %lf \n", xi);
			int i;
			for ( i=0; i<2; i++)
				printf("N[%ld]: %lf dNdXdY[%ld]: %lf, dNdY[%ld]: %lf\n", i, N[i], i, dNdSdN[i], i+4, dNdSdN[i+4]);
#endif	
 }

/*!
 return the shape functions for a triangle (in barycentric coordinates) given a polynomial approximation order
 */
void triangularBaryShapeFunctions(double* N, double* dNdSdN, double xi, double eta, double* coords, PolynomialOrderEnum order)
 {
	// page 91 kwon (finite element method in matlab)
	// the shape functions are different in bhatti fundamentals  pp. 361
	// but the matricies come out the same	
	
	double area2 = 1.0/(2.0*triArea(coords));
	
	switch (order){
		case cubic:	
		break;		
		case quadratic:
		break;
		case linear:
		default:
			double x1,x2,x3,y1,y2,y3;
			x1 = coords[0]; y1 = coords[3];
			x2 = coords[1]; y2 = coords[4];
			x3 = coords[2]; y3 = coords[5];
					
			N[0]=area2*((x2*y3-x3*y2)+(y2-y3)*xi+(x3-x2)*eta);
			N[1]=area2*((x3*y1-x1*y3)+(y3-y1)*xi+(x1-x3)*eta);
			N[2]=area2*((x1*y2-x2*y1)+(y1-y2)*xi+(x2-x1)*eta);

			dNdSdN[0]=area2*(y2-y3);
			dNdSdN[1]=area2*(y3-y1);
			dNdSdN[2]=area2*(y1-y2);  
			
			dNdSdN[3]=area2*(x3-x2);  
			dNdSdN[4]=area2*(x1-x3);  
			dNdSdN[5]=area2*(x2-x1);  
				
#ifdef DEBUG_PRINT
			printf("xi: %lf eta: %lf\n", xi, eta);
			int i;
			for ( i=0; i<3; i++)
				printf("N[%ld]: %lf dNdXdY[%ld]: %lf, dNdY[%ld]: %lf\n", i, N[i], i, dNdSdN[i], i+3, dNdSdN[i+3]);
#endif				
		break;
	}
 }
 
/*!
 return the shape functions for a triangle (in isoparametric coordinates) given a polynomial approximation order
 */
void triangularShapeFunctions(double* N, double* dNdSdN, double xi, double eta, PolynomialOrderEnum order)
 {	
 
  // p 171 kwon fea in matlab, get fundamentals book bhatti

	switch (order){
		case cubic:	
		break;		
		case quadratic:
		N[0] = (1 - xi - eta)*(1 - 2*xi-2*eta);
		N[1] = xi*(2*xi-1);
		N[2] = eta*(2*eta - 1);
		N[3] = 4*xi*(1 - xi - eta);
		N[4] = 4*xi*eta;
		N[5] = 4*eta*(1 - xi - eta);
		
		// cbnote finish
		
		break;
		case linear:	
		default:
		
		N[0] = 1.0 - xi - eta;
		N[1] = xi;
		N[2] = eta;
		
		dNdSdN[0] =  -1.0; dNdSdN[3] = -1.0;
		dNdSdN[1] =   1.0; dNdSdN[4] = 0.0;
		dNdSdN[2] =   0.0; dNdSdN[5] = 1.0;
				
#ifdef DEBUG_PRINT
			printf("xi: %lf eta: %lf\n", xi, eta);
			int i;
			for ( i=0; i<3; i++)
				printf("N[%ld]: %lf dNdXdY[%ld]: %lf, dNdY[%ld]: %lf\n", i, N[i], i, dNdSdN[i], i+3, dNdSdN[i+3]);
#endif				
		break;
	}
 }

/*!
 return the shape functions for a quadrilateral element (in isoparametric coordinates) given a polynomial approximation order
 */
inline  void quadrilateralShapeFunctions(double* N, double* dNdSdN, double xi, double eta, PolynomialOrderEnum order)
 {// page 534 reddy
	switch (order){
		case cubic:	
		break;		
		case quadratic:
			N[0]= 0.25*(pow(xi,2.0)-xi)*(pow(eta,2.0)-eta);
			N[1]= 0.5*(1.0-pow(xi,2.0))*(pow(eta,2.0)-eta);
			N[2]= 0.25*(pow(xi,2.0)+xi)*(pow(eta,2.0)-eta);
			N[3]= 0.5*(pow(xi,2.0)-xi)*(1.0-pow(eta,2.0));
			N[4]=     (1.0-pow(xi,2.0))*(1.0-pow(eta,2.0));
			N[5]= 0.5*(pow(xi,2.0)+xi)*(1.0-pow(eta,2.0));
			N[6]= 0.25*(pow(xi,2.0)+xi)*(1.0-pow(eta,2.0));
			N[7]= 0.5*(1.0-pow(xi,2.0))*(pow(eta,2.0)+eta);
			N[8]= 0.25*(pow(xi,2.0)+xi)*(pow(eta,2.0)+eta);
			
			dNdSdN[0]= (0.5*xi - 0.25)*(pow(eta,2.0)-eta);		dNdSdN[9]= (0.25*pow(xi,2.0) - 0.25*xi)*(2.0*eta - 1.0);
			dNdSdN[1]= -xi*(pow(eta,2.0)-eta);					dNdSdN[10]= (0.5 - 0.5*pow(xi,2.0))*(2.0*eta-1.0);
			dNdSdN[2]= (0.5*xi + 0.25)*(pow(eta,2.0)-eta);		dNdSdN[11]= (0.25*pow(xi,2.0) + 0.25*xi)*(2.0*eta - 1.0);
			dNdSdN[3]= (-0.5 + xi)*(1.0 - pow(eta,2.0));		dNdSdN[12]= -2.0*(0.5*pow(xi,2.0) - 0.5*xi)*eta;
			dNdSdN[4]= -2.0*xi*(1.0 - pow(eta,2.0));			dNdSdN[13]= -2.0*(.0-pow(xi,2.0))*eta;
			dNdSdN[5]= (0.5 + xi)*(1.0 - pow(eta,2.0));			dNdSdN[14]= -2.0*(0.5*pow(xi,2.0) + 0.5*xi)*eta;
			dNdSdN[6]= (0.5*xi + 0.25)*(1.0 - pow(eta,2.0));	dNdSdN[15]= -2.0*(0.25*pow(xi,2.0) + 0.25*xi)*eta;
			dNdSdN[7]= -xi*(pow(eta,2.0)+eta);					dNdSdN[16]= (0.5 - 0.5*pow(xi,2.0))*(2.0*eta+1.0);
			dNdSdN[8]= (0.5*xi + 0.25)*(pow(eta,2.0)+eta);		dNdSdN[17]= (0.25*pow(xi,2.0) + 0.25*xi)*(2.0*eta+1.0);														
		break;
		case linear:
		default:
			N[0]= 0.25*(1.0-xi)*(1.0-eta);
			N[1]= 0.25*(1.0+xi)*(1.0-eta);
			N[2]= 0.25*(1.0+xi)*(1.0+eta);
			N[3]= 0.25*(1.0-xi)*(1.0+eta);
			
			dNdSdN[0]= -0.25 + 0.25*eta; dNdSdN[4]= -0.25 + 0.25*xi;
			dNdSdN[1]=  0.25 - 0.25*eta; dNdSdN[5]= -0.25 - 0.25*xi;
			dNdSdN[2]=  0.25 + 0.25*eta; dNdSdN[6]=  0.25 + 0.25*xi;
			dNdSdN[3]= -0.25 - 0.25*eta; dNdSdN[7]=  0.25 - 0.25*xi;
		
#ifdef DEBUG_PRINT
			printf("xi: %lf eta: %lf\n", xi, eta);
			int i;
			for ( i=0; i<4; i++)
				printf("N[%ld]: %lf dNdXdY[%ld]: %lf, dNdY[%ld]: %lf\n", i, N[i], i, dNdSdN[i], i+4, dNdSdN[i+4]);
#endif	
			
		break;
	}
 }
 
 
/*!
 return the shape functions for a brick element (in isoparametric coordinates) given a polynomial approximation order
 */
inline  void brickShapeFunctions(double* N, double* dNdSdN, double xi, double eta, double phi, PolynomialOrderEnum order)
 {
	switch (order){
		case cubic:			
		break;
		case quadratic:					
		break;
		case linear:
		default:		
			N[0]=0.125*(1.0  - xi)*(1.0  - eta)*(1.0  - phi);
			N[1]=0.125*(1.0  + xi)*(1.0  - eta)*(1.0  - phi);
			N[2]=0.125*(1.0  + xi)*(1.0  + eta)*(1.0  - phi);
			N[3]=0.125*(1.0  - xi)*(1.0  + eta)*(1.0  - phi);
			 
			N[4]=0.125*(1.0  - xi)*(1.0  - eta)*(1.0  + phi);
			N[5]=0.125*(1.0  + xi)*(1.0  - eta)*(1.0  + phi);
			N[6]=0.125*(1.0  + xi)*(1.0  + eta)*(1.0  + phi);
			N[7]=0.125*(1.0  - xi)*(1.0  + eta)*(1.0  + phi);
			
			
			dNdSdN[0]=-0.125*(1.0  - eta)*(1.0  - phi);
			dNdSdN[1]=0.125*(1.0  - eta)*(1.0  - phi);
			dNdSdN[2]=0.125*(1.0  + eta)*(1.0  - phi);
			dNdSdN[3]=-0.125*(1.0  + eta)*(1.0  - phi);
			dNdSdN[4]=-0.125*(1.0  - eta)*(1.0  + phi);
			dNdSdN[5]=0.125*(1.0  - eta)*(1.0  + phi);
			dNdSdN[6]=0.125*(1.0  + eta)*(1.0  + phi);
			dNdSdN[7]=-0.125*(1.0  + eta)*(1.0  + phi);


			dNdSdN[8]=-0.125*(1.0  - xi)*(1.0  - phi);
			dNdSdN[9]=-0.125*(1.0  + xi)*(1.0  - phi);
			dNdSdN[10]=0.125*(1.0  + xi)*(1.0  - phi);
			dNdSdN[11]=0.125*(1.0  - xi)*(1.0  - phi);
			dNdSdN[12]=-0.125*(1.0  - xi)*(1.0  + phi);
			dNdSdN[13]=-0.125*(1.0  + xi)*(1.0  + phi);
			dNdSdN[14]=0.125*(1.0  + xi)*(1.0  + phi);
			dNdSdN[15]=0.125*(1.0  - xi)*(1.0  + phi);

			dNdSdN[16]=-0.125*(1.0  - xi)*(1.0  - eta);
			dNdSdN[17]=-0.125*(1.0  + xi)*(1.0  - eta);
			dNdSdN[18]=-0.125*(1.0  + xi)*(1.0  + eta);
			dNdSdN[19]=-0.125*(1.0  - xi)*(1.0  + eta);
			dNdSdN[20]=0.125*(1.0  - xi)*(1.0  - eta);
			dNdSdN[21]=0.125*(1.0  + xi)*(1.0  - eta);
			dNdSdN[22]=0.125*(1.0  + xi)*(1.0  + eta);
			dNdSdN[23]=0.125*(1.0  - xi)*(1.0  + eta);

		
#ifdef DEBUG_PRINT
			printf("xi: %lf eta: %lf phi: %lf\n", xi, eta,phi);
			int i;
			for ( i=0; i<8; i++)
				printf("N[%ld]: %lf dNdX[%ld]: %lf, dNdY[%ld]: %lf, dNdY[%ld]: %lf\n", i, N[i], i, dNdSdN[i], i+8, dNdSdN[i+8],i+16, dNdSdN[i+16]);
#endif	
			
		break;
	}
 }


#endif