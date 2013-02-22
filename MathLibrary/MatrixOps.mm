/*
 *  MatrixOps.mm
 *  ShapeSliders
 *
 *  Created by Cynthia Bruyns on 9/30/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "MatrixOps.h"

#include <math.h>

#pragma mark -
#pragma mark Double

#pragma mark ___2Vector___

double v3_squared(double* vec)
{
	return vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
}

double v3_length(double* vec) 
{
	return sqrt(v3_squared(vec));
}

void v3_cross(double* u, double* v, double* w)
{
	w[0] = u[1]*v[2] - u[2]*v[1];
	w[1] = u[2]*v[0] - u[0]*v[2];		
	w[2] = u[0]*v[1] - u[1]*v[0];
}

double v3_dot(double* u, double* v)
{
	return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
}

void v3_normalize(double*u)
{
	int i;
	double length = v3_length(u);
	for (i=0; i<3; i++)
		u[i]/=length;
}

void v3_sub(double* u, double* v, double* w)
{
	w[0] = u[0]-v[0];
	w[1] = u[1]-v[1];
	w[2] = u[2]-v[2];
}


#pragma mark ___VectorMatrixOps___

void v3m3_mult(double* vin, double* M, double* vout)
{
	vout[0] = vin[0]*M[0] + vin[1]*M[1] + vin[2]*M[2];
    vout[1] = vin[0]*M[4] + vin[1]*M[5] + vin[2]*M[6];
    vout[2] = vin[0]*M[8] + vin[1]*M[9] + vin[2]*M[10];
	
	//printf(	"vin: %lf %lf %lf\n vout: %lf %lf %lf\n", 
	//vin[0],vin[1],vin[2],
	//vout[0],vout[1],vout[2]);
	
	//http://www.gamedev.net/reference/articles/article695.asp
}


#pragma mark ___2Matrix___

// mat is a 2x2 matrix
// [ 0  2 ]
// [ 1  3 ]
//http://www.cvl.iis.u-tokyo.ac.jp/~miyazaki/tech/teche23.html
double m2_det(double* in)
{
	return in[0]*in[3] - in[1]*in[2]; 
}

// mat is a 2x2 matrix
// [ 0  2 ]
// [ 1  3 ]
int m2_inverse(double* out, double* in)
{
	double det = m2_det(in);
	if (det == 0.0) return 0;

	out[0] =  in[3]/det;
	out[1] = -in[1]/det;
	out[2] = -in[2]/det;
	out[3] =  in[0]/det;
	return (1);
}

#pragma mark ___3Matrix___
//http://mathworld.wolfram.com/images/equations/Determinant/equation1.gif
//http://mathworld.wolfram.com/Determinant.html
double m3_det( double* in )
{		
	double a1 = in[0];  double a2 = in[3];  double a3 = in[6];
	double b1 = in[1];  double b2 = in[4];  double b3 = in[7];
	double c1 = in[2];  double c2 = in[5];  double c3 = in[8];
	 	 	 
	double det = 	a1*b2*c3 + a2*b3*c1 + a3*b1*c2 -
					a1*b3*c2 - a2*b1*c3 - a3*b2*c1;
	
	return det;
}

int m3_inverse( double* out, double* in )
 {
	 double det = m3_det( in );
	 if (det == 0.0) return 0;

	 double a11 = in[0];
	 double a21 = in[1];
	 double a31 = in[2];
	 
	 double a12 = in[3];
	 double a22 = in[4];
	 double a32 = in[5];
	 
	 double a13 = in[6];
	 double a23 = in[7];
	 double a33 = in[8];
	
	 double invd = 1.0/det;
	 
	 out[0] =  ( a22*a33 - a23*a32 ) *invd;
	 out[1] =  ( a23*a31 - a21*a33 ) *invd;
	 out[2] =  ( a21*a32 - a22*a31 ) *invd;
	 
	 out[3] =  ( a13*a32 - a12*a33 ) *invd;
	 out[4] =  ( a11*a33 - a13*a31 ) *invd;
	 out[5] =  ( a12*a31 - a11*a32 ) *invd;
	 
	 out[6] =  ( a12*a23 - a13*a22 ) *invd;
	 out[7] =  ( a13*a21 - a11*a23 ) *invd;
	 out[8] =  ( a11*a22 - a12*a21 ) *invd;

	 return(1);
 }
	  
	  
void m3_mult(double* out, double* in, double* a)
{
	int i,j,k;
    for ( i=0;i<3;i++) {
        for ( j=0;j<3;j++) {
            out[i+j*3] = 0;
            for ( k=0;k<3;k++)
                out[i+j*3] += in[i+k*3]*a[k+j*3];
        }
    }
 }
	  
	  
#pragma mark ___4Matrix___
double m4_det( double* in )
{
	double m00 = in[0];
	double m10 = in[1];
	double m20 = in[2];
	double m30 = in[3];
	
	double m01 = in[4];
	double m11 = in[5];
	double m21 = in[6];
	double m31 = in[7];

	double m02 = in[8];
	double m12 = in[9];
	double m22 = in[10];
	double m32 = in[11];

	double m03 = in[12];
	double m13 = in[13];
	double m23 = in[14];
	double m33 = in[15];
		
	double det =
	m03 * m12 * m21 * m30-m02 * m13 * m21 * m30-m03 * m11 * m22 * m30+m01 * m13    * m22 * m30+
	m02 * m11 * m23 * m30-m01 * m12 * m23 * m30-m03 * m12 * m20 * m31+m02 * m13    * m20 * m31+
	m03 * m10 * m22 * m31-m00 * m13 * m22 * m31-m02 * m10 * m23 * m31+m00 * m12    * m23 * m31+
	m03 * m11 * m20 * m32-m01 * m13 * m20 * m32-m03 * m10 * m21 * m32+m00 * m13    * m21 * m32+
	m01 * m10 * m23 * m32-m00 * m11 * m23 * m32-m02 * m11 * m20 * m33+m01 * m12    * m20 * m33+
	m02 * m10 * m21 * m33-m00 * m12 * m21 * m33-m01 * m10 * m22 * m33+m00 * m11    * m22 * m33;

	return det;
}

	
int m4_inverse(double* out, double* in ) 
{
    double det = m4_det( in );	
	if (det == 0.0) return 0;
	
	double invd = 1.0/det;
	
	double m00 = in[0];
	double m10 = in[1];
	double m20 = in[2];
	double m30 = in[3];
	
	double m01 = in[4];
	double m11 = in[5];
	double m21 = in[6];
	double m31 = in[7];

	double m02 = in[8];
	double m12 = in[9];
	double m22 = in[10];
	double m32 = in[11];

	double m03 = in[12];
	double m13 = in[13];
	double m23 = in[14];
	double m33 = in[15];
	
	double o00 = m12*m23*m31 - m13*m22*m31 + m13*m21*m32 - m11*m23*m32 - m12*m21*m33 + m11*m22*m33;
	double o01 = m03*m22*m31 - m02*m23*m31 - m03*m21*m32 + m01*m23*m32 + m02*m21*m33 - m01*m22*m33;
	double o02 = m02*m13*m31 - m03*m12*m31 + m03*m11*m32 - m01*m13*m32 - m02*m11*m33 + m01*m12*m33;
	double o03 = m03*m12*m21 - m02*m13*m21 - m03*m11*m22 + m01*m13*m22 + m02*m11*m23 - m01*m12*m23;
	double o10 = m13*m22*m30 - m12*m23*m30 - m13*m20*m32 + m10*m23*m32 + m12*m20*m33 - m10*m22*m33;
	double o11 = m02*m23*m30 - m03*m22*m30 + m03*m20*m32 - m00*m23*m32 - m02*m20*m33 + m00*m22*m33;
	double o12 = m03*m12*m30 - m02*m13*m30 - m03*m10*m32 + m00*m13*m32 + m02*m10*m33 - m00*m12*m33;
	double o13 = m02*m13*m20 - m03*m12*m20 + m03*m10*m22 - m00*m13*m22 - m02*m10*m23 + m00*m12*m23;
	double o20 = m11*m23*m30 - m13*m21*m30 + m13*m20*m31 - m10*m23*m31 - m11*m20*m33 + m10*m21*m33;
	double o21 = m03*m21*m30 - m01*m23*m30 - m03*m20*m31 + m00*m23*m31 + m01*m20*m33 - m00*m21*m33;
	double o22 = m01*m13*m30 - m03*m11*m30 + m03*m10*m31 - m00*m13*m31 - m01*m10*m33 + m00*m11*m33;
	double o23 = m03*m11*m20 - m01*m13*m20 - m03*m10*m21 + m00*m13*m21 + m01*m10*m23 - m00*m11*m23;
	double o30 = m12*m21*m30 - m11*m22*m30 - m12*m20*m31 + m10*m22*m31 + m11*m20*m32 - m10*m21*m32;
	double o31 = m01*m22*m30 - m02*m21*m30 + m02*m20*m31 - m00*m22*m31 - m01*m20*m32 + m00*m21*m32;
	double o32 = m02*m11*m30 - m01*m12*m30 - m02*m10*m31 + m00*m12*m31 + m01*m10*m32 - m00*m11*m32;
	double o33 = m01*m12*m20 - m02*m11*m20 + m02*m10*m21 - m00*m12*m21 - m01*m10*m22 + m00*m11*m22;

	out[0] = invd*o00; out[4] = invd*o01; out[8]  = invd*o02; out[12] = invd*o03;
	out[1] = invd*o10; out[5] = invd*o11; out[9]  = invd*o12; out[13] = invd*o13;
	out[2] = invd*o20; out[6] = invd*o21; out[10] = invd*o22; out[14] = invd*o23;
	out[3] = invd*o30; out[7] = invd*o31; out[11] = invd*o32; out[15] = invd*o33;
	
	return 1;
}

double m4_sum(double* in)
{
	return (in[0]+in[1]+in[2]+in[3]);
}


#pragma mark -
#pragma mark Single



float v3_squared_f(float* vec)
{
	return vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
}

float v3_length_f(float* vec) 
{
	return sqrtf(v3_squared_f(vec));
}

void v3_cross_f(float* u, float* v, float* w)
{
	w[0] = u[1]*v[2] - u[2]*v[1];
	w[1] = u[2]*v[0] - u[0]*v[2];		
	w[2] = u[0]*v[1] - u[1]*v[0];
}

float v3_dot_f(float* u, float* v)
{
	return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
}

void v3_normalize_f(float*u)
{
	int i;
	float length = v3_length_f(u);
	for (i=0; i<3; i++)
		u[i]/=length;
}

void v3_sub_f(float* u, float* v, float* w)
{
	w[0] = u[0]-v[0];
	w[1] = u[1]-v[1];
	w[2] = u[2]-v[2];
}


#pragma mark ___VectorMatrixOps___

void v3m3_mult_f(float* vin, float* M, float* vout)
{
	vout[0] = vin[0]*M[0] + vin[1]*M[1] + vin[2]*M[2];
    vout[1] = vin[0]*M[4] + vin[1]*M[5] + vin[2]*M[6];
    vout[2] = vin[0]*M[8] + vin[1]*M[9] + vin[2]*M[10];
	
	//printf_f(	"vin: %lf %lf %lf\n vout: %lf %lf %lf\n", 
	//vin[0],vin[1],vin[2],
	//vout[0],vout[1],vout[2]);
	
	//http://www.gamedev.net/reference/articles/article695.asp
}


#pragma mark ___2Matrix___

// mat is a 2x2 matrix
// [ 0  2 ]
// [ 1  3 ]
//http://www.cvl.iis.u-tokyo.ac.jp/~miyazaki/tech/teche23.html
float m2_det_f(float* in)
{
	return in[0]*in[3] - in[1]*in[2]; 
}

// mat is a 2x2 matrix
// [ 0  2 ]
// [ 1  3 ]
int m2_inverse_f(float* out, float* in)
{
	float det = m2_det_f(in);
	if (det == 0.0) return 0;

	out[0] =  in[3]/det;
	out[1] = -in[1]/det;
	out[2] = -in[2]/det;
	out[3] =  in[0]/det;
	return (1);
}

#pragma mark ___3Matrix___

float m3_det_f( float* in )
{		
	 float a11 = in[0];
	 float a21 = in[1];
	 float a31 = in[2];
	 
	 float a12 = in[3];
	 float a22 = in[4];
	 float a32 = in[5];
	 
	 float a13 = in[6];
	 float a23 = in[7];
	 float a33 = in[8];

	float det = 	a11*a22*a33 + a21*a32*a13 +	a31*a12*a23 -
					a11*a32*a23 - a31*a22*a13 - a21*a12*a33;
	
	return det;
}

int m3_inverse_f( float* out, float* in )
 {
	 float det = m3_det_f( in );
	 if (det == 0.0) return 0;

	 float a11 = in[0];
	 float a21 = in[1];
	 float a31 = in[2];
	 
	 float a12 = in[3];
	 float a22 = in[4];
	 float a32 = in[5];
	 
	 float a13 = in[6];
	 float a23 = in[7];
	 float a33 = in[8];
	
	 float invd = 1.0/det;
	 
	 out[0] =  ( a22*a33 - a23*a32 ) *invd;
	 out[1] =  ( a23*a31 - a21*a33 ) *invd;
	 out[2] =  ( a21*a32 - a22*a31 ) *invd;
	 
	 out[3] =  ( a13*a32 - a12*a33 ) *invd;
	 out[4] =  ( a11*a33 - a13*a31 ) *invd;
	 out[5] =  ( a12*a31 - a11*a32 ) *invd;
	 
	 out[6] =  ( a12*a23 - a13*a22 ) *invd;
	 out[7] =  ( a13*a21 - a11*a23 ) *invd;
	 out[8] =  ( a11*a22 - a12*a21 ) *invd;
		   

	 return(1);
 }
	  
	  
void m3_mult_f(float* out, float* in, float* a)
{
	int i,j,k;
    for ( i=0;i<3;i++) {
        for ( j=0;j<3;j++) {
            out[i+j*3] = 0;
            for ( k=0;k<3;k++)
                out[i+j*3] += in[i+k*3]*a[k+j*3];
        }
    }
 }
	  
	  
#pragma mark ___4Matrix___
float m4_det_f( float* in )
{
	float a11 = in[0];
	float a21 = in[1];
	float a31 = in[2];
	float a41 = in[3];
	
	float a12 = in[4];
	float a22 = in[5];
	float a32 = in[6];
	float a42 = in[7];

	float a13 = in[8];
	float a23 = in[9];
	float a33 = in[10];
	float a43 = in[11];

	float a14 = in[12];
	float a24 = in[13];
	float a34 = in[14];
	float a44 = in[15];
	
	
	float det = 
	
	a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 + 
	
	a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41 + 
	
	a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 +
	 
	a14*a21*a33*a42 + a14*a22*a31*a43 + a14*a23*a32*a41 -
	
	a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42 - 
	
	a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 - 
	
	a13*a21*a34*a42 - a13*a22*a31*a44 - a13*a24*a32*a41 - 
	
	a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42;
	

    return det;
}

	
int m4_inverse_f(float* out, float* in ) 
{
    float det = m4_det_f( in );	
	if (det == 0.0) return 0;

	float a11 = in[0];
	float a21 = in[1];
	float a31 = in[2];
	float a41 = in[3];
	
	float a12 = in[4];
	float a22 = in[5];
	float a32 = in[6];
	float a42 = in[7];

	float a13 = in[8];
	float a23 = in[9];
	float a33 = in[10];
	float a43 = in[11];

	float a14 = in[12];
	float a24 = in[13];
	float a34 = in[14];
	float a44 = in[15];
	
	float invd = 1.0/det;

	float b11 = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42;
	float b12 = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43;
	float b13 = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42;
	float b14 = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33;
	
	float b21 = a21*a23*a43 + a23*a31*a44 + a24*a33*a41 - a21*a33*a44 - a23*a34*a41 - a24*a31*a43;
	float b22 = a11*a33*a44 + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 - a14*a33*a41;
	float b23 = a11*a24*a43 + a13*a21*a44 + a14*a23*a41 - a11*a23*a44 - a13*a24*a41 - a14*a21*a43;
	float b24 = a11*a23*a34 + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 - a14*a23*a31;
	
	float b31 = a21*a32*a44 + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 - a24*a32*a41;
	float b32 = a11*a34*a42 + a12*a31*a44 + a14*a32*a41 - a11*a32*a44 - a12*a34*a41 - a14*a31*a42;
	float b33 = a11*a22*a44 + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 - a14*a22*a41;
	float b34 = a11*a24*a32 + a12*a21*a34 + a14*a22*a31 - a11*a22*a34 - a12*a24*a31 - a14*a21*a32;
	
	float b41 = a21*a33*a42 + a22*a31*a43 + a23*a32*a41 - a21*a32*a43 - a22*a33*a41 - a23*a31*a42;
	float b42 = a11*a32*a43 + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 - a13*a32*a41;
	float b43 = a11*a23*a42 + a12*a21*a43 + a13*a22*a41 - a11*a22*a43 - a12*a23*a41 - a13*a21*a42;
	float b44 = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 - a13*a22*a31; 
	
	out[0] = invd*b11; out[4] = invd*b12; out[8]  = invd*b13; out[12] = invd*b14;
	out[1] = invd*b21; out[5] = invd*b22; out[9]  = invd*b23; out[13] = invd*b24;
	out[2] = invd*b31; out[6] = invd*b32; out[10] = invd*b33; out[14] = invd*b34;
	out[3] = invd*b41; out[7] = invd*b42; out[11] = invd*b43; out[15] = invd*b44;
	
	return 1;
}

float m4_sum_f(float* in)
{
	return (in[0]+in[1]+in[2]+in[3]);
}