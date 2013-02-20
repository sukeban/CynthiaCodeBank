#ifndef _MatrixOps_H_
#define _MatrixOps_H_

#pragma mark -
#pragma mark Double
#pragma mark ___2Vector___

double v3_squared(double* vec);
double v3_length(double* vec);
void v3_cross(double* u, double* v, double* w);
double v3_dot(double* u, double* v);
void v3_normalize(double*u);
void v3_sub(double* u, double* v, double* w);


#pragma mark ___VectorMatrixOps___

void v3m3_mult(double* vin, double* M, double* vout);

#pragma mark ___2Matrix___

double m2_det(double* in);
int m2_inverse(double* out, double* in);

#pragma mark ___3Matrix___

double m3_det( double* in );
int m3_inverse( double* out, double* in );	  	  
void m3_mult(double* out, double* in, double* a);	  
	  
#pragma mark ___4Matrix___
double m4_det( double* in );	
double m4_sum( double* in);

int m4_inverse(double* out, double* in );

	
#pragma mark -
#pragma mark Single
#pragma mark ___2Vector___

float v3_squared_f(float* vec);
float v3_length_f(float* vec);
void v3_cross_f(float* u, float* v, float* w);
float v3_dot_f(float* u, float* v);
void v3_normalize_f(float*u);
void v3_sub_f(float* u, float* v, float* w);


#pragma mark ___VectorMatrixOps___

void v3m3_mult_f(float* vin, float* M, float* vout);

#pragma mark ___2Matrix___

float m2_det_f(float* in);
int m2_inverse_f(float* out, float* in);

#pragma mark ___3Matrix___

float m3_det_f( float* in );
int m3_inverse_f( float* out, float* in );	  	  
void m3_mult_f(float* out, float* in, float* a);	  
	  
#pragma mark ___4Matrix___
float m4_det_f( float* in );	
float m4_sum_f( float* in);

int m4_inverse_f(float* out, float* in );


#endif