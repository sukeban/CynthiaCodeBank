/*
 *  AcceleratedResonator.cpp
 *  ResonatorTestApp
 *
 *  Created by Cynthia Bruyns on 10/14/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "AcceleratedResonator.h"

#include <stdlib.h>
#include <vecLib/cblas.h>



void	AcceleratedResonator::Initialize()
{
	if (location == -1) return;

	double itwopi = 1.0/(2.0*M_PI);	
	
	int i;
	t = 0.0;	
	dt = 1.0/samplerate;
		
	int dof = location*eigendata->numDOF + 2; // cbnote need a more general force here
//	force[dof] = thumpforce;
	
	// multiply V^T by the force matrix to get g
	double* g = (double*) calloc(num, sizeof(double));		
	//http://www.netlib.org/blas/dgemv.f	
	
//	cblas_dgemv(	CblasColMajor, CblasTrans,									
//					eigendata->numDOFS, eigendata->numFreq,	1.0,						
//					eigendata->EigenVectors, 
//					eigendata->numFreq, 
//					force,1,
//					0.0,
//					g,1); 

	for (i=0; i<num; i++){
		g[i] = thumpforce*eigendata->EigenVectors[dof + i*eigendata->numDOFS];	
	}
			
	double zeta2;
	double lambda;
	
	for (i=0; i< num; i++){
		// wn = sqrt(lambdai)/2*PI
		lambda = eigendata->EigenValues[i]*fscale;
		if (isnan(lambda) || lambda < 0.0)
			lambda = 0.0;
		//printf("lambda[%d]: %lf\n", i, lambda);
			
		wn[i] = sqrtf(lambda)*itwopi;
		printf("wn[%d]: %lf\n", i, wn[i]);
		
		// zeta = (alpha1*wn + alpha2);
		zeta[i] = (alpha1*wn[i] + alpha2)/(2*wn[i]);	
		//printf("b_real[%d]: %lf\n", i, zeta[i]);
		zeta2 = powf(zeta[i],2.0);
		
		if (zeta2 == 1.0){ // if critically damped
			// ep = exp(-wn*dt)
			ep1[i].mReal = ep2[i].mReal = expf(-wn[i]*dt);	
			//printf("ep1[%d]: %lf\n", i, ep1[i].mReal);
			//printf("critically damped!\n");
		}
		else {							// not critically damped
			if (zeta2 < 1.0){
				// w = sqrt( 1.0 - zeta^2)
				w[i].mImag = sqrt(1.0-zeta2);
				//printf("w_imag[%d]: %lf\n", i, w[i].mImag);								
			}
			else{
				 // w = sqrt(zeta^2 - 1.0)
				w[i].mReal = sqrt(zeta2-1.0);
				//printf("w_real[%d]: %lf\n", i, w[i].mReal);				
			}
			
			//w1 = wn*(-zeta/2 + w)
			w1[i].mReal = -wn[i]*zeta[i] + wn[i]*w[i].mReal;
			w1[i].mImag = wn[i]*w[i].mImag;
			//printf("w1_real[%d]: %lf\n", i, w1[i].mReal);
			//printf("w1_imag[%d]: %lf\n", i, w1[i].mImag);
			
			//w2 = wn*(-zeta/2 - w)
			w2[i].mReal = -wn[i]*zeta[i] - wn[i]*w[i].mReal;
			w2[i].mImag = -wn[i]*w[i].mImag;
			//printf("w2_real[%d]: %lf\n", i, w2[i].mReal);
			//printf("w2_imag[%d]: %lf\n", i, w2[i].mImag);
			
			//ep1 = exp(real(w1)*dt)*exp(imag(w1))*dt)			
			//ep1 = exp(real(w1)*dt)*(cos(imag(w1*dt)+i sin(imag(w1)*dt))
			ep1[i].mReal = expf(w1[i].mReal*dt)*cosf(w1[i].mImag*dt);
			ep1[i].mImag = expf(w1[i].mReal*dt)*sinf(w1[i].mImag*dt);
			//printf("ep1_real[%d]: %lf\n", i, ep1[i].mReal);
			//printf("ep1_imag[%d]: %lf\n", i, ep1[i].mImag);
			
			//ep2 = exp(real(w2)*dt)*(cos(imag(w2*dt)+i sin(imag(w2)*dt))
			ep2[i].mReal = expf(w2[i].mReal*dt)*cos(w2[i].mImag*dt);
			ep2[i].mImag = expf(w2[i].mReal*dt)*sin(w2[i].mImag*dt);
			//printf("ep2_real[%d]: %lf\n", i, ep2[i].mReal);
			//printf("ep2_imag[%d]: %lf\n", i, ep2[i].mImag);
				
			//printf("normally damped!\n");
	
		} // not critically damped
	
		// c1 = z0/2 + zeta*z0/2w + zdot0/2w
		// the force applied to the veloity  of c's
		c1[i].mReal = 0.0;
		c1[i].mImag = g[i]/(w2[i].mImag-w1[i].mImag);
		// c2 = z0/2 - zeta*z0/2w + zdot0/2w		
		c2[i].mReal = 0.0;
		c2[i].mImag = g[i]/(w1[i].mImag-w2[i].mImag);
		// cbnote if it was displacement then would get a different result
		
		//printf("c1_real[%d]: %lf\n", i, c1[i].mReal);
		//printf("c1_imag[%d]: %lf\n", i, c1[i].mImag);
		//printf("c2_real[%d]: %lf\n", i, c2[i].mReal);
		//printf("c2_imag[%d]: %lf\n", i, c2[i].mImag);
		
		z1[i] = c1[i];
		z2[i] = c2[i];		
		//printf("z1_real[%d]: %lf\n", i, z1[i].mReal);
		//printf("z1_imag[%d]: %lf\n", i, z1[i].mImag);
		//printf("z2_real[%d]: %lf\n", i, z2[i].mReal);
		//printf("z2_imag[%d]: %lf\n", i, z2[i].mImag);
		
	}	// num
	ready = true;
	free(g);
	
	//printf("\n\n");
}	// function



// when you thump you make a force, and set all the parameters etc.

void AcceleratedResonator::Render(UInt32 inNumberFrames, AudioBufferList& inBufferList ) 
{ 
	if (!eigendata || !ready) return; 
	
	// cbnote not sure if it is getting to render before all of the data is ready or if an old render process is still around	

	int numChans = inBufferList.mNumberBuffers;
	float *left, *right;
	left = (float*)inBufferList.mBuffers[0].mData;
	right = numChans == 2 ? (float*)inBufferList.mBuffers[1].mData : 0;
						
	int n = num;
	if (nummodes < n)
		n = nummodes;

	if (n < 1) return;	
		
#ifdef DEBUG_PRINT
	printf("using num modes: %d\n",n);			
#endif
	float ss = sscale; // /(float)n;

#ifdef DEBUG_PRINT
	printf("soundscale: %lf\n",ss);
#endif

	double z;
	Complex s;
	for (UInt32 index=0; index < inNumberFrames; index++){
		z = 0.0;
		for (int res = 0; res < n; res++){ 
		
			if (wn[res] <  MIN_AUDIBLE_FREQ ){
				continue;
			}
			else if (wn[res] > MAX_AUDIBLE_FREQ)
				break;
					
			if (powf(zeta[res],2.0) != 1.0){	
				// z = real(x1 + x2)+ imag(x1 + x2);  
				s = z2[res] + z1[res]; 
				z += (s.mReal + s.mImag);				
			}
			else {
				// z = (z1 + t*z2)
				z += (z1[res].mReal + t*z2[res].mReal);
			}
			//z1 *= e1p
			z1[res] = z1[res]*ep1[res];						            
			//printf("z1[%d]: %lf + %lfi\n", res, z1[res].mReal, z1[res].mImag);
					
			//z2 *= e2p
			z2[res] = z2[res]*ep2[res];	
			//printf("z2[%d]: %lf + %lfi\n", res, z2[res].mReal, z2[res].mImag);
			
		}
		
		// increase t
		t += dt;
			
		//printf("%lf\n", z);
		// cbnote u = V z----------------------------------------------------------------------------------------
		
		if (left) left[index] = ss*z;	
		if (right) right[index] = left[index];
		
#ifdef DEBUG_PRINT
	if (left[index] != 0.0)
		printf("left[%d]: %lf\n", index,left[index]);
#endif

	}

}

void AcceleratedResonator::Thump(float velocity, int lc, float f)
{	
	if (!eigendata) return;
	
	location = lc;
	thumpforce = velocity*f; // base force that the velocity scales

	Initialize();
}

