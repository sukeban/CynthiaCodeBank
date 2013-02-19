/*
 *  AcceleratedConvolutionResonator.cpp
 *  PlateReverb
 *
 *  Created by Cynthia Bruyns on 2/28/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "AcceleratedConvolutionResonator.h"

#include <stdlib.h>


//#include <Accelerate/Accelerate.h>


#include <vecLib/cblas.h>
#include <vecLib/vDSP.h>
#include <vecLib/vDSP_translate.h>
#include <vecLib/vectorOps.h>


void AcceleratedConvolutionResonator::Thump(int lc)
{	
	if (!eigendata) return;

	location[0] = lc;
	int dof = location[0]*eigendata->numDOF + 2;	 //  need a more general force here 
	// gn = V^T e_location
	
	int i;
	for ( i=0; i< nummodespow2; i++){
		gn[i] = eigendata->EigenVectors[dof + i*eigendata->numDOFS]; // cbnote times thumpforce
	}	

#ifdef DEBUG_PRINT
	for ( i=0; i< nummodespow2; i++){
		printf("gn[%d]: %lf\n", i, gn[i]);
	}
#endif


}

void AcceleratedConvolutionResonator::PickUp(int pc)
{	
	if (!eigendata) return;

	pickup[0] = pc;	
	int dof = pickup[0]*eigendata->numDOF + 2;	 //  need a more general pick up here 
	// fn = L V
	int i;
	for ( i=0; i< nummodespow2; i++){
		fn[i] =  eigendata->EigenVectors[dof + i*eigendata->numDOFS];
	}	
	
	for (i=0; i<nummodespow2; i++){
		if (fn[i] != 0.0)
			fn[i] = 1.0;
	}	

#ifdef BEBUG_PRINT
	for (int i=0; i< nummodespow2; i++){
		printf("fn[%d]: %lf\n", i, fn[i]);
	}
#endif
	
	ready = true;
}

void	AcceleratedConvolutionResonator::Initialize()
{
	if (!eigendata) return;
	
	//printf("initialize\n");

	float twopi = 2.0*M_PI;
	float itwopi = 1.0/twopi;	
	int i;
		
	dt = 1.0/samplerate;
			
	float lambda;
	float p_real, p_imag;
	float exp_p_imag_real, exp_p_imag_imag;
	
	memset(an,0, nummodespow2*sizeof(float));
	
	if (nummodes > nummodespow2)
		nummodes = nummodespow2;

	for (i=0; i< nummodes; i++){ // always less that power of two number
	
		// wn = sqrtf(lambdai)/2*PI
		lambda = eigendata->EigenValues[i];		
		if (isnan(lambda) || lambda < 0.0)
			lambda = 0.0;
			
		wn[i] = sqrtf(lambda)*itwopi*fscale;	

		if (wn[i] <  MIN_AUDIBLE_FREQ || wn[i] > MAX_AUDIBLE_FREQ){
			continue;				
		}	
		
		f[i] = wn[i]*itwopi;
		phi[i] = (alpha1*wn[i] + alpha2)/(2*wn[i]);	
		
		an[i] = 1;
		
		d[i] = f[i]*M_PI*tanf(phi[i]);
		
		w_real[i] = wn[i];
		w_imag[i] = d[i];
				
		// http://scholar.hw.ac.uk/site/maths/topic8.html
		// a	b		c		d
		// (0 + 1i) (w_real + w_imag)
		//(a + ib)(c + id) = (ac - bd) + i(bc + ad)
		
		p_real = -w_imag[i]*dt;
		p_imag = w_real[i]*dt;
		
		//http://mathforum.org/library/drmath/view/52251.html

		exp_p_imag_real = cosf(p_imag);
		exp_p_imag_imag = sinf(p_imag);
		
		//(a + ib)(c + id) = (ac - bd) + i(bc + ad)
		//(exp(p_real),0)(	exp_p_imag_real, exp_p_imag_imag);
		// a		b			c				d
		c_real[i] = expf(p_real)*exp_p_imag_real;
		c_imag[i] = expf(p_real)*exp_p_imag_imag;
		
			
	}	// num
	
	//  scale an by gn to ge rid of nonparticipating nodes
	vDSP_vmul (an, 1, gn,1, an,1, nummodespow2);
	
	for (i=0; i<nummodespow2; i++){
		if (an[i] != 0.0)
			an[i] = 1.0;
	}
 
#ifdef BEBUG_PRINT
	for (i=0; i< nummodespow2; i++){
		printf("an[%d]: %lf\n", i, an[i]);
	}
#endif
		
#ifdef BEBUG_PRINT
	for (i=0; i<nummodes; i++){
		printf("wn[%d]: %lf\t f[%d]: %lf\t phi[%d]: %lf \t d[%d]: %lf \t an[%d]: %lf \t c_real[%d]: %lf \t c_imag[%d]: %lf\n", 
		i, wn[i], i, f[i], i, phi[i], i, d[i], i, an[i], i, c_real[i], i, c_imag[i]);	
	}
#endif



}	



inline void AcceleratedConvolutionResonator::CheckOutputBuffer( Float32* output, UInt32 inNumberOfFrames)
{	
	for (UInt32 i = 0; i < inNumberOfFrames; ++i) {
		if (output[i] > 1.0 || 
			output[i] == NAN ||
			output[i] < -1.0) {
				printf("frame[%ld]: %lf\n",i, output[i]);
		}
		if (i > 0){
			Float32 diff;
			diff = output[i] - output[i-1];
			if (fabs(diff) > 0.3){
				printf("frame[%ld]: %lf jumped\n",i, output[i]);

				for (UInt32 j=0; j< inNumberOfFrames; j++){
					printf("output[%ld]: %lf\n", j,output[j]);
				
				}
			}
		}
		
		/*if (output[i] != 0.0) 
			printf("%lf\n", output[i]);*/
	}
}


void AcceleratedConvolutionResonator::Process(	const Float32		*inSourceP,
														Float32		*inDestP,
														UInt32		inFramesToProcess,
														UInt32		inNumChannels, // for version 2 AudioUnits inNumChannels is always 1
														bool		&ioSilence )
{ 
	if (!eigendata || !ready) return; 

	const Float32 *sourceP = inSourceP;
	Float32 *destP = inDestP;
										
#ifdef DEBUG_PRINT
	printf("using num modes: %d\n", nummodes);			
#endif
	float ss = sscale;

#ifdef DEBUG_PRINT
	printf("soundscale: %lf\n",ss);
#endif

	int nummodes_l = nummodes;
	int nummodespow2_l = nummodespow2;
	
	float* u_l = u;
	
	
	float* fn_l = fn;

	
	float* deltan_l = deltan;

	float* znew_real_l = znew_real;
	float* znew_imag_l = znew_imag;
	
	float* zprev_real_l = zprev_real;
	float* zprev_imag_l = zprev_imag;
		
	float* c_real_l = c_real;
	float* c_imag_l = c_imag;	 

	float delta;	
	float u_sum;
	
	int index; 
	int j;
	for ( index=0; index < inFramesToProcess; index++){
	
		delta = *sourceP;
		sourceP++;		
		
#ifdef DEBUG_PRINT
		if (delta != 0.0)
			printf(" %lf\n",delta);
#endif

		// e = ac-bd
		// c(i)*y(i-1)
		//Vector multiply, multiply, and subtract
		vDSP_vmmsb (	c_real_l,
						1,
						zprev_real_l,
						1,
						c_imag,
						1,
						zprev_imag_l,
						1,
						znew_real_l,
						1,
						nummodespow2_l);
		
		// E = bc + ad
		// c(i)*y(i-1)
		//Vector multiply, multiply, and add
		vDSP_vmma (		c_imag_l,
						1,
						zprev_real_l,
						1,
						c_real,
						1,
						zprev_imag_l,
						1,
						znew_imag_l,
						1,
						nummodespow2_l);
	
		vDSP_vsmul (	an, 1, &delta, deltan_l, 1, nummodespow2_l);
 						
		vDSP_vadd (		znew_real_l,
						1,
						deltan_l, 
						1,
						znew_real_l, 
						1,
						nummodespow2_l);
				
		// not fully transfomed though, still need Phi[j] premultiplied
		// u[j] = sum_i V[j][i] z[i] only at pickups and in direction
 		vDSP_vmul (fn_l, 1, znew_real_l, 1, u_l,1, nummodespow2);
 						
		vDSP_sve (u_l, 1, &u_sum, nummodespow2_l); 
		
		// cbnote the system can be made simpler if i use the plate elements too
					
		(*destP) = (float) (drymix * delta) + (float) (ss*wetmix*u_sum);	
		destP++;
		
		vDSP_vswap (znew_real_l,1, zprev_real_l,1, nummodespow2_l);
		vDSP_vswap (znew_imag_l,1, zprev_imag_l,1, nummodespow2_l);
  		
#ifdef DEBUG_PRINT
		if (destP[index] != 0.0)
			printf("%lf\n", destP[index]);
#endif
	}

}


