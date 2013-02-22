/*
 *  AcceleratedConvolutionResonator.h
 *  ResonatorTestApp
 *
 *  Created by Cynthia Bruyns on 10/14/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef AcceleratedConvolutionResonator_H
#define AcceleratedConvolutionResonator_H

#include <stdlib.h>
#import <AudioUnit/AudioUnit.h>
#import <AudioToolbox/AudioToolbox.h>

#include "AcceleratedConvolutionResonatorSharedData.h"

#include "ModelSharedData.h"
#include "ComplexNumber.h"

/*!
 a resonator class that contains methods for computing output values using convolution
 this class has capabilites for setting damping and scale parameters as well as wet/dry mix
 */
class AcceleratedConvolutionResonator {
public:

	AcceleratedConvolutionResonator()
	{
		eigendata = NULL;
		
		ready = false;
		
		u = NULL;
		
		wn = NULL;
		f = NULL;
		
		d = NULL;
		phi = NULL;
		
		fn = NULL;
		gn = NULL;
		
		an = NULL;	
		deltan = NULL;
		
		zprev_real = NULL;		zprev_imag = NULL;
		znew_real = NULL;		znew_imag = NULL;
		
		w_real = NULL;			w_imag = NULL;
		c_real = NULL;			c_imag = NULL;
		
		location[0] = -1;
		pickup[0] = -1;
	};
	
	
	void FreeData()
	{
		if (u) free(u);
	
		if (wn) free(wn);
		if (f) free(f);

		if (d) free(d);
		if (phi) free(phi);
		
		if (fn) free(fn);
		if (gn) free(gn);
	
		if (an) free(an);
		if (deltan) free(deltan);
			
		if (zprev_real) free(zprev_real);
		if (zprev_imag) free(zprev_imag);
		
		if (znew_real) free(znew_real);
		if (znew_imag) free(znew_imag);
		
		if (w_real) free(w_real);
		if (w_imag) free(w_imag);
		
		if (c_real) free(c_real);
		if (c_imag) free(c_imag);
		
		u = NULL;
		
		wn = NULL;
		f = NULL;
		
		d = NULL;
		phi = NULL;
		
		fn = NULL;
		gn = NULL;
		
		an = NULL;
		deltan = NULL;
		
		zprev_real = NULL; zprev_imag = NULL;
		znew_real = NULL; znew_imag = NULL;
		
		w_real = NULL; w_imag = NULL;
		c_real = NULL; c_imag = NULL; 
	};
		
	~AcceleratedConvolutionResonator(){	
		FreeData();
	};
		
	void				Initialize();
	
	void				CheckOutputBuffer( Float32* output, UInt32 inNumberOfFrames);

	// callled when you want to render
	 void				Process(	const Float32		*inSourceP,
                                          Float32		*inDestP,
                                          UInt32 		inFramesToProcess,
                                          UInt32		inNumChannels, // for version 2 AudioUnits inNumChannels is always 1
                                          bool			&ioSilence);

	// called from attack**********************************
	void				SetEigenData(EigenData* ed)
	{
						FreeData();

						eigendata = ed;	
					
						int floatsize = sizeof(float); 
												
						nummodes = eigendata->numFreq;
						nummodespow2 = nummodes;  // cbnote find next power of two, but then mult by the matrix would be off
						
						
						phi = (float*) calloc(nummodespow2, floatsize);		// zeta = (alpha1*wn + alpha2)/(2*wn);	
						d = (float*) calloc(nummodespow2, floatsize);		// wn = sqrt(lambdai)/2*PI	
						
						fn = (float*) calloc(nummodespow2, floatsize);		// the row of the tump dof of interest
						gn = (float*) calloc(nummodespow2, floatsize);		// the column of the pick dof of interest
	
						an = (float*) calloc(nummodespow2, floatsize);		// constant
						deltan = (float*) calloc(nummodespow2, floatsize);	// incoming

						wn = (float*) calloc(nummodespow2, floatsize);		// wn = sqrt(lambdai)/2*PI		
						f = (float*) calloc(nummodespow2, floatsize);		// f = wn/PI	

							
						zprev_real = (float*) calloc(nummodespow2, floatsize);	
						zprev_imag = (float*) calloc(nummodespow2, floatsize);	
						
						znew_real = (float*) calloc(nummodespow2, floatsize);	
						znew_imag = (float*) calloc(nummodespow2, floatsize);	
						
																
						w_real = (float*) calloc(nummodespow2, floatsize);	// w = sqrt(b^2 - 4 wn)/2
						w_imag = (float*) calloc(nummodespow2, floatsize);	// w = sqrt(b^2 - 4 wn)/2

						c_real = (float*) calloc(nummodespow2, floatsize);	// w1 = -b + w
						c_imag = (float*) calloc(nummodespow2, floatsize);	// w1 = -b + w
						
						
						u = (float*) calloc(nummodespow2, floatsize);							
								
	};
	
	void				SetSampleRate(float sr)
	{ 
						samplerate = sr;
						dt = 1.0/samplerate;
	};
	
	void				SetSoundScale(float ss)
	{
						sscale = ss; 
	};	
	
	void				SetMix(float mix)
	{
						wetmix = mix;
						drymix = sqrtf(1.0- wetmix);
//						printf(" sscale: %lf \t wetmix: %lf ", sscale, wetmix); 
	};
	
	void				SetRenderParameters(int nm, float fs, float da1, float da2)
	{
						nummodes = nm;
						
						fscale = fs;
						
						alpha1 =  da1;
						alpha2 =  da2;						

	//					printf("nummodes: %d \t fscale: %lf \t alpha1: %lf \t alpha2: %lf\n", nummodes, fscale, alpha1, alpha2);
	};
	
	void				Thump(int hit); 	
	void				PickUp(int hit);
		
private: 
	
	float *u;			// the transformed values

	float *wn;			// wn = sqrt(lambdai)/2*PI
	float *f;			// d = wn*/PI

	float *phi;			// zeta = (alpha1*wn + alpha2);
	float *d;			// d = wn*PI*tan(phi)
	
	float *fn;			// row of pickup dof
	float *gn;			// column of thump location

	float *an;			// constant before the force
	float *deltan;		// scaled force
	
	float *w_real;			// w = wn + d i
	float *w_imag;
	
	float *c_real;			// c = exp((i * w)/fs)
	float *c_imag;			// c = exp((i * w)/fs)
	
	float* zprev_real;		// the previous value for each resoonator
	float* zprev_imag;		// the previous value for each resoonator
	
	float* znew_real;		// the current value for each resoonator
	float* znew_imag;		// the current value for each resoonator
	
	int					nummodes;			//  number of modes to use
	int					nummodespow2;		// alighed to power of two
	
	// vars**********************************
	EigenData*			eigendata;
		
	// params********************************
	int					location[MAXNUMTHUMPERS];
	
	// cbnote pickup locations
	// if no pickup locations are found, using rayleigh summation for sound synthesis
	
	int					pickup[MAXNUMPICKUPS];
	

	float				dt;
	float				samplerate;
	
	float				fscale; 			
	float				alpha1; 
	float				alpha2;
	float				sscale;
	
	float				wetmix;
	float				drymix;
	
	bool				ready;
	

};




#endif