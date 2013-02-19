/*
 *  AcceleratedResonator.h
 *  ResonatorTestApp
 *
 *  Created by Cynthia Bruyns on 10/14/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef AcceleratedResonator_H
#define AcceleratedResonator_H

#include <stdlib.h>
#import <AudioUnit/AudioUnit.h>
#import <AudioToolbox/AudioToolbox.h>

#include "AcceleratedResonatorSharedData.h"

#include "ModelSharedData.h"
#include "ComplexNumber.h"

class AcceleratedResonator {
public:

	AcceleratedResonator()
	{
		eigendata = NULL;
		
		wn = NULL;
		zeta = NULL;
		force = NULL;
		w = NULL;
		w1 = NULL; w2 = NULL;
		c1 = NULL; c2 = NULL;
		ep1 = NULL; ep2 = NULL;
		z1 = NULL; z2 = NULL;
		
		location = -1;
	};
	
	
	void FreeData()
	{
		if (wn) free(wn);
		if (zeta) free(zeta);
		if (force) free(force);
		if (w) free(w);
		if (w1) free(w1);
		if (w2) free(w2);
		if (c1) free(c1);
		if (c2) free(c2);
		if (ep1) free(ep1);
		if (ep2) free(ep2);
		if (z1) free(z1);
		if (z2) free(z2);

		wn = NULL;
		zeta = NULL;
		force = NULL;
		w = NULL;
		w1 = NULL; w2 = NULL;
		c1 = NULL; c2 = NULL;
		ep1 = NULL; ep2 = NULL;
		z1 = NULL; z2 = NULL;	
	};
		
	~AcceleratedResonator(){	
		FreeData();
	};
		
	void				Initialize();
	
	// callled when you want to render
	void				Render(UInt32 inNumFrames, AudioBufferList& inBufferList);

	// called from attack**********************************
	void				SetEigenData(EigenData* ed)
	{
						eigendata = ed;	
						ready = false;						
					
						FreeData();

						int doublesize = sizeof(double); 
						int complexsize = sizeof(Complex);
												
						num = eigendata->numFreq;					// only audible ones are included, i.e. pruning has happened
						
						wn = (double*) calloc(num, doublesize);		// wn = sqrt(lambdai)/2*PI		
						zeta = (double*) calloc(num, doublesize);	// zeta = (alpha1*wn + alpha2)/(2*wn);	
						force = (double*) calloc(num, sizeof(double));
							
						w = (Complex*) calloc(num, complexsize);	// w = sqrt(b^2 - 4 wn)/2

						w1 = (Complex*) calloc(num, complexsize);	// w1 = -b + w
						w2 = (Complex*) calloc(num, complexsize);	// w2 = -b - w
						
						c1 = (Complex*) calloc(num, complexsize);	// c1 = z0/2 + b*z0/w + zdot0/w
						c2 = (Complex*) calloc(num, complexsize);	// c2 = z0/2 - b*z0/w + zdot0/w

						ep1 = (Complex*) calloc(num, complexsize);	// e1p = exp(w1*dt)*(cos(phase1*dt) + i sin(phase1*dt))
						ep2 = (Complex*) calloc(num, complexsize);	// e2p = exp(w1*dt)*(cos(phase2*dt) + i sin(phase2*dt))

						z1 = (Complex*) calloc(num, complexsize);	// rolling values
						z2 = (Complex*) calloc(num, complexsize);	// rolling values
									
	};
	
	void				SetSampleRate(float sr)
	{ 
						samplerate = sr;
						dt = 1.0/samplerate;
	};
	
	void				SetRenderParameters(int nm, float fs, float da1, float da2, float ss)
	{
						nummodes = nm;
						fscale = fs;
						alpha1 =  da1;
						alpha2 =  da2;
						sscale = ss;
	};
	
	void				Thump(float velocity, int hit, float f);
		
private: 
	
	double *wn;			// wn = sqrt(lambdai)/2*PI
	
	double *zeta;		// zeta = (alpha1*wn + alpha2);
	double *force;		// forces
	
	Complex *w;			// w = sqrt(b^2 - 4 wn)/2

	Complex *w1;		// w1 = -b + w
	Complex *w2;		// w2 = -b - w
	
	Complex *c1;		// c1 = z0/2 + b*z0/w + 2*zdot0/w
	Complex *c2;		// c2 = z0/2 - b*z0/w + 2*zdot0/w

	Complex *ep1;		// e1p = exp(w1*dt)*(cos(phase1*dt) + i sin(phase1*dt))
	Complex *ep2;		// e2p = exp(w1*dt)*(cos(phase2*dt) + i sin(phase2*dt))

	Complex *z1;		// rolling values
	Complex *z2;		// rolling values
	
	
	// vars**********************************
	EigenData*			eigendata;
		
	// params********************************
	int					location;
	double				thumpforce;
	
	int					nummodes;	// number to render
	int					num;		// number given 
	
	int					lowest;
	int					highest; 
	
	double				dt;
	double				t;
	double				samplerate;
	
	double				fscale; 			
	double				alpha1; 
	double				alpha2;
	double				sscale;

	bool				ready;

};




#endif