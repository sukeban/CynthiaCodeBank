/*
 *  ResonatorSharedData.h
 *  VISynth
 *
 *  Created by Cynthia Bruyns on 9/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
 
#import <AudioUnit/AudioUnit.h>
#import <AudioToolbox/AudioToolbox.h>

#ifndef ResonatorSharedData_H
#define ResonatorSharedData_H


#define LOW_KEY				41			//  these will change depending on the keyboard
#define NUM_KEYS			10

#define MAXNUMTHUMPERS		10


#define MIN_AUDIBLE_FREQ	20.0
#define MAX_AUDIBLE_FREQ	6000.0 //cbnote


// parameters
enum {
	kResonatorParam_NumberModes = 100,
	
	kResonatorParam_Resolution,
	kResonatorParam_Material,
	
	kResonatorParam_KeyMode, 
	
	kResonatorParam_DeltaFreqScale, 	
	kResonatorParam_DeltaAlpha1,
	kResonatorParam_DeltaAlpha2,
	
	kResonatorParam_DeltaSoundScale,
	
	kResonatorParam_StrikeRadius,
	kResonatorParam_ThumpBase
};

enum {
	kParam_NoteOnEvent  = 1000,
	kParam_NoteOffEvent = 1001
};


enum{
		
	kResonatorProperty_ModelName	= 100002,
	kResonatorProperty_NumRes		= 100003,
	kResonatorProperty_NumMat		= 100004,
	kResonatorProperty_NumModes		= 100000,
	
	kResonatorProperty_Note			= 100001,
	kResonatorProperty_NumThump		= 100005,
	kResonatorProperty_LastThump	= 100006,
	
	kResonatorProperty_GeomData		= 100007,
	kResonatorProperty_EigenData	= 100008
};

class SimulationData { // cbnote this is seperated out in case the view wants it
	public:

		float		
					*omega_plus, *omega_minus, 
					
					*c1_real, *c2_real, 
					*c1_imag, *c2_imag, 
					
					*increment_plus_real, *increment_minus_real, 
					*increment_plus_imag, *increment_minus_imag, 
					
					*initial_phase_plus, *initial_phase_minus, 
					
					*frequency;
		
		bool		*underdamped;
		
		float		*forces;	
			
		int			*audible;
		int			numaudible;
		
		
		
		
		
		SimulationData(int numFreq, int numDOFS){
		
			underdamped = (bool*) calloc(numFreq, sizeof(bool));
			
			frequency = (float*) calloc (numFreq,sizeof(float));		
			
			omega_plus = (float*) calloc (numFreq,sizeof(float));		
			omega_minus = (float*) calloc (numFreq,sizeof(float));	
						
			c1_real = (float*) calloc (numFreq,sizeof(float));			// the coefficients
			c2_real = (float*) calloc (numFreq,sizeof(float));		
			
			c1_imag = (float*) calloc (numFreq,sizeof(float));			// only for underdamped case	
			c2_imag = (float*) calloc (numFreq,sizeof(float));		

			initial_phase_plus = (float*) calloc (numFreq,sizeof(float));	// only for underdamped case
			initial_phase_minus = (float*) calloc (numFreq,sizeof(float));
			
			increment_plus_real = (float*) calloc (numFreq,sizeof(float));	
			increment_minus_real = (float*) calloc (numFreq,sizeof(float));	
			
			increment_plus_imag = (float*) calloc (numFreq,sizeof(float));	
			increment_minus_imag = (float*) calloc (numFreq,sizeof(float));	

			forces = (float*) calloc(numDOFS,sizeof(float));  
		
			audible = (int*) calloc(numFreq,sizeof(int));
			numaudible = 0;
					
		};
		
		~SimulationData(){
		
			if (underdamped) free(underdamped);
		
			if (frequency) free(frequency);
			
			if (omega_plus) free(omega_plus);
			if (omega_minus) free(omega_minus);
			
			if (initial_phase_plus) free(initial_phase_plus);
			if (initial_phase_minus) free(initial_phase_minus);
			
			if (c1_real) free(c1_real);
			if (c2_real) free(c2_real);

			if (c1_imag) free(c1_imag);
			if (c2_imag) free(c2_imag);

			if (increment_plus_real) free(increment_plus_real);
			if (increment_minus_real) free(increment_minus_real);
			
			if (increment_plus_imag) free(increment_plus_imag);
			if (increment_minus_imag) free(increment_minus_imag);

			if (forces) free(forces);

			if (audible) free(audible);	
		};
	
};

#endif