/*
 *  AcceleratedConvolutionResonatorSharedData.h
 *  PlateReverb
 *
 *  Created by Cynthia Bruyns on 3/21/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
 
#import <AudioUnit/AudioUnit.h>
#import <AudioToolbox/AudioToolbox.h>

#ifndef AccelaratedConvolutionResonatorSharedData_H
#define AccelaratedConvolutionResonatorSharedData_H


#define MAXNUMTHUMPERS		1
#define MAXNUMPICKUPS		2

#define MIN_AUDIBLE_FREQ	20.0
#define MAX_AUDIBLE_FREQ	20000.0


static CFStringRef kNumberModesParameter = CFSTR("Number Modes");
static int kDefaultValue_NumberModes = 10;

static CFStringRef kDeltaFreqScaleParameter = CFSTR("Freq Scale"); 
static float kDefaultValue_DeltaFreqScale = 1.0;

static CFStringRef kThumpBaseParameter = CFSTR("Thump Base "); 
static float kDefaultValue_ThumpBase = 1.0;

static CFStringRef kDeltaAlpha1Parameter = CFSTR("Alpha1"); 
static float kDefaultValue_DeltaAlpha1 = 0.1;

static CFStringRef kDeltaAlpha2Parameter = CFSTR("Alpha2"); 
static float kDefaultValue_DeltaAlpha2 = 10.0;

static CFStringRef kDeltaSoundScaleParameter = CFSTR("Sound Scale"); 
static float kDefaultValue_DeltaSoundScale = 0.0; 

static CFStringRef kWetDryParameter = CFSTR("Wet Dry Mix"); 
static float kDefaultValue_WetDry = 50.0; 


// parameters
enum {
	kResonatorParam_NumberModes = 100,
	
	kResonatorParam_Resolution,
	kResonatorParam_Material,
		
	kResonatorParam_DeltaFreqScale, 	
	kResonatorParam_DeltaAlpha1,
	kResonatorParam_DeltaAlpha2,
	
	kResonatorParam_DeltaSoundScale,
	
	kResonatorParam_WetDryMix,
	
};


enum{
		
	kResonatorProperty_ModelName	= 100000,
		
	kResonatorProperty_NumModes		= 100001,
	
	kResonatorProperty_Note			= 100002,
	
	kResonatorProperty_NumThump		= 100003,
	kResonatorProperty_LastThump	= 100004,
	kResonatorProperty_MoveThump	= 100005,
	
	kResonatorProperty_NumPickUps	= 100006,
	kResonatorProperty_LastPickUp	= 100007,
	kResonatorProperty_MovePickUp	= 100008,
	
	kResonatorProperty_GeomData		= 100009,
	kResonatorProperty_EigenData	= 100010
};



	
#endif