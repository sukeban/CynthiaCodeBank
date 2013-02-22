/*
 *  AccelaratedResonatorSharedData.h
 *  ModalComparison
 *
 *  Created by Cynthia Bruyns on 10/16/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

 
#import <AudioUnit/AudioUnit.h>
#import <AudioToolbox/AudioToolbox.h>

#ifndef AccelaratedResonatorSharedData_H
#define AccelaratedResonatorSharedData_H

/*!
 a collection of parameter and property definions for the synthesizers using the accelerated resonator objects
 */

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
		
	kResonatorParam_DeltaFreqScale, 	
	kResonatorParam_DeltaAlpha1,
	kResonatorParam_DeltaAlpha2,
	
	kResonatorParam_DeltaSoundScale,
	
	kResonatorParam_ThumpBase,
	
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

AudioUnitParameter parameter[] = {	{ 0, kResonatorParam_NumberModes, kAudioUnitScope_Global, 0 },		//0

									{ 0, kResonatorParam_Resolution, kAudioUnitScope_Global, 0 },		//1
									{ 0, kResonatorParam_Material, kAudioUnitScope_Global, 0 },			//2

                                    { 0, kResonatorParam_DeltaFreqScale, kAudioUnitScope_Global, 0 },	//3
									{ 0, kResonatorParam_DeltaAlpha1, kAudioUnitScope_Global, 0 },		//4
									{ 0, kResonatorParam_DeltaAlpha2, kAudioUnitScope_Global, 0 },		//5

									{ 0, kResonatorParam_DeltaSoundScale, kAudioUnitScope_Global, 0 },	//6
									
									{ 0, kResonatorParam_ThumpBase, kAudioUnitScope_Global, 0 },		//7
									
									{ 0, kParam_NoteOnEvent, kAudioUnitScope_Global, 0 },				//8
									{ 0, kParam_NoteOffEvent, kAudioUnitScope_Global, 0 },
									
}; // for easy listener creation in the view

	
#endif