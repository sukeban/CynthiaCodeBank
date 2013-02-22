/*
 *  ResonatorInstrumentBase.h
 *  VISynth
 *
 *  Created by Cynthia Bruyns on 9/7/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef ResonatorInstrumentBase_H_
#define ResonatorInstrumentBase_H_

#include "ResonatorNote.h"

static const UInt32 kInitOutputBuses = 1;
static const UInt32 kMaxNotes = 4;				
static const AUChannelInfo sChannels = { 0, 2 };


#pragma mark ___Params___

static CFStringRef kNumberModesParameter = CFSTR("Number Modes");
static int kDefaultValue_NumberModes = 20;

static CFStringRef kStrikeRadiusParameter = CFSTR("Strike Radius");  // cbnote not used right now
static int kDefaultValue_StrikeRadius = 0; 

static CFStringRef kDeltaFreqScaleParameter = CFSTR("Freq Scale"); 
static float kDefaultValue_DeltaFreqScale = 1000.0;

static CFStringRef kThumpBaseParameter = CFSTR("Thump Base "); 
static float kDefaultValue_ThumpBase = 1.0;

static CFStringRef kResolutionParameter = CFSTR("Resolution"); 
static int	 kDefaultValue_Resolution = 0; 

static CFStringRef kMaterialParameter = CFSTR("Material"); 
static int	 kDefaultValue_Material = 0; 

static CFStringRef kKeyModeParameter = CFSTR("Key mode"); // cbnote keymode would raise and lower the pitch depending on the key, not we have freq scale for that now
static bool	 kDefaultValue_KeyMode = false;


static CFStringRef kDeltaAlpha1Parameter = CFSTR("Alpha1"); 
static float kDefaultValue_DeltaAlpha1 = 0.0;

static CFStringRef kDeltaAlpha2Parameter = CFSTR("Alpha2"); 
static float kDefaultValue_DeltaAlpha2 = 0.0;

static CFStringRef kDeltaSoundScaleParameter = CFSTR("Sound Scale "); 
static float kDefaultValue_DeltaSoundScale = 1.0; // cbnote change these from delta

// cbnote parameters to add here: 
//  the direction of the strike
//  the duration of the strike


/*!
 an AUMonotimbralInstrumentBase subclass used to render notes
 this class is used to register a synthesizer to be used in Audio Unit hosts
 */
class ResonatorInstrumentBase : public AUMonotimbralInstrumentBase {
private:

public: 
	ResonatorNote mMyNotes[kMaxNotes]; 

	std::vector<GeomData*>	geomdata;	// for different resolutions of the same object
	std::vector<EigenData*>	eigendata;	// also for different materials of the same geometry

	// public for the note
	UInt32	res;
	UInt32	mat;
	
	UInt32	numRes;
	UInt32	numMat;
	
	UInt32	noteHit;						
	UInt32	numThump;						
	int*	thumpers;
	
			
public: 

	void	ClearEigenData();
	void	ClearGeomData();
	
	ResonatorInstrumentBase(ComponentInstance inComponentInstance);
	~ResonatorInstrumentBase();	
	
	virtual ComponentResult		Initialize();
//	virtual void				Cleanup();
//	virtual ComponentResult		Version();
	

	virtual ComponentResult		GetParameterInfo(	AudioUnitScope			inScope,
													AudioUnitParameterID	inParameterID,
													AudioUnitParameterInfo	&outParameterInfo );

	/*! @method GetPropertyInfo */
	virtual ComponentResult		GetPropertyInfo(	AudioUnitPropertyID		inID,
													AudioUnitScope			inScope,
													AudioUnitElement		inElement,
													UInt32 &				outDataSize,
													Boolean &				outWritable);

	virtual ComponentResult		SetProperty(		AudioUnitPropertyID 	inID,
													AudioUnitScope 			inScope,
													AudioUnitElement 		inElement,
													const void *			inData,
													UInt32 					inDataSize);
	
	virtual ComponentResult		GetProperty(		AudioUnitPropertyID 	inID,
													AudioUnitScope 			inScope,
													AudioUnitElement		inElement,
													void *					outData);

		
	virtual OSStatus			HandleMidiEvent(	UInt8 	inStatus,
													UInt8 	inChannel,
													UInt8 	inData1,
													UInt8 	inData2,
													long 	inStartFrame);

	virtual	bool				SupportsTail () { return false; }
			
};

#endif
