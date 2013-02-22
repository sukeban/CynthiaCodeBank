/*
 *  ResonatorInstrumentBase.cpp
 *  ResonatorInstrumentBase
 *
 *  Created by Cynthia Bruyns on 9/7/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "ResonatorInstrumentBase.h"


#pragma mark ___Init Functions___

ComponentResult ResonatorInstrumentBase::Initialize()
{		
	
	AUMonotimbralInstrumentBase::Initialize();
	SetNotes(kMaxNotes, kMaxNotes, mMyNotes, sizeof(ResonatorNote));	
	
	return noErr;
}

/*
AUMonotimbralInstrumentBase::AUMonotimbralInstrumentBase(
							ComponentInstance				inInstance, 
							UInt32							numInputs,
							UInt32							numOutputs,
							UInt32							numGroups,
							UInt32							numParts)
*/
							
ResonatorInstrumentBase::ResonatorInstrumentBase(ComponentInstance inComponentInstance)
: AUMonotimbralInstrumentBase(inComponentInstance, 0, 1, 1)
{
	CreateElements();
	
	res = 0; mat = 0;	
	numRes = 0;	numMat = 0;	
	
	noteHit = 0; numThump = 0;	
		
	thumpers = (int*) calloc(MAXNUMTHUMPERS, sizeof(int));
	
//	Globals()->UseIndexedParameters (kNumberOfParameters);
	Globals()->SetParameter(kResonatorParam_NumberModes,		kDefaultValue_NumberModes);
	Globals()->SetParameter(kResonatorParam_StrikeRadius,		kDefaultValue_StrikeRadius);
	Globals()->SetParameter(kResonatorParam_DeltaFreqScale,		kDefaultValue_DeltaFreqScale);
	Globals()->SetParameter(kResonatorParam_Resolution,			kDefaultValue_Resolution);
	Globals()->SetParameter(kResonatorParam_Material,			kDefaultValue_Material);
	Globals()->SetParameter(kResonatorParam_KeyMode,			kDefaultValue_KeyMode);
	Globals()->SetParameter(kResonatorParam_DeltaAlpha1,		kDefaultValue_DeltaAlpha1);
	Globals()->SetParameter(kResonatorParam_DeltaAlpha2,		kDefaultValue_DeltaAlpha2);
	Globals()->SetParameter(kResonatorParam_DeltaSoundScale,	kDefaultValue_DeltaSoundScale);
	Globals()->SetParameter(kResonatorParam_ThumpBase,			kDefaultValue_ThumpBase);
		
}

void ResonatorInstrumentBase::ClearGeomData()
{
  // start from the beginning of the array
  std::vector<GeomData*>::iterator itPos = geomdata.begin();
  // clear all elements from the array
  for(; itPos < geomdata.end(); itPos++)
    delete *itPos;    // free the element from memory
	
	geomdata.clear(); 
}

void ResonatorInstrumentBase::ClearEigenData()
{
  // start from the beginning of the array
  std::vector<EigenData*>::iterator itPos = eigendata.begin();
  // clear all elements from the array
  for(; itPos < eigendata.end(); itPos++)
    delete *itPos;    // free the element from memory
	
	eigendata.clear(); 
}

ResonatorInstrumentBase::~ResonatorInstrumentBase()
{

	ClearGeomData();
	ClearEigenData();
	
	if (thumpers) free(thumpers);
	
}

#pragma mark ___Parameter Functions___

ComponentResult ResonatorInstrumentBase::GetParameterInfo(	AudioUnitScope			inScope,
														AudioUnitParameterID	inParameterID,
														AudioUnitParameterInfo	&outParameterInfo )
{
	ComponentResult result = noErr;
	
	outParameterInfo.flags = 	kAudioUnitParameterFlag_IsWritable +	kAudioUnitParameterFlag_IsReadable;
	
	if (inScope == kAudioUnitScope_Global) {
		switch(inParameterID)  {
			case kResonatorParam_NumberModes:
				AUBase::FillInParameterName (outParameterInfo, kNumberModesParameter, false);
				outParameterInfo.unit = kAudioUnitParameterUnit_Indexed;
				outParameterInfo.minValue = 0;
				outParameterInfo.maxValue = 500;
				outParameterInfo.defaultValue = kDefaultValue_NumberModes;
				break;
			case kResonatorParam_StrikeRadius:
				AUBase::FillInParameterName (outParameterInfo, kStrikeRadiusParameter, false);
				outParameterInfo.unit = kAudioUnitParameterUnit_Indexed;
				outParameterInfo.minValue = 0;
				outParameterInfo.maxValue = 20;
				outParameterInfo.defaultValue = kDefaultValue_StrikeRadius;
				break;
			case kResonatorParam_DeltaFreqScale:
				AUBase::FillInParameterName (outParameterInfo, kDeltaFreqScaleParameter, false);
				outParameterInfo.unit = kAudioUnitParameterUnit_Generic;
				outParameterInfo.minValue = 0.0;
				outParameterInfo.maxValue = 20000.0;
				outParameterInfo.defaultValue = kDefaultValue_DeltaFreqScale;
				break;
			case kResonatorParam_Resolution:
				AUBase::FillInParameterName (outParameterInfo, kResolutionParameter, false);
				outParameterInfo.unit = kAudioUnitParameterUnit_Indexed;
				outParameterInfo.minValue = 0;
				outParameterInfo.maxValue = 2;
				outParameterInfo.defaultValue = 0;
				break;			
			case kResonatorParam_Material:
				AUBase::FillInParameterName (outParameterInfo, kMaterialParameter, false);
				outParameterInfo.unit = kAudioUnitParameterUnit_Indexed;
				outParameterInfo.minValue = 0;
				outParameterInfo.maxValue = 2;
				outParameterInfo.defaultValue = 0;
				break;			
			case kResonatorParam_KeyMode:
				AUBase::FillInParameterName (outParameterInfo, kKeyModeParameter, false);
				outParameterInfo.unit = kAudioUnitParameterUnit_Boolean;
				outParameterInfo.defaultValue = kDefaultValue_KeyMode;
				break;

			case kResonatorParam_DeltaAlpha1:
				AUBase::FillInParameterName (outParameterInfo, kDeltaAlpha1Parameter, false);
				outParameterInfo.unit = kAudioUnitParameterUnit_Generic;
				outParameterInfo.minValue = 0;
				outParameterInfo.maxValue = 100.0;
				outParameterInfo.defaultValue = kDefaultValue_DeltaAlpha1;
				break;
			case kResonatorParam_DeltaAlpha2:
				AUBase::FillInParameterName (outParameterInfo, kDeltaAlpha2Parameter, false);
				outParameterInfo.unit = kAudioUnitParameterUnit_Generic;
				outParameterInfo.minValue = 0;
				outParameterInfo.maxValue = 100.0;
				outParameterInfo.defaultValue = kDefaultValue_DeltaAlpha2;
				break;
			case kResonatorParam_DeltaSoundScale:
				AUBase::FillInParameterName (outParameterInfo, kDeltaSoundScaleParameter, false);
				outParameterInfo.unit = kAudioUnitParameterUnit_Generic;
				outParameterInfo.minValue = 0.0;
				outParameterInfo.maxValue = 1.0;
				outParameterInfo.defaultValue = kDefaultValue_DeltaSoundScale;
				break;
			case kResonatorParam_ThumpBase:
				AUBase::FillInParameterName (outParameterInfo, kThumpBaseParameter, false);
				outParameterInfo.unit = kAudioUnitParameterUnit_Generic;
				outParameterInfo.minValue = 0.0;
				outParameterInfo.maxValue = 100.0;
				outParameterInfo.defaultValue = kDefaultValue_ThumpBase;
				break;	
				
			default:
				result = kAudioUnitErr_InvalidParameter;
				break;
		}
	} else {
		result = kAudioUnitErr_InvalidParameter;
	}
	
	return result;
}

#pragma mark ___Property Functions___

ComponentResult		ResonatorInstrumentBase::SetProperty(	AudioUnitPropertyID		inID,
														AudioUnitScope 			inScope,
														AudioUnitElement 		inElement,
														const void *			inData,
														UInt32 					inDataSize)
{	
	switch (inID) 
	{			
		case kResonatorProperty_NumThump:
			numThump = 0;
			//cerr << "\n\nclearing thumpers " << endl;
			return noErr;
			break;
			
		case kResonatorProperty_LastThump:{
		if (numThump >= MAXNUMTHUMPERS) return noErr;
			UInt32* pt = (UInt32*) inData;
			thumpers[numThump] = *pt;
			printf("\n\naddingToThump: %d", thumpers[numThump]);
			numThump++;
			printf("making number of thumpers: %ld ",  numThump);
			return noErr;
			break;
		}
				
		default:
			return 
			AUInstrumentBase::SetProperty (inID, inScope, inElement, inData, inDataSize);
	}
	return noErr;
}

ComponentResult		ResonatorInstrumentBase::GetPropertyInfo (	AudioUnitPropertyID		inID,
															AudioUnitScope			inScope,
															AudioUnitElement		inElement,
															UInt32 &				outDataSize,
															Boolean &				outWritable)
{
	if (inScope == kAudioUnitScope_Global) {
		switch (inID) {
		
			case kAudioUnitProperty_ParameterValueFromString:
				outWritable = false;
				outDataSize = sizeof (AudioUnitParameterValueFromString);
				return noErr;
				
			case kAudioUnitProperty_ParameterValueName:
				outWritable = false;
				outDataSize = sizeof (AudioUnitParameterValueName);
				return noErr;
		
			case kResonatorProperty_NumRes:
				outWritable = false;
				outDataSize = sizeof(UInt32);
				return noErr;
				
			case kResonatorProperty_NumMat:
				outWritable = false;
				outDataSize = sizeof(UInt32);
				return noErr;
									
			case kResonatorProperty_NumModes:
				outWritable = false;
				outDataSize = sizeof(UInt32);
				return noErr;


			case	kResonatorProperty_GeomData:
				outWritable = false;
				outDataSize = sizeof (GeomData**);
				return noErr;	
				
			case	kResonatorProperty_EigenData:
				outWritable = false;
				outDataSize = sizeof (EigenData**);
				return noErr;	
			

				
			case kResonatorProperty_Note:
				outWritable = false;
				outDataSize = sizeof(UInt32);
				return noErr;
				
			case kResonatorProperty_NumThump:
				outWritable = false;
				outDataSize = sizeof(UInt32);
				return noErr;
				
			case kResonatorProperty_LastThump:
				outWritable = false;
				outDataSize = sizeof(UInt32);
				return noErr;
				
		}
	}
	return AUInstrumentBase::GetPropertyInfo (inID, inScope, inElement, outDataSize, outWritable);
}


ComponentResult		ResonatorInstrumentBase::GetProperty (	AudioUnitPropertyID			inID,
														AudioUnitScope 				inScope,
														AudioUnitElement			inElement,
														void *						outData)
{
	if (inScope == kAudioUnitScope_Global) {
		switch (inID) {            
			case kAudioUnitProperty_ParameterValueFromString:
			{
				OSStatus retVal = kAudioUnitErr_InvalidPropertyValue;
				return retVal; 
			}
		
			
			
			case	kResonatorProperty_GeomData:
			
				res = (UInt32) Globals()->GetParameter(kResonatorParam_Resolution);
				mat = (UInt32) Globals()->GetParameter(kResonatorParam_Material);
				
				//cerr << "size geomdata: " << geomdata.size() << endl;
				
				if (geomdata.size() == 0) 
					*(static_cast<GeomData**>(outData)) = NULL;
				else
					*(static_cast<GeomData**>(outData)) = geomdata[res];
				return noErr;

			case	kResonatorProperty_EigenData:
			
				res = (UInt32) Globals()->GetParameter(kResonatorParam_Resolution);
				mat = (UInt32) Globals()->GetParameter(kResonatorParam_Material);
				
				if (eigendata.size() ==0) 
					*(static_cast<EigenData**>(outData)) = NULL;
				else
					*(static_cast<EigenData**>(outData)) = eigendata[res*numRes + mat];
				return noErr;
			
			case kResonatorProperty_NumModes:
				if (eigendata.size() == 0) {
						*(static_cast<UInt32*>(outData)) = 0;
						return noErr;
				} 
				
				res = (UInt32) Globals()->GetParameter(kResonatorParam_Resolution);
				mat = (UInt32) Globals()->GetParameter(kResonatorParam_Material);
				
				*(static_cast<UInt32*>(outData)) = eigendata[res*numRes + mat]->numFreq;
				
				return noErr;
				
			case kResonatorProperty_Note:
				*(static_cast<UInt32*>(outData)) = noteHit;
				return noErr;
				
				
			case kResonatorProperty_NumRes:
				*(static_cast<UInt32*>(outData)) = numRes;
				return noErr;
				
			case kResonatorProperty_NumMat:
				*(static_cast<UInt32*>(outData)) = numMat;
				return noErr;
				
			case kResonatorProperty_NumThump:
				*(static_cast<UInt32*>(outData)) = numThump;
				return noErr;
				
			case kResonatorProperty_LastThump:
				*(static_cast<UInt32*>(outData)) = thumpers[numThump-1];
				return noErr;
				
		}
	}
	return AUInstrumentBase::GetProperty (inID, inScope, inElement, outData);
}



#pragma mark ___AUBase Functions___

OSStatus ResonatorInstrumentBase::HandleMidiEvent(	UInt8 	inStatus,
													UInt8 	inChannel,
													UInt8 	inData1,
													UInt8 	inData2,
													long	inStartFrame)
{
	return AUMIDIBase::HandleMidiEvent(inStatus, inChannel, inData1, inData2, inStartFrame);
}

#pragma mark ___Note Functions___

