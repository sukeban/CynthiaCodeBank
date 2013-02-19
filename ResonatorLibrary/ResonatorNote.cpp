/*
 *  ResonatorNote.cpp
 *  VISynth
 *
 *  Created by Cynthia Bruyns on 9/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "ResonatorNote.h"
#include "ResonatorInstrumentBase.h"

void ResonatorNote::Release(UInt32 inFrame)
{
	SynthNote::Release(inFrame);	
}

void ResonatorNote::FastRelease(UInt32 inFrame) 
{
	SynthNote::Release(inFrame);
}

void ResonatorNote::Kill(UInt32 inFrame) // voice is being stolen.
{
	SynthNote::Kill(inFrame);
}

void ResonatorNote::SendNoteNotification (AudioUnitEventType inEventType, UInt32 inParamID)
{
	AudioUnitEvent event;
	event.mEventType = inEventType;
	event.mArgument.mParameter.mAudioUnit = GetAudioUnit()->GetComponentInstance();
	event.mArgument.mParameter.mParameterID = inParamID;
	event.mArgument.mParameter.mScope = kAudioUnitScope_Global;
	event.mArgument.mParameter.mElement = 0;
	
	AUEventListenerNotify(NULL, NULL, &event);
}

#pragma mark ___ActualFunctions___

ComponentResult ResonatorNote::Render(UInt32 inNumberFrames, AudioBufferList& inBufferList ) 
{ 
	
	switch (GetState()){
	
		case kNoteState_FastReleased:
		case kNoteState_Released:{
		
			UInt32 endFrame = 0xFFFFFFFF;
			for (UInt32 i=0; i<inNumberFrames; i++)	{
				if (endFrame == 0xFFFFFFFF) endFrame = i;				
			}
			if (endFrame != 0xFFFFFFFF) {
				NoteEnded(endFrame);
			}
			break;
		}
	
		case kNoteState_Attacked :
		case kNoteState_Sostenutoed :
		case kNoteState_ReleasedButSostenutoed :
		case kNoteState_ReleasedButSustained :{									
			resonator->Render(inNumberFrames, inBufferList);						
			break;
		}

		default :
			break;
	} //end switch
	
	return noErr;
}

void ResonatorNote::Attack(const MusicDeviceNoteParams &inParams)
{	
	ResonatorInstrumentBase* synth = (ResonatorInstrumentBase*) GetAudioUnit();
	
	synth->res = (int)synth->Globals()->GetParameter(kResonatorParam_Resolution);
	synth->mat = (int)synth->Globals()->GetParameter(kResonatorParam_Material);
	synth->noteHit = (abs((int) inParams.mPitch - LOW_KEY))%synth->numThump; 	
											
#ifdef DEBUG_PRINT
	printf("pitch: %d noteHit: %ld\n,",(int) inParams.mPitch,synth->noteHit);
	printf("computing with res: %ld using material: %ld\n",synth->res,synth->mat);
#endif

	resonator->SetEigenData(synth->eigendata[synth->res*synth->numRes + synth->mat]); // cbnote all of the inaudible ones removed
	resonator->SetSampleRate(SampleRate());	
	resonator->SetRenderParameters(
			(int) synth->Globals()->GetParameter (kResonatorParam_NumberModes),
			synth->Globals()->GetParameter(kResonatorParam_DeltaFreqScale),
			synth->Globals()->GetParameter (kResonatorParam_DeltaAlpha1),
			synth->Globals()->GetParameter (kResonatorParam_DeltaAlpha2),
			synth->Globals()->GetParameter(kResonatorParam_DeltaSoundScale)
	);

	resonator->Thump(inParams.mVelocity/128.0, synth->thumpers[synth->noteHit], synth->Globals()->GetParameter (kResonatorParam_ThumpBase)); // resonator only needs to know where
	SendNoteNotification (kAudioUnitEvent_BeginParameterChangeGesture, kParam_NoteOnEvent);
}

