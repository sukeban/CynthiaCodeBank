/*
 *  ResonatorNote.h
 *  VISynth
 *
 *  Created by Cynthia Bruyns on 9/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef ResonatorNote_H_
#define ResonatorNote_H_

#include "AUInstrumentBase.h"
#include "AcceleratedResonator.h"

class ResonatorInstrumentBase;

/*!
 an Audio Unit wrapper around a SynthNode object
 this wrapper is used to render the resonator values
 */
class ResonatorNote : public SynthNote{
public:

	ResonatorNote(){
		resonator = new AcceleratedResonator;
	};
	
	~ResonatorNote(){
		if (resonator) delete(resonator);
	}
	
	void		Attack(const MusicDeviceNoteParams &inParams);
	
	void		Kill(UInt32 inFrame);
	void		Release(UInt32 inFrame);
	void		FastRelease(UInt32 inFrame);
	
	Float32		Amplitude() { return amp; } 
	OSStatus	Render(UInt32 inNumFrames, AudioBufferList& inBufferList);
	
	void		SendNoteNotification (AudioUnitEventType inEventType, UInt32 inParamID);
			
	// funcs**********************************
	ResonatorInstrumentBase*		GetSynth() { return (ResonatorInstrumentBase*)GetAudioUnit(); }

	//dumb note vars*************************
	float							amp;

	//the actual data*************************
	AcceleratedResonator*			resonator; // the guy who makes sound

};

#endif