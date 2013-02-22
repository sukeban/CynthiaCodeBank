//
//  SoundRenderingController.h
//  ModelViewer
//
//  Created by Cynthia Bruyns on 3/22/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>

/*!
 an object for controlling the sound rendering parameters
 */
@interface SoundRenderingController : NSObject {

	IBOutlet NSObject		*myController; // for telling the au about the picked points

//#pragma mark ____ParameterHandling___

	IBOutlet NSSlider		* lambda1Slider;			// damping
	IBOutlet NSTextField	* lambda1TextField;
	IBOutlet NSSlider		* lambda2Slider;			// damping
	IBOutlet NSTextField	* lambda2TextField;
	
	IBOutlet NSSlider		* numModeSlider;			// how many modes to compute sound with
	IBOutlet NSTextField	* numModeStaticTextField;	// how many possible
	IBOutlet NSTextField	* numModeTextField;
	
	IBOutlet NSSlider		* thumpBaseSlider;			// how hard to hit it
	IBOutlet NSTextField	* thumpBaseTextField;
	
	IBOutlet NSSlider		* soundScaleSlider;			// scale sound up
	IBOutlet NSTextField	* soundScaleTextField;
	
	IBOutlet NSSlider		* freqScaleSlider;			// scale sound up
	IBOutlet NSTextField	* freqScaleTextField;

	IBOutlet NSSlider		* wetDrySlider;				// wet dry mix
	IBOutlet NSTextField	* wetDryTextField;

}

// actions
- (IBAction) lambda1Changed: (id) sender;
- (IBAction) lambda2Changed: (id) sender;
- (IBAction) numModesChanged: (id) sender;
- (IBAction) thumpForceChanged: (id) sender;
- (IBAction) soundScaleChanged: (id) sender;
- (IBAction) freqScaleChanged: (id) sender;
- (IBAction) wetDryChanged: (id) sender;


// accessors for ui
- (NSSlider*) lambda1Slider;			
- (NSTextField*) lambda1TextField;

- (NSSlider*)lambda2Slider;			
- (NSTextField*) lambda2TextField;
	
- (NSSlider*) numModeSlider;			
- (NSTextField*) numModeStaticTextField;	
- (NSTextField*) numModeTextField;
	
- (NSSlider*) thumpBaseSlider;			
- (NSTextField*) thumpBaseTextField;
	
- (NSSlider*) soundScaleSlider;			
- (NSTextField*) soundScaleTextField;
	
- (NSSlider*) freqScaleSlider;			
- (NSTextField*) freqScaleTextField;

- (NSSlider*) wetDrySlider;				
- (NSTextField*) wetDryTextField;


@end
