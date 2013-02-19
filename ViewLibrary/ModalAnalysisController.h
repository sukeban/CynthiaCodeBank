//
//  ModalAnalysisController.h
//  ModelViewer
//
//  Created by Cynthia Bruyns on 3/19/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>



@interface ModalAnalysisController : NSObject {
	
	IBOutlet NSObject		*myController; // for telling the au about the picked points

//#pragma mark ____ParameterHandling___

	IBOutlet NSSlider		* stiffnessSlider;		// for changing elastic modulus
	IBOutlet NSTextField	* stiffnessTextField;

	IBOutlet NSSlider		* densitySlider;		// for changing density
	IBOutlet NSTextField	* densityTextField;	

	IBOutlet NSSlider		* poissonsSlider;		// for changing Poisson's ratio
	IBOutlet NSTextField	* poissonsTextField;	

	IBOutlet NSSlider		* thicknessSlider;		// for changing thickness
	IBOutlet NSTextField	* thicknessTextField;	

	IBOutlet NSSlider		* objectScaleSlider;	// for changing size of the object
	IBOutlet NSTextField	* objectScaleTextField;	
	
	IBOutlet NSPopUpButton	* lumpedMassButton;

}

// IBActions from user interfaces

- (IBAction) stiffnessChanged: (id) sender;
- (IBAction) densityChanged: (id) sender;
- (IBAction) poissonsChanged: (id) sender;
- (IBAction) thicknessChanged: (id) sender;
- (IBAction) objectScaleChanged: (id) sender;
- (IBAction) lumpedMassChanged: (id) sender;

- (IBAction) redoAnalysis: (id) sender;


// accessors for ui 
- (NSSlider*) stiffnessSlider;
- (NSTextField*) stiffnessTextField;

- (NSSlider*) densitySlider;
- (NSTextField*) densityTextField;

- (NSSlider*) poissonsSlider;
- (NSTextField*) poissonsTextField;

- (NSSlider*) thicknessSlider;
- (NSTextField*) thicknessTextField;

- (NSSlider*) objectScaleSlider;
- (NSTextField*) objectScaleTextField;

- (NSPopUpButton*) lumpedMassButton;

@end
