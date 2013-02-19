//
//  ModalAnalysisController.mm
//  ModelViewer
//
//  Created by Cynthia Bruyns on 3/19/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "ModalAnalysisController.h"

@implementation ModalAnalysisController


#pragma mark ____CleanUp
-(void) dealloc
{	
	[super dealloc];
}


#pragma mark ___AnalaysisParams___

- (IBAction) thicknessChanged: (id) sender
{
	if (sender == thicknessSlider){
		[thicknessTextField setFloatValue: [sender floatValue]];
	}
	else{
		[thicknessSlider setFloatValue: [sender floatValue]];	
	}
	[myController thicknessChanged: sender];
}

- (IBAction) stiffnessChanged: (id) sender
{
	if (sender == stiffnessSlider){
		[stiffnessTextField setFloatValue: [sender floatValue]*1e9];
	}
	else{
		[stiffnessSlider setFloatValue: [sender floatValue]/1e9];	
	}
	[myController stiffnessChanged: sender];
}

- (IBAction) poissonsChanged: (id) sender
{
	if (sender == poissonsSlider){
		[poissonsTextField setFloatValue: [sender floatValue]];
	}
	else{
		[poissonsSlider setFloatValue: [sender floatValue]];	
	}
	[myController poissonsChanged: sender];
}

- (IBAction) objectScaleChanged: (id) sender
{
	if (sender == objectScaleSlider){
		[objectScaleTextField setFloatValue: [sender floatValue]];
	}
	else{
		[objectScaleSlider setFloatValue: [sender floatValue]];	
	}
	[myController objectScaleChanged: sender];
}

- (IBAction) densityChanged: (id) sender
{
	if (sender == densitySlider){
		[densityTextField setFloatValue: [sender floatValue]];
	}
	else{
		[densitySlider setFloatValue: [sender floatValue]];	
	}
	[myController densityChanged: sender];
}

- (IBAction) lumpedMassChanged: (id) sender
{
	[myController lumpedMassChanged: sender];
}



- (IBAction) redoAnalysis: (id) sender
{

	printf("going to compute eigenvalues\n");
		
	[myController redoAnalysis: sender];
	
	printf("done computing eigenvalues\n");
	
}

#pragma mark ___Accessors___


- (NSSlider*) stiffnessSlider
{	
	return stiffnessSlider;
}

- (NSTextField*) stiffnessTextField
{
	return stiffnessTextField;
}


- (NSSlider*) densitySlider
{
	return densitySlider;
}
- (NSTextField*) densityTextField
{
	return densityTextField;
}


- (NSSlider*) poissonsSlider
{
	return poissonsSlider;
}

- (NSTextField*) poissonsTextField
{
	return poissonsTextField;
}



- (NSSlider*) thicknessSlider
{
	return thicknessSlider;
}

- (NSTextField*) thicknessTextField
{
	return thicknessTextField;
}


- (NSSlider*) objectScaleSlider
{
	return objectScaleSlider;
}

- (NSTextField*) objectScaleTextField
{
	return objectScaleTextField;
}


- (NSPopUpButton*) lumpedMassButton
{
	return lumpedMassButton;
}

@end
