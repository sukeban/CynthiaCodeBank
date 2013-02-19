//
//  SoundRenderingController.mm
//  ModelViewer
//
//  Created by Cynthia Bruyns on 3/22/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "SoundRenderingController.h"



@implementation SoundRenderingController

#pragma mark ___SoundParams___
- (IBAction) numModesChanged: (id) sender
{	
	if (sender == numModeSlider){
		[numModeTextField setIntValue: [sender intValue]];
	}
	else{
		[numModeSlider setIntValue: [sender intValue]];	
	}	
	
	[myController numModesChanged: sender];
}

- (void) makeEmodPoissons
{
	double lambda = [lambda1TextField doubleValue]*1e-7;
	double nu = [lambda2TextField doubleValue]*1e-3;
	
	double emod = nu*(3.0*lambda+2.0*nu)/(lambda+nu);
	double poissons = lambda/(lambda+nu)/2.0;

}

- (IBAction) lambda1Changed: (id) sender
{
	if (sender == lambda1Slider){
		[lambda1TextField setFloatValue: [sender doubleValue]];
	}
	else{
		[lambda1Slider setFloatValue: [sender doubleValue]];	
	}
	
	//[self makeEmodPoissons];
	[myController lambda1Changed: sender];
}

- (IBAction) lambda2Changed: (id) sender
{
	if (sender == lambda2Slider){
		[lambda2TextField setFloatValue: [sender doubleValue]];
	}
	else{
		[lambda2Slider setFloatValue: [sender doubleValue]];	
	}
	
	//[self makeEmodPoissons];
	[myController lambda2Changed: sender];
}

- (IBAction) thumpForceChanged: (id) sender
{
	if (sender == thumpBaseSlider){
		[thumpBaseTextField setFloatValue: [sender doubleValue]];
	}
	else{
		[thumpBaseSlider setFloatValue: [sender doubleValue]];	
	}
	[myController thumpForceChanged: sender];
}

- (IBAction) soundScaleChanged: (id) sender
{
	if (sender == soundScaleSlider){
		[soundScaleTextField setFloatValue: [sender doubleValue]];
	}
	else{
		[soundScaleSlider setFloatValue: [sender doubleValue]];	
	}
	[myController soundScaleChanged: sender];
}

- (IBAction) freqScaleChanged: (id) sender
{
	if (sender == freqScaleSlider){
		[freqScaleTextField setFloatValue: [sender doubleValue]];
	}
	else{
		[freqScaleSlider setFloatValue: [sender doubleValue]];	
	}
	[myController freqScaleChanged: sender];
}

- (IBAction) wetDryChanged: (id) sender
{
	if (sender == wetDrySlider){
		[wetDryTextField setFloatValue: [sender doubleValue]];
	}
	else{
		[wetDrySlider setFloatValue: [sender doubleValue]];	
	}
	[myController wetDryChanged: sender];
}


#pragma mark ___Accessors___
- (NSSlider*) lambda1Slider
{
	return lambda1Slider;
}
		
- (NSTextField*) lambda1TextField
{
	return lambda1TextField;
}

- (NSSlider*)lambda2Slider
{
	return lambda2Slider;
}
			
- (NSTextField*) lambda2TextField
{
	return lambda2TextField;
}
	
- (NSSlider*) numModeSlider
{
	return numModeSlider;
}

- (NSTextField*) numModeStaticTextField
{
	return numModeStaticTextField;
}

		
- (NSTextField*) numModeTextField
{
	return numModeTextField;
}
	
- (NSSlider*) thumpBaseSlider
{
	return thumpBaseSlider;
}
		
- (NSTextField*) thumpBaseTextField
{
	return thumpBaseTextField;
}
	
- (NSSlider*) soundScaleSlider
{
	return soundScaleSlider;
}

- (NSTextField*) soundScaleTextField
{
	return soundScaleTextField;
}
	
- (NSSlider*) freqScaleSlider
{
	return freqScaleSlider;
}
		
- (NSTextField*) freqScaleTextField
{
	return freqScaleTextField;
}

- (NSSlider*) wetDrySlider
{
	return wetDrySlider;
}
			
- (NSTextField*) wetDryTextField
{
	return wetDryTextField;
}


@end
