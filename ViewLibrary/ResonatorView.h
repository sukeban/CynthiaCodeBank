//
//  ResonatorView.h
//  VISynth
//
//  Created by Cynthia Bruyns on 9/8/06.
//  Copyright 2006 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <Carbon/Carbon.h>

#import <OpenGL/OpenGL.h>
#import <OpenGL/gl.h>
#import <OpenGL/glext.h>
#import <OpenGL/glu.h>

#import "AUOpenGLTrackballView.h"

#import "ModelSharedData.h"

#define NUMBUFFERS     2
#define MAXNUMTHUMPERS 10
#define NUM_KEYS       31

static GLfloat default_back_mat_ambient[] =
{ 1.0f, 1.0f, 1.0f, 1.0f };
static GLfloat default_back_mat_diffuse[] = 
{ 1.0f, 1.0f, 1.0f, 1.0f };
static GLfloat default_back_mat_specular[] = 
{ 1.0f, 1.0f, 1.0f, 1.0f };
static GLfloat default_back_mat_shininess[] = 
{ 0.0 };

static GLfloat default_mat_ambient[] = 
{ 1.0f, 0.9f, 0.0f, 1.0f };
static GLfloat default_mat_diffuse[] = 
{ 1.0f, 0.9f, 0.0f, 1.0f };
static GLfloat default_mat_specular[] = 
{ 1.0f, 0.9f, 0.0f, 1.0f };
static GLfloat default_mat_shininess[] = 
{ 0.0 };


static GLfloat  default_light_pos[] = 
{ 0.0f, 0.0f, 0.0f };

static GLfloat  default_ambientIntensity[] = 
{ 1.0f,1.0f,0.0f,1.0f };
static GLfloat  default_diffuseIntensity[] = 
{ 1.0f,1.0f,0.0f,1.0f };
static GLfloat  default_specularIntensity[] = 
{ 1.0f,1.0f,0.0f,1.0f };
static GLfloat  default_spotExponent[] = 
{0.0f};

/*!
 a trackball view subclass with more parameters for controlling a 3d model for rendering
*/
@interface ResonatorView : AUOpenGLTrackballView {
	
	IBOutlet NSObject		*myController; // for telling the au about the picked points

// view stuff
	IBOutlet NSButton		*resetButton;
	
	IBOutlet NSSlider		*curModeSlider;
	IBOutlet NSTextField	*curModeStaticText;	// for max number
	IBOutlet NSTextField	*curModeTextField;
			
	IBOutlet NSTextField	*hertzTextField; // a static display of the current resonant frequency
	
	IBOutlet NSSlider		*deformScaleSlider;
	IBOutlet NSTextField	*deformScaleTextField;
	
	IBOutlet NSButton		*drawWireframeButton;
	IBOutlet NSSlider		*wireframeOffsetSlider;
	IBOutlet NSTextField	*wireframeOffsetTextField;
	
	IBOutlet NSButton		*drawNormalsButton;
	IBOutlet NSSlider		*normalSizeSlider;
	IBOutlet NSTextField	*normalSizeTextField;
	
	IBOutlet NSButton		*drawThumpersButton;	
	IBOutlet NSSlider		*thumperSizeSlider;
	IBOutlet NSTextField	*thumperSizeTextField;
	
	IBOutlet NSButton		*drawFacesButton;
	IBOutlet NSButton		*drawAxesButton;
	IBOutlet NSButton		*drawPickingRayButton; // pick size?	
	
// open gl stuff	
	BOOL		fGLInit;
	
	NSTimer		*timer;

	double		zoomvalue;
	
	GLfloat		mat_ambient[4];
	GLfloat		mat_diffuse[4];
	GLfloat		mat_specular[4];
	GLfloat		mat_shininess[1];

	GLfloat		back_mat_ambient[4];
	GLfloat		back_mat_diffuse[4];
	GLfloat		back_mat_specular[4];
	GLfloat		back_mat_shininess[1];
	
	GLfloat clear[3];


	double		pickOrigin[3];
	double		pickDirection[3];
	
	double		thumpers[MAXNUMTHUMPERS*3];
	int			numthumpers;
	
	int			picked_face;
	
	int			curThump;
	
	bool		justthumped;
	int			numtimesteps;
			
	GeomData*	geomdata;
	EigenData*	eigendata;

}


// IBActions
- (IBAction) redraw: (id) sender;
- (IBAction) resetView:(id)sender;

- (IBAction) setWireframeOffset:(id)sender;
- (IBAction) setNormalScale:(id)sender;
- (IBAction) setThumperScale:(id)sender;

- (IBAction) setDeformScale: (id) sender;
- (IBAction) setCurMode:(id)sender;

- (IBAction) resetThumpers: (id) sender;

// for drawing that it became active
- (int) getIntersectionPoint: (double*) intPt withWindowPoint: (NSPoint) point;

- (void) Thump: (int) intValue;
- (void) updateGeometry:(NSTimer *)t;
- (void) setNumThumpers: (int) num;

// from the controller
- (void) setMaxFreq: (UInt32) maxF;

- (void) setGeomData: (GeomData*) data;
- (void) swapInGeomData: (GeomData*) data;

- (void) setEigenData: (EigenData*) data;
- (void) swapInEigenData: (EigenData*) data;

// mouse handling
- (void)thumperMouseDown:(NSEvent *)theEvent;
- (void) thumperMouseDragged:(NSEvent *)theEvent;
//drawing
- (void) drawRect:(NSRect)rect;
- (void) drawThumpers;
- (void) drawGeometry;
- (void) drawWireframe;
- (void) drawPickingRay;

@end
