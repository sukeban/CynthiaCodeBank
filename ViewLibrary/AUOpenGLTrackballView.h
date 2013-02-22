//
//  AUOpenGLTrackballView.h
//  SDIFSynth
//
//  Created by Cynthia Bruyns on 8/11/05.
//  Copyright 2005 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>

#include "OpenGLSharedData.h"

/*!
 a NSOpenGLView subclass used to move the viewing camera around in 3d space
 */

@interface AUOpenGLTrackballView : NSOpenGLView {

	GLfloat			light_pos[4];
	GLfloat			ambientIntensity[4];
	GLfloat			diffuseIntensity[4];
	GLfloat			specularIntensity[4];
	GLfloat			spotExponent[1];

	recCamera		camera;

	GLint			gDollyPanStartPoint[2];
	GLfloat			gTrackBallRotation [4];
	GLboolean		gDolly;
	GLboolean		gPan;
	GLboolean		gTrackball;
	
	GLfloat			worldRotation[4];
	GLfloat			shapeSize;
	
	NSOpenGLView	*gTrackingViewInfo;
}

+ (NSOpenGLPixelFormat*) basicPixelFormat;

- (void) updateProjection;
- (void) updateModelView;
- (void) resizeGL;

- (void) mouseDown:(NSEvent *)theEvent;
- (void) rightMouseDown:(NSEvent *)theEvent;
- (void) otherMouseDown:(NSEvent *)theEvent;
- (void) mouseUp:(NSEvent *)theEvent;
- (void) rightMouseUp:(NSEvent *)theEvent;
- (void) otherMouseUp:(NSEvent *)theEvent;
- (void) mouseDragged:(NSEvent *)theEvent;
- (void) scrollWheel:(NSEvent *)theEvent;
- (void) rightMouseDragged:(NSEvent *)theEvent;
- (void) otherMouseDragged:(NSEvent *)theEvent;

- (BOOL) acceptsFirstResponder;
- (BOOL) becomeFirstResponder;
- (BOOL) resignFirstResponder;

@end
