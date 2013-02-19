//
//  AUOpenGLTrackballView.mm
//  SDIFSynth
//
//  Created by Cynthia Bruyns on 8/11/05.
//  Copyright 2005 __MyCompanyName__. All rights reserved.
//

#import "AUOpenGLTrackballView.h"
#import "trackball.h"

// directional lighting is fastest (last value on pos is 0 versus 1)
static GLfloat  default_light_pos[] = { 1.0, 1.0, 2.0, 0.0 };

static GLfloat  default_ambientIntensity[] = { 1.0, 1.0, 1.0, 1.0 };
static GLfloat  default_diffuseIntensity[] = { 1.0, 1.0, 1.0, 1.0 };
static GLfloat  default_specularIntensity[] = { 1.0, 1.0, 1.0, 1.0 };
static GLfloat  default_spotExponent[] = {128.0};


@implementation AUOpenGLTrackballView

#pragma mark ___basicPixelFormat___

+ (NSOpenGLPixelFormat*) basicPixelFormat
{
    NSOpenGLPixelFormatAttribute attributes [] = {
        NSOpenGLPFAWindow,
        NSOpenGLPFADoubleBuffer,	// double buffered
        NSOpenGLPFADepthSize, 
		(NSOpenGLPixelFormatAttribute) 16, // 16 bit depth buffer
        (NSOpenGLPixelFormatAttribute) nil
    };
    return [[[NSOpenGLPixelFormat alloc] initWithAttributes:attributes] autorelease];
}


#pragma mark ___init___

- (void) initValues
{	
	memcpy(light_pos,			default_light_pos,			4*sizeof(GLfloat));
	memcpy(ambientIntensity,	default_ambientIntensity,	4*sizeof(GLfloat));
	memcpy(diffuseIntensity,	default_diffuseIntensity,	4*sizeof(GLfloat));
	memcpy(specularIntensity,	default_specularIntensity,	4*sizeof(GLfloat));
	memcpy(spotExponent, default_spotExponent, 1*sizeof(GLfloat));

	gDollyPanStartPoint[0] = gDollyPanStartPoint[1] = 0;
	gTrackBallRotation[0] =  gTrackBallRotation[1] = gTrackBallRotation[2] = gTrackBallRotation[3] = 0.0f;

	gDolly = GL_FALSE;
	gPan = GL_FALSE;
	gTrackball = GL_FALSE;
	
	gTrackingViewInfo = NULL;
}

- (id)initWithFrame:(NSRect)frame pixelFormat: (NSOpenGLPixelFormat*) pf
{
    self = [super initWithFrame:frame pixelFormat: pf];
    if (self) {
       [self initValues];
    }
    return self;
}

#pragma mark ___MouseFunctions___

// move camera in z axis
-(void)mouseDolly: (NSPoint) location
{
	GLfloat dolly = (gDollyPanStartPoint[1] -location.y) * -camera.viewPos.z / 300.0f;
	camera.viewPos.z += dolly;
	if (camera.viewPos.z == 0.0) // do not let z = 0.0
		camera.viewPos.z = 0.0001;
	gDollyPanStartPoint[0] = (GLint) location.x;
	gDollyPanStartPoint[1] = (GLint) location.y;
}
		
// move camera in x/y plane
- (void)mousePan: (NSPoint) location
{
	GLfloat panX = (gDollyPanStartPoint[0] - location.x) / (900.0f / -camera.viewPos.z);
	GLfloat panY = (gDollyPanStartPoint[1] - location.y) / (900.0f / -camera.viewPos.z);
	camera.viewPos.x -= panX;
	camera.viewPos.y -= panY;
	gDollyPanStartPoint[0] = (GLint) location.x;
	gDollyPanStartPoint[1] = (GLint) location.y;
}

- (void) mouseDown:(NSEvent *)theEvent // trackball
{
	if ([theEvent modifierFlags] & NSControlKeyMask) // send to pan
		[self rightMouseDown:theEvent];
	else if ([theEvent modifierFlags] & NSAlternateKeyMask) // send to dolly
		[self otherMouseDown:theEvent];
	else {
		NSPoint location = [self convertPoint:[theEvent locationInWindow] fromView:nil];
		location.y = camera.viewHeight - location.y;
		gDolly = GL_FALSE; // no dolly
		gPan = GL_FALSE; // no pan
		gTrackball = GL_TRUE;
		startTrackball ((long int) location.x, (long int) location.y, 0, 0, camera.viewWidth, camera.viewHeight);
		gTrackingViewInfo = self;
	}
}

- (void)rightMouseDown:(NSEvent *)theEvent // pan
{
	NSPoint location = [self convertPoint:[theEvent locationInWindow] fromView:nil];
	location.y = camera.viewHeight - location.y;
	if (gTrackball) { // if we are currently tracking, end trackball
		if (gTrackBallRotation[0] != 0.0)
			addToRotationTrackball (gTrackBallRotation, worldRotation);
		
		gTrackBallRotation [0] = gTrackBallRotation [1] = gTrackBallRotation [2] = gTrackBallRotation [3] = 0.0f;
	}
	gDolly = GL_FALSE; // no dolly
	gPan = GL_TRUE; 
	gTrackball = GL_FALSE; // no trackball
	gDollyPanStartPoint[0] = (GLint) location.x;
	gDollyPanStartPoint[1] = (GLint) location.y;
	gTrackingViewInfo = self;
}

- (void)otherMouseDown:(NSEvent *)theEvent //dolly
{
	NSPoint location = [self convertPoint:[theEvent locationInWindow] fromView:nil];
	location.y = camera.viewHeight - location.y;
	if (gTrackball) { // if we are currently tracking, end trackball
		if (gTrackBallRotation[0] != 0.0)
			addToRotationTrackball (gTrackBallRotation, worldRotation);
		gTrackBallRotation [0] = gTrackBallRotation [1] = gTrackBallRotation [2] = gTrackBallRotation [3] = 0.0f;
	}
	gDolly = GL_TRUE;
	gPan = GL_FALSE; // no pan
	gTrackball = GL_FALSE; // no trackball
	gDollyPanStartPoint[0] = (GLint) location.x;
	gDollyPanStartPoint[1] = (GLint) location.y;
	gTrackingViewInfo = self;
}

- (void)mouseUp:(NSEvent *)theEvent
{
	if (gDolly) { // end dolly
		gDolly = GL_FALSE;
	} else if (gPan) { // end pan
		gPan = GL_FALSE;
	} else if (gTrackball) { // end trackball
		gTrackball = GL_FALSE;
		if (gTrackBallRotation[0] != 0.0)
			addToRotationTrackball (gTrackBallRotation, worldRotation);
		gTrackBallRotation [0] = gTrackBallRotation [1] = gTrackBallRotation [2] = gTrackBallRotation [3] = 0.0f;
	} 
	gTrackingViewInfo = NULL;
}

- (void)rightMouseUp:(NSEvent *)theEvent
{
	[self mouseUp:theEvent];
}

- (void)otherMouseUp:(NSEvent *)theEvent
{
	[self mouseUp:theEvent];
}

- (void)mouseDragged:(NSEvent *)theEvent
{
	NSPoint location = [self convertPoint:[theEvent locationInWindow] fromView:nil];
	location.y = camera.viewHeight - location.y;
	if (gTrackball) {
		rollToTrackball ((long int) location.x, (long int) location.y, gTrackBallRotation);
		[self setNeedsDisplay: YES];
	} else if (gDolly) {
		[self mouseDolly: location];
		[self updateProjection];  // update projection matrix (not normally done on draw)
		[self setNeedsDisplay: YES];
	} else if (gPan) {
		[self mousePan: location];
		[self setNeedsDisplay: YES];
	}
}

- (void)scrollWheel:(NSEvent *)theEvent
{
	float wheelDelta = [theEvent deltaX] +[theEvent deltaY] + [theEvent deltaZ];
	if (wheelDelta)
	{
		GLfloat deltaAperture = wheelDelta * -camera.aperture / 200.0f;
		camera.aperture += deltaAperture;
		if (camera.aperture < 0.1) // do not let aperture <= 0.1
			camera.aperture = 0.1;
		if (camera.aperture > 179.9) // do not let aperture >= 180
			camera.aperture = 179.9;
		[self updateProjection]; // update projection matrix
		[self setNeedsDisplay: YES];
		
		//printf("camera.aperture: %lf\n", camera.aperture);
	}
}

- (void)rightMouseDragged:(NSEvent *)theEvent
{
	[self mouseDragged: theEvent];
}

- (void)otherMouseDragged:(NSEvent *)theEvent
{
	[self mouseDragged: theEvent];
}

- (void) keyPress:(NSEvent*)theEvent
{
	printf("key pressed");	
}


#pragma mark ___Update___

- (void) updateProjection
{
	GLdouble ratio, radians, wd2;
	GLdouble left, right, top, bottom, near, far;

	[[self openGLContext] makeCurrentContext];

	// set projection
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	near = -camera.viewPos.z - shapeSize * 0.5;
	if (near < 0.00001)
		near = 0.00001;
	far = -camera.viewPos.z + shapeSize * 0.5;
	if (near < 1.0)
		near = 1.0;
	radians = 0.0174532925 * camera.aperture / 2.; 
	// half aperture degrees to radians 
	wd2 = near * tan(radians);
	ratio = camera.viewWidth / (Float32) camera.viewHeight;
	if (ratio >= 1.0) {
		left  = -ratio * wd2;
		right = ratio * wd2;
		top = wd2;
		bottom = -wd2;	
	} else {
		left  = -wd2;
		right = wd2;
		top = wd2 / ratio;
		bottom = -wd2 / ratio;	
	}
	glFrustum (left, right, bottom, top, near, far);
	
}

- (void) updateModelView
{
    [[self openGLContext] makeCurrentContext];
	
	// move view
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity ();
	gluLookAt (camera.viewPos.x, camera.viewPos.y, camera.viewPos.z,
			   camera.viewPos.x + camera.viewDir.x,
			   camera.viewPos.y + camera.viewDir.y,
			   camera.viewPos.z + camera.viewDir.z,
			   camera.viewUp.x, camera.viewUp.y ,camera.viewUp.z);
			
	// if we have trackball rotation to map (this IS the test I want as it can be explicitly 0.0f)
	if ((gTrackingViewInfo == self) && gTrackBallRotation[0] != 0.0f) 
		glRotatef (gTrackBallRotation[0], gTrackBallRotation[1], gTrackBallRotation[2], gTrackBallRotation[3]);
	else {
	}
	// accumlated world rotation via trackball
	glRotatef (worldRotation[0], worldRotation[1], worldRotation[2], worldRotation[3]);

	
#ifdef DEBUG_PRINT	GLdouble m[16];
	GLdouble p[16];
	
	glGetDoublev( GL_MODELVIEW_MATRIX, m );
	glGetDoublev( GL_PROJECTION_MATRIX, p );

	printf("Projection Matrix\n");
	printf("%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\n",
		p[0],p[4],p[8],p[12],p[1],p[5],p[9],p[13],p[2],p[6],p[10],p[14],p[3],p[7],p[11],p[15]);
	
	printf("ModelView Matrix\n");
	printf("%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\n",
		m[0],m[4],m[8],m[12],m[1],m[5],m[9],m[13],m[2],m[6],m[10],m[14],m[3],m[7],m[11],m[15]);
#endif	

}

- (void) resizeGL
{
	NSRect rectView = [self bounds];
	
	// ensure camera knows size changed
	if ((camera.viewHeight != rectView.size.height) ||
	    (camera.viewWidth != rectView.size.width)) {
		camera.viewHeight = (GLint) rectView.size.height;
		camera.viewWidth = (GLint) rectView.size.width;
		
		glViewport (0, 0, camera.viewWidth, camera.viewHeight);
		[self updateProjection]; 
	}
}


- (void)update  // moved or resized
{
	NSRect rect;
	
	[super update];
	
	[[self openGLContext] makeCurrentContext];
	[[self openGLContext] update];
	
	rect = [self bounds];
		
	[self setNeedsDisplay:true];
}

- (void)reshape	// scrolled, moved or resized
{
	NSRect rect;
	
	[super reshape];
	
	[[self openGLContext] makeCurrentContext];
	[[self openGLContext] update];
	
	rect = [self bounds];
		
	[self setNeedsDisplay:true];
}

#pragma mark ___ResponderFuncs___

- (BOOL)acceptsFirstResponder
{
  return YES;
}

- (BOOL)becomeFirstResponder
{
  return  YES;
}

- (BOOL)resignFirstResponder
{
  return YES;
}


@end
