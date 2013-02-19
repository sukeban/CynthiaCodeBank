/*
 *  OpenGLSharedData.h
 *  VISynth
 *
 *  Created by Cynthia Bruyns on 9/8/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef OpenGLSharedData_H_
#define OpenGLSharedData_H_

#import <OpenGL/gl.h>
#import <OpenGL/glext.h>
#import <OpenGL/glu.h>


typedef struct {
	GLdouble x,y,z;
} recVec;

typedef struct {
	recVec		viewPos;				// View position
	recVec		viewDir;				// View direction vector
	recVec		viewUp;					// View up direction
	recVec		rotPoint;				// Point to rotate about
	GLdouble	aperture;				// pContextInfo->camera aperture
	GLint		viewWidth, viewHeight;	// current window/screen height and width
} recCamera;


#endif