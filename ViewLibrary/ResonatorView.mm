//
//  ResonatorView.mm
//  VISynth
//
//  Created by Cynthia Bruyns on 9/8/06.
//  Copyright 2006 __MyCompanyName__. All rights reserved.
//

#import "ResonatorView.h"

#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iostream>

#include <GLUT/glut.h>

#import "MatrixOps.h"


@implementation ResonatorView

#pragma mark ___Thump___

- (void) setNumThumpers: (int) num
{
	 numthumpers = num;
}

- (void) updateGeometry:(NSTimer *)t
{	
	numtimesteps++;
//	printf("update: %d\n", numtimesteps);
	if (numtimesteps % 10 == 0) {
		justthumped = false;
//		printf("clearing thump: %d\n", numtimesteps);
		numtimesteps = 0;
		
		// cbnote stop the timer
	}
	
	[self setNeedsDisplay: YES];
		
}

- (NSTimer *) timer {
	return [[timer retain] autorelease];
}

- (void) setTimer: (NSTimer *) value {
	if ( timer != value ) {
		[timer release];
		timer = [value retain];
	}
}

-(void) Thump: (int) intValue
{
//	curThump = intValue%geomdata->thumpnum;
	curThump = intValue%numthumpers;
	//printf("curThump: %d\n", curThump);
	
	justthumped = true;
		
/*	[self setTimer: [NSTimer scheduledTimerWithTimeInterval: (1.0/160.0)
											 target: self
											 selector: @selector(updateGeometry:)
											 userInfo: nil
											repeats: YES]];	
*/	[self setNeedsDisplay: YES];
}

- (void) modulate:(int) idx
{
	if (!geomdata || !eigendata) return;
		
	//cerr << "hertzTextField["<< idx << "]:" << eigendata->EigenValues[idx] << endl;
	//cerr << " numverts: " << geomdata->numverts << endl;
	//cerr << "numdof: " << eigendata->numDOF << endl;
	//cerr << "numdofs: " << eigendata->numDOFS << endl;
	
	int numverts = geomdata->numverts;
	
	int numDOFS = eigendata->numDOFS;
	int numDOF = eigendata->numDOF;
	
	double scale = geomdata->radius/10.0; 
	if (deformScaleSlider)
		scale = [deformScaleSlider doubleValue];
		
	double off[3];
	int i;
	for (i=0; i< numverts; i++) {
		
		off[0] = eigendata->EigenVectors[idx * numDOFS + i*numDOF    ];
		off[1] = eigendata->EigenVectors[idx * numDOFS + i*numDOF + 1];
		off[2] = eigendata->EigenVectors[idx * numDOFS + i*numDOF + 2];
		
		geomdata->vertex_pos[i + 0			] = geomdata->vertex_pos_init[i + 0		    ] + scale*off[0];
		geomdata->vertex_pos[i + 1*numverts	] = geomdata->vertex_pos_init[i + 1*numverts] + scale*off[1];
		geomdata->vertex_pos[i + 2*numverts	] = geomdata->vertex_pos_init[i + 2*numverts] + scale*off[2];	

																										
	}
	[self setNeedsDisplay: YES];
}

#pragma mark ___Reset___

- (void) resetObject
{	
	picked_face = -1;
	if (!geomdata) return;
	geomdata->reset();
	[self setNeedsDisplay: YES];
}

- (void) resetCamera
{
	if (geomdata)
		shapeSize = 10*geomdata->radius; 
	else
		shapeSize = 0.0;
	
	camera.aperture = 6.0;  

	camera.rotPoint.x = 0.0;
	camera.rotPoint.y = 0.0;
	camera.rotPoint.z = 0.0;

	camera.viewPos.x = 0.0;
	camera.viewPos.y = 0.0;
	camera.viewPos.z = -2*shapeSize;
	
	camera.viewDir.x = -camera.viewPos.x; 
	camera.viewDir.y = -camera.viewPos.y; 
	camera.viewDir.z = -camera.viewPos.z;

	camera.viewUp.x = 0;  
	camera.viewUp.y = 1; 
	camera.viewUp.z = 0;

	light_pos[0] = 0.0;
	light_pos[1] = 0.0;
	light_pos[2] = 0.0;
	
	if (!geomdata) return;

	light_pos[0] = geomdata->radius;
	light_pos[1] = geomdata->radius;
	light_pos[2] = geomdata->radius;

}


#pragma mark ___UICallbacks__

- (IBAction) redraw: (id) sender
{
	[self setNeedsDisplay: YES];
}

- (IBAction)  resetView: (id) sender
{
	[self resetCamera];
	
	[self updateProjection];
	
	[self resetObject];

	[self prepareOpenGL];
	[self setNeedsDisplay: YES];	
}

- (IBAction) setWireframeOffset:(id)sender
{
	UInt32 index = [sender intValue];
	if (sender == wireframeOffsetSlider){
		[wireframeOffsetTextField setIntValue: index];
	}
	else{
		[wireframeOffsetSlider setIntValue: index];
	}
		
		
	[self redraw: self];
}
- (IBAction) setNormalScale:(id)sender
{
	UInt32 index = [sender intValue];
	if (sender == normalSizeSlider){
		[normalSizeTextField setIntValue: index];
	}
	else{
		[normalSizeSlider setIntValue: index];
	}
		
		
	[self redraw: self];
}
- (IBAction) setThumperScale:(id)sender
{
	double index = [sender doubleValue];
	if (sender == thumperSizeSlider){
		[thumperSizeTextField setDoubleValue: index];
	}
	else{
		[thumperSizeSlider setDoubleValue: index];
	}
		
		
	[self redraw: self];
}


- (IBAction) setCurMode:(id) sender
{
	if (!eigendata) return;
	
	UInt32 index = [sender intValue];
	if (sender == curModeSlider){
		[curModeTextField setIntValue: index];
	}
	else{
		[curModeSlider setIntValue: index];
	}
		
	if (!eigendata->EigenValues) return;
	if (hertzTextField) {
		float v = sqrtf(fabs(eigendata->EigenValues[index])/2/M_PI);
		NSNumber* value = [[NSNumber alloc] initWithFloat: v]; 
		[hertzTextField setFloatValue: [value doubleValue]];
	}
	
	[self modulate: index];

//	printf("setting cur mode to: %ld\n", index);
}


- (IBAction) resetThumpers: (id) sender
{
	numthumpers = 0;
	[myController resetThumpers];
	[self setNeedsDisplay: YES];	   
}

- (IBAction) setDeformScale: (id) sender
{
	double scale = [sender doubleValue];
	if (sender == deformScaleSlider){
		[deformScaleTextField setFloatValue: scale];
	}
	else{
		[deformScaleSlider setFloatValue: scale];
	}
	
	[self modulate: [curModeSlider intValue]];

	[self setNeedsDisplay: YES];
}

#pragma mark ___Intersection___

- (void) setPickDirection: (NSPoint) point
{
	
	GLint viewport[4];
	GLdouble m[16];
	GLdouble p[16];
	GLfloat winX, winY, winZ;
	GLdouble posX, posY, posZ;
	
	glGetDoublev( GL_MODELVIEW_MATRIX, m );
	glGetDoublev( GL_PROJECTION_MATRIX, p );
	glGetIntegerv( GL_VIEWPORT, viewport );
	
#ifdef DEBUG_PRINT
	printf("HIT TEST\n");
	printf("Projection Matrix\n");
	printf("%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\n",
		p[0],p[4],p[8],p[12],p[1],p[5],p[9],p[13],p[2],p[6],p[10],p[14],p[3],p[7],p[11],p[15]);
	
	printf("ModelView Matrix\n");
	printf("%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\n%f\t%f\t%f\t%f\n",
		m[0],m[4],m[8],m[12],m[1],m[5],m[9],m[13],m[2],m[6],m[10],m[14],m[3],m[7],m[11],m[15]);
	
	printf("View Port: %d %d %d %d\n", viewport[0], viewport[1], viewport[2], viewport[3]);
#endif	

	winX = point.x;
	winY = point.y;  // NSOpenGLView does not need to do this
	glReadPixels( (int) point.x, (int) winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
	gluUnProject( winX, winY, winZ, m, p, viewport, &posX, &posY, &posZ); // this is a point in 3d space
	
	pickOrigin[0] = camera.viewPos.x;
	pickOrigin[1] = camera.viewPos.y;
	pickOrigin[2] = -camera.viewPos.z;
	
	double po[3];
	po[0] = pickOrigin[0];
	po[1] = pickOrigin[1];
	po[2] = pickOrigin[2];
			
	double oRot[3];
	v3m3_mult(po, m, oRot);
	pickOrigin[0] = oRot[0];
	pickOrigin[1] = oRot[1];
	pickOrigin[2] = oRot[2];

	double vo[3]; 
	vo[0] = posX;
	vo[1] = posY;
	vo[2] = posZ; 
	
	double v[3];
	v3_sub(vo, pickOrigin,v);
	//v3_normalize(v);
	
	pickDirection[0] = v[0]*10;
	pickDirection[1] = v[1]*10;
	pickDirection[2] = v[2]*10;
			
}

- (int) getIntersectionPoint: (double*) intPt withWindowPoint: (NSPoint) point
{
	[self setPickDirection: point];

	int face = geomdata->pickedOnSurface(pickOrigin, pickDirection, intPt);
	return face;
}

#pragma mark ___Mouse/Keyboard Events___

- (bool) madeThumper: (NSPoint) point
{
	bool hit = false;
	double intPt[3];
	picked_face = [self getIntersectionPoint: intPt withWindowPoint: point];

	if (picked_face > -1){
		printf("you hit a face\n" );
		
		if (numthumpers < MAXNUMTHUMPERS) {

			thumpers[numthumpers					] = intPt[0]; 
			thumpers[numthumpers+ MAXNUMTHUMPERS	] = intPt[1];
			thumpers[numthumpers+2*MAXNUMTHUMPERS	] = intPt[2]; 
			
			printf("at %lf %lf %lf\n", intPt[0], intPt[1], intPt[2]);

			[myController addThumper: geomdata->faces[picked_face]]; // cbnote the first node in the face, not an interpolated location
			numthumpers++;
		
		}
		
		hit = true;
	}
	return hit;
}

- (bool) moveThumper: (NSPoint) point
{
	bool hit = false;
	double intPt[3];
	picked_face = [self getIntersectionPoint: intPt withWindowPoint: point];
	int num  = numthumpers-1;
	if (picked_face > -1){
		printf("you hit a face\n" );
		
		thumpers[num				] = intPt[0]; 
		thumpers[num+ MAXNUMTHUMPERS] = intPt[1];
		thumpers[num+2*MAXNUMTHUMPERS] = intPt[2]; 

		[myController moveThumper: geomdata->faces[picked_face]];	 // cbnote the first node in the face, not an interpolated location	
				
		hit = true;
	}
	return hit;
}

- (void)thumperMouseDown:(NSEvent *)theEvent 
{
	// convert event location to our coordinate system
	NSPoint point = [self convertPoint:[theEvent locationInWindow] fromView:nil];
	//printf("point: %lf %lf\n", point.x, point.y);
	
	// only take notice if the user clicked in the shape
	if ( [self  madeThumper: point] ) {
		[self setNeedsDisplay: YES];
	}	
}


- (void)thumperMouseDragged:(NSEvent *)theEvent 
{
	// convert event location to our coordinate system
	NSPoint point = [self convertPoint:[theEvent locationInWindow] fromView:nil];
	//printf("point: %lf %lf\n", point.x, point.y);
	
	// only take notice if the user clicked in the shape
	if ( [self  moveThumper: point] ) {
		[self setNeedsDisplay: YES];
	}	
}


#pragma mark -
#pragma mark ___Drawing___


- (void) turnOnLight
{		
	glFrontFace(GL_CCW);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 

	glEnable(GL_LIGHTING);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, default_ambientIntensity);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, default_diffuseIntensity);
	glLightfv(GL_LIGHT0, GL_SPECULAR, default_specularIntensity);
	glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, default_spotExponent[0]);
	glLightfv(GL_LIGHT0, GL_POSITION, default_light_pos);
	
	// we provide normals ourselves for everything
	glEnable(GL_RESCALE_NORMAL);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glShadeModel(GL_SMOOTH);    

	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

	glColor4f(mat_ambient[0],mat_ambient[1],mat_ambient[2],mat_ambient[3]) ; //cbnote needs colorwell

	glMaterialfv(GL_BACK, GL_AMBIENT, back_mat_ambient);
	glMaterialfv(GL_BACK, GL_DIFFUSE, back_mat_diffuse);
	glMaterialfv(GL_BACK, GL_SPECULAR, back_mat_specular);
	glMaterialfv(GL_BACK, GL_SHININESS, back_mat_shininess);
					
	//glEnable(GL_LINE_SMOOTH);

}



- (void) clearBuffer
{
	glClearColor( clear[0], clear[1], clear[2], 1 ) ;	
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);	
}


- (void) drawAxes
{	
	double xlen = geomdata->radius*2;
	double ylen = xlen;
	double zlen = ylen;
	
	glDisable(GL_LIGHTING);

	//-------------x	
	glColor3f(1,0,0);
	glBegin(GL_LINES);
	glVertex3d(xlen, 0.0, 0.0);
	glVertex3d(0.0, 0.0, 0.0);
	glEnd();
	
	
	//-------------y-
	glColor3f(0,1,0);
	glBegin(GL_LINES);
	glVertex3d(0.0, ylen, 0.0);
	glVertex3d(0.0, 0.0, 0.0);
	glEnd();
	

	//-------------z
	glColor3f(0,0,1);
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0,0.0);
	glVertex3d(0.0, 0.0,zlen);
	glEnd();	
	
	
	glEnable(GL_LIGHTING);

}

-(void) drawPickingRay
{
	//if (!justthumped) return;

	glDisable(GL_LIGHTING);
	glColor3f(0,1,1);  
	glBegin(GL_LINES);

	glVertex3d(pickOrigin[0],pickOrigin[1],pickOrigin[2]);
	glVertex3d(pickDirection[0], pickDirection[1], pickDirection[2]);
#ifdef DEBUG_PRINT	
	printf("pick pickOrigin: %lf %lf %lf\n", pickOrigin[0], pickOrigin[1], pickOrigin[2]);
	printf("pick pickDirection: %lf %lf %lf\n", pickDirection[0], pickDirection[1], pickDirection[2]);
#endif		
	glEnd();
	glEnable(GL_LIGHTING);
	
}

-(void) drawThumpers
{
#ifdef DEBUG_PRINT
	printf("numthumpers: %d\n", numthumpers);
#endif	
	int i;
	for ( i = 0; i < numthumpers; i++) {
		if (justthumped && (i == curThump)) {
			glColor3f(1,0,0); 
		}
		else {
			glColor3f(0,1,0);  
		}
		double inter_rad = geomdata->radius/50.0;
		if (thumperSizeSlider) inter_rad = [thumperSizeSlider doubleValue];			
		glMatrixMode(GL_MODELVIEW);	
		glPushMatrix();
		
#ifdef DEBUG_RINT
		printf("drawing at %lf %lf %lf\n", thumpers[i],thumpers[i+MAXNUMTHUMPERS], thumpers[i+2*MAXNUMTHUMPERS]);
#endif
		glTranslatef(thumpers[i],thumpers[i+MAXNUMTHUMPERS], thumpers[i+2*MAXNUMTHUMPERS]);
		glutSolidSphere(GLdouble(inter_rad), 16, 16);
		glPopMatrix();
	}	
}

- (void) drawNormals
{
	double center[3];
	double normal[3];
	double pt[3];
	
	double scale =  geomdata->radius/100.0;
	if (normalSizeSlider) scale = [normalSizeSlider floatValue];
	
	glDisable(GL_LIGHTING);
	
	int i;
	
	for ( i=0; i<geomdata->numfaces; i++){
	
		if (!geomdata->outter_face[i])
			continue; 

		geomdata->ComputeFaceCenter(i,center);
		geomdata->getFaceNormal(i,normal);
		
		pt[0] = center[0] + scale*normal[0];
		pt[1] = center[1] + scale*normal[1];
		pt[2] = center[2] + scale*normal[2];
		
		//-------------x	
		glColor3f(1,0,1);
		glBegin(GL_LINES);
		glVertex3d(center[0], center[1], center[2]);
		glVertex3d(	pt[0], 
					pt[1], 
					pt[2]);
		
#ifdef DEBUG_PRINT
	printf("drawing line from: %lf %lf %lf to %lf %lf %lf\n",
		center[0], center[1], center[2],
		pt[0], pt[1], pt[2]);
#endif
		
		glEnd();
	
	}
	
	glEnable(GL_LIGHTING);

}


-(void) drawWireframe
{	
	int numverts = geomdata->numverts;
	int numfaces = geomdata->numfaces;
	int vertexperface = geomdata->vertexperface;
	
	glDisable(GL_LIGHTING);
	
		
	double scale = 1.0*(geomdata->radius/100.0);
	if (wireframeOffsetTextField) [wireframeOffsetTextField doubleValue];	
	int i;
	for (i=0; i<geomdata->numfaces; i++){
	
		if (i == picked_face){
			glColor3f(1,0,1);  
			glLineWidth(2.0);
		}
		else{
			glColor3f(0,0,0);		
			glLineWidth(1.0);
		}
	
		if (!geomdata->outter_face[i])
			continue; 

		glBegin(GL_LINE_STRIP);
		glVertex3d(
				    scale*geomdata->vertex_normals[geomdata->faces[i]				] + geomdata->vertex_pos[geomdata->faces[i]	              ],
				    scale*geomdata->vertex_normals[geomdata->faces[i]+numverts	] + geomdata->vertex_pos[geomdata->faces[i] +	  numverts],
				    scale*geomdata->vertex_normals[geomdata->faces[i]+2*numverts	] + geomdata->vertex_pos[geomdata->faces[i] +	2*numverts]);
		
		glVertex3d(
				   scale*geomdata->vertex_normals[geomdata->faces[i + numfaces]			]+ geomdata->vertex_pos[geomdata->faces[i + numfaces]				],
				   scale*geomdata->vertex_normals[geomdata->faces[i + numfaces]+numverts	]+ geomdata->vertex_pos[geomdata->faces[i + numfaces] +	numverts	],
				   scale*geomdata->vertex_normals[geomdata->faces[i + numfaces]+2*numverts]+ geomdata->vertex_pos[geomdata->faces[i + numfaces] + 2*numverts]);
		
		glVertex3d(
				    scale*geomdata->vertex_normals[geomdata->faces[i + 2*numfaces]			]+ geomdata->vertex_pos[geomdata->faces[i + 2*numfaces]				],
				    scale*geomdata->vertex_normals[geomdata->faces[i + 2*numfaces]+numverts	]+ geomdata->vertex_pos[geomdata->faces[i + 2*numfaces] +	numverts],
				    scale*geomdata->vertex_normals[geomdata->faces[i + 2*numfaces]+2*numverts	]+ geomdata->vertex_pos[geomdata->faces[i + 2*numfaces] + 2*numverts]);
					   
		if (vertexperface >3 ){
				glVertex3d(
				    scale*geomdata->vertex_normals[geomdata->faces[i + 3*numfaces]			]+ geomdata->vertex_pos[geomdata->faces[i + 3*numfaces]				],
				    scale*geomdata->vertex_normals[geomdata->faces[i + 3*numfaces]+numverts	]+ geomdata->vertex_pos[geomdata->faces[i + 3*numfaces] + 1*numverts],
				    scale*geomdata->vertex_normals[geomdata->faces[i + 3*numfaces]+2*numverts	]+ geomdata->vertex_pos[geomdata->faces[i + 3*numfaces] + 2*numverts]);
		}
		
		glVertex3d(
				    scale*geomdata->vertex_normals[geomdata->faces[i]				] + geomdata->vertex_pos[geomdata->faces[i]	              ],
				    scale*geomdata->vertex_normals[geomdata->faces[i]+numverts	] + geomdata->vertex_pos[geomdata->faces[i] +	  numverts],
				    scale*geomdata->vertex_normals[geomdata->faces[i]+2*numverts	] + geomdata->vertex_pos[geomdata->faces[i] +	2*numverts]);
					   
		glEnd();
	}
	
	glEnable(GL_LIGHTING);	
}



-(void) drawGeometry
{ 	 
	[self turnOnLight];
	
	int numverts = geomdata->numverts;
	int numfaces = geomdata->numfaces;

	//printf("geom.numfaces = %d\n", geomdata->numfaces);
	int vertexperface = geomdata->vertexperface;
	int i;
	 for (i=0; i<geomdata->numfaces; i++){

		if (!geomdata->outter_face[i])
			continue; // dont add in interior faces

		glBegin(GL_POLYGON);

		glNormal3f(
			geomdata->vertex_normals[geomdata->faces[i]				],
			geomdata->vertex_normals[geomdata->faces[i] + numverts	],
			geomdata->vertex_normals[geomdata->faces[i] + 2*numverts]);
		glVertex3d(
			geomdata->vertex_pos[geomdata->faces[i]				],
			geomdata->vertex_pos[geomdata->faces[i] + numverts	],
			geomdata->vertex_pos[geomdata->faces[i] + 2*numverts]);

		glNormal3f(
			geomdata->vertex_normals[geomdata->faces[i + numfaces]				],
			geomdata->vertex_normals[geomdata->faces[i + numfaces] + numverts	],
			geomdata->vertex_normals[geomdata->faces[i + numfaces] + 2*numverts	]);
		glVertex3d(
			geomdata->vertex_pos[geomdata->faces[i + numfaces]				],
			geomdata->vertex_pos[geomdata->faces[i + numfaces] + numverts	],
			geomdata->vertex_pos[geomdata->faces[i + numfaces] + 2*numverts	]);

		glNormal3f(
			geomdata->vertex_normals[geomdata->faces[i + 2*numfaces]			 ],
			geomdata->vertex_normals[geomdata->faces[i + 2*numfaces] + numverts	 ],
			geomdata->vertex_normals[geomdata->faces[i + 2*numfaces] + 2*numverts]);					
		glVertex3d(
			geomdata->vertex_pos[geomdata->faces[i + 2*numfaces]			 ],
			geomdata->vertex_pos[geomdata->faces[i + 2*numfaces] + numverts	 ],
			geomdata->vertex_pos[geomdata->faces[i + 2*numfaces] + 2*numverts]);		
		
		if (vertexperface >3){
			glNormal3f(
				geomdata->vertex_normals[geomdata->faces[i + 3*numfaces]			 ],
				geomdata->vertex_normals[geomdata->faces[i + 3*numfaces] + numverts	 ],
				geomdata->vertex_normals[geomdata->faces[i + 3*numfaces] + 2*numverts]);	
			glVertex3d(
				geomdata->vertex_pos[geomdata->faces[i + 3*numfaces]			 ],
				geomdata->vertex_pos[geomdata->faces[i + 3*numfaces] + numverts	 ],
				geomdata->vertex_pos[geomdata->faces[i + 3*numfaces] + 2*numverts]);	
		}
		glEnd();
		
/*		printf("face[%d]: %lf %lf %lf\n %lf %lf %lf \n %lf %lf %lf\n %lf %lf %lf\n\n ", 
		i,
		geomdata->vertex_pos[geomdata->faces[i		     ] ],geomdata->vertex_pos[geomdata->faces[i             ] + numverts],geomdata->vertex_pos[geomdata->faces[i		    ] + 2*numverts],
		geomdata->vertex_pos[geomdata->faces[i + numfaces  ] ],geomdata->vertex_pos[geomdata->faces[i + numfaces  ] + numverts],geomdata->vertex_pos[geomdata->faces[i + numfaces ] + 2*numverts],
		geomdata->vertex_pos[geomdata->faces[i + 2*numfaces] ],geomdata->vertex_pos[geomdata->faces[i + 2*numfaces] + numverts],geomdata->vertex_pos[geomdata->faces[i+ 2*numfaces] + 2*numverts],		
		geomdata->vertex_pos[geomdata->faces[i + 3*numfaces] ],geomdata->vertex_pos[geomdata->faces[i + 3*numfaces] + numverts],geomdata->vertex_pos[geomdata->faces[i+ 3*numfaces] + 2*numverts]);		
*/		
/*		printf("face[%d]: %lf %lf %lf\n %lf %lf %lf \n %lf %lf %lf\n \n ", 
		i,
		geomdata->vertex_pos[geomdata->faces[i		     ] ],geomdata->vertex_pos[geomdata->faces[i             ] + numverts],geomdata->vertex_pos[geomdata->faces[i		    ] + 2*numverts],
		geomdata->vertex_pos[geomdata->faces[i + numfaces  ] ],geomdata->vertex_pos[geomdata->faces[i + numfaces  ] + numverts],geomdata->vertex_pos[geomdata->faces[i + numfaces ] + 2*numverts],
		geomdata->vertex_pos[geomdata->faces[i + 2*numfaces] ],geomdata->vertex_pos[geomdata->faces[i + 2*numfaces] + numverts],geomdata->vertex_pos[geomdata->faces[i+ 2*numfaces] + 2*numverts]		
		);		
*/		
	}
}



- (void) drawRect:(NSRect)rect
{	
	if (NO == fGLInit) {
		[self prepareOpenGL];
	}
	fGLInit = YES;
	
	[self resizeGL];			// forces projection matrix update (does test for size changes)
	[self updateModelView];		// update model view matrix for object

	[self clearBuffer];

	if (geomdata){
		
		if (!drawFacesButton || [drawFacesButton state] == NSOnState)
			[self drawGeometry];
			
		if (!drawWireframeButton || [drawWireframeButton state] == NSOnState)
			[self drawWireframe];
			
		if (drawNormalsButton && [drawNormalsButton state] == NSOnState)
			[self drawNormals];
				
		if (drawThumpersButton && [drawThumpersButton state] == NSOnState)
			[self drawThumpers];

		if (drawPickingRayButton && [drawPickingRayButton state] == NSOnState)
			[self drawPickingRay];
		
		if (drawAxesButton && [drawAxesButton state] == NSOnState)
			[self drawAxes];
	}
	
	if ([self inLiveResize])
		glFlush ();
	else
		[[self openGLContext] flushBuffer];
}

#pragma mark ___VIEWSTUFF___


- (void) prepareOpenGLContext
{
	 long swapInt = 1;

    [[self openGLContext] setValues:&swapInt forParameter:NSOpenGLCPSwapInterval]; // set to vbl sync
	
	[self resetCamera];
}

- (void) prepareOpenGL
{
	[self prepareOpenGLContext];

	justthumped = false;	
	numtimesteps = 0;
		
	numthumpers = 0;
	
	memcpy(mat_ambient,default_mat_ambient,4*sizeof(GLfloat));
	memcpy(mat_diffuse,default_mat_diffuse,4*sizeof(GLfloat));
	memcpy(mat_specular,default_mat_specular,4*sizeof(GLfloat));
	memcpy(mat_shininess,default_mat_shininess,1*sizeof(GLfloat));

	memcpy(back_mat_ambient,default_back_mat_ambient,4*sizeof(GLfloat));
	memcpy(back_mat_diffuse,default_back_mat_diffuse,4*sizeof(GLfloat));
	memcpy(back_mat_specular,default_back_mat_specular,4*sizeof(GLfloat));
	memcpy(back_mat_shininess,default_back_mat_shininess,1*sizeof(GLfloat));

	pickDirection[0] = pickDirection[1] = pickDirection[2] = 0.0;
	pickOrigin[0] = 0.0; pickOrigin[1] = 0.0; pickOrigin[2] = camera.viewPos.z;
	
	clear[0] = clear[1] = clear[2] = 1.0f;

}

#pragma mark ___init___
/*
- (id)initWithFrame:(NSRect)frame {
    self = [super initWithFrame:frame];
    if (self) {
        // Initialization code here.
    }
	
    return self;
}
*/

-(id) initWithFrame: (NSRect) frameRect
{
	 NSOpenGLPixelFormatAttribute attrs [] = {
        NSOpenGLPFAWindow,
        NSOpenGLPFADoubleBuffer,	
        NSOpenGLPFADepthSize, (NSOpenGLPixelFormatAttribute)16, 
        (NSOpenGLPixelFormatAttribute)nil
    };	
	
	NSOpenGLPixelFormat* pixFmt = [[NSOpenGLPixelFormat alloc] initWithAttributes:attrs];
	if(!pixFmt)
	{
		NSLog(@"No pixel format -- exiting");
		exit(1);
	}
	
	self = [super initWithFrame: frameRect pixelFormat: pixFmt];
	if (self){
		//geomdata = NULL; eigendata = NULL;
		picked_face = -1;
	}

	return self;
}


- (void) awakeFromNib
{	

}

#pragma mark ___dealloc___

- (void)dealloc 
{
	[[self timer] invalidate];
	[self setTimer: nil];
		
	[super dealloc];
}

#pragma mark ___Communicating With Controller___

- (void) setMaxFreq: (UInt32) maxF
{
	[curModeSlider setMaxValue: maxF];
	[curModeStaticText setIntValue: maxF];
	//printf("max mode is: %ld\n", maxF);
}

- (void) setGeomData: (GeomData*) data 
{
	geomdata = data;
	
	geomdata->RemoveStrandedVerticies();

	geomdata->ComputeOutterFaces();
	if (geomdata->numtetras){
		geomdata->CleanTetraModel();	// this order always
	}
	geomdata->ComputeFaceNormals();
	geomdata->ComputeVertexNormals();


	[self resetView: self];		
	[self setNeedsDisplay: YES];
}

- (void) swapInGeomData: (GeomData*) data
{
	geomdata = data;
		
	[self setNeedsDisplay: YES];
}

- (void) setEigenData: (EigenData*) data
{
	numthumpers = 0;
	eigendata = data;
	[self setMaxFreq: (eigendata->numFreq-1)];
}

- (void) swapInEigenData: (EigenData*) data
{
	eigendata = data;
		
	[self setNeedsDisplay: YES];
}

@end
