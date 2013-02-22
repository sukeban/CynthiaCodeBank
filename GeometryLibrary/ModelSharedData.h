/*
 *  ModelSharedData.h
 *  VISynth
 *
 *  Created by Cynthia Bruyns on 9/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ModelSharedData_H
#define ModelSharedData_H

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*!
 an object that holds the eigendata for an object
 these eigenvalues represent the resonant frequencies and these eignvectors represent the mode shapes for the object
 */
class EigenData {

public:
		
	double* 	EigenVectors; // changed to doubles because parsing old files required it, some numbers were out of range
	
	double* 	EigenValues;	
	
	int			numFreq;
	int			numDOFS;
	int			numDOF;
						
	EigenData(int numF, int numDS, int numD){
	
		numFreq = numF;		// how many modes computed
		numDOFS = numDS;	// how many degrees of freedom in the system
		numDOF = numD;		// how many per node
		
		EigenValues = (double*) calloc(numFreq,sizeof(double));		
		
		EigenVectors = (double*) calloc(numDOFS*numFreq,sizeof(double));
		
	};
	
	void Copy(EigenData* data) {
		
		memcpy(EigenValues, data->EigenValues, numFreq*sizeof(double));
		memcpy(EigenVectors, data->EigenVectors, numDOFS*numFreq*sizeof(double));
	
	};
	
	void SetEigenVectors(double* data){		
		memcpy(EigenVectors, data, numDOFS*numFreq*sizeof(double));			
#ifdef DEBUG_PRINT
		printf("copying %d x %d eigenvectors\n", numDOFS, numFreq);	
#endif
	};
	
	~EigenData(){
		if (EigenValues) free(EigenValues);	
		if (EigenVectors) free(EigenVectors);		
	};
		
};

/*!
 an object that stores the geometric information about an object
 that is the points, faces and normals for use in rendering
 
 this object also defines a set of useful graphics operations like picking on a surface and testing for intersection
 */
class GeomData {
	public:

		double*	vertex_pos_init;	// initial position, good for resetting
		double*	vertex_pos;			// current position, for deforming		
		double*	vertex_normals;		// for wireframe drawing, and other functions		
		int		numverts;			// how many points
				
		int		vertexperface;		// number of points per face
		int*	faces;				// the faces
		double*	face_normals;		// the normal of each face
		bool*	outter_face;		// for the tetra models, i.e. does it face outwards
		int		numfaces;			// how many faces do we have

		int*	tetras;				// the tetrahedra defined by the given verticies
		int		numtetras;			// how many tetrahedra
		
		double	radius;				// the bounding sphere for this object
		double	max[3];				// mins and max dimensions
		double	min[3];
		
		
		GeomData(int numv, int numf, int numt, int nel);
		~GeomData();
		
#pragma mark ___Utility___	
	
		void reset();
		
		void Copy(GeomData* other)
		{			
			memcpy(vertex_pos_init		, other->vertex_pos_init		, numverts*3*sizeof(double));
			memcpy(vertex_pos			, other->vertex_pos				, numverts*3*sizeof(double));
			memcpy(vertex_normals		, other->vertex_normals			, numverts*3*sizeof(double));
			
			memcpy(faces		, other->faces			, numfaces*vertexperface*sizeof(int));
			memcpy(face_normals	, other->face_normals	, numfaces*3*sizeof(double));
			memcpy(outter_face	, other->outter_face	, numfaces*sizeof(bool));
			
			if (numtetras)
				memcpy(tetras, other->tetras, numtetras*4*sizeof(int));
			
			radius = other->radius;
			memcpy(max, other->max, 3*sizeof(double));
			memcpy(min, other->min, 3*sizeof(double));
		
		}
		
#pragma mark ___Faces___		
			
		void getFaceNormal(int i, double* normal);				
		int  getFaceIndex(int node0, int node1, int node2);
		
		int MatchingFace(int i);		
		void FlipFace(int i);
		void AlignFaces();
					
		void ComputeFaceCenter(int i, double* center);
		void ComputeFaceNormal(int i, double* normal);		

		void ComputeFaceNormals();
		void ComputeOutterFaces();
				
#pragma mark ___Verticies___

		void ComputeVertexNormals();	
		
		void RemoveStrandedVerticies();
				
#pragma mark ___Tetras___	
		bool insideTetrahedra(int i, double* intPt);

		void ComputeTetraCenter(int i, double* center);		
		bool TetraFaceIsPointingInwards(int ti, int fi);
			
		void CleanTetraModel();
		
#pragma mark ___Picking___
		bool insideTriangle(int i, double* intPt);
		void determinePlaneEquation(int i, double &a, double &b, double &c, double &d);

		bool checkFace(int i, double* P1, double* P2, double* intPt);
		bool checkClosest(double* intPt, double* cIntPt, double* P1);
		
		int pickedOnSurface(double* P1, double* P2, double* intPt);

};



#endif