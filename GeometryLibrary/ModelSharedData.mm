/*
 *  ModelSharedData.mm
 *  ShapeSliders
 *
 *  Created by Cynthia Bruyns on 9/30/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "ModelSharedData.h"
#include "MatrixOps.h"

#include "SharedFEData.h"
#include "GeneralFE.h"


#define EPSILON				1e-6

GeomData::GeomData(int numv, int numf, int numt, int nel)
{
	numverts = numv;
	numfaces = numf;
	numtetras = numt;
	
	vertexperface = nel;
	
	radius = 
	min[0] = min[1] = min[2] = 
	max[0] = max[1] = max[2] = 0.0;

	vertex_pos_init = (double*) calloc(3*numverts,sizeof(double)); 
	vertex_pos = (double*) calloc(3*numverts,sizeof(double)); 
	vertex_normals = (double*) calloc(3*numverts,sizeof(double));
	
	faces = (int*) calloc(vertexperface*numfaces,sizeof(int)); 
	face_normals = (double*) calloc(3*numfaces, sizeof(double));
	outter_face = (bool*) calloc(numfaces, sizeof(bool));
	
	tetras = NULL;
	if (numtetras){
		tetras = (int*) calloc(4*numtetras, sizeof(int));
	}														
}

GeomData::~GeomData(){

	if (vertex_pos_init) free(vertex_pos_init);
	if (vertex_pos) free(vertex_pos);
	if (vertex_normals) free(vertex_normals);
	
	if (faces) free (faces);	
	if (face_normals) free(face_normals);
	if (outter_face) free(outter_face);
	
	if (tetras) free(tetras);
}

#pragma mark ___Utility___		

void GeomData::reset()
{
	memcpy(vertex_pos, vertex_pos_init, numverts*3*sizeof(double));
}

#pragma mark ___Faces___		

int GeomData::MatchingFace(int i)
{
	// search through all the faces and see if the same set of 3 verticies appears
	
	int node0 = faces[i];
	int node1 = faces[i+numfaces];
	int node2 = faces[i+numfaces*2];
	
	int index = getFaceIndex(node0, node1, node2);
	if (index == i) 
		return -1;
	else 
		return index;
	
}

void GeomData::FlipFace(int i)
{
	/*
			0
			
		1      2	
	*/
	
	/*
	  0     3
	  
	  1     2
	
	*/

	// swap nodes 1 and 2
	int node1 = faces[i+1*numfaces];
	int node2 = faces[i+2*numfaces];
	
	faces[i+1*numfaces] = node2;
	faces[i+2*numfaces] = node1;
	
}

void GeomData::ComputeFaceCenter(int i, double* center)
{	
	int j;
	center[0] = center[1] = center[2] = 0.0;
	for ( j=0; j< vertexperface; j++)
	{
		center[0] += vertex_pos[faces[i + j*numfaces]             ];
		center[1] += vertex_pos[faces[i + j*numfaces] + numverts  ];
		center[2] += vertex_pos[faces[i + j*numfaces] + 2*numverts];
	}
	center[0]/=vertexperface;
	center[1]/=vertexperface;
	center[2]/=vertexperface;
	
#ifdef DEBUG_PRINT
	printf("face[%d]: center: %lf %lf %lf\n", i,center[0],center[1],center[2]);
#endif
}

void GeomData::getFaceNormal(int i, double* normal)
{
	normal[0] = face_normals[i           ];
	normal[1] = face_normals[i+numfaces  ];
	normal[2] = face_normals[i+2*numfaces];
}
		
void GeomData::ComputeFaceNormal(int i, double* normal)
{
	double u[3], v[3];

	u[0] = vertex_pos[faces[i + numfaces]				] - 
	   vertex_pos[faces[i]								];

	u[1] = vertex_pos[faces[i + numfaces] + numverts	] - 
	   vertex_pos[faces[i]			      + numverts	];

	u[2] = vertex_pos[faces[i + numfaces] + 2*numverts	] - 
	   vertex_pos[faces[i]			      + 2*numverts	]; // 1 - 0

	v[0] = vertex_pos[faces[i + 2*numfaces]				] - 
	   vertex_pos[faces[i]								];

	v[1] = vertex_pos[faces[i + 2*numfaces] + numverts	] - 
	   vertex_pos[faces[i]				    + numverts	];

	v[2] =  vertex_pos[faces[i + 2*numfaces] + 2*numverts] - 
		vertex_pos[faces[i]				     + 2*numverts]; // 2 - 0
		
	v3_cross(u,v,normal);
	v3_normalize(normal);
}

void GeomData::ComputeFaceNormals()
{
	double normal[3];
	
	int i;
	for (i=0; i< numfaces; i++){
		
		ComputeFaceNormal(i, normal);

		face_normals[i             ] = normal[0];
		face_normals[i+    numfaces] = normal[1];
		face_normals[i + 2*numfaces] = normal[2];

#ifdef DEBUG_PRINT
		printf("normal[%d]: %lf %lf %lf\n", i, face_normals[i], face_normals[i + numfaces], face_normals[i + 2*numfaces]);
#endif
	
	}
}

void GeomData::AlignFaces()
{
	// for plane objects with zero z dimension
		
	double normal[3];
	
	int i;
	for (i=0; i< numfaces; i++){		
		ComputeFaceNormal(i, normal);

		if (normal[2] < 0){
			FlipFace(i);
		}
	}
	
	ComputeFaceNormals();
}

void GeomData::ComputeOutterFaces()
{	
	int i;
	for ( i=0; i<numfaces; i++){
		if (MatchingFace(i) < 0 || numtetras == 0){
			outter_face[i] = 1;
#ifdef DEBUG_PRINT
			printf("face[%d] is an outter face\n", i);
#endif
		}
		else{
			printf("face[%d] is not an outter face\n", i);
			outter_face[i] = 0;
		}
	}
}

#pragma mark ___Verticies___
void GeomData::ComputeVertexNormals()
{		
	// expensive but one once			
	int i,j,k;
	for (i=0; i<numverts; i++){

		double vn[3];
		vn[0] =vn[1] = vn[2] = 0.0;
		int num =0;

		for (j=0; j<numfaces; j++){
		
			if (!outter_face[j])
				continue; // dont add in interior faces, for the tetras on the inside
		
			for (k=0; k< vertexperface; k++){
			
				if (faces[j + k*numfaces] == i){
				
					// add in that face's normal
					vn[0] += face_normals[j           ];
					vn[1] += face_normals[j+numfaces  ];
					vn[2] += face_normals[j+2*numfaces];
					num++;
				
				}					
			}				
		}
		
#ifdef DEBUG_PRINT
		printf("num[%d]: %d\n", i, num);
#endif			
		vn[0]/=num;
		vn[1]/=num;
		vn[2]/=num;
		
		v3_normalize(vn);
		
		vertex_normals[i             ] = vn[0];
		vertex_normals[i +   numverts] = vn[1];
		vertex_normals[i + 2*numverts] = vn[2];
#ifdef DEBUG_PRINT
		printf("vnormal[%d]: %lf %lf %lf\n", i, vertex_normals[i], vertex_normals[i + numverts], vertex_normals[i + 2*numverts]);
#endif	
	}
}	

void GeomData::RemoveStrandedVerticies()
{
	// for each vertex, check to see if it is rerenced in the faces
	// if it is, copy it over to the new vertex list
	// otherwise do not copy it over

/*	
	for (int i=0; i< numtetras; i++){
		printf("tetra[%d] %d %d %d %d\n", i, tetras[i], tetras[i+numtetras], tetras[i+2*numtetras], tetras[i+3*numtetras]);
	}
	
	for (int i=0; i< numfaces; i++){
		printf("faces[%d] %d %d %d\n", i, faces[i], faces[i+numfaces], faces[i+2*numfaces]);
	}
*/	
	printf("checking %d verts for being stranded\n", numverts);		
	printf("going from\n");	
	bool* found = (bool*) calloc(numverts, sizeof(bool));
	
	int numfound = 0;
	for (int i=0; i<numverts; i++){

		found[i] = 0;
		for (int j=0; j<numfaces; j++){		
			for (int k=0; k< vertexperface; k++){			
				if (faces[j+k*numfaces] == i){
					found[i] = 1;
				}			
			}	
		}
		
		/*
		found[i] = 0;
		for (int j=0; j<numtetras; j++){		
			for (int k=0; k< 4; k++){			
				if (tetras[j+k*numtetras] == i){
					found[i] = 1;
					break;
				}			
			}	
		}*/
		
		if (found[i])
			numfound++;
		
	}

	printf("to\n");
	double* newverts = (double*) calloc(numfound*3, sizeof(double));
	
	memset(min,0,3*sizeof(double));
	memset(max,0,3*sizeof(double));

	int num = 0;
	for (int i=0; i<numverts; i++){
	
		//printf("verts[%d]: %lf %lf %lf - %d\n", i, 
		//vertex_pos_init[i], vertex_pos_init[i+numverts], vertex_pos_init[i+2*numverts], found[i]);
	
		if (found[i]){
			newverts[num			 ] = vertex_pos_init[i	      	   ];
			newverts[num +   numfound] = vertex_pos_init[i +   numverts];
			newverts[num + 2*numfound] = vertex_pos_init[i + 2*numverts];
		
			//printf("newverts[%d]: %lf %lf %lf\n", num, newverts[num ], newverts[num+numfound], newverts[num+ 2*numfound]);
			
			if (newverts[num		    ]< min[0]) min[0] = newverts[num		    ];
			if (newverts[num+ numfound  ]< min[1]) min[1] = newverts[num+ numfound  ];
			if (newverts[num+ 2*numfound]< min[2]) min[2] = newverts[num+ 2*numfound];

			if (newverts[num		    ]> max[0]) max[0] = newverts[num		    ];
			if (newverts[num+ numfound  ]> max[1]) max[1] = newverts[num+ numfound  ];
			if (newverts[num+ 2*numfound]> max[2]) max[2] = newverts[num+ 2*numfound];

			num++;			
		}
		else{
			// for every face that refrenced a vertex higher than that, subtract one
			for (int j=0; j<numfaces; j++){
				for (int k=0; k< vertexperface; k++){
					if (faces[j+k*numfaces] >= i){
						faces[j+k*numfaces] = faces[j+k*numfaces]-1;
					}
				}	
			}
			for (int j=0; j<numtetras; j++){
				for (int k=0; k< 4; k++){
					if (tetras[j+k*numtetras] >= i){
						tetras[j+k*numtetras] = tetras[j+k*numtetras]-1;
					}				
				}	
			}
		}
	}

	// copy the new vertex list over
	numverts = num;
	printf("wound up with %d verts\n", numverts);
//#ifdef DEBUG_PRINT
	printf("data min: %lf %lf %lf max: %lf %lf %lf\n", 
	min[0], min[1], min[2], 
	max[0], max[1], max[2]);
//#endif		
	
/*	for (int i=0; i< numtetras; i++){
		printf("tetra[%d] %d %d %d %d\n", i, tetras[i] , tetras[i+numtetras] ,tetras[i+2*numtetras] ,tetras[i+3*numtetras]);
	}
	
	for (int i=0; i< numfaces; i++){
		printf("faces[%d] %d %d %d\n", i, faces[i] , faces[i+numfaces] ,faces[i+2*numfaces]);
	}
*/
	free(vertex_pos_init);
	free(vertex_pos);
	vertex_pos_init = (double*) calloc(num*3, sizeof(double));
	vertex_pos = (double*) calloc(num*3, sizeof(double));
	
	memcpy(vertex_pos_init, newverts, num*3*sizeof(double));
	memcpy(vertex_pos, newverts, num*3*sizeof(double));
	
	free(vertex_normals);
	vertex_normals = (double*) calloc(3*num,sizeof(double));

	free(found);
	free(newverts);
}

#pragma mark ___Tetras___	
void GeomData::ComputeTetraCenter(int i, double* center)
{
	int j;
	
	center[0] = center[1] = center[2] = 0.0;
	for ( j=0; j< 4; j++)
	{
		center[0] += vertex_pos[tetras[i + j*numtetras]             ];
		center[1] += vertex_pos[tetras[i + j*numtetras] + numverts  ];
		center[2] += vertex_pos[tetras[i + j*numtetras] + 2*numverts];
	}
	center[0]/=4;
	center[1]/=4;
	center[2]/=4;
	
#ifdef DEBUG_PRINT
	printf("tetra[%d]: center: %lf %lf %lf\n", i,center[0],center[1],center[2]);
#endif		
}

bool GeomData::TetraFaceIsPointingInwards(int ti, int fi)
{			
	double normal[3];
	ComputeFaceNormal(fi, normal);
		
	double fcenter[3];
	ComputeFaceCenter(fi, fcenter);
	
	double tcenter[3];
	ComputeTetraCenter(ti, tcenter);
		
	double w[3];			
	v3_sub(fcenter,tcenter, w);
	v3_normalize(w);
	
	double dist = v3_dot(normal,w);

	if (dist < 0.0) return 1;
	
	dist = v3_dot(w,normal);
	if (dist < 0.0) return 1;
	
	return 0;		
}

int GeomData::getFaceIndex(int node0, int node1, int node2)
{
	int j,k;
	
	for (j=0; j<numfaces; j++){
			
		bool has[3];
		has[0] = has[1] = has[2] = 0;
	
		for (k=0; k<3; k++){ //cbnote for triangles only
		
			int n0 = faces[j+k*numfaces];
			if (n0 == node0){
				has[k] = 1;
			}
			int n1 = faces[j+k*numfaces];
			 if (n1 == node1){
				has[k] = 1;
			}
			int n2 = faces[j+k*numfaces];
			if ( n2 == node2){
				has[k] = 1;
			}
			
		}
		if (has[0] && has[1] && has[2])
			return j;					
	}
	
	return -1;

}

void GeomData::CleanTetraModel()
{
	int i;
	int face_index;
	int i1, i2, i3, i4;
	
	for (i=0; i<numtetras; i++){
		
		i1 = tetras[i	         ];
		i2 = tetras[i+1*numtetras];
		i3 = tetras[i+2*numtetras];
		i4 = tetras[i+3*numtetras];
								
		face_index = getFaceIndex(i1, i2, i3);		
		if (face_index < 0) return;
					
		if (TetraFaceIsPointingInwards(i,face_index)){
			FlipFace(face_index);							
		}
		
		face_index = getFaceIndex(i2, i3, i4);
		if (face_index < 0) return;
					
		if (TetraFaceIsPointingInwards(i,face_index)){
			FlipFace(face_index);							
		}
		
		face_index = getFaceIndex(i3, i4, i1);
		if (face_index < 0) return;
					
		if (TetraFaceIsPointingInwards(i,face_index)){
			FlipFace(face_index);							
		}
		
		face_index = getFaceIndex(i4, i1, i2);
		if (face_index < 0) return;
				
		if (TetraFaceIsPointingInwards(i,face_index)){
			FlipFace(face_index);							
		}
	}
}

#pragma mark ___Picking___
// cbnote need a way to find out if the point is inside a tetra or not
bool GeomData::insideTetrahedra(int i, double* intPt)
{
	//http://steve.hollasch.net/cgindex/geometry/ptintet.html
		
	double x = intPt[0];
	double y = intPt[1];
	double z = intPt[2];
	
	int sign0, sign1, sign2, sign3, sign4;
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
	
	double coords[12];
	
	GeomInfo ginfo;
	ginfo.dimension = 3;
	ginfo.nodedofs = 3;
	ginfo.numelements = numtetras;
	ginfo.numelnodes = 4;
	ginfo.numnodes = numverts;
	ginfo.scale = 1.0;
	
	get3dCoords(coords, i, vertex_pos, tetras, ginfo);
	//checkTetraElement(coords); 
	
	double close = 1e-6;
	double cutoff = 1e-6;
	// EPSILON = 1e-6
	
	x1 =  coords[0];
	y1 =  coords[4];
	z1 =  coords[8];	
	
	if ( dabs(x-x1)< close && dabs(y-y1)< close && dabs(z-z1)< close )
			return true;
	
	x2 =  coords[1];
    y2 =  coords[5];
	z2 =  coords[9];
	
	if ( dabs(x-x2)< close && dabs(y-y2)< close && dabs(z-z2)< close )
			return true;

	x3 =  coords[2];
    y3 =  coords[6];
	z3 =  coords[10];	

	if ( dabs(x-x3)< close  && dabs(y-y3)< close  && dabs(z-z3)< close  )
			return true;

	x4 =  coords[3];
    y4 =  coords[7];
	z4 =  coords[11];	
	
	if ( dabs(x-x4)< close && dabs(y-y4)< close && dabs(z-z4)< close )
			return true;

	double d0[16];
	d0[0] = x1; d0[4] = y1; d0[8] =  z1; d0[12] = 1.;
	d0[1] = x2; d0[5] = y2; d0[9] =  z2; d0[13] = 1.;
	d0[2] = x3; d0[6] = y3; d0[10] = z3; d0[14] = 1.;
	d0[3] = x4; d0[7] = y4; d0[11] = z4; d0[15] = 1.;
	
	double s0 = m4_det(d0);
	if (s0  == 0.0){
		printf("degenerate tetrahedra\n");
		exit(1);
	}
	else if (s0 < 0){
		sign0 = -1;
	}
	else {
		sign0 = 1;
	}
	
	double d1[16];
	d1[0] =  x; d1[4] =  y; d1[8] =   z; d1[12] = 1.;
	d1[1] = x2; d1[5] = y2; d1[9] =  z2; d1[13] = 1.;
	d1[2] = x3; d1[6] = y3; d1[10] = z3; d1[14] = 1.;
	d1[3] = x4; d1[7] = y4; d1[11] = z4; d1[15] = 1.;
	
	double s1 = m4_det(d1);
	if (dabs(s1)  < cutoff){
		sign1 = sign0;
	}
	else if (s1 < 0){
		sign1 = -1;
	}
	else {
		sign1 = 1;
	}

	double d2[16];
	d2[0] = x1; d2[4] = y1; d2[8] =  z1; d2[12] = 1.;
	d2[1] =  x; d2[5] =  y; d2[9] =   z; d2[13] = 1.;
	d2[2] = x3; d2[6] = y3; d2[10] = z3; d2[14] = 1.;
	d2[3] = x4; d2[7] = y4; d2[11] = z4; d2[15] = 1.;
	
	double s2 = m4_det(d2);
	if (dabs(s2)  < cutoff){
		sign2 = sign0;
	}
	else if (s2 < 0){
		sign2 = -1;
	}
	else {
		sign2 = 1;
	}
		
	double d3[16];
	d3[0] = x1; d3[4] = y1; d3[8] =  z1; d3[12] = 1.;
	d3[1] = x2; d3[5] = y2; d3[9] =  z2; d3[13] = 1.;
	d3[2] =  x; d3[6] =  y; d3[10] =  z; d3[14] = 1.;
	d3[3] = x4; d3[7] = y4; d3[11] = z4; d3[15] = 1.;
	
	double s3 = m4_det(d3);
	if (dabs(s3)  <  cutoff){
		sign3 = sign0;
	}
	else if (s3 < 0){
		sign3 = -1;
	}
	else {
		sign3 = 1;
	}
		
	double d4[16];
	d4[0] = x1; d4[4] = y1; d4[8] =  z1; d4[12] = 1.;
	d4[1] = x2; d4[5] = y2; d4[9] =  z2; d4[13] = 1.;
	d4[2] = x3; d4[6] = y3; d4[10] = z3; d4[14] = 1.;
	d4[3] =  x; d4[7] =  y; d4[11] =  z; d4[15] = 1.;
	
	double s4 = m4_det(d4);
	if (dabs(s4) < cutoff){
		sign4 = sign0;
	}
	else if (s4 < 0){
		sign4 = -1;
	}
	else {
		sign4 = 1;
	}
	
	int sum = sign0 + sign1 + sign2 + sign3 + sign4;
	//printf("sum: %d\n", sum);
	
	for (int i=0; i<16; i++){
		if ( dabs(s0 -(s1 + s2 + s3 + s4)) > EPSILON ){
			printf("bunk tetrahedra\n");
			exit(1);
		}
	}
	
	
	if ( (sum != -5) && (sum != 5) ){
		return false;
	}


	return true;

}

/*bool GeomData::InsideQuadrilateral(int i, double* intPt)
{
	Point3DBase<T> verts[4];
	verts[0] = p0;
	verts[1] = p1;
	verts[2] = p2;
	verts[3] = p3; 
	
	Point3DBase<T> p4,p5;
  
	double m1,m2;
	double anglesum=0,costheta;

	for (int i=0;i<4;i++) {

  	p4.x = verts[i].x - x;
  	p4.y = verts[i].y - y;
  	p4.z = verts[i].z - z;
  	p5.x = verts[(i+1)%4].x - x;
  	p5.y = verts[(i+1)%4].y - y;
  	p5.z = verts[(i+1)%4].z - z;

  	m1 = p4.Length();
  	m2 = p5.Length();
  	if (m1*m2 <= EPSILON)
    	 return(1); // We are on a node, consider this inside 
  	else
    	 costheta = (p4.x*p5.x + p4.y*p5.y + p4.z*p5.z) / (m1*m2);

		anglesum += acos(costheta);
	}
	
	double diff = anglesum - TWOPI;
	if (sqrtf(diff*diff) <= EPSILON) return (1);
	return (0);

}
*/


bool GeomData::insideTriangle(int i, double* intPt)
{
	
#ifdef DEBUG_PRINT
	printf("face: %d intPt: %lf %lf %lf\n", i , intPt[0], intPt[1], intPt[2]);
#endif	

	double verts[3][3];
	int j;
	
	verts[0][0] =  vertex_pos[faces[i			]			   ];
    verts[0][1] =  vertex_pos[faces[i			] + numverts   ];
	verts[0][2] =  vertex_pos[faces[i			] + 2*numverts ];		
	
	verts[1][0] =  vertex_pos[faces[i + numfaces]			  ];
    verts[1][1] =  vertex_pos[faces[i + numfaces] + numverts  ];
	verts[1][2] =  vertex_pos[faces[i + numfaces] + 2*numverts];
	
	verts[2][0] =  vertex_pos[faces[i + 2*numfaces]			    ];
    verts[2][1] =  vertex_pos[faces[i + 2*numfaces] + numverts  ];
	verts[2][2] =  vertex_pos[faces[i + 2*numfaces] + 2*numverts];
	
#ifdef DEBUG_PRINT
	int k;
	printf("face[%d]\n", i);
	for (j=0; j<3; j++){
		for (k=0; k<3; k++)
			printf("v[%d]: %lf \t",j, verts[j][k]);
			
		printf("\n");
	}
#endif
	
	double p4[3];
	double p5[3]; 
	
	double m1,m2;
	double anglesum=0.0;
	double costheta;

	for ( j=0;j<3;j++) { 
		 
		p4[0] = verts[j][0] - intPt[0];
		p4[1] = verts[j][1] - intPt[1];
		p4[2] = verts[j][2] - intPt[2];
		
		p5[0] = verts[(j+1)%3][0] - intPt[0];
		p5[1] = verts[(j+1)%3][1] - intPt[1];
		p5[2] = verts[(j+1)%3][2] - intPt[2];
		
#ifdef DEBUG_PRINT
		printf("p4: %lf %lf %lf \n p5: %lf %lf %lf\n", p4[0],p4[1],p4[2], p5[0],p5[1],p5[2]);
#endif		
		m1 = v3_length(p4);
		m2 = v3_length(p5);
		
#ifdef DEBUG_PRINT
		printf("length p4: %lf length p5: %lf\n", m1, m2);
#endif		
		if (m1*m2 <= EPSILON)
			return(true); /* We are on a node, consider this inside */
		else
			costheta = v3_dot(p4,p5)/(m1*m2);
			
		anglesum += acos(costheta); 
		
	}
	
	double diff = anglesum - 2.0*M_PI;
#ifdef DEBUG_PRINT
	printf("anglesum: %lf\n", anglesum);
	printf("diff: %lf\n", diff);
#endif		
	if (sqrt(diff*diff) <= EPSILON) {
#ifdef DEBUG_PRINT	
		printf("inside the triangle[%d]\n", i);
#endif
		return (true);
	}
		
	return (false);
	
}	

void GeomData::determinePlaneEquation(int i, double &a, double &b, double &c, double &d)
{
	double n0[3];
	n0[0] = vertex_pos[faces[i	  ]			    ];
	n0[1] = vertex_pos[faces[i	  ] + numverts  ];
	n0[2] = vertex_pos[faces[i    ] + 2*numverts];
			
	double normal[3];
	normal[0] =     face_normals[i			   ];
	normal[1] =     face_normals[i +  numfaces ];
	normal[2] =     face_normals[i + 2*numfaces];
			
	a = normal[0];
	b = normal[1];
	c = normal[2];
	d = - (a*n0[0] + b*n0[1] + c*n0[2]);
}

bool GeomData::checkFace(int i, double* P1, double* P2, double* intPt)
{

	double a,b,c,d;
	determinePlaneEquation(i, a,b,c,d);

	double n0[3];
	n0[0] = vertex_pos[faces[i	  ]			    ];
	n0[1] = vertex_pos[faces[i	  ] + numverts  ];
	n0[2] = vertex_pos[faces[i    ] + 2*numverts];

	double alpha, t;		
	alpha = a * (P1[0] - P2[0]) + b * (P1[1] - P2[1]) + c * (P1[2] - P2[2]); 
	if(alpha == 0) return 0;
	
	t = -(a * P2[0] + b * P2[1] + c * P2[2] + d) / alpha;
	if ((t < 0) || (t > 1)) return 0;
	
	intPt[0] = (P1[0] - P2[0]) * t + P2[0];
	intPt[1] = (P1[1] - P2[1]) * t + P2[1];
	intPt[2] = (P1[2] - P2[2]) * t + P2[2];
	
#ifdef DEBUG_PRINT
	printf("crosses plane of face[%d]\n",i);
#endif		
	
	if (insideTriangle(i, intPt))
		return 1;	
		
	return 0;
}

bool GeomData::checkClosest(double* intPt, double* cIntPt, double* P1)
{
	double u[3];
	u[0] = intPt[0] - P1[0];
	u[1] = intPt[1] - P1[1];
	u[2] = intPt[2] - P1[2]; // vector from camera to int pt
	
	v3_sub(intPt,P1,u);
	
	double v[3];
	v[0] = cIntPt[0] - P1[0];
	v[1] = cIntPt[1] - P1[1];
	v[2] = cIntPt[2] - P1[2]; // vector from camera to current cloest
	
	v3_sub(cIntPt,P1,v);
	
	if (v3_length(u) < v3_length(v)){
		cIntPt[0] = intPt[0];
		cIntPt[1] = intPt[1];
		cIntPt[2] = intPt[2];	
		return true;
	}
	return false;
}

int GeomData::pickedOnSurface(double* P1, double* P2, double* intPt)
{	
#ifdef DEBUG_PRINT
	printf("pickOrigin: %lf %lf %lf\n",    P1[0], P1[1], P1[2]);
	printf("pickDirection: %lf %lf %lf\n", P2[0], P2[1], P2[2]);
#endif
	
	double cIntPt[3];
	memset(cIntPt, 1E6, 3*sizeof(double));
		
	int cFace = -1;
	int i;
	for( i=0; i < numfaces; i++) { // the slowest way possible	
								
		if ( checkFace(i, P1, P2, intPt)){
#ifdef DEBUG_PRINT	
			printf("intPt: %lf %lf %lf\n", intPt[0], intPt[1], intPt[2]);
#endif
			if (checkClosest(intPt, cIntPt, P1))
				cFace = i;
#ifdef DEBUG_PRINT	
			printf("cintPt: %lf %lf %lf\n", cIntPt[0], cIntPt[1], cIntPt[2]);				
#endif
			
		}									
	}
	if (cFace > -1){
	
		intPt[0]= cIntPt[0];
		intPt[1]= cIntPt[1];
		intPt[2]= cIntPt[2];
		
#ifdef DEBUG_PRINT	
		printf("closest face: %d\n", cFace);
#endif
		return cFace; 
	}
	
	return -1;
}

