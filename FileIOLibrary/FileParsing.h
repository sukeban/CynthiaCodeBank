/*
 *  FileParsing.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 9/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef _FileParsing_H_
#define _FileParsing_H_

#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#define MAXLINE 1024

#include "ModelSharedData.h"

#include <algorithm>

/*!
 wrapper around fgets to throw out comment lines (lines starting with #)
 and possibly in the future check for other special lines in data files 
*/
char *myfgets(char *s, int n, FILE *stream)
{
    char *t;                     // return string
    do {
        t = fgets(s, n, stream);
        if (t == NULL)			// no more to read
            return t;
        else {                             // is non-null string
            // kill newline at end of string
            if (t[strlen(t)-1] == '\n') t[strlen(t)-1] = '\0'; 
            
            // eliminate leading white space
            while (isspace(t[0])) t++;        
            
            // return if not a comment or blank line
            if (t[0] != '#' && t[0] != '\0')
                return t;                      
        }
    } while (1);
}

#pragma mark ___SMF___

/*!
    opens a SMF (http://gicl.cs.drexel.edu/wiki/SMF) file and returns the number of nodes, faces and verticies per element
 */
int  ReadSMFNumbers(char* filename, int &num_nodes, int &num_faces, int &verts_per_face) 
{	
	FILE *file;
    
    // open the node data file
    if ((file = fopen(filename, "r")) == NULL) {
        printf("Could not open data file\n");
        return 0;
	}
	
    // read file once to get the number of verts
 	num_nodes = 0;
    while (!feof(file)) {
		char s[MAXLINE];
        fgets(s, MAXLINE, file);
        switch (s[0]) {
        case 'v' : num_nodes++; break;// nodes
        case '#' : // comment
        default  : // who knows?
            break;
        }
    }
    
    fseek(file, 0, SEEK_SET);

	num_faces = 0;
	while (!feof(file)) {
		char s[MAXLINE];
		fgets(s, MAXLINE, file);
        if (feof(file)) break;
        switch (s[0]) {
		
        case 'f' : {
			num_faces++; 
			int i1, i2, i3, i4;
			verts_per_face = sscanf(s, "%*c %d %d %d %d", &i1, &i2, &i3, &i4);
			break;
		}
        case 'e' : //edge
        case 'v' : // vertex 
        case '#' : // comment
        default  : // who knows?
            break;
        }
    }
    
    // reset input file pointer
    printf("%d nodes %d faces verts_per_face %d\n", num_nodes, num_faces, verts_per_face);

	fclose(file);
	return 1;
}

/*!
 opens a SMF (http://gicl.cs.drexel.edu/wiki/SMF) file and puts the data in an pre-allocated GeomData object
 */
int ReadSMFData(char* filename, GeomData* data)
{
	FILE *file;
    
    // open the node data file
    if ((file = fopen(filename, "r")) == NULL) {
        printf("Could not open data file\n");
        return 0;
	}

	double x, y, z;
	int i;
	
	int num_nodes = data->numverts;
	int num_faces = data->numfaces;
	
#ifdef DEBUG_PRINT
	printf("verts=[\n");
#endif 
	// loop through and get data for each node
    for ( i = 0; i < num_nodes; i++) {
		char s[MAXLINE];
        do {			
           if (myfgets(s, MAXLINE, file) == NULL) {
                printf("Not enough nodes in data file\n");
                goto close;
            }
        } while (s[0] != 'v'); // advance past any comments
        
        // initialize node with current data line
           
        sscanf(s, "%*c %lf %lf %lf", &x, &y, &z);
        data->vertex_pos_init[i] = x;
		data->vertex_pos_init[i+  num_nodes] = y;
		data->vertex_pos_init[i+2*num_nodes] = z;
		
		if (x< data->min[0]) data->min[0] = x;
		if (y< data->min[1]) data->min[1] = y;
		if (z< data->min[2]) data->min[2] = z;
		
		if (x> data->max[0]) data->max[0] = x;
		if (y> data->max[1]) data->max[1] = y;
		if (z> data->max[2]) data->max[2] = z;

#ifdef DEBUG_PRINT
		printf("%lf %lf %lf;\n", x,y,z);  
#endif     
	
    }
#ifdef DEBUG_PRINT
	printf("];\n");
#endif 
	data->reset();
	data->radius = std::max(data->max[0] - data->min[0], data->max[1] - data->min[1]);
	data->radius = std::max(data->max[2] - data->min[2], data->radius);
	
#ifdef DEBUG_PRINT
	printf("data min: %lf %lf %lf max: %lf %lf %lf\n", 
	data->min[0], data->min[1], data->min[2], 
	data->max[0], data->max[1], data->max[2]);
#endif	
	// loop through and get data for each face
    // reset input file pointer
    fseek(file, 0, SEEK_SET);
 
#ifdef DEBUG_PRINT
	printf("faces=[\n");
#endif	
	for (i=0; i<num_faces; i++) {
		char s[MAXLINE];
		do {
            if (myfgets(s, MAXLINE, file) == NULL) {
                printf("Not enough faces in data file\n");
                goto close;
            }
        } while (s[0] != 'f'); // advance past any comments
		   
		// initialize face with current data line
		int i1, i2, i3, i4;
		int numread = sscanf(s, "%*c %d %d %d %d", &i1, &i2, &i3, &i4);
		if (numread < 3) {
			printf("Bad read on face %d- '%d %d %d %d' - IGNORED\n", i, i1,i2,i3,i4);
			continue;
		}
            
         // SMF is 1-based, we are 0-index-based
         i1--; i2--; i3--; i4--;
            
		data->faces[i			  ] = i1;
		data->faces[i +num_faces  ] = i2;
		data->faces[i +2*num_faces] = i3;
		
		if (i4> -1){
			data->faces[i +3*num_faces] = i4;
		}
		
#ifdef DEBUG_PRINT
		printf("%d %d %d %d;\n",i1,i2,i3,i4);
#endif
	}
#ifdef DEBUG_PRINT
	printf("];\n");
#endif
	
close:
	fclose(file);
	
	printf("data min: %lf %lf %lf data max: %lf %lf %lf\n", 
	data->min[0], data->min[1], data->min[2], data->max[0], data->max[1], data->max[2]);
	return 1;

}

#pragma mark ___Amira___

/*!
 opens an AMIRA (http://www.mpi-inf.mpg.de/~weinkauf/notes/amiramesh.html)file and returns the number of nodes and number of tetrahedra in the file
 */
int  ReadAmiraNumbers(char* filename, int &num_nodes, int &num_tetras)
{	
	FILE *file;
    
    // open the node data file
    if ((file = fopen(filename, "r")) == NULL) {
        printf("Could not open data file\n");
        return 0;
	}
	
	char s[MAXLINE];
	myfgets(s, MAXLINE, file);
	int num_args = sscanf(s, "%d %d", &num_nodes, &num_tetras);
	if (num_args < 2) {
	 	printf("bad file\n");
		return 0;
	}

	fclose(file);
	return 1;
}

/*!
 opens an AMIRA (http://gicl.cs.drexel.edu/wiki/SMF) file and puts the data in an pre-allocated GeomData object
 */
int ReadAmiraData(char* filename, GeomData* data)
{
	FILE *file;
    
    // open the node data file
    if ((file = fopen(filename, "r")) == NULL) {
        printf("Could not open data file\n");
        return 0;
	}

	int num_nodes = 0;
	int num_tetras = 0;
	
	char s[MAXLINE];
	myfgets(s, MAXLINE, file);
	int num_args = sscanf(s, "%d %d", &num_nodes, &num_tetras);
	if (num_args < 2) {
	 	printf("bad file\n");
		return 0;
	}
	
	int num_faces = 4*num_tetras;
	int i;
	for ( i = 0; i < num_nodes; i++) {
		if (myfgets(s, MAXLINE, file) == NULL) {
			printf("Not enough nodes in data file\n");
			return(0);
		}
					
		double x, y, z;
		int node_num;
		int numread = sscanf(s, "%d %lf %lf %lf", &node_num, &x, &y, &z);
		if (numread < 3) {
			printf("Bad read on point%d- '%lf %lf %lf' - IGNORED\n", i, x,y,z);
			continue;
		}
		data->vertex_pos_init[i] = x;
		data->vertex_pos_init[i+  num_nodes] = y;
		data->vertex_pos_init[i+2*num_nodes] = z;
		
		if (x< data->min[0]) data->min[0] = x;
		if (y< data->min[1]) data->min[1] = y;
		if (z< data->min[2]) data->min[2] = z;
		
		if (x> data->max[0]) data->max[0] = x;
		if (y> data->max[1]) data->max[1] = y;
		if (z> data->max[2]) data->max[2] = z;
		
		data->radius = std::max(data->max[0] - data->min[0], data->max[1] - data->min[1]);
		data->radius = std::max(data->max[2] - data->min[2], data->radius);
		
		
#ifdef DEBUG_PRINT
		printf("node: %lf %lf %lf\n",x,y,z);  
#endif
    }
	
	data->reset();
	
#ifdef DEBUG_PRINT
	printf("data min: %lf %lf %lf max: %lf %lf %lf\n", 
	data->min[0], data->min[1], data->min[2], 
	data->max[0], data->max[1], data->max[2]);
#endif	
	// loop through and get data for each face

	for (i=0; i<num_tetras; i++) {
		if (myfgets(s, MAXLINE, file) == NULL) {
			printf("Not enough tetras in data file\n");
			return(0);
		}
    
		// initialize face with current data line
		int num, dum, i1, i2, i3, i4;
		char ptype[4];
		int numread = sscanf(s, "%d %d %s %d %d %d %d", 
							 &num, &dum, ptype, &i1, &i2, &i3, &i4);

		if (numread < 4) {
			printf("Bad read on tetra %d- '%d %d %d %d' - IGNORED\n", i, i1,i2,i3,i4);
			continue;
		}
		
		data->tetras[i				] = i1;
		data->tetras[i+  num_tetras	] = i2;
		data->tetras[i+2*num_tetras	] = i3;
		data->tetras[i+3*num_tetras	] = i4;

#ifdef DEBUG_PRINT
		printf("i1: %d %d %d %d\n", i1,i2,i3,i4);
#endif		
		data->faces[i*4				    ] = i1;
		data->faces[i*4    +num_faces	] = i3;
		data->faces[i*4    +2*num_faces	] = i2;		
#ifdef DEBUG_PRINT
		printf("face 1: %d %d %d \n", i1,i3,i2);
#endif		
		data->faces[i*4 + 1				] = i2;
		data->faces[i*4 + 1+num_faces	] = i3;
		data->faces[i*4 + 1+2*num_faces	] = i4;

#ifdef DEBUG_PRINT
		printf("face 2: %d %d %d \n", i2,i3,i4);
#endif
		data->faces[i*4 + 2				] = i1;
		data->faces[i*4 + 2+num_faces	] = i2;
		data->faces[i*4 + 2+2*num_faces	] = i4;

#ifdef DEBUG_PRINT
		printf("face 3: %d %d %d \n", i1,i2,i4);
#endif
		data->faces[i*4 + 3				] = i1;
		data->faces[i*4 + 3+num_faces	] = i4;
		data->faces[i*4 + 3+2*num_faces	] = i3;
		
#ifdef DEBUG_PRINT
		printf("face 4: %d %d %d \n", i1,i4,i3);
#endif
	}
	
	fclose(file);
	
	printf("data min: %lf %lf %lf data max: %lf %lf %lf\n", data->min[0], data->min[1], data->min[2], data->max[0], data->max[1], data->max[2]);
	return 1;

}

#pragma mark ___Mesh___
/*!
 opens a MESH (http://wordwood.merseine.us/TitanQuest-MeshFileFormat/) file and returns the number of nodes, faces and tetrahedra
 */
int  ReadMeshNumbers(char* filename, int &num_nodes, int &num_faces, int &num_tetras)
{	
	FILE *file;
    
    // open the node data file
    if ((file = fopen(filename, "r")) == NULL) {
        printf("Could not open data file\n");
        return 0;
	}
	
	num_nodes = 0;
	num_faces = 0;
	num_tetras = 0;
	
	char v[MAXLINE];
	sprintf(v, "Vertices");
	
	char t[MAXLINE];
	sprintf(t, "Triangles");
	
	char h[MAXLINE];
	sprintf(h, "Tetrahedra");

	char s[MAXLINE];     
	fgets(s, MAXLINE, file);
	
	while (!feof(file)) {
	  
		char attr_name[80]; sscanf(s, "%s", attr_name);
		
		// check for the word Vericies, Triangles, or Tetrahedra
		// then get the next line and assign it to the variable
		if (!strcmp(attr_name, v)){
			fgets(s, MAXLINE, file);
			sscanf(s, "%d", &num_nodes); 
		}
		else if (!strcmp(attr_name, t)){
			fgets(s, MAXLINE, file);
			sscanf(s, "%d", &num_faces);
		}
		else if (!strcmp(attr_name, h)){
			fgets(s, MAXLINE, file);
			sscanf(s, "%d", &num_tetras);
		}
		fgets(s, MAXLINE, file);
		
	}

	fclose(file);
	return 1;
}
      
/*!
 opens a MESH (http://wordwood.merseine.us/TitanQuest-MeshFileFormat/) file and puts the data in a pre-allocated GeomData object
 */
int ReadMeshData(char* filename, GeomData* data)
{
	FILE *file;
    
    // open the node data file
    if ((file = fopen(filename, "r")) == NULL) {
        printf("Could not open data file\n");
        return 0;
	}	
	
	int i;
	int numread;	
	double x, y, z;		
	int  i1, i2, i3, i4;

	
	char s[MAXLINE];
	myfgets(s, MAXLINE, file); // skip 4 lines
	myfgets(s, MAXLINE, file);
	myfgets(s, MAXLINE, file);
	myfgets(s, MAXLINE, file);
	
	int num_nodes = 0;
	myfgets(s, MAXLINE, file);
	int num_args = sscanf(s, "%d", &num_nodes); // and read off the number (n) of verticies to expect
	if (num_args < 1) {
	 	printf("bad file\n");
		return 0;
	}	
			
#ifdef DEBUG_PRINT
	printf("verts=[\n");
#endif
	for ( i = 0; i < num_nodes; i++) { // parse the next (n) lines as verticies
		if (myfgets(s, MAXLINE, file) == NULL) {
			printf("Not enough nodes in data file\n");
			return(0);
		}					
		
		numread = sscanf(s, "%lf %lf %lf", &x, &y, &z);
		if (numread < 3) {
			printf("Bad read on point%d '%lf %lf %lf' - IGNORED\n",i,x,y,z);
			continue;
		}
		data->vertex_pos_init[i				] = x;
		data->vertex_pos_init[i+  num_nodes	] = y;
		data->vertex_pos_init[i+2*num_nodes	] = z;
		
		if (x< data->min[0]) data->min[0] = x;
		if (y< data->min[1]) data->min[1] = y;
		if (z< data->min[2]) data->min[2] = z;
		
		if (x> data->max[0]) data->max[0] = x;
		if (y> data->max[1]) data->max[1] = y;
		if (z> data->max[2]) data->max[2] = z;
		
		data->radius = std::max(data->max[0] - data->min[0], data->max[1] - data->min[1]);
		data->radius = std::max(data->max[2] - data->min[2], data->radius);
		
		
#ifdef DEBUG_PRINT
		printf("%lf %lf %lf;\n",x,y,z);  
#endif
    }
#ifdef DEBUG_PRINT
	printf("];\n");	
#endif	
	data->reset();
	
#ifdef DEBUG_PRINT
	printf("data min: %lf %lf %lf max: %lf %lf %lf\n", 
	data->min[0], data->min[1], data->min[2], 
	data->max[0], data->max[1], data->max[2]);
#endif	

	if (data->numfaces){

		myfgets(s, MAXLINE, file); // then skip a line

		int num_faces = 0;
		myfgets(s, MAXLINE, file);
		num_args = sscanf(s, "%d", &num_faces); // and read the number (m) of faces to expect
		if (num_args < 1) {
			printf("bad file\n");
			return 0;
		}	
		
	#ifdef DEBUG_PRINT
		printf("faces=[\n");
	#endif	
		for (i=0; i<num_faces; i++) { 	// read the next (m) lines as faces

			if (myfgets(s, MAXLINE, file) == NULL) {
				printf("Not enough faces in data file\n");
				return(0);
			}
			numread = sscanf(s, "%d %d %d", 
								 &i1, &i2, &i3);

			if (numread <  3) {
				printf("Bad read on face %d- '%d %d %d' - IGNORED\n", i, i1,i2,i3);
				continue;
			}
			
			i1--; i2--; i3--;   
			
			data->faces[i				    ] = i1;
			data->faces[i    +num_faces		] = i2;
			data->faces[i    +2*num_faces	] = i3;
			
	#ifdef DEBUG_PRINT
			printf("%d %d %d;\n", i1, i2, i3);
	#endif	
		}
	#ifdef DEBUG_PRINT
		printf("];\n");
	#endif	
	}

	if (data->numtetras){
	
		// read the next (v) lines as tetra
		myfgets(s, MAXLINE, file); // then skip a line

		int num_tetras = 0;
		myfgets(s, MAXLINE, file);
		num_args = sscanf(s, "%d", &num_tetras); // and read the number (v) of tetra to expect
		//if (num_args < 1) {
		 //	printf("bad file\n");
		//	return 0;
		//}	

		int num_faces = 4*num_tetras;
		data->numfaces = num_faces;
		
		// if there werent any faces allocate them now
		if (data->faces) free(data->faces);
		data->faces = (int*) calloc(num_faces*3, sizeof(int));		
		
		if(data->face_normals) free(data->face_normals);
		data->face_normals = (double*) calloc(num_faces*3, sizeof(double));
		
		if (data->outter_face) free(data->outter_face);
		data->outter_face = (bool*) calloc(num_faces, sizeof(bool));
		
		
	#ifdef DEBUG_PRINT
		printf("tetras=[\n");
	#endif
		for (i=0; i<num_tetras; i++) { 	// read the next (m) lines as tetras

			if (myfgets(s, MAXLINE, file) == NULL) {
				printf("Not enough faces in data file\n");
				return(0);
			}
			numread = sscanf(s, "%d %d %d %d", 
								 &i1, &i2, &i3, &i4);

			if (numread < 4) {
				printf("Bad read on tetra %d- '%d %d %d' - IGNORED\n", i, i1,i2,i3,i4);
				continue;
			}
			
			 i1--; i2--; i3--; i4--;
   
			
			data->tetras[i				] = i1;
			data->tetras[i+  num_tetras	] = i2;
			data->tetras[i+2*num_tetras	] = i3;
			data->tetras[i+3*num_tetras	] = i4; 
	#ifdef DEBUG_PRINT
			printf("%d %d %d %d;\n", i1, i2, i3, i4);
	#endif	
	
			data->faces[i*4				    ] = i1;
			data->faces[i*4    +num_faces	] = i3;
			data->faces[i*4    +2*num_faces	] = i2;		
	#ifdef DEBUG_PRINT
			printf("face 1: %d %d %d \n", i1,i3,i2);
	#endif		
			data->faces[i*4 + 1				] = i2;
			data->faces[i*4 + 1+num_faces	] = i3;
			data->faces[i*4 + 1+2*num_faces	] = i4;

	#ifdef DEBUG_PRINT
			printf("face 2: %d %d %d \n", i2,i3,i4);
	#endif
			data->faces[i*4 + 2				] = i1;
			data->faces[i*4 + 2+num_faces	] = i2;
			data->faces[i*4 + 2+2*num_faces	] = i4;

	#ifdef DEBUG_PRINT
			printf("face 3: %d %d %d \n", i1,i2,i4);
	#endif
			data->faces[i*4 + 3				] = i1;
			data->faces[i*4 + 3+num_faces	] = i4;
			data->faces[i*4 + 3+2*num_faces	] = i3;
			
	#ifdef DEBUG_PRINT
			printf("face 4: %d %d %d \n", i1,i4,i3);
	#endif
		
		}
		
	#ifdef DEBUG_PRINT
		printf("];\n");
	#endif
		
	}
	
	fclose(file);
	
#ifdef DEBUG_PRINT
	printf("data min: %lf %lf %lf data max: %lf %lf %lf\n", data->min[0], data->min[1], data->min[2], data->max[0], data->max[1], data->max[2]);
#endif
	return 1;

}


#endif