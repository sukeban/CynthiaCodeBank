/*
 *  SharedFEData.h //cbnote rename to FESharedData.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 1/9/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _SharedData_H
#define _SharedData_H

#include <mach/mach.h>
#include <mach/mach_time.h>

/*!
 constans used in CBLAS
 */
char	L		= 'L';
char	N		= 'N';
char	V		= 'V';
char	SIDE	= 'L';
char	TRANS	= 'T';
char	I		= 'I';
char	U		= 'U';
char	S		= 'S';

/*!
 Structure for Stress State
 */
enum StressEnum{
	planestress,
	planestrain,
	pure_axissymmetric,// makes a 3d problem a 2d one
	solid,
	acoustic
};

/*!
 Structure for Material Properties
 */
struct MaterialInfo {
	double		elasticModulus;
	double		poisson;
	double		density;
	double		soundspeed;
	StressEnum	state;
};

struct PlateInfo{
	double thickness;
};

/*!
 Structure for Shell Properties
 */
struct ShellInfo{
	double	thickness;
	int		numgauss_s;
};

struct BeamInfo{
	double crossSectionalArea;
	double momentOfInertia;
};
 
/*!
 Structure for Geometry Properties
 */
struct GeomInfo{
	int dimension;		// dimension of the space	
	
	int numnodes;		// total number of nodes
	int numelements;	// total number of elements
	int numelnodes;		// number of nodes per element
	
	int nodedofs;		// number of dofs per node
	
	double scale;
};

/*!
 Structure for Fluid-Object Coupling
 */
struct CoupleInfo{
	
	int numbeamverts;		// total number of beam nodes
	int numfluidverts;	    // total number of fluid nodes

	GeomInfo geominfo;
};

/*!
 Enumeration for the Approximating Polynomial Type
 */
enum PolynomialOrderEnum {
	linear, 
	quadratic, 
	cubic
};
 
/*!
 Structure conatining the Approximating Polynomial order and the type of Integration used for the Mass Matrix
 */
struct IntegrationInfo{
	PolynomialOrderEnum order;
	int				    numgauss;
	bool				lumpedmass;
};

/*!
 Structure with Constraint Information
 */
struct ConstraintInfo{
	int*				constraints;	// the constrained degrees of freedom
	int					numconstraints; // how many constrained degrees of freedom
};

/*!
 Structure with all of the Analysis Information
 */
struct AnalysisInfo{
	GeomInfo			geominfo;
	MaterialInfo		materialinfo;
	ConstraintInfo		constinfo;
	IntegrationInfo		intinfo;
};


#pragma mark ___Misc___
/*! 
 Returns the current time in seconds
*/
static __inline__ double
currentTime(void)
{
    static double scale = 0.0;

    if (scale == 0.0) {
        mach_timebase_info_data_t info;
        mach_timebase_info(&info);
        scale = info.numer / info.denom * 1.0e-9;
    }

    return mach_absolute_time() * scale;
}

#endif
