/*
 *  NewmarkSolve.h
 *  FiniteElement
 *
 *  Created by Cynthia Bruyns on 1/10/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef NewmarkSolve_H_
#define NewmarkSolve_H_

#include "CLPKRealMatrixOpsCM.h"

void newmarkSteup(	double* Kc, double* Mc, 
					double* f0,  double* u0, 
					int cdofs, 
					int* bcs, int numbcs, 
					int numsteps	) //cbnote not much use if u_t doesnt come back out
{
	double dt = 1.0/120.0;

	double delta = 1.0/2.0;
	double alpha = 1.0/4.0;

	double a0 = 1.0/alpha/powf(dt,2.0);	
	double a1 = delta/alpha/dt;
	double a2 = 1.0/alpha/dt;
	double a3 = (1.0/2.0/alpha) - 1.0;
	double a4 = (delta/alpha) -1.0;
	double a5 = dt/2.0 * ((delta/alpha) - 2.0);
	double a6 = dt*(1.0-delta);
	double a7 = delta*dt;
	
	int dsize = sizeof(double);
	int bsize = cdofs*dsize;

	double* Keff = (double*) calloc(cdofs*cdofs, dsize);
	double* temp = (double*) calloc(cdofs*cdofs, dsize);
	realDoubleMatrixScale(temp, Mc, a0, cdofs);
	realDoubleAdd2Matricies(Keff, Kc, temp, cdofs);//Keff = Kc+a0*Mc;

	double* Keffcopy = (double*) calloc(cdofs*cdofs, dsize);
	
#ifdef DEBUG_PRINT
	for (int i=0; i<cdofs; i++){	
		for (int j=0; j<cdofs; j++){		
			if ( Keff[i + j*cdofs] != 0.0) printf("Keff[%d][%d]: %lf\n", i,j,Keff[i + j*cdofs]);					
		}		
		printf("\n");
	}
#endif

	double* u_t			= (double*) calloc(cdofs,dsize);
	double* u_t_dt		= (double*) calloc(cdofs,dsize);
	
	double* udot_t		= (double*) calloc(cdofs,dsize);
	double* udot_t_dt	= (double*) calloc(cdofs,dsize);

	double* uddot_t		= (double*) calloc(cdofs,dsize);
	double* uddot_t_dt	= (double*) calloc(cdofs,dsize);

	double* temp1 = (double*) calloc(cdofs, dsize);	
	double* temp2 = (double*) calloc(cdofs, dsize);	
	double* temp3 = (double*) calloc(cdofs, dsize);	
	double* temp4 = (double*) calloc(cdofs, dsize);	 
	double* temp5 = (double*) calloc(cdofs, dsize);	
	double* temp6 = (double*) calloc(cdofs, dsize);	
	double* temp7 = (double*) calloc(cdofs, dsize);	
	double* temp8 = (double*) calloc(cdofs, dsize);	
	double* temp9 = (double*) calloc(cdofs, dsize);	
	double* temp10 = (double*) calloc(cdofs, dsize);
	double* temp11 = (double*) calloc(cdofs, dsize);
	double* temp12 = (double*) calloc(cdofs, dsize);
	
	double* Feff_t = (double*) calloc(cdofs, dsize);
	
	for (int i=0; i<numsteps; i++){
	
	double time = currentTime();
	
	
		realDoubleVectorScale(temp1, u_t_dt, a0, cdofs);			
	#ifdef DEBUG_PRINT
		printf("\n");
		printf("temp1\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n", temp1[j]);
	#endif		

		realDoubleVectorScale(temp2, udot_t_dt, a2, cdofs);		
	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("temp2\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n", temp2[j]);
	#endif

		realDoubleVectorScale(temp3, uddot_t_dt, a3, cdofs);	
	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("temp3\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n",  temp3[j]);
	#endif		

		realDoubleVectorAdd(temp4, temp1, temp2, cdofs);
	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("temp4\n");
		for (int j=0; j<cdofs; j++)
			printf("%lf\n", temp4[j]);
	#endif		

		realDoubleVectorAdd(temp5, temp4, temp3, cdofs);
	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("temp5[%d]:\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n", temp5[j]);
	#endif		

		realDoubleMatrixVectorMult(temp6, Mc, temp5, cdofs); // M*();		
	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("temp6\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n", temp6[j]);
	#endif		
		realDoubleVectorAdd(Feff_t, f0, temp6, cdofs);//Feff_t = f0 + M*();

		for (int j=0; j<ctinfo.numconstraints; j++)
			Feff_t[ctinfo.constraints[j]] = 0;

	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("Feff_t\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n", Feff_t[j]);
	#endif

		//	realDoubleMatrixVectorMult(u_t, Keff_inv, Feff_t, cdofs); //u_t = Keff\Feff_t;
		//	solve Keff*u_t = Feff_t
		memcpy(u_t, Feff_t, bsize);	
		memcpy(Keffcopy, Keff, bsize*cdofs);
		realDoubleSolve(Keffcopy, u_t, cdofs); //u_t is Feff_t copy

	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("u_t\n");
		for (int j=0; j<cdofs; j++)
			printf(" %e\n", u_t[j]);
	#endif
			
		realDoubleVectorScale(temp7, u_t, a0, cdofs);		

	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("temp7\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n", temp7[j]);
	#endif			

		realDoubleVectorSubtract(temp8, temp7, temp1, cdofs);	//a0*(u_t - u_t_dt)

	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("temp8\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n", temp8[j]);
	#endif

		realDoubleVectorSubtract(temp9, temp8, temp2, cdofs);	//- a2*udot_t_dt

	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("temp9\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n", temp9[j]);
	#endif

		realDoubleVectorSubtract(uddot_t, temp9, temp3, cdofs);//- a3*uddot_t_dt

	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("uddot_t\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n", uddot_t[j]);
		//uddot_t = a0*(u_t - u_t_dt) - a2*udot_t_dt - a3*uddot_t_dt;
	#endif

		realDoubleVectorScale(temp10, uddot_t_dt, a6, cdofs);	

	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("temp10\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n", temp10[j]);
	#endif

		realDoubleVectorScale(temp11, uddot_t, a7, cdofs);	

	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("temp11\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n", temp11[j]);
	#endif

		realDoubleVectorAdd(temp12, temp10, temp11, cdofs);// a6*uddot_t_dt + a7*uddot_t;

	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("temp12\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n", temp12[j]);
	#endif

		realDoubleVectorAdd(udot_t, udot_t_dt, temp12, cdofs);//udot_t  = udot_t_dt + a6*uddot_t_dt + a7*uddot_t;

	#ifdef DEBUG_PRINT	
		printf("\n");
		printf("udot_t\n");
		for (int j=0; j<cdofs; j++)
			printf(" %lf\n", udot_t[j]);
	#endif
			
		memcpy(u_t_dt,u_t,bsize);
		memcpy(udot_t_dt,udot_t,bsize);
		memcpy(uddot_t_dt,uddot_t,bsize);

		time = currentTime() - time;
		printf("total newmark method to solve %d equations =  %lf ms\n", cdofs, time*1E3);

		memset(f0, 0, bsize);	
		
	
	}

	free(temp);
	
	free(u_t), free(u_t_dt), free(udot_t), free(udot_t_dt), free(uddot_t), free(uddot_t_dt);

	free(Keff), free(Keffcopy);
	free(Feff_t);
	free(temp1), free(temp2), free(temp3), free(temp4), free(temp5), free(temp6), free(temp7), free(temp8), free(temp9), free(temp10), free(temp11), free(temp12);

}

	
	
	
	

#endif
