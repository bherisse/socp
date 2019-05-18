/*
 * testGoddard.cpp
 *
 *  Created on: February 02, 2017
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */
/*
 * The algorithm is presented in the paper "Singular Arcs in the Generalized Goddard’s Problem", F. Bonnans, P. Martinon, E. Trélat (J Optim Theory Appl (2008) 139: 439–461)
 */
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "../src/socp/shooting.hpp"
#include "../src/models/goddard/goddard.hpp"

// main function
int main(int argc, char** argv) {

	/********************************************************/
	/***************** OCP initialization  ******************/
	/********************************************************/
	std::cout << "Initialize the OCP: " << std::endl;
	
	// new model object
	std::string my_fileTrace("../../trace/goddard/trace.dat");		// file full path for trace
	goddard my_goddard(my_fileTrace);								// The OCP is defined in the goddard class
	int goddardDim = my_goddard.GetDim();							// GetDim return the dimension of the dynamic model

	// new shooting object
	int nMulti = 6;													// multiple shooting parameter (>0)
	shooting my_shooting(my_goddard, nMulti, 1);					// model / number for multiple shooting (1 if simple shooting) / number of threads (1 for one thread)
	my_shooting.SetPrecision(1e-6);									// relative tolerance of the solver

	// vectors to read the trajectory
	std::vector<real> vt(nMulti+1);									// time line (t0 .. t(nMulti+1))
	std::vector<model::mstate> vX(nMulti+1);						// state and costate at times vt

	// set mode for final time and state
	int mode_tf = 1;												// 0 for fixed time, 1 for free optimized time
	std::vector<int> mode_Xf(goddardDim,0);							// 0 for fixed state, 1 for free optimized state
	mode_Xf[3] = 1;	// vxf is free
	mode_Xf[4] = 1;	// vyf is free
	mode_Xf[5] = 1;	// vzf is free
	mode_Xf[6] = 1;	// mf is free
	my_shooting.SetMode(mode_tf, mode_Xf);							// setting mode to the solver						

	// set initial guess for the shooting (trivial guess)
	real ti = 0;					// initial time
	model::mstate Xi(2*goddardDim);	// with the adjoint vector, the full initial state dimension is twice this dimension
	Xi[0] =		0.999949994;		// x
	Xi[1] =		0.0001;				// y
	Xi[2] =		0.01;				// z
	Xi[3] =		1e-10;				// vx			// !=0 to avoid dividing by 0
	Xi[4] =		1e-10;				// vy
	Xi[5] =		1e-10;				// vz	
	Xi[6] =		1.0;				// mass	
	Xi[7] =		0.1;				// p_x			// the initial costate vector is the guess
	Xi[8] =		0.1;				// p_y				
	Xi[9] =		0.1;				// p_z
	Xi[10] =	0.1;				// p_vx
	Xi[11] =	0.1;				// p_vy
	Xi[12] =	0.1;				// p_vz
	Xi[13] =	0.1;				// p_mass
	real tf = 0.1;					// initial guess for tf
	model::mstate Xf(2*goddardDim);	// no need to specify the costate vector for Xf
	Xf[0] =		1.01;				// x
	Xf[1] =		0.0;				// y	
	Xf[2] =		0.0;				// z
	Xf[3] =		0.0;				// vx (free)
	Xf[4] =		0.0;				// vy (free)
	Xf[5] =		0.0;				// vz (free)
	Xf[6] =		0.0;				// mass (free)
	my_shooting.InitShooting(ti, Xi, tf, Xf);		// init shooting with this guess
	
	std::cout << "	done." << std::endl;
	
	/********************************************************/
	/********************* solve OCP ************************/
	/********************************************************/
	std::cout << "Solve the OCP: " << std::endl;
	
	// to return error flag
	int info = 1;
	
	// to return computing time
	long int time1, time2;
	double time;
	time1 = clock();
	
	// start with no drag to compute a better guess from the trivial guess
	my_goddard.GetParameterData().KD = 0.0;

	// solve OCP
	std::cout << "	Solve OCP without drag (KD = 0): ";
	if (info==1)	
		info = my_shooting.SolveOCP(0.0);			// solve OCP without continuation solve OCP without continuation (continuationstep = 0)
	std::cout << "OK = " << info << std::endl;

	// then, start a continuation on KD
	std::cout << "	Continuation on parameter KD (KD = 0 to 310): ";
	if (info==1)	
		info = my_shooting.SolveOCP(1.0, my_goddard.GetParameterData().KD, 310);		// solve OCP with continuation step 1
	std::cout << "OK = " << info << std::endl;

	// then, start a continuation on muL2 
	std::cout << "	Continuation on parameter muL2 (muL2 = 1 to 0.2): ";
	if (info == 1)
		info = my_shooting.SolveOCP(1.0, my_goddard.GetParameterData().muL2, 0.2); // solve OCP with continuation step 1
	std::cout << "OK = " << info << std::endl;

	// to set muL2 = 0, we need to change the structure of the solution (Bang-Singular-Off)
	my_shooting.GetSolution(vt,vX);								// guess from previous solution
	tf = vt[nMulti];											// tf from last computed solution
	real SwitchingTimes_1 = 0.0227; 							// Define switching times for singular arcs (estimated from observation of the previous solution)
	real SwitchingTimes_2 = 0.08;						
	vt[0] = ti; 												// Define a new timeline according to switching times
	vt[1] = SwitchingTimes_1/2; 
	vt[2] = SwitchingTimes_1; 
	vt[3] = (SwitchingTimes_2+SwitchingTimes_1)/2; 
	vt[4] = SwitchingTimes_2; 
	vt[5] = (SwitchingTimes_2+tf)/2;
	vt[6] = tf;
	for(int i=0;i<nMulti+1;i++){
		vX[i] = my_shooting.Move(vt[i]);						// Compute the guess on this new timeline
	}

	// update mode (define the new structure of the trajectory)
	std::vector<int> mode_t(nMulti+1,2);						// define time modes on the timeline		
	mode_t[0] = 0;	// ti (fixed)
	mode_t[2] = 1;	// SwitchingTimes_1 (free)
	mode_t[4] = 1;	// SwitchingTimes_2 (free)
	mode_t[nMulti] = mode_tf;	// tf
	std::vector< std::vector<int> > mode_X(nMulti+1);
	mode_X[0] = std::vector<int>(goddardDim,0); 				// define state modes on the timeline
	mode_X[nMulti] = mode_Xf;									// Xi is fixed
	for (int i=1; i<nMulti; i++){
		mode_X[i] = std::vector<int>(goddardDim, 2);
	}
	my_shooting.SetMode(mode_t, mode_X);						// setting mode to the solver
	
	// init
	my_shooting.InitShooting(vt, vX);							// init shooting with the guess

	// change cost of the model (L1 norm only to maximize the mass)
	my_goddard.GetParameterData().muL2 = 0;
	// set singular control 
	//my_goddard->GetParameterData()->singularControl = 0.6;	// constant approximation for the singular control
	my_goddard.GetParameterData().singularControl = -1;			// if singularControl < 0, then true explicit singular control

	// solve OCP
	std::cout << "	Solve OCP with singular control (muL2 = 0): ";
	if (info==1)	
		info = my_shooting.SolveOCP(0.0);						// solve OCP with continuation step 0 (no continuation possible here)
	std::cout << "OK = " << info << std::endl;

	// print computing time
	time2 = clock();
	time = (double)(time2 - time1) / CLOCKS_PER_SEC;
	std::cout << "	Computing time: " << time << " secondes." << std::endl;

	/********************************************************/
	/****************** trace solution **********************/
	/********************************************************/
	std::cout << "Trace the solution: " << std::endl;
	std::cout << "	Writting trace/goddard/trace.dat" << std::endl;
	my_shooting.Trace(); 	// trace 

	std::cout << "\nPress enter to exit !" << std::endl;

	return 0;
}

