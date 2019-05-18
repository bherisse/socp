/*
* testInterceptor.cpp
*
*  Created on: June 30, 2016
*      Author:  Riccardo BONALLI & Bruno HERISSE (ONERA/DTIS)
*/
/*
* The algorithm is presented in the paper "Analytical Initialization of a Continuation-Based Indirect
Method for Optimal Control of Endo-Atmospheric Launch Vehicle Systems", R. Bonalli, B. Hérissé, E. Trélat (IFAC WC 2017)
*/
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "../src/socp/shooting.hpp"
#include "../src/models/interceptor/interceptor.hpp"

// to initialize the costate
int initState(real & ti, model::mstate & Xi, real & tf, model::mstate & Xf);
// to compute the solution 
int solve(real & ti, model::mstate & Xi, real & tf, model::mstate & Xf, std::string & fileTrace);

// main function
int main(int argc, char** argv) {
	int info;						// to return error flag

	model::mstate Xi(2 * 6);		// initial state		
	model::mstate Xf(2 * 6);		// final state									
	real ti, tf;
	std::string fileTrace;

	/********************************************************/
	/**** compute optimal control  for Scenario 1 ***********/
	/********************************************************/
	ti = 0;
	Xi[0] = 3000;					//h
	Xi[1] = 1000;					//v
	Xi[2] = -M_PI / 6;				//gamma
	Xi[3] = 0.0;					//chi
	Xi[4] = 5454661 / 6378145.0;	//L			
	Xi[5] = 46086 / 6378145.0;		//l
	tf = 20;
	Xf[0] = 12000;					//h
	Xf[1] = 1000;					//v	
	Xf[2] = 0.0;					//gamma
	Xf[3] = M_PI / 8;				//chi
	Xf[4] = 5475000 / 6378145.0;	//L
	Xf[5] = 42000 / 6378145.0;		//l

	// solve OCP
	std::cout << "Solve S1... ";
	fileTrace = std::string("../../trace/interceptor/trace_S1.dat");
	info = solve(ti, Xi, tf, Xf, fileTrace);
	std::cout << "OK = " << info << ", done" << std::endl;

	/********************************************************/
	/**** compute optimal control  for Scenario 2 ***********/
	/********************************************************/
	ti = 0;
	Xi[0] = 3000;					//h
	Xi[1] = 1000;					//v
	Xi[2] = M_PI / 4;				//gamma
	Xi[3] = 0.0;					//chi
	Xi[4] = 5454661 / 6378145.0;	//L			
	Xi[5] = 46086 / 6378145.0;		//l
	tf = 20;
	Xf[0] = 12000;					//h
	Xf[1] = 1000;					//v	
	Xf[2] = -M_PI / 4;					//gamma
	Xf[3] = -M_PI / 2;				//chi
	Xf[4] = 5485000 / 6378145.0;	//L
	Xf[5] = 36178 / 6378145.0;		//l

	// solve OCP
	std::cout << "Solve S2... ";
	fileTrace = std::string("../../trace/interceptor/trace_S2.dat");
	info = solve(ti, Xi, tf, Xf, fileTrace);
	std::cout << "OK = " << info << ", done" << std::endl;

	/********************************************************/
	/**** compute optimal control  for Scenario 3 ***********/
	/********************************************************/
	ti = 0;
	Xi[0] = 3000;					//h
	Xi[1] = 1000;					//v
	Xi[2] = 0.0;					//gamma
	Xi[3] = 0.0;					//chi
	Xi[4] = 5454661 / 6378145.0;	//L			
	Xi[5] = 46086 / 6378145.0;		//l
	tf = 20;
	Xf[0] = 3000;					//h
	Xf[1] = 1000;					//v	
	Xf[2] = 0.0;					//gamma
	Xf[3] = 0.0;					//chi
	Xf[4] = 5485000 / 6378145.0;	//L
	Xf[5] = 46086 / 6378145.0;		//l

	// solve OCP
	std::cout << "Solve S3... ";
	fileTrace = std::string("../../trace/interceptor/trace_S3.dat");
	info = solve(ti, Xi, tf, Xf, fileTrace);
	std::cout << "OK = " << info << ", done" << std::endl;

	return 0;
}

int solve(real & ti, model::mstate & Xi, real & tf, model::mstate & Xf, std::string & fileTrace) {
	int info = 1;

	/********************************************************/
	/***************** OCP initialization  ******************/
	/********************************************************/
	// new model object
	interceptor my_interceptor(fileTrace);							// The OCP is defined in the interceptor class

	// new shooting object
	shooting my_shooting(my_interceptor, 1, 1);						// model / number for multiple shooting (1 if simple shooting) / number of threads (1 for one thread)
	my_shooting.SetPrecision(1e-8);

	int mode_tf = 1;												// 0 for fixed time, 1 for free optimized time
	std::vector<int> mode_X(my_interceptor.GetDim(), 0);			// 0 for fixed state, 1 for free optimized state
	mode_X[1] = 1;	// vf is free
	my_shooting.SetMode(mode_tf, mode_X);							// setting mode to the solver	

	// set initial guess for the shooting 
	model::mstate X0(2 * my_interceptor.GetDim());		
	model::mstate X1(2 * my_interceptor.GetDim());
	real t0, t1;

	info = initState(t0, X0, t1, X1);								// compute first guess from explicit solution and continuation on the model
	my_shooting.InitShooting(t0, X0, t1, X1);						// init shooting with this guess

	/********************************************************/
	/*************** compute optimal control ****************/
	/********************************************************/
	my_shooting.SetDesiredState(ti, Xi, tf, Xf);					// set the desired states

	if (info==1) info = my_shooting.SolveOCP(0.1);					// solve OCP with continuation step 0.1

	/********************************************************/
	/****************** trace solution **********************/
	/********************************************************/
	my_shooting.Trace(); 	// trace

	return info;

}

int initState(real & ti, model::mstate & Xi, real & tf, model::mstate & Xf) {
	// initial and final states for initialization
	ti = 0;
	Xi[0] = 1000;								//h
	Xi[1] = 1000;								//v
	Xi[2] = M_PI / 4;							//gamma
	Xi[3] = 0.0;								//chi
	Xi[4] = 5454661 / 6378145.0;				//L			
	Xi[5] = 46086 / 6378145.0;					//l
	tf = 10;
	Xf[0] = 6000;								//h
	Xf[1] = 1000;								//v	
	Xf[2] = 0.01*M_PI;							//gamma
	Xf[3] = 0.01*M_PI;							//chi
	Xf[4] = (5454661 + 27829.0) / 6378145.0;	//L
	Xf[5] = 46086 / 6378145.0;					//l

	interceptor ini_interceptor(std::string(""));
	int nMulti = 1;
	shooting ini_shooting(ini_interceptor, nMulti, 1);	
	ini_shooting.SetPrecision(1e-8);

	int mode_tf = 1;									
	std::vector<int> mode_X(ini_interceptor.GetDim(), 0);
	mode_X[1] = 1;	// vf is free
	ini_shooting.SetMode(mode_tf, mode_X);							
																	
	model::mstate X0(2 * ini_interceptor.GetDim());
	model::mstate X1(2 * ini_interceptor.GetDim());
	real t0, t1;
	t0 = ti;
	X0 = Xi;
	t1 = tf;
	X1 = Xf;
	ini_interceptor.GetParameterData().mu_gft = 0;			// set the model without considering gravity and propulsion
	ini_interceptor.InitAnalytical(t0, X0, t1, X1);			// compute analytical guess for costate
	ini_shooting.InitShooting(t0, X0, t1, X1);				// init shooting with this guess

	// solve OCP with analytical guess
	int info = ini_shooting.SolveOCP(0.0);
	
	// continuation on parameter mu_gft from 0 to 1
	if (info == 1)	info = ini_shooting.SolveOCP(0.1, ini_interceptor.GetParameterData().mu_gft, 1);

	std::vector<real> vt(nMulti + 1);								// Timeline
	std::vector<model::mstate> vX(nMulti + 1);						// States on the timeline
	ini_shooting.GetSolution(vt, vX);								// Get solution
	ti = vt[0];
	tf = vt[nMulti];
	Xi = vX[0];
	Xf = vX[nMulti];

	return info;

}