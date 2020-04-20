/*
 * testCovid19.cpp
 *
 *  Created on: April 19, 2020
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */
/*
 * Using a SEIR model for the covid-19 epidemic, the problem considered here is the following:
 *		- Control lockdown with minimum effort (in a quadratic way). This is motivated by the fact that lockdown effort has a great impact on economy.
 *		- One year after end of lockdown (2020, May 11 in France), the total removed (recovered + isolated + fatalities) population achieves 80% (which corresponds approximatively to natural herd immunity).
 *
 * Even for this "worst case" scenario (80% herd immunity will induce many fatalities), results confirm that end of lockdown sould be very progressive to restrain the number of infections actively circulating along the next year. 
 *
 * Note that initial covid-19 parameters should be refined.
 *
 * Future work: control the reproduction number to minimize lockdown effort while constraining the number of infected people.
 * 
 */
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "../src/socp/shooting.hpp"
#include "../src/models/covid19/covid19.hpp"

// main function
int main(int argc, char** argv) {

	/********************************************************/
	/***************** OCP initialization  ******************/
	/********************************************************/
	std::cout << "Initialize the Optimal Control Problem: " << std::endl;
	
	// new model object
	std::string my_fileTrace("../../trace/covid19/trace.dat");		// file full path for trace
	covid19 my_covid19(my_fileTrace);								// The OCP is defined in the covid19 class
	int covid19Dim = my_covid19.GetDim();							// GetDim return the dimension of the dynamic model

	// parameters of SEIR model for covid-19
	my_covid19.GetParameterData().R0 = 4;
	my_covid19.GetParameterData().Tinf = 10;
	my_covid19.GetParameterData().Tinc = 5;

	// new shooting object
	int nMulti = 20;												// multiple shooting parameter (>0) to tackle sensitivity of the problem
	shooting my_shooting(my_covid19, nMulti, 4);					// model / number for multiple shooting (1 if simple shooting) / number of threads (1 for one thread)
	my_shooting.SetPrecision(1e-8);	

	// vectors to read the trajectory
	std::vector<real> vt(nMulti+1);									// time line (t0 .. t(nMulti+1))
	std::vector<model::mstate> vX(nMulti+1);						// state and costate at times vt

	// set mode for final time and state
	int mode_tf = 0;												// 0 for fixed time, 1 for free optimized time
	std::vector<int> mode_Xf(covid19Dim,0);							// 0 for fixed state, 1 for free optimized state
	mode_Xf[0] = 1;	// S is free
	mode_Xf[1] = 1;	// E is free
	mode_Xf[2] = 1;	// I is free
	mode_Xf[3] = 0;	// R is fixed
	my_shooting.SetMode(mode_tf, mode_Xf);							// setting mode to the solver						

	// set initial guess for the shooting (trivial guess)
	real ti = 0;					// initial time (corresponding to May 11)
	model::mstate Xi(2* covid19Dim);	// with the adjoint vector, the full initial state dimension is twice this dimension
	Xi[0] =		0.9;				// S (estimate on May 11)
	Xi[1] =		0.003;				// E (estimate on May 11)
	Xi[2] =		0.01;				// I (estimate on May 11)
	Xi[3] =		0.087;				// R (estimate on May 11)
	Xi[4] =		-0.001;				// pS (parameter guess for optimization)
	Xi[5] =		0.001;				// pE (parameter guess for optimization)
	Xi[6] =		0*0.01;				// pI (parameter guess for optimization)
	Xi[7] =		0;					// pR (parameter guess for optimization)	
	real tf = 30;					// final simulation time (start with 30 days to make it work)
	model::mstate Xf(2* covid19Dim);	// no need to specify the costate vector for Xf
	Xf[3] =		0.6;					// target for R (start with 60% to make it work)
	my_shooting.InitShooting(ti, Xi, tf, Xf);		// init shooting with this guess
	
	std::cout << "	done." << std::endl;
	
	/********************************************************/
	/********************* solve OCP ************************/
	/********************************************************/
	std::cout << "Solve the OCP: " << std::endl;
	
	// to return error flag
	int info = 1;

	// solve OCP with initial guess
	std::cout << "	Solve OCP with initial guess: ";
	if (info==1)	
		info = my_shooting.SolveOCP(0.0);			// solve OCP without continuation solve OCP without continuation (continuationstep = 0)
	std::cout << "OK = " << info << std::endl;

	// solve OCP to achieve 80% of removed population after 30 days
	Xf[3] = 0.8;
	std::cout << "	Solve OCP to achieve 80% of removed population after 30 days: ";
	my_shooting.SetDesiredState(ti, Xi, tf, Xf);
	info = my_shooting.SolveOCP(0.1);
	std::cout << "OK = " << info << std::endl;

	// solve OCP to achieve 80% of removed population after 1 year
	tf = 365;
	std::cout << "	Solve OCP to achieve 80% of removed population after 1 year: ";
	my_shooting.SetDesiredState(ti, Xi, tf, Xf);
	info = my_shooting.SolveOCP(0.01);
	std::cout << "OK = " << info << std::endl;

	/********************************************************/
	/****************** trace solution **********************/
	/********************************************************/
	std::cout << "Trace the solution: " << std::endl;
	std::cout << "	Writting trace/covid19/trace.dat" << std::endl;
	my_shooting.Trace(); 	// trace 

	std::cout << "\nPress enter to exit !" << std::endl;

	return 0;
}

