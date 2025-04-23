/*
 * testDoubleIntegrator.cpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "../src/socp/shooting.hpp"
#include "../src/models/doubleIntegrator/doubleIntegrator.hpp"

// main function
int main(int argc, char** argv) {

	/********************************************************/
	/***************** OCP initialization  ******************/
	/********************************************************/
	std::cout << "Initialize the OCP: " << std::endl;

	/// create a new model object (a doubleIntegrator here)
	std::string my_fileTrace("../../trace/doubleIntegrator/trace.dat");		// path to the log file
	doubleIntegrator my_doubleIntegrator(1, my_fileTrace);	
	my_doubleIntegrator.SetStepNumber(20);									// set number of integration steps

	/// create a new shooting object for solving the OCP problem
	int nMulti = 1;
	shooting my_shooting(my_doubleIntegrator, nMulti, 1);	// model / number for multiple shooting (1 if simple shooting) / number of threads for multithreading
	my_shooting.SetPrecision(1e-8);							// relative tolerance
	my_shooting.SetContinuationMinStep(1e-12);				// minimum continuation step

	/// set your final constraints
	int mode_tf = 1;											// 0 for fixed final time, 1 for free final time
	my_shooting.SetMode(mode_tf,std::vector<int>(my_doubleIntegrator.GetDim(),0));							// free final time / static final state (default)

	/// set initial guess for the shooting 
	real ti = 0;											// initial time
	doubleIntegrator::mstate Xi(2*my_doubleIntegrator.GetDim());		// GetDim return the dimension of the dynamic model. With the adjoint vector, the full state is twice this dimension
	Xi[0] =		0.0;				//x
	Xi[1] =		0.0;				//y
	Xi[2] =		0.0;				//z
	Xi[3] =		0.0;				//vx
	Xi[4] =		0.0;				//vy
	Xi[5] =		0.0;				//vz		
	Xi[6] =		0.01;				//p_x					// the adjoint vector is the guess (rough guess here)
	Xi[7] =		0.01;				//p_y				
	Xi[8] =		0.01;				//p_z
	Xi[9] =		0.01;				//p_vx
	Xi[10] =	0.01;				//p_vy
	Xi[11] =	0.01;				//p_vz
	real tf = 10;											// initial guess for tf (tf is free)
	doubleIntegrator::mstate Xf(2*my_doubleIntegrator.GetDim());		// No need to specify the adjoint vector for Xf
	Xf[0] =		10.0;				//x
	Xf[1] =		15.0;				//y	
	Xf[2] =		0.0;				//z
	Xf[3] =		0.0;				//vx
	Xf[4] =		0.0;				//vy
	Xf[5] =		0.0;				//vz
	my_shooting.InitShooting(ti, Xi, tf, Xf);				// init shooting with this guess
	std::cout << "	done" << std::endl;

	/********************************************************/
	/********************* solve OCP ************************/
	/********************************************************/
	std::cout << "Solve the OCP: " << std::endl;

	/// to return computing time
	long int time1, time2;
	double time;

	/// print initial parameters for shooting (=tf+adjoint)
	std::vector<real> vt(nMulti+1);							// initial + final time (nMulti = 1)
	std::vector<doubleIntegrator::mstate> vX(nMulti+1);		// initial + final state (nMulti = 1)
	for(int i=0;i<nMulti+1;i++) vX[i] = doubleIntegrator::mstate(2*my_doubleIntegrator.GetDim());
	my_shooting.GetSolution(vt, vX);
	std::cout << "	INIT : ";
	std::cout	<< "tf = " << vt[1] << ", "
				<< "p_x = " << vX[0][my_doubleIntegrator.GetDim()+0] << ", "
				<< "p_y = " << vX[0][my_doubleIntegrator.GetDim()+1] << ", "
				<< "p_z = " << vX[0][my_doubleIntegrator.GetDim()+2] << ", "
				<< "p_vx = " << vX[0][my_doubleIntegrator.GetDim()+3] << ", "
				<< "p_vy = " << vX[0][my_doubleIntegrator.GetDim()+4] << ", "
				<< "p_vz = " << vX[0][my_doubleIntegrator.GetDim()+5] << std::endl;

	/// solve OCP
	std::cout << "	Solve OCP... ";
	time1 = clock();
	int info = my_shooting.SolveOCP(0.0);		// solve OCP with continuation step 0 (no continuation)
	time2 = clock();
	std::cout << "	done";
	time = (double)(time2 - time1) / CLOCKS_PER_SEC;
	std::cout <<  ",	Algo returned " << info << ",";	// the algorithm converged if info = 1
	std::cout <<  "	Computing time : " << time << " secondes." << std::endl;

	// print solution
	my_shooting.GetSolution(vt, vX);
	std::vector<int> nCallFunc = my_shooting.GetCallNumber();
	std::cout  << "	SOL :  (nCallFunc = " << nCallFunc[0] << ", nCallJac = " << nCallFunc[1] <<  ") ";
	std::cout	<< "tf = " << vt[1] << ", "
				<< "p_x = " << vX[0][my_doubleIntegrator.GetDim()+0] << ", "
				<< "p_y = " << vX[0][my_doubleIntegrator.GetDim()+1] << ", "
				<< "p_z = " << vX[0][my_doubleIntegrator.GetDim()+2] << ", "
				<< "p_vx = " << vX[0][my_doubleIntegrator.GetDim()+3] << ", "
				<< "p_vy = " << vX[0][my_doubleIntegrator.GetDim()+4] << ", "
				<< "p_vz = " << vX[0][my_doubleIntegrator.GetDim()+5] << std::endl;

	

	/*************** change initial and/or final time/state optimal control ****************/
	/// if you want to use a different desired state from the previous guess
	Xf[1] = 20;
	my_shooting.SetDesiredState(ti, Xi, tf, Xf);	// change boundaries

	/// solve OCP
	std::cout << "	Change final state... ";
	time1 = clock();
	info = my_shooting.SolveOCP(1.0);		// solve OCP with continuation step 1 
	time2 = clock();
	std::cout << "	done";
	time = (double)(time2 - time1) / CLOCKS_PER_SEC;
	std::cout <<  ",	Algo returned " << info << ",";	// the algorithm converged if info = 1
	std::cout << "	Computing time : " << time << " secondes." << std::endl;

	// print solution
	my_shooting.GetSolution(vt, vX);
	std::cout << "	SOL :  ";
	std::cout	<< "tf = " << vt[1] << ", "
				<< "p_x = " << vX[0][my_doubleIntegrator.GetDim()+0] << ", "
				<< "p_y = " << vX[0][my_doubleIntegrator.GetDim()+1] << ", "
				<< "p_z = " << vX[0][my_doubleIntegrator.GetDim()+2] << ", "
				<< "p_vx = " << vX[0][my_doubleIntegrator.GetDim()+3] << ", "
				<< "p_vy = " << vX[0][my_doubleIntegrator.GetDim()+4] << ", "
				<< "p_vz = " << vX[0][my_doubleIntegrator.GetDim()+5] << std::endl;

	/*********** continuation on a model parameter **********/
	// model parameters can be changed using continuation 
	std::cout  << "	Continuation on parameter muT... ";		// increase muT to obtain minimum time solution

	/// solve OCP
	time1 = clock();
	if (info==1)	info = my_shooting.SolveOCP(1.0, my_doubleIntegrator.GetParameterData().muT, 0.02);		// solve OCP with continuation step 1 (no continuation)
	time2 = clock();
	std::cout << "	done";
	time = (double)(time2 - time1) / CLOCKS_PER_SEC;
	std::cout << ",	Algo returned " << info << ",";	// the algorithm converged if info = 1
	std::cout << "	Computing time : " << time << " secondes." << std::endl;

	// print solution
	my_shooting.GetSolution(vt, vX);
	std::cout << "	SOL :  ";
	std::cout	<< "tf = " << vt[1] << ", "
				<< "p_x = " << vX[0][my_doubleIntegrator.GetDim()+0] << ", "
				<< "p_y = " << vX[0][my_doubleIntegrator.GetDim()+1] << ", "
				<< "p_z = " << vX[0][my_doubleIntegrator.GetDim()+2] << ", "
				<< "p_vx = " << vX[0][my_doubleIntegrator.GetDim()+3] << ", "
				<< "p_vy = " << vX[0][my_doubleIntegrator.GetDim()+4] << ", "
				<< "p_vz = " << vX[0][my_doubleIntegrator.GetDim()+5] << std::endl;

	/********************************************************/
	/****************** trace solution **********************/
	/********************************************************/
	std::cout << "Trace the solution: " << std::endl;
	std::cout << "	Writting trace/doubleIntegrator/trace.dat" << std::endl;
	my_shooting.Trace(); 	// trace
	
	std::cout << "\nPress enter to exit !" << std::endl;

	return 0;
}