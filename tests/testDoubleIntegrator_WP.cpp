/*
 * testDoubleIntegrator_WP.cpp
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

	/// new map object (environment parameters, obstacles, etc.)
	//map* my_map = new map();
	/// new model object
	std::string my_fileTrace("../../trace/doubleIntegrator/trace.dat");
	doubleIntegrator my_doubleIntegrator = doubleIntegrator(1, my_fileTrace);
	// new shooting object
	int nMulti = 2;
	shooting my_shooting = shooting(my_doubleIntegrator, nMulti, 2);	// model / number for multiple shooting (1 if simple shooting)
	int mode_tf = 1;													// 0 for fixed final time, 1 for free final time
	
	// Set mode
	std::vector<int> mode_t(nMulti+1);
	for(int i=0;i<nMulti+1;i++) mode_t[i] = 1;
	mode_t[0] = 0;	// t0 fixed
	mode_t[nMulti] = 1;
	std::vector< std::vector<int> > mode_X(nMulti+1);
	mode_X[0] = std::vector<int>(my_doubleIntegrator.GetDim(),0);		// fixed
	//mode_X[0][0] = 1; // to test free initial state
	mode_X[nMulti] = std::vector<int>(my_doubleIntegrator.GetDim(),0); // fixed
	//mode_X[nMulti][3] = 1;	// vf is free
	//mode_X[nMulti][4] = 1;	// vf is free
	//mode_X[nMulti][5] = 1;	// vf is free
	for (int i=1; i<nMulti; i++){ 
		mode_X[i] = std::vector<int>(my_doubleIntegrator.GetDim(),0);	// fixed
		mode_X[i][3] = 2;	// v is free
		mode_X[i][4] = 2;	// v is free
		mode_X[i][5] = 2;	// v is free
	}
	my_shooting.SetMode(mode_t, mode_X);
	
	// set initial guess for the shooting 
	real t0 = 0;											// initial guess for ti
	model::mstate X0(2*my_doubleIntegrator.GetDim());		// GetDim return the dimension of the dynamic model. With the adjoint vector, the full state is twice this dimension
	X0[0] =		0.0;				//x
	X0[1] =		0.0;				//y
	X0[2] =		0.0;				//z
	X0[3] =		0.0;				//vx
	X0[4] =		0.0;				//vy
	X0[5] =		0.0;				//vz		
	X0[6] =		0.001;				//p_x					// the adjoint vector is the guess
	X0[7] =		0.001;				//p_y				
	X0[8] =		0.001;				//p_z
	X0[9] =		0.001;				//p_vx
	X0[10] =	0.001;				//p_vy
	X0[11] =	0.001;				//p_vz
	real t1 = 30;											// initial guess for tf
	model::mstate X1(2*my_doubleIntegrator.GetDim());		// No need to specify the adjoint vector for Xf
	X1[0] =		10.0;				//x
	X1[1] =		0.0;				//y	
	X1[2] =		0.0;				//z
	X1[3] =		0.0;				//vx
	X1[4] =		0.0;				//vy
	X1[5] =		0.0;				//vz
	X1[6] =		0.001;				//p_x					// the adjoint vector is the guess
	X1[7] =		0.001;				//p_y				
	X1[8] =		0.001;				//p_z
	X1[9] =		0.001;				//p_vx
	X1[10] =	0.001;				//p_vy
	X1[11] =	0.001;				//p_vz
	real t2 = 60;
	model::mstate X2(2*my_doubleIntegrator.GetDim());		// No need to specify the adjoint vector for Xf
	X2[0] =		20;				//x
	X2[1] =		0.0;				//y	
	X2[2] =		0.0;				//z
	X2[3] =		0.0;				//vx
	X2[4] =		0.0;				//vy
	X2[5] =		0.0;				//vz
	
	std::vector<real> vt(nMulti+1);								// Timeline
	vt[0] = t0; vt[1] = t1; vt[2] = t2;
	//for(int i=0;i<nMulti+1;i++) vt[i] = i*tf/nMulti;
	std::vector<model::mstate> vX(nMulti+1);					// States on the timeline
	vX[0] = X0; vX[1] = X1; vX[2] = X2;
	//for(int i=1;i<nMulti;i++){
	//	vX[i] = my_shooting.Move(vt[0], X0, vt[i], 0);
	//}
	my_shooting.InitShooting(vt, vX);					// init shooting with this guess
	std::cout << "	done" << std::endl;

	/********************************************************/
	/********************* solve OCP ************************/
	/********************************************************/
	std::cout << "Solve the OCP: " << std::endl;

	// to return computing time
	long int time1, time2;
	double time;
	
	// if you want to use a different desired state from the intial guess
	//my_shooting.SetDesiredState(ti, Xi, tf, Xf);	

	// solve OCP
	std::cout << "	Solve OCP... ";
	time1 = clock();
	int info = my_shooting.SolveOCP(0.0);		// solve OCP with continuation step 0

	// change mission objective
	vX[1][1] =	15.0;
	//vX[1][2] =	15.0;
	//vX[2][2] =	-15.0;
	vX[2][1] = 5.0;
	vX[2][2] = 10.0;
	my_shooting.SetDesiredState(vt, vX);
	info = my_shooting.SolveOCP(1.0);		// solve OCP with continuation step 1

	time2 = clock();
	time = (double)(time2 - time1) / CLOCKS_PER_SEC;
	std::cout << "	Algo returned " << info << ", ";
	std::cout << "	Computing time : " << time << " secondes." << std::endl;

	// print solution
	my_shooting.GetSolution(vt, vX);
	std::vector<int> nCallFunc = my_shooting.GetCallNumber();
	std::cout << "	SOL :  (nCallFunc = " << nCallFunc[0] << ", nCallJac = " << nCallFunc[1]<< ") ";
	std::cout << "tf = " << vt[1] << ", "
		<< "p_x = " << vX[0][my_doubleIntegrator.GetDim() + 0] << ", "
		<< "p_y = " << vX[0][my_doubleIntegrator.GetDim() + 1] << ", "
		<< "p_z = " << vX[0][my_doubleIntegrator.GetDim() + 2] << ", "
		<< "p_vx = " << vX[0][my_doubleIntegrator.GetDim() + 3] << ", "
		<< "p_vy = " << vX[0][my_doubleIntegrator.GetDim() + 4] << ", "
		<< "p_vz = " << vX[0][my_doubleIntegrator.GetDim() + 5] << std::endl;

	/********************************************************/
	/*********** continuation on a model parameter **********/
	/********************************************************/
	// model parameters can be changed using continuation 
	std::cout << "	Continuation on parameter muT... ";
	if (info==1)	info = my_shooting.SolveOCP(1.0, my_doubleIntegrator.GetParameterData().muT, 0.02);
	std::cout << "OK = " << info << ", done" << std::endl;

	/********************************************************/
	/****************** trace solution **********************/
	/********************************************************/
	std::cout << "Trace the solution: " << std::endl;
	std::cout << "	Writting trace/doubleIntegrator/trace.dat" << std::endl;
	my_shooting.Trace(); 	// trace

	std::cout << "\nPress enter to exit !" << std::endl;

	return 0;
}

