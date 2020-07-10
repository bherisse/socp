/*
 * testVtolUAV.cpp
 *
 *  Created on: June 30, 2019
 *      Author: Bruno HERISSE
 */

#include <iostream>     // std::cout
#include <stdlib.h>
#include <fstream>      // std::ifstream
#include <sstream>
#include <string>
#include <time.h>
#include <math.h>

#include "../src/socp/shooting.hpp"
#include "../src/models/vtolUAV/vtolUAV.hpp"
#include "../src/maps/obstacle/obstacle.hpp"

// modeMPP
int modeMPP = 0; // 0: (MPP)_{1,sigma} ; 1: (MPP)_{1,0} ; 2: (MPP)_{1,inf}

// data
std::string my_fileWP("../../data/vtolUAV/waypoints");
std::string my_fileObstacles("../../data/vtolUAV/obstacles");
std::string my_fileTrace("../../trace/vtolUAV/trace.dat");

// functions
void ReadInput(std::string const& fileMCWP, std::vector<std::vector<real>> & pathStateWP, std::vector<real> & pathDistanceWP);
int PathContinuation(vtolUAV & _votlUAV, shooting & _shooting_wp, std::vector<std::vector<real>> const& pathStateWP, std::vector<long int> & timeWP);

// main function
int main(int argc, char** argv){

	/********************************************************/
	/***************** Read args ****************************/
	/********************************************************/
	real sigma=2, mu=1;
	if (argc > 1 && argc <= 4) {
		modeMPP = atoi(argv[1]);		//	modeMPP
		sigma = atof(argv[2]);			// sigma
		if (sigma == 0) {
			std::cout << "sigma must be greater than 0" << std::endl;
			return 0;
		}
		mu = atof(argv[3]);				// mu	
	}
	else {
		std::cout << "Number of arguments must be equal to 3, ";
		std::cout << "running with default values." << std::endl;
	}

	/********************************************************/
	/***************** OCP initialization  ******************/
	/********************************************************/
	std::cout << "Initialization... ";
	/// new map object (environment parameters, obstacles, etc.)
	obstacle my_obstacle = obstacle(my_fileObstacles, my_fileWP);
	/// new model object
	vtolUAV my_vtolUAV(my_obstacle, my_fileTrace);
	shooting my_shooting_wp(my_vtolUAV, 1, 1);
	my_shooting_wp.SetPrecision(1e-4);
	my_vtolUAV.SetODEIntPrecision(1e-5);

	// load path as waypoints
	std::vector<std::vector<real>> pathStateWP;				// waypoints position
	std::vector<real> pathDistanceWP;						// distance between wp
	ReadInput(my_fileWP, pathStateWP, pathDistanceWP);
	std::cout << " done" << " (modeMPP = " << modeMPP << ",	sigma = " << sigma << ",	mu = " << mu << ")" << std::endl;

	// to return computing time
	std::vector<long int> timeWP(pathStateWP.size(), 0);
	long int time1, time2;
	double computingTime;
	time1 = clock();
	timeWP[0] = time1;

	int info = 1;

	/********************************************************/
	/*********** Continuation along the loaded path *********/
	/********************************************************/
	// launch trajectory planning
	std::cout << "Continuation along the path... ";
	if (info == 1)
		info = PathContinuation(my_vtolUAV, my_shooting_wp, my_obstacle.GetPath(), timeWP);
	std::cout << "OK = " << info << ", done" << std::endl;

	/********************************************************/
	/*********** continuation on a model parameter **********/
	/********************************************************/
	std::cout << "Continuation on parameter invSigma... ";
	if (info == 1) info = my_shooting_wp.SolveOCP(0.1, my_vtolUAV.GetParameterData().invSigmaXwp, 1/ sigma);
	std::cout << "OK = " << info << ", done" << std::endl;

	std::cout << "Continuation on parameter muObs... "; 
	if (info == 1) info = my_shooting_wp.SolveOCP(0.1, my_obstacle.GetParameterData().muObs, mu);
	std::cout << "OK = " << info << ", done" << std::endl;

	//std::cout << "Continuation on parameter u_max... ";
	if (info == 1) info = my_shooting_wp.SolveOCP(1.0, my_vtolUAV.GetParameterData().u_max, 1.0);
	std::cout << "OK = " << info << ", done" << std::endl;

	// print computing time
	time2 = clock();
	computingTime = (double)(time2 - time1) / CLOCKS_PER_SEC;
	std::cout << "Computing time : " << computingTime << " secondes." << std::endl;

	/********************************************************/
	/*************** trace results **************************/
	/********************************************************/
	// print the trajectory in my_fileTrace from the last obtained solution
	my_shooting_wp.Trace(); 	// trace

	return 0;
}

/**
* Read parameters in file
*/
void ReadInput(std::string const& fileMCWP, std::vector<std::vector<real>> & pathStateWP, std::vector<real> & pathDistanceWP) {
	std::ifstream is(fileMCWP);
	int nWP;
	if (is) {
		std::string line;
		// get n
		std::getline(is, line); // get "n_wp:"
		std::getline(is, line); // get the integer n
		std::istringstream isn(line);
		if (!(isn >> nWP)) { return; } // error

		// get positions
		pathStateWP = std::vector<std::vector<real>>(nWP);
		std::getline(is, line); // get "position_wp:"
		for(int i=0;i<nWP;i++){
			pathStateWP[i] = std::vector<real>(6,0);
			std::getline(is, line);
			std::istringstream ispos(line);
			if (!(ispos >> pathStateWP[i][0] >> pathStateWP[i][1] >> pathStateWP[i][2])) { return; } // error
		}

		// get distances and directions
		pathDistanceWP = std::vector<real>(nWP,0);
		for (int i = 0; i<nWP-1; i++) {
			std::vector<real> deltaX(3);
			deltaX[0] = pathStateWP[i + 1][0] - pathStateWP[i][0];
			deltaX[1] = pathStateWP[i + 1][1] - pathStateWP[i][1];
			deltaX[2] = pathStateWP[i + 1][2] - pathStateWP[i][2];
			pathDistanceWP[i+1] = pathDistanceWP[i] + sqrt(deltaX[0] * deltaX[0] + deltaX[1] * deltaX[1] + deltaX[2] * deltaX[2]);
			pathStateWP[i][3] = deltaX[0] / pathDistanceWP[i];
			pathStateWP[i][4] = deltaX[1] / pathDistanceWP[i];
			pathStateWP[i][5] = deltaX[2] / pathDistanceWP[i];
		}
	}
}

/**
* Path Continuation
*/
int PathContinuation(vtolUAV & _votlUAV, shooting & _shooting_wp, std::vector<std::vector<real>> const& pathStateWP, std::vector<long int> & timeWP) {
	real modelDim = _shooting_wp.GetModel().GetDim();

	// shooting object
	int nMulti = 1; //n_wp-1; // number of points for multi-shooting
	int nThread = 4; // number of threads for parallel computing
	_shooting_wp.Resize(nMulti, nThread);

	// Setting modeMPP
	std::vector<int> mode_t(nMulti+1);
	mode_t[0] = 0;
	for(int i=1;i<nMulti;i++) mode_t[i] = 1 + modeMPP;
	mode_t[nMulti] = 1;	// tf free
	std::vector< std::vector<int> > mode_X(nMulti+1);
	mode_X[0] = std::vector<int>(modelDim,0);
	//mode_X[0][3] = 1; mode_X[0][4] = 1; mode_X[0][5] = 1; //v0 free
	mode_X[nMulti] = std::vector<int>(modelDim,0);
	mode_X[nMulti][3] = 1; mode_X[nMulti][4] = 1; mode_X[nMulti][5] = 1; //vf free
	for (int i=1; i<nMulti; i++){
		mode_X[i] = std::vector<int>(modelDim,2);
		if (modeMPP==0){ mode_X[i][0] = 0; mode_X[i][1] = 0; mode_X[i][2] = 0;}
	}
	_shooting_wp.SetMode(mode_t, mode_X);
	// guess 
	std::vector<real> vt(nMulti+1);								// Timeline
	vt[0] = 0;
	std::vector<model::mstate> vX(nMulti + 1);
	for(int i=0;i<nMulti+1;i++){
		vX[i] = std::vector<real>(12,0);
		vX[i][0] = pathStateWP[i][0];
		vX[i][1] = pathStateWP[i][1];
		vX[i][2] = pathStateWP[i][2];
		vX[i][3] = 0.001*pathStateWP[i][3];	// to avoid dividing by 0, should be != 0
		vX[i][4] = 0.001*pathStateWP[i][4];
		vX[i][5] = 0.001*pathStateWP[i][5];
	}
	// guess first solution
	std::vector<real> dX(3);
	dX[0] = (vX[nMulti][0] - vX[0][0]);
	dX[1] = (vX[nMulti][1] - vX[0][1]);
	dX[2] = (vX[nMulti][2] - vX[0][2]);
	real normDX = sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
	vt[nMulti] = pow(4.5*normDX*normDX / _votlUAV.GetParameterData().alphaT, 0.25);
	vX[0][6] = -3 * dX[0] / vt[nMulti] / vt[nMulti] / vt[nMulti];
	vX[0][7] = -3 * dX[1] / vt[nMulti] / vt[nMulti] / vt[nMulti];
	vX[0][8] = -3 * dX[2] / vt[nMulti] / vt[nMulti] / vt[nMulti];
	vX[0][9] = -3 * dX[0] / vt[nMulti] / vt[nMulti];
	vX[0][10] = -3 * dX[1] / vt[nMulti] / vt[nMulti];
	vX[0][11] = -3 * dX[2] / vt[nMulti] / vt[nMulti];
	//std::cout << vX[0][6] << "\t" << vX[0][7] << "\t" << vX[0][8] << "\t" << vX[0][9] << "\t" << vX[0][10] << "\t" << vX[0][11] << "\t" << std::endl;
	_shooting_wp.InitShooting(vt, vX);				// init shooting with the guess
	//_shooting_wp.SetDesiredState(0, vX[0], vt[nMulti], vX[nMulti]);

	/// loop
	int info = 1; 
	int Nwp = nMulti;
	int end = pathStateWP.size();
	double lambda = 0;
	while(info==1 && Nwp<end){
		//std::cout << "nMulti =  " << nMulti << std::endl;

		//lambda += pathDistanceWP[Nwp -1];
		//std::cout << "lambda =  " << lambda << std::endl;

		/// Compute the traj with multi shooting from x_wp[0] (start) to x_wp[iter]
		info = _shooting_wp.SolveOCP(1.0);		// solve OCP 
		if (nMulti == 1 && info) {
			info = _shooting_wp.SolveOCP(1, _votlUAV.GetParameterData().ca, 0.05);
			//info = _shooting_wp.SolveOCP(1, _votlUAV.GetParameterData().alphaV, 0.00);
		}
		if (info!=1){
			std::cout << "failed iter " << nMulti << std::endl;
			break;
		}
		_shooting_wp.GetSolution(vt,vX);

		// get time
		timeWP[Nwp] = clock();

		/// update for next iteration
		nMulti = nMulti +1;		// to comment for no multi
		Nwp += 1;
		if (Nwp < end) {
			vt.resize(nMulti + 1);					// Timeline
			vX.resize(nMulti + 1);					// States on the timeline
			mode_t.resize(nMulti + 1);
			mode_X.resize(nMulti + 1);

			if (modeMPP > 1) { // (MPP)_{1,inf}
				real tf = vt[nMulti - 1] + 0.5 * 10;		// to comment for no multi
				for (int i = 0; i < nMulti + 1; i++) vt[i] = i*tf / nMulti;		// to comment for no multi
				for (int i = 0; i < nMulti + 1; i++) {
					vX[i] = _shooting_wp.Move(vt[i]);
				}
				// Setting modeMPP
				mode_t[0] = 0;
				for (int i = 1; i < nMulti; i++) mode_t[i] = 2;
				mode_t[nMulti] = 1;	// tf free
				mode_X[0] = std::vector<int>(modelDim, 0);
				//mode_X[0][3] = 1; mode_X[0][4] = 1; mode_X[0][5] = 1; //v0 free
				mode_X[nMulti] = std::vector<int>(modelDim, 0);
				if (nMulti < end - 1) {
					mode_X[nMulti][3] = 1; mode_X[nMulti][4] = 1; mode_X[nMulti][5] = 1; //vf free
				}
				for (int i = 1; i < nMulti; i++) {
					mode_X[i] = std::vector<int>(modelDim, 2);
				}
			}else{
				vt[nMulti] = vt[nMulti - 1] + 0.3 * 10;
				vX[nMulti] = _shooting_wp.Move(vt[nMulti]);
				for (int i = 1; i < nMulti; i++) mode_t[i] = 1; // free times on the timeline
				mode_t[nMulti] = 1;	// tf free
				mode_X[0] = std::vector<int>(modelDim, 0);
				//mode_X[0][3] = 1; mode_X[0][4] = 1; mode_X[0][5] = 1; //v0 free
				mode_X[nMulti] = std::vector<int>(modelDim, 0);
				if (nMulti < end - 1) {
					mode_X[nMulti][3] = 1; mode_X[nMulti][4] = 1; mode_X[nMulti][5] = 1; //vf free
				}
				for (int i = 1; i < nMulti; i++) {
					mode_X[i] = std::vector<int>(modelDim, 2);
					mode_X[i][0] = 0; mode_X[i][1] = 0; mode_X[i][2] = 0;
					if (modeMPP == 0) {
						mode_X[i][0] = 1; mode_X[i][1] = 1; mode_X[i][2] = 1;		// to use cost at intermediate waypoints
					}
				}
			}

			// Resize shooting
			_shooting_wp.Resize(nMulti, nThread);
			// Mode
			_shooting_wp.SetMode(mode_t, mode_X);
			// Guess 
			_shooting_wp.InitShooting(vt, vX);				// init shooting with the guess
			// New desired final state
			for (int i = 0; i < Nwp+1; i++) {
				vX[i][0] = pathStateWP[i][0];
				vX[i][1] = pathStateWP[i][1];
				vX[i][2] = pathStateWP[i][2];
				vX[i][3] = 0.001*pathStateWP[i][3];
				vX[i][4] = 0.001*pathStateWP[i][4];
				vX[i][5] = 0.001*pathStateWP[i][5];
			}
			_shooting_wp.SetDesiredState(vt, vX);

		}

	}

	return info;
}

