/*
 * obstacle.cpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE
 */
#include <math.h>
#ifdef WIN32 /* Windows */
	#include <float.h>
	#define isnan _isnan 
#endif

#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <sstream>
#include <string>

#include "obstacle.hpp"

/**
* Map data
*/
struct obstacle::data_struct{
	int obstacleNbr;									///< obstacle number
	std::vector<std::vector<real>> obstaclePosition;	///< obstacle position
	std::vector<std::vector<real>> obstacleRadius;		///< obstacle radius
	int WPNbr;											///< waypoint number
	std::vector<std::vector<real>> WPState;				///< waypoint position
	std::vector<real> type;								///< type of obstacle (0 = circle, 1 = rectangle)
	int  modelChoice;									///< what model is used 
	parameters_struct parameters;						///< parameters can be used for continuation
	real penalizationRange;								///< range of the penalization around obstacles or waypoints
	std::string fileObstacles;							///< obstacle file
	std::string fileWP;									///< obstacle file
};

/**
* Constructor
*/
obstacle::obstacle(std::string the_fileObstacles, std::string the_fileWP){

	// Map data
	data = new data_struct;
	data->obstacleNbr = 0;
	data->WPNbr = 0;
		data->parameters.phiObs = 1;				// weight of obstacle function
		data->parameters.psiWP = 0.03;				// weight of waypoint function
		data->parameters.muObs = 1;					// diff between real radius and applied radius in grad_obs
		data->parameters.sigmaWP = 2.5;				// diff between real radius and applied radius in grad_obs
	data->penalizationRange = 10;
	data->fileObstacles = the_fileObstacles;		// obstacle file
	data->fileWP = the_fileWP;						// obstacle file

	ReadObstacleInput();
	ReadWPInput();

};

/**
* Destructor
*/
obstacle::~obstacle(){
	delete(data);
};

/**
* Read parameters in obstacle file
*/
void obstacle::ReadObstacleInput() {
	std::ifstream is(data->fileObstacles);

	if (is) {
		std::string line;
		// get obstacleNbr
		std::getline(is, line); // get "obstacleNbr:"
		std::getline(is, line); // get the integer obstacleNbr
		std::istringstream isn(line);
		if (!(isn >> data->obstacleNbr)) { return; } // error

		// get type
		data->type = std::vector<real>(data->obstacleNbr);
		std::getline(is, line); // get "type:"
		for(int i=0;i<data->obstacleNbr;i++){
			std::getline(is, line);
			std::istringstream ispos(line);
			if (!(ispos >> data->type[i])) { return; } // error
		}

		// get position
		data->obstaclePosition = std::vector<std::vector<real>>(data->obstacleNbr);
		std::getline(is, line); // get "position:"
		for(int i=0;i<data->obstacleNbr;i++){
			data->obstaclePosition[i] = std::vector<real>(3);
			std::getline(is, line);
			std::istringstream ispos(line);
			if (!(ispos >> data->obstaclePosition[i][0] >> data->obstaclePosition[i][1] >> data->obstaclePosition[i][2])) { return; } // error
		}

		// get radius
		data->obstacleRadius = std::vector<std::vector<real>>(data->obstacleNbr);
		std::getline(is, line); // get "radius:"
		for(int i=0;i<data->obstacleNbr;i++){
			data->obstacleRadius[i] = std::vector<real>(3);
			std::getline(is, line);
			std::istringstream israd(line);
			if (!(israd >> data->obstacleRadius[i][0] >> data->obstacleRadius[i][1] >> data->obstacleRadius[i][2])) { return; } // error
		}
	}
}

/**
* Read parameters in waypoint file
*/
void obstacle::ReadWPInput() {
	std::ifstream is(data->fileWP);

	if (is) {
		std::string line;
		// get n
		std::getline(is, line); // get "n_wp:"
		std::getline(is, line); // get the integer n
		std::istringstream isn(line);
		if (!(isn >> data->WPNbr)) { return; } // error

		// get positions
		data->WPState = std::vector<std::vector<real>>(data->WPNbr);
		std::getline(is, line); // get "position_wp:"
		for (int i = 0; i<data->WPNbr; i++) {
			data->WPState[i] = std::vector<real>(6, 0);
			std::getline(is, line);
			std::istringstream ispos(line);
			if (!(ispos >> data->WPState[i][0] >> data->WPState[i][1] >> data->WPState[i][2])) { return; } // error
		}

		// get directions
		for (int i = 0; i<data->WPNbr - 1; i++) {
			std::vector<real> deltaX(3);
			deltaX[0] = data->WPState[i + 1][0] - data->WPState[i][0];
			deltaX[1] = data->WPState[i + 1][1] - data->WPState[i][1];
			deltaX[2] = data->WPState[i + 1][2] - data->WPState[i][2];
			real distance = sqrt(deltaX[0] * deltaX[0] + deltaX[1] * deltaX[1] + deltaX[2] * deltaX[2]);
			data->WPState[i][3] = deltaX[0] / distance;
			data->WPState[i][4] = deltaX[1] / distance;
			data->WPState[i][5] = deltaX[2] / distance;
		}
		data->WPState[data->WPNbr - 1][3] = data->WPState[data->WPNbr - 2][3];
		data->WPState[data->WPNbr - 1][4] = data->WPState[data->WPNbr - 2][4];
		data->WPState[data->WPNbr - 1][5] = data->WPState[data->WPNbr - 2][5];
	}
}

/**
* Total function
*/
void obstacle::Function(std::vector<real> const& position, real & funcTot) const {
	real funcObs = 0;
	real funcWP = 0;

	ObstaclePenalizationFunction(position, funcObs);
	//WPPenalizationFunction(position, funcWP);

	funcTot = data->parameters.phiObs*funcObs + data->parameters.psiWP*funcWP;
};

/**
* Total gradient
*/
void obstacle::Gradient(std::vector<real> const& position, std::vector<real> & gradTot) const {
	std::vector<real> gradObs(3, 0);
	std::vector<real> gradWP(3, 0);

	ObstaclePenalizationGradient(position, gradObs);
	//WPPenalizationGradient(position, gradWP);

	gradTot[0] = data->parameters.phiObs*gradObs[0] + data->parameters.psiWP*gradWP[0];
	gradTot[1] = data->parameters.phiObs*gradObs[1] + data->parameters.psiWP*gradWP[1];
	gradTot[2] = data->parameters.phiObs*gradObs[2] + data->parameters.psiWP*gradWP[2];
};

/**
* Obstacle function
*/
void obstacle::ObstaclePenalizationFunction(std::vector<real> const& position, real & funcObs) const {
	funcObs = 0;

	real radx, rady, radz, rad;
	real x, y, z, d;
	real hx, hy, hz, h;

	for (int i = 0; i<data->obstacleNbr; i++) {

		x = data->obstaclePosition[i][0];
		y = data->obstaclePosition[i][1];
		z = data->obstaclePosition[i][2];
		radx = data->obstacleRadius[i][0];
		rady = data->obstacleRadius[i][1];
		radz = data->obstacleRadius[i][2];

		if (data->type[i] == 0) {
			// obstacles circulaires
			hx = (position[0] - x);
			hy = (position[1] - y);
			hz = (position[2] - z);
			d = sqrt(hx*hx + hy*hy + hz*hz);
			rad = d / sqrt(hx*hx / radx / radx + hy*hy / rady / rady + hz*hz / radz / radz);
			h = (d - rad) / data->parameters.muObs;
			funcObs = funcObs + (1 - tanh(h)) / 2;
		}
		else if (data->type[i] == 1) {
			// obstacles rectangulaires
			hx = (fabs(position[0] - x) - radx) / data->parameters.muObs;
			hy = (fabs(position[1] - y) - rady) / data->parameters.muObs;
			hz = (fabs(position[2] - z) - radz) / data->parameters.muObs;

			//if (hx*data->parameters.muObs <= data->penalizationRange &&
			//	hy*data->parameters.muObs <= data->penalizationRange &&
			//	hz*data->parameters.muObs <= data->penalizationRange)

			funcObs = funcObs + (1 - tanh(hx))*(1 - tanh(hy))*(1 - tanh(hz)) / 8;

			//h = fmin(hx, fmin(hy, hz));
			//funcObs = funcObs + (1 - tanh(h)) / 2;
		}
		else {
			funcObs = funcObs;
		}

	}

	if (isnan(funcObs)) funcObs = 0.0;
}

/**
* Obstacle gradient
*/
void obstacle::ObstaclePenalizationGradient(std::vector<real> const& position, std::vector<real> & gradObs) const {
	gradObs[0] = 0;
	gradObs[1] = 0;
	gradObs[2] = 0;

	real radx, rady, radz, rad;
	real x, y, z, d;
	real hx, hy, hz, h;

	for (int i = 0; i<data->obstacleNbr; i++) {

		x = data->obstaclePosition[i][0];
		y = data->obstaclePosition[i][1];
		z = data->obstaclePosition[i][2];
		radx = data->obstacleRadius[i][0];
		rady = data->obstacleRadius[i][1];
		radz = data->obstacleRadius[i][2];

		if (data->type[i] == 0) {
			// obstacles circulaires
			hx = (position[0] - x);
			hy = (position[1] - y);
			hz = (position[2] - z);
			d = sqrt(hx*hx + hy*hy + hz*hz);
			rad = d / sqrt(hx*hx / radx / radx + hy*hy / rady / rady);
			h = (d - rad) / data->parameters.muObs;
			real rho2 = (radx*radx - rady*rady) / (radx*radx*rady*rady) / sqrt(hx*hx / radx / radx + hy*hy / rady / rady) / sqrt(hx*hx / radx / radx + hy*hy / rady / rady) / sqrt(hx*hx / radx / radx + hy*hy / rady / rady);
			gradObs[0] = gradObs[0] - hx / d*(1 - hy*hy*rho2) / data->parameters.muObs*(1 - tanh(h)*tanh(h)) / 2;
			gradObs[1] = gradObs[1] - hy / d*(1 + hx*hx*rho2) / data->parameters.muObs*(1 - tanh(h)*tanh(h)) / 2;
			gradObs[2] = gradObs[2] - 0; //TODO
		}
		else if (data->type[i] == 1) {
			// obstacles rectangulaires
			hx = (fabs(position[0] - x) - radx) / data->parameters.muObs;
			hy = (fabs(position[1] - y) - rady) / data->parameters.muObs;
			hz = (fabs(position[2] - z) - radz) / data->parameters.muObs;

			//if (hx*data->parameters.muObs <= data->penalizationRange &&
			//	hy*data->parameters.muObs <= data->penalizationRange &&
			//	hz*data->parameters.muObs <= data->penalizationRange)
			{
				gradObs[0] = gradObs[0] - (position[0] - x) / fabs(position[0] - x) / data->parameters.muObs*(1 - tanh(hx)*tanh(hx))*(1 - tanh(hy))*(1 - tanh(hz)) / 8;
				gradObs[1] = gradObs[1] - (position[1] - y) / fabs(position[1] - y) / data->parameters.muObs*(1 - tanh(hy)*tanh(hy))*(1 - tanh(hx))*(1 - tanh(hz)) / 8;
				gradObs[2] = gradObs[2] - (position[2] - z) / fabs(position[2] - z) / data->parameters.muObs*(1 - tanh(hz)*tanh(hz))*(1 - tanh(hx))*(1 - tanh(hy)) / 8;
			}

			//real gradx, grady, gradz;
			//h = hx;
			//gradx = -(position[0] - x) / fabs(position[0] - x) / data->parameters.muObs*(1 - tanh(hx)*tanh(hx)) / 2;
			//grady = 0;
			//gradz = 0;
			//if (h > hy) {
			//	h = hy;
			//	gradx = 0;
			//	grady = -(position[1] - y) / fabs(position[1] - y) / data->parameters.muObs*(1 - tanh(hy)*tanh(hy)) / 2;
			//	gradz = 0;
			//}
			//if (h > hz) {
			//	h = hz;
			//	gradx = 0;
			//	grady = 0;
			//	gradz = -(position[2] - z) / fabs(position[2] - z) / data->parameters.muObs*(1 - tanh(hz)*tanh(hz)) / 2;
			//}
			//gradObs[0] = gradObs[0] + gradx;
			//gradObs[1] = gradObs[1] + grady;
			//gradObs[2] = gradObs[2] + gradz;


		}
		else {
			gradObs[0] = gradObs[0];
			gradObs[1] = gradObs[1];
			gradObs[2] = gradObs[2];
		}

	}

	if (isnan(gradObs[0])) gradObs[0] = 0.0;
	if (isnan(gradObs[1])) gradObs[1] = 0.0;
	if (isnan(gradObs[2])) gradObs[2] = 0.0;
}

/**
* WayPoints function
*/
void obstacle::WPPenalizationFunction(std::vector<real> const& position, real & funcWP) const {
	funcWP = 0;

	real x, y, z, d;
	real hx, hy, hz;

	for (int i = 0; i<data->WPNbr; i++) {

		x = data->WPState[i][0];
		y = data->WPState[i][1];
		z = data->WPState[i][2];
		hx = (position[0] - x) / data->parameters.sigmaWP;
		hy = (position[1] - y) / data->parameters.sigmaWP;
		hz = (position[2] - z) / data->parameters.sigmaWP;
		d = sqrt(hx*hx + hy*hy + hz*hz);

		if (d*data->parameters.sigmaWP <= data->penalizationRange)
			funcWP = funcWP - exp(- d * d / 2.0);
	}

	if (isnan(funcWP)) funcWP = 0.0;
}

/**
* WayPoints gradient
*/
void obstacle::WPPenalizationGradient(std::vector<real> const& position, std::vector<real> & gradWP) const {
	gradWP[0] = 0;
	gradWP[1] = 0;
	gradWP[2] = 0;

	real x, y, z, d;
	real hx, hy, hz;

	for (int i = 0; i<data->WPNbr; i++) {

		x = data->WPState[i][0];
		y = data->WPState[i][1];
		z = data->WPState[i][2];
		hx = (position[0] - x) / data->parameters.sigmaWP;
		hy = (position[1] - y) / data->parameters.sigmaWP;
		hz = (position[2] - z) / data->parameters.sigmaWP;
		d = sqrt(hx*hx + hy*hy + hz*hz);

		if (d*data->parameters.sigmaWP <= data->penalizationRange)
		{
			gradWP[0] = gradWP[0] + hx / data->parameters.sigmaWP * exp(-d * d / 2.0);
			gradWP[1] = gradWP[1] + hy / data->parameters.sigmaWP * exp(-d * d / 2.0);
			gradWP[2] = gradWP[2] + hz / data->parameters.sigmaWP * exp(-d * d / 2.0);
		}
	}

	if (isnan(gradWP[0])) gradWP[0] = 0.0;
	if (isnan(gradWP[1])) gradWP[1] = 0.0;
	if (isnan(gradWP[2])) gradWP[2] = 0.0;
}

/**
* Get parameters structure
*/
obstacle::parameters_struct & obstacle::GetParameterData(){
	return data->parameters;
}

/**
* Get path as a vector
*/
const std::vector<std::vector<real>> & obstacle::GetPath() {
	return data->WPState;
}
