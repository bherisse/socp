/*
 * vtolUAV.cpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */
#include <fstream>

#include <math.h>

#include "vtolUAV.hpp"

 ////////////////////////////////////////////////////////////////////////////////////////////
struct vtolUAV::data_struct{
	int n;									///< state dimension
	parameters_struct parameters;			///< those parameters can be used for continuation
	std::vector<real> switchingTimes;		///< switching times for intermediate points
	int  stepNbr;							///< step number for ModelInt
	std::string strFileTrace;				///< trace file
};

////////////////////////////////////////////////////////////////////////////////////////////
vtolUAV::vtolUAV(map & map, std::string the_fileTrace) : model(6), myMap(map){

	// vehicle data
	data = new data_struct;
	data->n = dim;									// state dimension
		data->parameters.u_max = 10;				// max normalized control
		data->parameters.a_max = 0.3;				// max acceleration
		data->parameters.alphaT = 0.05;				// weight for time cost
		data->parameters.alphaV = 0*0.05;			// weight for Vd
		data->parameters.invSigmaXwp = 0.5;			// weight for cost at intermediate points 
		data->parameters.Vd = 1;					// desired velocity
		data->parameters.ca = 0*0.05;//0.05			// drag coefficient
	data->switchingTimes = std::vector<real>();
	data->stepNbr = 10;//5;							// step number for ModelInt
	data->strFileTrace = the_fileTrace;				// trace file

	// trace file
	std::ofstream fileTrace;
  	fileTrace.open(data->strFileTrace.c_str(), std::ios::trunc);	// erase file
	fileTrace.close();
};

////////////////////////////////////////////////////////////////////////////////////////////
vtolUAV::~vtolUAV(){
	delete(data);
};

////////////////////////////////////////////////////////////////////////////////////////////
map & vtolUAV::GetMap() const {
	return myMap;
}

////////////////////////////////////////////////////////////////////////////////////////////
vtolUAV::mstate vtolUAV::Model(real const& t, mstate const& X) const{
	mstate Xdot(X.size());
	
	// current state value
	real x = 		X[0];
	real y = 		X[1];
	real z = 		X[2];
	real vx = 		X[3];
	real vy = 		X[4];
	real vz = 		X[5];
	real p_x = 		X[6];
    real p_y =		X[7];
	real p_z = 		X[8];
	real p_vx = 	X[9];
	real p_vy =		X[10];
	real p_vz = 	X[11];

	real normV = sqrt(vx*vx+vy*vy+vz*vz);

	// temporary variables
	real a_max = data->parameters.a_max;

    // control computation
	mcontrol u = Control(t, X);

	// get obstacle gradients
	std::vector<real> gradObsValue(3);
	std::vector<real> position(3);
	position[0] = x;
	position[1] = y;
	position[2] = z;
	myMap.Gradient(position, gradObsValue);

	// state and costate equations
    Xdot[0] = vx;
	Xdot[1] = vy;
	Xdot[2] = vz;
	Xdot[3] = a_max*u[0] - data->parameters.ca*vx*normV;
	Xdot[4] = a_max*u[1] - data->parameters.ca*vy*normV;
	Xdot[5] = a_max*u[2] - data->parameters.ca*vz*normV;
	Xdot[6] = 0 - gradObsValue[0];
	Xdot[7] = 0 - gradObsValue[1];
	Xdot[8] = 0 - gradObsValue[2];
	Xdot[9] = -p_x	 + data->parameters.ca*( p_vx*(normV+vx*vx/normV) + p_vy*(vy*vx/normV) + p_vz*(vz*vx/normV) ) - data->parameters.alphaV*vx/normV*(normV-data->parameters.Vd);
	Xdot[10] = -p_y	 + data->parameters.ca*( p_vy*(normV+vy*vy/normV) + p_vx*(vx*vy/normV) + p_vz*(vz*vy/normV) ) - data->parameters.alphaV*vy/normV*(normV-data->parameters.Vd);
	Xdot[11] = -p_z	 + data->parameters.ca*( p_vz*(normV+vz*vz/normV) + p_vx*(vx*vz/normV) + p_vy*(vy*vz/normV) ) - data->parameters.alphaV*vz/normV*(normV-data->parameters.Vd);

	return Xdot;

};

////////////////////////////////////////////////////////////////////////////////////////////
vtolUAV::mcontrol vtolUAV::Control(real const& t, mstate const& X) const{
	// current state value
	real x = 		X[0];
	real y = 		X[1];
	real z = 		X[2];
	real vx = 		X[3];
	real vy = 		X[4];
	real vz = 		X[5];
	real p_x = 		X[6];
    real p_y =		X[7];
	real p_z = 		X[8];
	real p_vx = 	X[9];
	real p_vy =		X[10];
	real p_vz = 	X[11];

	real norm_v = sqrt(vx*vx+vy*vy+vz*vz);

	// temporary variables
	real a_max = data->parameters.a_max;

    // control computation
	real u[3] = {-p_vx/a_max, -p_vy/a_max, -p_vz/a_max}; 

	real norm_u = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
	if (norm_u>data->parameters.u_max){
		u[0] = u[0]/norm_u*data->parameters.u_max;
		u[1] = u[1]/norm_u*data->parameters.u_max;
		u[2] = u[2]/norm_u*data->parameters.u_max;
		norm_u = data->parameters.u_max;
	}

	mcontrol control(3);
	control[0] = u[0];
	control[1] = u[1];
	control[2] = u[2];

	return control;

}

////////////////////////////////////////////////////////////////////////////////////////////
real vtolUAV::Hamiltonian(real const& t, mstate const& X) const{
	// current state value
	real x = 		X[0];
	real y = 		X[1];
	real z = 		X[2];
	real vx = 		X[3];
	real vy = 		X[4];
	real vz = 		X[5];
	real p_x = 		X[6];
    real p_y =		X[7];
	real p_z = 		X[8];
	real p_vx = 	X[9];
	real p_vy =		X[10];
	real p_vz = 	X[11];

	real normV = sqrt(vx*vx+vy*vy+vz*vz);

	// temporary variables
	real a_max = data->parameters.a_max;

	// control computation
	mcontrol u = Control(t, X);
	real norm_u = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);

	// get obstacle function
	real ObsValue = 0;
	std::vector<real> position(3);
	position[0] = x;
	position[1] = y;
	position[2] = z;
	myMap.Function(position, ObsValue);

	// H is computed
	real H = data->parameters.alphaT*1 
			+ data->parameters.alphaV/2*(normV-data->parameters.Vd)*(normV-data->parameters.Vd) 
			+ ObsValue
			+ a_max*a_max*norm_u*norm_u/2 
			+ p_x*vx + p_y*vy + p_z*vz 
			+ (p_vx*(a_max*u[0] - data->parameters.ca*vx*normV) + p_vy*(a_max*u[1] - data->parameters.ca*vy*normV) + p_vz*(a_max*u[2] - data->parameters.ca*vz*normV));

	return H;
}

////////////////////////////////////////////////////////////////////////////////////////////
vtolUAV::mstate vtolUAV::ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace){
	std::stringstream ss;
	
	real t = t0;								// time
	real dt = (tf - t0) / data->stepNbr;		// time step
	mstate Xs = X;

	if (isTrace) {
		integrate(modelStruct(this), Xs, t0, tf, dt, observerStruct(this, ss));
		// write in trace file
		std::ofstream fileTrace;
		fileTrace.open(data->strFileTrace.c_str(), std::ios::app);
		fileTrace << ss.str();
		fileTrace.close();
	}
	else {
		integrate(modelStruct(this), Xs, t0, tf, dt);
	}

	return Xs;
};

////////////////////////////////////////////////////////////////////////////////////////////
vtolUAV::parameters_struct & vtolUAV::GetParameterData(){
	return data->parameters;
}

////////////////////////////////////////////////////////////////////////////////////////////
void vtolUAV::FinalFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const {
	// Compute the function to solve
	for (int j = 0; j<dim; j++) {
		if (mode_X[j] == 1) {
			// continuity => X(k+n) = 0 (transversality condition)
			fvec[j] = X_tf[j + dim];
		}
		else {
			// continuity => X(k) = Xf(k)
			fvec[j] = X_tf[j] - Xf[j];
		}
	}
};

////////////////////////////////////////////////////////////////////////////////////////////
void vtolUAV::FinalHFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const {
	// Compute the function to solve
	for (int j = 0; j<dim; j++) {
		if (mode_X[j] == 1) {
			// continuity => X(k+n) = 0 (transversality condition)
			fvec[j] = X_tf[j + dim];
		}
		else {
			// continuity => X(k) = Xf(k)
			fvec[j] = X_tf[j] - Xf[j];
		}
	}
	fvec[dim] = Hamiltonian(tf, X_tf);
};

////////////////////////////////////////////////////////////////////////////////////////////
void vtolUAV::SwitchingTimesUpdate(std::vector<real> const& switchingTimes) {
	// Switching times
	data->switchingTimes.resize(switchingTimes.size());
	for (int i = 0; i<switchingTimes.size(); i++) data->switchingTimes[i] = switchingTimes[i];
}

////////////////////////////////////////////////////////////////////////////////////////////
void vtolUAV::SwitchingStateFunction(real const& t, int const& stateID, mstate const& X, mstate const& Xp, mstate const& Xd, mstate & fvec) const {
	// function to implement if mode_X = 1
	fvec[stateID] = (X[stateID] - Xp[stateID]);
	fvec[stateID + data->n] = (X[stateID + data->n] - Xp[stateID + data->n]);
	if (stateID < 3) {
		fvec[stateID + data->n] = (X[stateID + data->n] - Xp[stateID + data->n]) - data->parameters.invSigmaXwp*(X[stateID] - Xd[stateID]);
	}
};