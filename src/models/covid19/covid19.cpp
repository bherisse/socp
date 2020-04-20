/*
* covid19.cpp
*
*  Created on: April 19, 2020
*      Author: Bruno HERISSE (ONERA/DTIS)
*/
/*
* The model used here is a SEIR model with estimated covid19 parameters
*/

#include <fstream>
#include <sstream>
#include <math.h>

#include "covid19.hpp"

////////////////////////////////////////////////////////////////////////////////////////////
struct covid19::data_struct {
	int n;										///< state dimension
	parameters_struct parameters;				///< those parameters can be used for continuation
	int  stepNbr;								///< step number for ModelInt
	std::string strFileTrace;					///< trace file
};

////////////////////////////////////////////////////////////////////////////////////////////
covid19::covid19(std::string the_fileTrace) : model(4) {
	// vehicle data
	data = new data_struct;
	data->n = dim;									// state dimension
		data->parameters.R0 = 4;					// the number of secondary infections each infected individual produces
		data->parameters.Tinf = 10;					// duration patient is infectious
		data->parameters.Tinc = 5;				// incubation period
		data->parameters.N = 1;						// size of population (= 1 if normalized))
	data->stepNbr = 1000;							// step number for ModelInt
	data->strFileTrace = the_fileTrace;				// trace file

	// trace file
	std::ofstream fileTrace;
	fileTrace.open(data->strFileTrace.c_str(), std::ios::trunc);	// erase file
	fileTrace.close();
};

////////////////////////////////////////////////////////////////////////////////////////////
covid19::~covid19() {
	delete(data);
};

////////////////////////////////////////////////////////////////////////////////////////////
covid19::mstate covid19::Model(real const& t, mstate const& X) const {
	mstate Xdot(X.size());

	// current state and costate value
	real S = X[0];
	real E = X[1];
	real I = X[2];
	real R = X[3];
	real pS = X[4];
	real pE = X[5];
	real pI = X[6];
	real pR = X[7];

	real R0 = data->parameters.R0;
	real Tinf = data->parameters.Tinf;
	real Tinc = data->parameters.Tinc;
	real N = data->parameters.N;

	// control computation
	mcontrol u = Control(t, X);
	real Rt = R0 * (1 - u[0]);

	// state and costate equations from SEIR model
	Xdot[0] = -Rt / Tinf / N*S*I;
	Xdot[1] = Rt / Tinf / N*S*I - E / Tinc;
	Xdot[2] = E / Tinc - I / Tinf;
	Xdot[3] = I / Tinf;
	Xdot[4] = (pS - pE)*R*I / Tinf / N;
	Xdot[5] = (pE - pI) / Tinc;
	Xdot[6] = (pS - pE)*R*S / Tinf / N + (pI - pR) / Tinf;
	Xdot[7] = 0;

	return Xdot;

};

////////////////////////////////////////////////////////////////////////////////////////////
covid19::mcontrol covid19::Control(real const& t, mstate const& X) const {
	// current state and costate value
	real S = X[0];
	real E = X[1];
	real I = X[2];
	real R = X[3];
	real pS = X[4];
	real pE = X[5];
	real pI = X[6];
	real pR = X[7];

	real R0 = data->parameters.R0;
	real Tinf = data->parameters.Tinf;
	real Tinc = data->parameters.Tinc;
	real N = data->parameters.N;

	mcontrol control(1);

	// optimal control
	control[0] = (pE - pS)*S*I / Tinf / N * R0; 

	return control;

}

////////////////////////////////////////////////////////////////////////////////////////////
real covid19::Hamiltonian(real const& t, mstate const& X) const {
	// current state and costate value
	real S = X[0];
	real E = X[1];
	real I = X[2];
	real R = X[3];
	real pS = X[4];
	real pE = X[5];
	real pI = X[6];
	real pR = X[7];

	real R0 = data->parameters.R0;
	real Tinf = data->parameters.Tinf;
	real Tinc = data->parameters.Tinc;
	real N = data->parameters.N;

	// control computation
	mcontrol u = Control(t, X);
	real Rt = R0 * (1 - u[0]);

	// H is computed
	real H = u[0]*u[0] / 2
		+ pS * (-Rt / Tinf / N*S*I)
		+ pE * (Rt / Tinf / N*S*I - E / Tinc)
		+ pI * (E / Tinc - I / Tinf)
		+ pR * (I / Tinf);

	return H;
}

////////////////////////////////////////////////////////////////////////////////////////////
covid19::mstate covid19::ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace) {
	std::stringstream ss;

	real t = t0;							// time
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
covid19::parameters_struct & covid19::GetParameterData() {
	return data->parameters;
}