/*
* goddard.cpp
*
*  Created on: February 02, 2017
*      Author: Bruno HERISSE (ONERA/DTIS)
*/
/*
* The model used here is presented in the paper "Singular Arcs in the Generalized Goddard's Problem", F. Bonnans, P. Martinon, E. Trélat (J Optim Theory Appl (2008) 139: 439-461)
*/

#include <fstream>
#include <sstream>
#include <math.h>

#include "goddard.hpp"

////////////////////////////////////////////////////////////////////////////////////////////
struct goddard::data_struct {
	int n;										///< state dimension
	parameters_struct parameters;				///< those parameters can be used for continuation
	std::vector<real> switchingTimes;			///< switching times for singular control
	int  stepNbr;								///< step number for ModelInt
	std::string strFileTrace;					///< trace file
};

////////////////////////////////////////////////////////////////////////////////////////////
goddard::goddard(std::string the_fileTrace) : model(7) {
	// vehicle data
	data = new data_struct;
	data->n = dim;									// state dimension
		data->parameters.C = 3.5;					// coefficient for thrust
		data->parameters.b = 7;						// coefficient for mass flow rate
		data->parameters.KD = 310;					// coefficient for drag
		data->parameters.kr = 500;					// coefficient for density of air
		data->parameters.u_max = 1;					// max normalized control
		data->parameters.mu1 = 1;					// weight for the cost on the norm of the control in [0,1]
		data->parameters.mu2 = 0;					// weight for the quadratic cost on the control in [0,1]
		data->parameters.singularControl = -1;		// a constant approximation of the singular arc
	data->switchingTimes = std::vector<real>(2);
	data->switchingTimes[0] = 0.0227;				// first switch time for singular control
	data->switchingTimes[1] = 0.08;					// second switch time for singular control
	data->stepNbr = 10;								// step number for ModelInt
	data->strFileTrace = the_fileTrace;				// trace file

	// trace file
	std::ofstream fileTrace;
	fileTrace.open(data->strFileTrace.c_str(), std::ios::trunc);	// erase file
	fileTrace.close();
};

////////////////////////////////////////////////////////////////////////////////////////////
goddard::~goddard() {
	delete(data);
};

////////////////////////////////////////////////////////////////////////////////////////////
goddard::mstate goddard::Model(real const& t, mstate const& X) const {
	mstate Xdot(X.size());

	// current state and costate value
	real x = X[0];
	real y = X[1];
	real z = X[2];
	real vx = X[3];
	real vy = X[4];
	real vz = X[5];
	real mass = X[6];
	real p_x = X[7];
	real p_y = X[8];
	real p_z = X[9];
	real p_vx = X[10];
	real p_vy = X[11];
	real p_vz = X[12];
	real p_mass = X[13];

	real r = sqrt(x*x + y*y + z*z);
	real v = sqrt(vx*vx + vy*vy + vz*vz);
	real pvdotv = p_vx*vx + p_vy*vy + p_vz*vz;
	real b = data->parameters.b;
	real C = data->parameters.C;
	real KD = data->parameters.KD;
	real kr = data->parameters.kr;
	real g = 1 / r / r;											// normalized gravity (g=1 for r=1 that corresponds to Earth radius)
	real norm_pv = sqrt(p_vx*p_vx + p_vy*p_vy + p_vz*p_vz);

	// control computation
	mcontrol u = Control(t, X);
	real norm_u = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
	real pvdotu = p_vx*u[0] + p_vy*u[1] + p_vz*u[2];

	// state and costate equations
	Xdot[0] = vx;
	Xdot[1] = vy;
	Xdot[2] = vz;
	Xdot[3] = -KD*v*vx*exp(-kr*(r - 1)) / mass - g*x / r + C*u[0] / mass;
	Xdot[4] = -KD*v*vy*exp(-kr*(r - 1)) / mass - g*y / r + C*u[1] / mass;
	Xdot[5] = -KD*v*vz*exp(-kr*(r - 1)) / mass - g*z / r + C*u[2] / mass;
	Xdot[6] = -b*norm_u;
	Xdot[7] = -kr*KD / mass*v*exp(-kr*(r - 1))*x / r*pvdotv + g*(p_vx*(1 - 3 * x*x / r / r) / r - p_vy * 3 * x*y / r / r / r - p_vz * 3 * x*z / r / r / r);
	Xdot[8] = -kr*KD / mass*v*exp(-kr*(r - 1))*y / r*pvdotv + g*(-p_vx * 3 * y*x / r / r / r + p_vy*(1 - 3 * y*y / r / r) / r - p_vz * 3 * y*z / r / r / r);
	Xdot[9] = -kr*KD / mass*v*exp(-kr*(r - 1))*z / r*pvdotv + g*(-p_vx * 3 * z*x / r / r / r - p_vy * 3 * z*y / r / r / r + p_vz*(1 - 3 * z*z / r / r) / r);
	Xdot[10] = -p_x + KD / mass*exp(-kr*(r - 1))*(pvdotv*vx / v + p_vx*v);
	Xdot[11] = -p_y + KD / mass*exp(-kr*(r - 1))*(pvdotv*vy / v + p_vy*v);
	Xdot[12] = -p_z + KD / mass*exp(-kr*(r - 1))*(pvdotv*vz / v + p_vz*v);
	Xdot[13] = -KD*exp(-kr*(r - 1)) / mass / mass*v*pvdotv + C / mass / mass*pvdotu;

	return Xdot;

};

////////////////////////////////////////////////////////////////////////////////////////////
goddard::mcontrol goddard::Control(real const& t, mstate const& X) const {
	// current state value
	real x = X[0];
	real y = X[1];
	real z = X[2];
	real vx = X[3];
	real vy = X[4];
	real vz = X[5];
	real mass = X[6];
	real p_x = X[7];
	real p_y = X[8];
	real p_z = X[9];
	real p_vx = X[10];
	real p_vy = X[11];
	real p_vz = X[12];
	real p_mass = X[13];

	real r = sqrt(x*x + y*y + z*z);
	real v = sqrt(vx*vx + vy*vy + vz*vz);
	real pvdotv = p_vx*vx + p_vy*vy + p_vz*vz;
	real b = data->parameters.b;
	real C = data->parameters.C;
	real KD = data->parameters.KD;
	real kr = data->parameters.kr;
	real g = 1 / r / r;											// normalized gravity (g=1 for r=1 that corresponds to R_Earth)
	real norm_pv = sqrt(p_vx*p_vx + p_vy*p_vy + p_vz*p_vz);

	// control computation (minimisation of the Hamiltoninan)
	real Switch;
	real alpha_u = 0;
	real u[3];
	Switch = data->parameters.mu1 - b*p_mass - C / mass*norm_pv;		// switching function

	if (data->parameters.mu2 > 0) {
		// if a quadratic cost is used 
		if (Switch<0) {
			alpha_u = -Switch / 2 / data->parameters.mu2;
		}
		else {			//(Switch>=0)
			alpha_u = 0;
		}
	}
	else {
		// the strcture of the control is imposed : Bang-Singular-Off
		if (t <= (data->switchingTimes[0])) {
			alpha_u = 1.0;
		}
		else if (t>data->switchingTimes[0] && t <= data->switchingTimes[1]) {
			if (data->parameters.singularControl < 0) {
				alpha_u = GetSingularControl(t, X);		// true singular control 
			}
			else {
				alpha_u = data->parameters.singularControl;		// approximation of the singular control by a constant parameter
			}
		}
		else {
			alpha_u = 0;
		}
	}
	u[0] = -p_vx*alpha_u / norm_pv;
	u[1] = -p_vy*alpha_u / norm_pv;
	u[2] = -p_vz*alpha_u / norm_pv;

	real norm_u = fabs(alpha_u);

	// saturation
	if (norm_u>data->parameters.u_max) {
		u[0] = u[0] / norm_u*data->parameters.u_max;
		u[1] = u[1] / norm_u*data->parameters.u_max;
		u[2] = u[2] / norm_u*data->parameters.u_max;
		norm_u = data->parameters.u_max;
	}

	mcontrol control(3);
	control[0] = u[0];
	control[1] = u[1];
	control[2] = u[2];

	return control;

}

////////////////////////////////////////////////////////////////////////////////////////////
real goddard::GetSingularControl(real t, mstate const& X) const {
	// current state and costate value
	real x = X[0];
	real y = X[1];
	real z = X[2];
	real vx = X[3];
	real vy = X[4];
	real vz = X[5];
	real mass = X[6];
	real p_x = X[7];
	real p_y = X[8];
	real p_z = X[9];
	real p_vx = X[10];
	real p_vy = X[11];
	real p_vz = X[12];
	real p_mass = X[13];

	real r = sqrt(x*x + y*y + z*z);
	real v = sqrt(vx*vx + vy*vy + vz*vz);
	real rdotv = x*vx + y*vy + z*vz;
	real pvdotv = p_vx*vx + p_vy*vy + p_vz*vz;
	real b = data->parameters.b;
	real C = data->parameters.C;
	real KD = data->parameters.KD;
	real kr = data->parameters.kr;
	real g = 1 / r / r;										// normalized gravity (g=1 for r=1 that corresponds to Earth radius)
	real norm_pv = sqrt(p_vx*p_vx + p_vy*p_vy + p_vz*p_vz);
	real D = KD*exp(-kr*(r - 1));

	real p_xdot = -kr*KD / mass*v*exp(-kr*(r - 1))*x / r*pvdotv + g*(p_vx*(1 - 3 * x*x / r / r) / r - p_vy * 3 * x*y / r / r / r - p_vz * 3 * x*z / r / r / r);
	real p_ydot = -kr*KD / mass*v*exp(-kr*(r - 1))*y / r*pvdotv + g*(-p_vx * 3 * y*x / r / r / r + p_vy*(1 - 3 * y*y / r / r) / r - p_vz * 3 * y*z / r / r / r);
	real p_zdot = -kr*KD / mass*v*exp(-kr*(r - 1))*z / r*pvdotv + g*(-p_vx * 3 * z*x / r / r / r - p_vy * 3 * z*y / r / r / r + p_vz*(1 - 3 * z*z / r / r) / r);
	real p_vxdot = -p_x + KD / mass*exp(-kr*(r - 1))*(pvdotv*vx / v + p_vx*v);
	real p_vydot = -p_y + KD / mass*exp(-kr*(r - 1))*(pvdotv*vy / v + p_vy*v);
	real p_vzdot = -p_z + KD / mass*exp(-kr*(r - 1))*(pvdotv*vz / v + p_vz*v);

	real prdotdotv = p_xdot*vx + p_ydot*vy + p_zdot*vz;
	real prdotdotpv = p_xdot*p_vx + p_ydot*p_vy + p_zdot*p_vz;
	real prdotpvdot = p_x*p_vxdot + p_y*p_vydot + p_z*p_vzdot;
	real prdotpv = p_x*p_vx + p_y*p_vy + p_z*p_vz;
	real pvdotdotg = p_vxdot*g*x / r + p_vydot*g*y / r + p_vzdot*g*z / r;
	real pvdotdotv = p_vxdot*vx + p_vydot*vy + p_vzdot*vz;
	real pvdotdotpv = p_vxdot*p_vx + p_vydot*p_vy + p_vzdot*p_vz;
	real vdotg = vx*g*x / r + vy*g*y / r + vz*g*z / r;
	real pvdotg = p_vx*g*x / r + p_vy*g*y / r + p_vz*g*z / r;
	real prdotv = p_x*vx + p_y*vy + p_z*vz;
	real prdotg = p_x*g*x / r + p_y*g*y / r + p_z*g*z / r;

	// singular control computation (minimisation of the Hamiltoninan)
	real alpha_u = 0;

	real au = 2 * norm_pv*C / mass*pvdotv
		+ 2 * pvdotv*C / mass*norm_pv
		- b / mass*(2 * pvdotv*pvdotv + norm_pv*norm_pv*v*v)
		- b / D*v*prdotpv - C / D*prdotpv / v*pvdotv / norm_pv;

	real bu = -2 * norm_pv*norm_pv*(vdotg + D / mass*v*v*v) + 2 * v*v*pvdotdotpv
		- 2 * pvdotv*(pvdotg + D / mass*v*pvdotv - pvdotdotv)
		+ b / C*(2 * norm_pv*pvdotv*(vdotg + D / mass*v*v*v) + norm_pv*v*v*(pvdotg + D / mass*v*pvdotv - pvdotdotv) - v*v*pvdotv / norm_pv*pvdotdotpv)
		- mass / D*kr*rdotv / r*v*prdotpv + mass / D*prdotpv / v*(vdotg + D / mass*v*v*v) - mass / D*v*(prdotdotpv + prdotpvdot);


	alpha_u = bu / au;

	return alpha_u;
}

////////////////////////////////////////////////////////////////////////////////////////////
real goddard::Hamiltonian(real const& t, mstate const& X) const {
	// current state value
	real x = X[0];
	real y = X[1];
	real z = X[2];
	real vx = X[3];
	real vy = X[4];
	real vz = X[5];
	real mass = X[6];
	real p_x = X[7];
	real p_y = X[8];
	real p_z = X[9];
	real p_vx = X[10];
	real p_vy = X[11];
	real p_vz = X[12];
	real p_mass = X[13];

	real r = sqrt(x*x + y*y + z*z);
	real v = sqrt(vx*vx + vy*vy + vz*vz);
	real pvdotv = p_vx*vx + p_vy*vy + p_vz*vz;
	real b = data->parameters.b;
	real C = data->parameters.C;
	real KD = data->parameters.KD;
	real kr = data->parameters.kr;
	real g = 1 / r / r;								// normalized gravity (g=1 for r=1 that corresponds to Earth radius)

													// control computation
	mstate u = Control(t, X);
	real norm_u = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);

	// H is computed
	real H = data->parameters.mu1*norm_u + data->parameters.mu2*norm_u*norm_u
		+ p_x*vx + p_y*vy + p_z*vz
		+ p_vx*(-KD*v*vx*exp(-kr*(r - 1)) / mass - g*x / r + C*u[0] / mass)
		+ p_vy*(-KD*v*vy*exp(-kr*(r - 1)) / mass - g*y / r + C*u[1] / mass)
		+ p_vz*(-KD*v*vz*exp(-kr*(r - 1)) / mass - g*z / r + C*u[2] / mass)
		- p_mass*b*norm_u;

	return H;
}

////////////////////////////////////////////////////////////////////////////////////////////
goddard::mstate goddard::ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace) {
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
void goddard::Trace(real const& t, mstate const& X, std::stringstream & file) const {
	// control computation
	mstate control = Control(t, X);

	// H is computed
	real H = Hamiltonian(t, X);

	// write in trace file
	file << t << "\t";
	for (int k = 0; k<2 * dim; k++) {
		file << X[k] << "\t";
	}
	for (int k = 0; k<control.size(); k++) {
		file << control[k] << "\t";
	}
	file << H << "\t";

	// additional trace
	real Switch = data->parameters.mu1 - data->parameters.b*X[13] - data->parameters.C / X[6] * sqrt(X[10] * X[10] + X[11] * X[11] + X[12] * X[12]);
	file << Switch << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////
goddard::parameters_struct & goddard::GetParameterData() {
	return data->parameters;
}

////////////////////////////////////////////////////////////////////////////////////////////
void goddard::SwitchingTimesFunction(real const& t, mstate const& X, real & fvec) const {
	// current state value
	real x = X[0];
	real y = X[1];
	real z = X[2];
	real vx = X[3];
	real vy = X[4];
	real vz = X[5];
	real mass = X[6];
	real p_x = X[7];
	real p_y = X[8];
	real p_z = X[9];
	real p_vx = X[10];
	real p_vy = X[11];
	real p_vz = X[12];
	real p_mass = X[13];

	real b = data->parameters.b;
	real C = data->parameters.C;
	real norm_pv = sqrt(p_vx*p_vx + p_vy*p_vy + p_vz*p_vz);

	// Switching function
	//real Switch = data->parameters.mu1 - b*p_mass - C/mass*norm_pv;

	// Using Switch or Hamiltonian is equivalent in this case (free final time)
	//fvec = Switch;
	fvec = Hamiltonian(t, X);
}

////////////////////////////////////////////////////////////////////////////////////////////
void goddard::SwitchingTimesUpdate(std::vector<real> const& switchingTimes) {
	// Switching times
	data->switchingTimes.resize(switchingTimes.size());
	for (int i = 0; i<switchingTimes.size(); i++) data->switchingTimes[i] = switchingTimes[i];
}

////////////////////////////////////////////////////////////////////////////////////////////
void goddard::SetParameterDataName(std::string name, real value) {
	if (name == std::string("KD"))
		data->parameters.KD = value;
	if (name == std::string("C"))
		data->parameters.C = value;
	if (name == std::string("b"))
		data->parameters.b = value;
	if (name == std::string("kr"))
		data->parameters.kr = value;
	if (name == std::string("mu2"))
		data->parameters.mu2 = value;
	if (name == std::string("mu1"))
		data->parameters.mu1 = value;
	if (name == std::string("singularControl"))
		data->parameters.singularControl = value;
	if (name == std::string("u_max"))
		data->parameters.u_max = value;
}

////////////////////////////////////////////////////////////////////////////////////////////
real & goddard::GetParameterDataName(std::string name) {
	if (name == std::string("KD"))
		return data->parameters.KD;
	if (name == std::string("C"))
		return data->parameters.C;
	if (name == std::string("b"))
		return data->parameters.b;
	if (name == std::string("kr"))
		return data->parameters.kr;
	if (name == std::string("mu2"))
		return data->parameters.mu2;
	if (name == std::string("mu1"))
		return data->parameters.mu1;
	if (name == std::string("singularControl"))
		return data->parameters.singularControl;
	if (name == std::string("u_max"))
		return data->parameters.u_max;

	real res = 0.0;
	return res;
}