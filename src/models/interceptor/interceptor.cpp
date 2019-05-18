/*
 * interceptor.cpp
 *
 *  Created on: June 30, 2016
 *      Author: Riccardo BONALLI & Bruno HERISSE (ONERA/DTIS)
 */
/*
 * The algorithm is presented in the paper "Analytical Initialization of a Continuation-Based Indirect
Method for Optimal Control of Endo-Atmospheric Launch Vehicle Systems", R. Bonalli, B. Hérissé, E. Trélat (IFAC WC 2017)
 */

#include <fstream>
#include <math.h>

#include "Eigen/Dense"

#include "interceptor.hpp"

 ////////////////////////////////////////////////////////////////////////////////////////////
struct interceptor::data_struct{
	int n;							///< state dimension
	parameters_struct parameters;	///< those parameters can be used for continuation
	real R_Earth;					///< Earth radius (m)
	real mu0; 						///< gravity constant = 3.986e14
	int  stepNbr;					///< step number for ModelInt
	std::string strFileTrace;		///< trace file
	int currentChart;				///< current chart used for the state (1 = default NED frame)
	real chartLimit;				///< parameter for charts
	int stageMode;					///< current mode: mode=1 (propulsion), mode=0 (no propulsion)
};

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::interceptor(std::string the_fileTrace) : model(6){
	// vehicle data
	data = new data_struct;
	data->n = dim;								// state dimension
		data->parameters.c0 = 0.00075;			// max curvature at ground level (1/m)
		data->parameters.hr = 7500;				// reference altitude (m)
		data->parameters.d0 = 0.00005;			// drag at ground level (1/m)	
		data->parameters.eta = 0.442;			// coeff of efficiency
		data->parameters.propellant_mass = 200;	// propellant mass (kg)
		data->parameters.empty_mass = 200; 		// empty mass (kg)
		data->parameters.q = 10; 				// mass flow rate (kg/s)
		data->parameters.ve = 1500; 			// gas speed (m/s)
		data->parameters.alpha_max = M_PI/6; 	// max angle of attack (rd)
		data->parameters.u_max = 1;				// max of the normalized control (for saturation)
		data->parameters.a_max = 1500;			// max acceleration allowed (for saturation)
		data->parameters.mu_gft = 1;			// parameter for considering gravity and propulsion
		data->parameters.muT = 0;				// weight for time cost
		data->parameters.muV = 1;				// weight for velocity cost
		data->parameters.muC = 0;				// weight for quadratic control cost
	data->R_Earth = 6378145;					// Earth radius (m)
	data->mu0 = 3.986e14; 						// gravity constant = 3.986e14
	data->stepNbr = 50;							// step number for ModelInt
	data->strFileTrace = the_fileTrace;			// trace file
	data->currentChart = 1;						// current chart used for the state (1 = default NED frame)
	data->chartLimit = 0.1;						// parameter for charts
	data->stageMode = 0;						// current mode: mode=1 (propulsion), mode=0 (no propulsion)

	// trace file
	std::ofstream fileTrace;
	fileTrace.open(data->strFileTrace.c_str(), std::ios::trunc);	// erase file
	fileTrace.close();
};

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::~interceptor(){
	delete(data);
};

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::mstate interceptor::Model(real const& t, mstate const& X) const{
	if (data->currentChart==1){
		return Model_1(t, X);
	}else{
		return Model_2(t, X);
	}
};

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::mstate interceptor::static_Model(real const& t, interceptor::mstate const& X, void *model){ // static model to be passed as an argument of the ODE solver
	return static_cast<interceptor*>(model)->Model(t, X);
} 

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::mcontrol interceptor::Control(real const& t, mstate const& X) const{
	if (data->currentChart==1){
		return Control_1(t, X);
	}else{
		return Control_2(t, X);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
real interceptor::Hamiltonian(real const& t, mstate const& X) const{
	if (data->currentChart==1){
		return Hamiltonian_1(t, X);
	}else{
		return Hamiltonian_2(t, X);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::mstate interceptor::ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace){
	std::stringstream ss;
	
	real t = t0;							// time
	real dt = (tf-t0)/data->stepNbr;		// time step
	mstate Xs = X;

	if (isTrace) Trace(t, Xs, ss);
  	for (int i = 0; i < data->stepNbr; i++) {
		// determine which chart to be used
		Xs = SetChart(t, Xs);

		// Solve ODE
		Xs = odeTools::RK4(t, Xs, dt, &static_Model, (void*) this);
		t += dt;

		if (isTrace) Trace(t, Xs, ss);
	}
	if (isTrace) {
		std::ofstream fileTrace;
		fileTrace.open(data->strFileTrace.c_str(), std::ios::app);
		fileTrace << ss.str();
		fileTrace.close();
	}

	return Xs;
};

////////////////////////////////////////////////////////////////////////////////////////////
void interceptor::Trace(real const& t, mstate const& X, std::stringstream & file) const{
	// control computation
	mstate control = Control(t,X);

	// H is computed
	real H = Hamiltonian(t, X);

	// write in trace file
	file << t << "\t";
	for (int k=0;k<2*dim; k++){
		file << X[k] << "\t";
	}
	for (int k=0;k<control.size(); k++){
		file << control[k] << "\t";
	}
	file << H << "\t";

	// additional trace
	file << data->currentChart << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////
int interceptor::GetMode(real const& t, mstate const& X) const{
	return data->stageMode;
};

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::parameters_struct & interceptor::GetParameterData(){
	return data->parameters;
}

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::mstate interceptor::ComputeTraj(real const& t0, mstate const& X0, real const& tf, int isTrace){
	// initialize chart (X0 is in chart 1 by default)
	data->currentChart=1;

	// time
	real t1 = data->parameters.propellant_mass/data->parameters.q;

	// check number of stages
	int n_stages = 1;
	if (t0 < t1) n_stages = 2;

	mstate X_t0(X0.size()), X1(X0.size()), X_tf(X0.size());
	X_t0 = X0;
	switch (n_stages) {
		case 2 :
			//************************************************************//
			//************ Stage 1 : powered vehicle *********************//
			//************************************************************//
			// integrate from t0 to t1
			//mode = 1;				// mode=1 (propulsion), mode=0 (no propulsion)
			data->stageMode = 1;
			if (tf>t1){
				X1 = ModelInt(t0,X_t0,t1,isTrace);
			}else{
				X_tf = ModelInt(t0,X_t0,tf,isTrace);
				break;
			}

			//****************************************************************//
			//************ Stage 2 : non powered vehicle *********************//
			//****************************************************************//
			// integrate from t1 to tf
			//mode = 0;
			data->stageMode = 0;
			X_tf = ModelInt(t1,X1,tf,isTrace);
			break;
		case 1 :
			//****************************************************************//
			//************ Stage 2 : non powered vehicle *********************//
			//****************************************************************//
    		// integrate from t0 to tf
			//mode = 0;
			data->stageMode = 0;
			X_tf = ModelInt(t0,X_t0,tf,isTrace);
			break;

		//default : /* Optional */
	}

	return X_tf;

};

////////////////////////////////////////////////////////////////////////////////////////////
void interceptor::FinalFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const{
	mstate Xc = X_tf;
	mstate Xfc = Xf;

	// chart handling
	if (data->currentChart==2) 
		Xc = ConversionState21(X_tf);

	// Compute the function to solve
	for (int j=0;j<data->n;j++){
		if (mode_X[j]==1){
			// X(k+n) (transversality condition)
			fvec[j] = Xc[j+data->n];
			if (j==1){
				fvec[j] = Xc[j+data->n] + data->parameters.muV;	// optimality of the final velocity
			}
		}else{
			// continuity => X(k) = Xf(k)
			fvec[j] = Xc[j] - Xfc[j];
			if (j == 0) {
				fvec[j] = fvec[j] / data->parameters.hr; // scaling
			}
			if (j == 3 && fabs(cos(Xfc[2]))<1e-5) {
				fvec[j] = Xc[j+data->n]; // if cos(gamma)=0, chi is free, then p_chi = 0 
			}
		}
	}

};

////////////////////////////////////////////////////////////////////////////////////////////
void interceptor::FinalHFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const{
	mstate Xc = X_tf;
	mstate Xfc = Xf;

	// chart handling
	if (data->currentChart == 2) 
		Xc = ConversionState21(X_tf);

	// Compute the function to solve
	for (int j = 0; j<data->n; j++) {
		if (mode_X[j] == 1) {
			// X(k+n) (transversality condition)
			fvec[j] = Xc[j + data->n];
			if (j == 1) {
				fvec[j] = Xc[j + data->n] + data->parameters.muV;	// optimality of the final velocity
			}
		}
		else {
			// continuity => X(k) = Xf(k)
			fvec[j] = Xc[j] - Xfc[j];
			if (j == 0) {
				fvec[j] = fvec[j] / data->parameters.hr; // scaling
			}
			if (j == 3 && fabs(cos(Xfc[2]))<1e-5) {
				fvec[j] = Xc[j + data->n]; // if cos(gamma)=0, chi is free, then p_chi = 0 
			}
		}
	}

	fvec[data->n] = Hamiltonian(tf,X_tf) + data->parameters.muT; // X here since data->currentChart has not changed !

};

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::mstate interceptor::Model_1(real const& t, mstate const& X) const{
	mstate Xdot(X.size());

   	// current state value (state + costate)
	real h = 		X[0];		//h
	real v = 		X[1];		//v
	real gamma = 	X[2];		//gamma
	real chi = 		X[3];		//chi
	real L = 		X[4];		//L
	real l = 		X[5];		//l
	real p_h = 		X[6];		//p_h
    real p_v =		X[7];		//p_v
	real p_gamma =	X[8];		//p_gamma
	real p_chi = 	X[9];		//p_chi
	real p_L = 		X[10];		//p_L
	real p_l = 		X[11];		//p_l

	// if mode = 1 : propulsion, if mode = 0 : no propulsion
	int mode = GetMode(t, X);

	// temporary variables
	real qm = mode*data->parameters.q*data->parameters.mu_gft;			// mass flow rate
	real mass = ComputeMass(t, X);										//mass
	real c_max = data->parameters.c0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;	// max curvature at altitude h
	real d = data->parameters.d0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;		// max drag at altitude h
	real r = h + data->R_Earth;											// distance to the Earth center
	real g = data->mu0/r/r*data->parameters.mu_gft;						// gravity
    real ft = data->parameters.ve*qm;									// propulsion
    real eta = data->parameters.eta;
    real hr = data->parameters.hr;
    real alpha_max = data->parameters.alpha_max;
    real R_Earth = data->R_Earth;
    real u_max = data->parameters.u_max;
    real a_max = data->parameters.a_max;

    // control computation (linear approximation, no singular control considered here)
	mcontrol control = Control_1(t, X);
	real u = control[0];
	real beta = control[1];
	real alpha = alpha_max*u;	// angle of attack

	// state and costate equations
    Xdot[0]		=	v*sin(gamma);
	Xdot[1]		=	-(d+eta*c_max*u*u)*v*v - g*sin(gamma) + ft*cos(alpha)/mass;
	Xdot[2]		=	v*c_max*u*cos(beta) - g/v*cos(gamma) + ft*sin(alpha)*cos(beta)/mass/v + v*cos(gamma)/r;
	Xdot[3]		=	v*c_max*u*sin(beta)/cos(gamma) + ft*sin(alpha)*sin(beta)/cos(gamma)/mass/v + v*cos(gamma)*tan(L)*sin(chi)/r;
	Xdot[4]		=	v*cos(gamma)*cos(chi)/r;
	Xdot[5]		=	v*cos(gamma)*sin(chi)/cos(L)/r;
	Xdot[6]		=	-p_v/hr*(d+eta*c_max*u*u)*v*v - 2*g/r*(p_gamma/v*cos(gamma) + p_v*sin(gamma))
				+ p_L*v*cos(gamma)*cos(chi)/r/r 			+ p_gamma*v*cos(gamma)/r/r + p_gamma*v*c_max*u*cos(beta)/hr
				+ p_l*v*cos(gamma)*sin(chi)/cos(L)/r/r 	+ p_chi*v*cos(gamma)*tan(L)*sin(chi)/r/r + p_chi*v*c_max*u*sin(beta)/cos(gamma)/hr;
	Xdot[7]		=	- ( p_L*cos(gamma)*cos(chi)/r + p_l*cos(gamma)*sin(chi)/cos(L)/r + p_h*sin(gamma)
				+ p_gamma*(c_max*u*cos(beta) + g/v/v*cos(gamma) - ft*sin(alpha)*cos(beta)/mass/v/v + cos(gamma)/r)
				+ p_chi*(c_max*u*sin(beta)/cos(gamma) - ft*sin(alpha)*sin(beta)/cos(gamma)/mass/v/v + cos(gamma)*tan(L)*sin(chi)/r)
				- p_v*2*(d+eta*c_max*u*u)*v );
	Xdot[8]		=	v * ( p_L*sin(gamma)*cos(chi)/r + p_l*sin(gamma)*sin(chi)/cos(L)/r - p_h*cos(gamma) )
				- g * ( p_gamma/v*sin(gamma) - p_v*cos(gamma) )
				+ p_gamma*v*sin(gamma)/r + p_chi*v*sin(gamma)*tan(L)*sin(chi)/r
				- p_chi*(v*c_max*u*sin(beta) + ft*sin(alpha)*sin(beta)/mass/v)*sin(gamma)/cos(gamma)/cos(gamma);
	Xdot[9]	=	v * ( p_L*cos(gamma)*sin(chi)/r - p_l*cos(gamma)*cos(chi)/cos(L)/r - p_chi*cos(gamma)*tan(L)*cos(chi)/r);
	Xdot[10]	=	- p_l*v*cos(gamma)*sin(chi)*sin(L)/cos(L)/cos(L)/r - p_chi*v*cos(gamma)*(1+tan(L)*tan(L))*sin(chi)/r;
	Xdot[11]	=	0.0;

    return Xdot;
};

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::mcontrol interceptor::Control_1(real const& t, mstate const& X) const{
	// current state value (state + costate)
	real h =		X[0];		//h
	real v =		X[1];		//v
	real gamma =	X[2];		//gamma
	real chi =		X[3];		//chi
	real L =		X[4];		//L
	real l =		X[5];		//l
	real p_h =		X[6];		//p_h
	real p_v =		X[7];		//p_v
	real p_gamma =	X[8];		//p_gamma
	real p_chi =	X[9];		//p_chi
	real p_L =		X[10];		//p_L
	real p_l =		X[11];		//p_l

	// if mode = 1 : propulsion, if mode = 0 : no propulsion
	int mode = GetMode(t, X);

	// temporary variables
	real qm = mode*data->parameters.q*data->parameters.mu_gft;			// mass flow rate
	real mass = ComputeMass(t, X);										//mass
	real c_max = data->parameters.c0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;	// max curvature at alatitude h
	real d = data->parameters.d0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;		// max drag at altitude h
	real r = h + data->R_Earth;											// distance to the Earth center
	real g = data->mu0/r/r*data->parameters.mu_gft;						// gravity
	real ft = data->parameters.ve*qm;									// propulsion
    real eta = data->parameters.eta;
    real hr = data->parameters.hr;
    real alpha_max = data->parameters.alpha_max;
    real R_Earth = data->R_Earth;
    real u_max = data->parameters.u_max;
    real a_max = data->parameters.a_max;

    // control computation (linear approximation, no singular control considered here)
	real beta = atan2(p_chi,p_gamma*cos(gamma)); //atan(p_chi/p_gamma/cos(gamma));
	real u = ( 	p_gamma*(v*c_max*cos(beta)+ft*cos(beta)*alpha_max/mass/v)
					+ p_chi*(v*c_max*sin(beta)/cos(gamma)+ft*sin(beta)/cos(gamma)*alpha_max/mass/v)
					)/( p_v*(2*eta*c_max*v*v+ft*alpha_max*alpha_max/mass) - data->parameters.muC);
	if (fabs(u)>u_max){
		u = u_max * u/fabs(u);
	}

	mcontrol control(2);
	control[0] = u;
	control[1] = beta; 

	return control;
}

////////////////////////////////////////////////////////////////////////////////////////////
real interceptor::Hamiltonian_1(real const& t, mstate const& X) const{
	// current state value (state + costate)
	real h =		X[0];		//h
	real v =		X[1];		//v
	real gamma =	X[2];		//gamma
	real chi =		X[3];		//chi
	real L =		X[4];		//L
	real l =		X[5];		//l
	real p_h =		X[6];		//p_h
	real p_v =		X[7];		//p_v
	real p_gamma =	X[8];		//p_gamma
	real p_chi =	X[9];		//p_chi
	real p_L =		X[10];		//p_L
	real p_l =		X[11];		//p_l

	// if mode = 1 : propulsion, if mode = 0 : no propulsion
	int mode = GetMode(t, X);

	// temporary variables
	real qm = mode*data->parameters.q*data->parameters.mu_gft;			// mass flow rate
	real mass = ComputeMass(t, X);										//mass
	real c_max = data->parameters.c0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;	// max curvature at alatitude h
	real d = data->parameters.d0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;		// max drag at altitude h
	real r = h + data->R_Earth;											// distance to the Earth center
	real g = data->mu0/r/r*data->parameters.mu_gft;						// gravity
    real ft = data->parameters.ve*qm;									// propulsion
    real eta = data->parameters.eta;
    real hr = data->parameters.hr;
    real alpha_max = data->parameters.alpha_max;
    real R_Earth = data->R_Earth;
    real u_max = data->parameters.u_max;
    real a_max = data->parameters.a_max;

    // control computation (linear approximation, no singular control considered here)
	mcontrol control = Control_1(t, X);
	real u = control[0];
	real beta = control[1];
	real alpha = alpha_max*u;	// angle of attack

	// H is computed
	real H = 	p_L*v*cos(gamma)*cos(chi)/r
				+p_l*v*cos(gamma)*sin(chi)/cos(L)/r
				+p_h*v*sin(gamma)
				+p_gamma*(v*c_max*u*cos(beta) - g/v*cos(gamma) + ft*sin(alpha)*cos(beta)/mass/v + v*cos(gamma)/r)
				+p_chi*(v*c_max*u*sin(beta)/cos(gamma) + ft*sin(alpha)*sin(beta)/cos(gamma)/mass/v + v*cos(gamma)*tan(L)*sin(chi)/r)
				-p_v*((d+eta*c_max*u*u)*v*v + g*sin(gamma) - ft*cos(alpha)/mass)
				+ data->parameters.muC*u*u/2;

	return H;

}

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::mstate interceptor::Model_2(real const& t, mstate const& X) const{
	mstate Xdot(X.size());

   	// current state value (state + costate)
	real h = 		X[0];		//h
	real v = 		X[1];		//v
	real theta = 	X[2];		//theta
	real phi = 		X[3];		//phi
	real L = 		X[4];		//L
	real l = 		X[5];		//l
	real p_h = 		X[6];		//p_h
    real p_v =		X[7];		//p_v
	real p_theta =	X[8];		//p_theta
	real p_phi = 	X[9];		//p_phi
	real p_L = 		X[10];		//p_L
	real p_l = 		X[11];		//p_l

	// if mode = 1 : propulsion, if mode = 0 : no propulsion
	int mode = GetMode(t, X);

	// temporary variables
	real qm = mode*data->parameters.q*data->parameters.mu_gft;			// mass flow rate
	real mass = ComputeMass(t, X);										//mass
	real c_max = data->parameters.c0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;	// max curvature at alatitude h
	real d = data->parameters.d0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;		// max drag at altitude h
	real r = h + data->R_Earth;											// distance to the Earth center
	real g = data->mu0/r/r*data->parameters.mu_gft;						// gravity
    real ft = data->parameters.ve*qm;									// propulsion
    real eta = data->parameters.eta;
    real hr = data->parameters.hr;
    real alpha_max = data->parameters.alpha_max;
    real R_Earth = data->R_Earth;
    real u_max = data->parameters.u_max;
    real a_max = data->parameters.a_max;

    // control computation (linear approximation, no singular control considered here)
	mcontrol control = Control_2(t, X);
	real u = control[0];
	real beta = control[1];
	real alpha = alpha_max*u;	// angle of attack

	// state and costate equations
    Xdot[0]		=	-v*cos(theta)*cos(phi);
	Xdot[1]		=	-(d+eta*c_max*u*u)*v*v + g*cos(theta)*cos(phi) + ft*cos(alpha)/mass;
	Xdot[2]		=	v*c_max*u*cos(beta) + v*sin(theta)*(cos(phi) + sin(phi)*tan(L))/r
					+ (ft*sin(alpha)*cos(beta)/(mass*v) - g*sin(theta)*cos(phi)/v);
	Xdot[3]		=	-v*c_max*u*sin(beta)/cos(theta)
					+ v*cos(theta)*(sin(phi) + tan(theta)*tan(theta)*(sin(phi) - tan(L)*cos(phi)))/r
					- (ft*sin(alpha)*sin(beta)/(mass*v*cos(theta)) + g*sin(phi)/(v*cos(theta)));
	Xdot[4]		=	v*cos(theta)*sin(phi)/r;
	Xdot[5]		=	v*sin(theta)/(r*cos(L));
	Xdot[6]		=	-p_v/hr*(d+eta*c_max*u*u)*v*v - 2*g/r*(p_theta*sin(theta)*cos(phi)/v + p_phi*sin(phi)/cos(theta)/v - p_v*cos(theta)*cos(phi))
				+ p_L*v*cos(theta)*sin(phi)/r/r + v*p_theta*sin(theta)*(cos(phi) + sin(phi)*tan(L))/r/r + p_theta*v*c_max*u*cos(beta)/hr
				+ p_l*v*sin(theta)/cos(L)/r/r 	+ v*p_phi*cos(theta)*(sin(phi) + tan(theta)*tan(theta)*(sin(phi) - tan(L)*cos(phi)))/r/r - p_phi*v*c_max*u*sin(beta)/cos(theta)/hr;
	Xdot[7]		=	- ( p_L*cos(theta)*sin(phi)/r + p_l*sin(theta)/(r*cos(L)) - p_h*cos(theta)*cos(phi)
				+ p_theta*(c_max*u*cos(beta) + g/v/v*sin(theta)*cos(phi) - ft*sin(alpha)*cos(beta)/mass/v/v + sin(theta)*(cos(phi) + sin(phi)*tan(L))/r)
				+ p_phi*(-c_max*u*sin(beta)/cos(theta) + g/v/v*sin(phi)/cos(theta) + ft*sin(alpha)*sin(beta)/cos(theta)/mass/v/v + cos(theta)*(sin(phi) + tan(theta)*tan(theta)*(sin(phi) - tan(L)*cos(phi)))/r)
				- p_v*2*(d+eta*c_max*u*u)*v );
	Xdot[8]		=	- v * (- p_L*sin(theta)*sin(phi)/r + p_l*cos(theta)/(r*cos(L)) + p_h*sin(theta)*cos(phi) )
				- g * ( - p_theta*cos(theta)*cos(phi)/v - p_phi*sin(phi)*tan(theta)/(v*cos(theta)) - p_v*sin(theta)*cos(phi) )
				- p_theta*v*cos(theta)*(cos(phi) + sin(phi)*tan(L))/r + p_phi*v*sin(theta)*(sin(phi) + tan(theta)*tan(theta)*(sin(phi) - tan(L)*cos(phi)))/r
				- p_phi*v*cos(theta)*(2*tan(theta)*(1+tan(theta)*tan(theta))*(sin(phi) - tan(L)*cos(phi)))/r
				- p_phi*(- v*c_max*u*sin(beta) - ft*sin(alpha)*sin(beta)/mass/v)*tan(theta)/cos(theta);
	Xdot[9]	=	-v * ( p_h*cos(theta)*sin(phi) + p_L*cos(theta)*cos(phi)/r)
					-g * ( p_theta*sin(theta)*sin(phi)/v - p_phi*cos(phi)/(v*cos(theta)) - p_v*cos(theta)*sin(phi) )
					-p_theta* (v*sin(theta)*(-sin(phi) + cos(phi)*tan(L))/r)
					-p_phi* v*cos(theta)*(cos(phi) + tan(theta)*tan(theta)*(cos(phi) + tan(L)*sin(phi)))/r;
	Xdot[10]	=	- p_l*v*sin(theta)*tan(L)/cos(L)/r - v*(1 + tan(L)*tan(L))*(p_theta*sin(theta)*sin(phi) - p_phi*cos(theta)*cos(phi)*tan(theta)*tan(theta))/r;
	Xdot[11]	=	0.0;

    return Xdot;
};

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::mcontrol interceptor::Control_2(real const& t, mstate const& X) const{
	// current state value (state + costate)
	real h =		X[0];		//h
	real v =		X[1];		//v
	real theta =	X[2];		//theta
	real phi =		X[3];		//phi
	real L =		X[4];		//L
	real l =		X[5];		//l
	real p_h =		X[6];		//p_h
	real p_v =		X[7];		//p_v
	real p_theta =	X[8];		//p_theta
	real p_phi =	X[9];		//p_phi
	real p_L =		X[10];		//p_L
	real p_l =		X[11];		//p_l

	// if mode = 1 : propulsion, if mode = 0 : no propulsion
	int mode = GetMode(t, X);

	// temporary variables
	real qm = mode*data->parameters.q*data->parameters.mu_gft;			// mass flow rate
	real mass = ComputeMass(t, X);										//mass
	real c_max = data->parameters.c0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;	// max curvature at alatitude h
	real d = data->parameters.d0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;		// max drag at altitude h
	real r = h + data->R_Earth;											// distance to the Earth center
	real g = data->mu0/r/r*data->parameters.mu_gft;						// gravity
    real ft = data->parameters.ve*qm;									// propulsion
    real eta = data->parameters.eta;
    real hr = data->parameters.hr;
    real alpha_max = data->parameters.alpha_max;
    real R_Earth = data->R_Earth;
    real u_max = data->parameters.u_max;
    real a_max = data->parameters.a_max;

    // control computation (linear approximation, no singular control considered here)
	real beta = atan2(-p_phi,p_theta*cos(theta)); 
	real u = ( 	p_theta*(v*c_max*cos(beta)+ft*cos(beta)*alpha_max/mass/v)
					- p_phi*(v*c_max*sin(beta)/cos(theta)+ft*sin(beta)/cos(theta)*alpha_max/mass/v)
					)/( p_v*(2*eta*c_max*v*v+ft*alpha_max*alpha_max/mass) - data->parameters.muC);
	if (fabs(u)>u_max){
		u = u_max * u/fabs(u);
	}

	mcontrol control(2);
	control[0] = u;
	control[1] = beta; 

	return control;
}

////////////////////////////////////////////////////////////////////////////////////////////
real interceptor::Hamiltonian_2(real const& t, mstate const& X) const{
	// current state value (state + costate)
	real h =		X[0];		//h
	real v =		X[1];		//v
	real theta =	X[2];		//theta
	real phi =		X[3];		//phi
	real L =		X[4];		//L
	real l =		X[5];		//l
	real p_h =		X[6];		//p_h
	real p_v =		X[7];		//p_v
	real p_theta =	X[8];		//p_theta
	real p_phi =	X[9];		//p_phi
	real p_L =		X[10];		//p_L
	real p_l =		X[11];		//p_l

	// if mode = 1 : propulsion, if mode = 0 : no propulsion
	int mode = GetMode(t, X);

	// temporary variables
	real qm = mode*data->parameters.q*data->parameters.mu_gft;			// mass flow rate
	real mass = ComputeMass(t, X);										//mass
	real c_max = data->parameters.c0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;	// max curvature at alatitude h
	real d = data->parameters.d0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;		// max drag at altitude h
	real r = h + data->R_Earth;											// distance to the Earth center
	real g = data->mu0/r/r*data->parameters.mu_gft;						// gravity
    real ft = data->parameters.ve*qm;									// propulsion
    real eta = data->parameters.eta;
    real hr = data->parameters.hr;
    real alpha_max = data->parameters.alpha_max;
    real R_Earth = data->R_Earth;
    real u_max = data->parameters.u_max;
    real a_max = data->parameters.a_max;

    // control computation (linear approximation, no singular control considered here)
	mcontrol control = Control_2(t, X);
	real u = control[0];
	real beta = control[1];
	real alpha = alpha_max*u;	// angle of attack

	// H is computed
	real H = 	p_L*v*cos(theta)*sin(phi)/r
				+p_l*v*sin(theta)/(r*cos(L))
				-p_h*v*cos(theta)*cos(phi)
				+p_theta*(v*c_max*u*cos(beta) + v*sin(theta)*(cos(phi) + sin(phi)*tan(L))/r	+ (ft*sin(alpha)*cos(beta)/(mass*v) - g*sin(theta)*cos(phi)/v))
				+p_phi*(-v*c_max*u*sin(beta)/cos(theta)	+ v*cos(theta)*(sin(phi) + tan(theta)*tan(theta)*(sin(phi) - tan(L)*cos(phi)))/r - (ft*sin(alpha)*sin(beta)/(mass*v*cos(theta)) + g*sin(phi)/(v*cos(theta))))
				-p_v*((d+eta*c_max*u*u)*v*v - g*cos(theta)*cos(phi) - ft*cos(alpha)/mass)
				+ data->parameters.muC*u*u/2;

	return H;

}

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::mstate  interceptor::ConversionState12(mstate const& X1) const{
	mstate X2 = X1;

	const real eps = 1e-18;
	// current state value (state + costate)
	real h = 		X1[0];		//h
	real v = 		X1[1];		//v
	real gamma = 	X1[2];		//gamma
	real chi = 		X1[3];		//chi
	real L = 		X1[4];		//L
	real l = 		X1[5];		//l

	real r = h + data->R_Earth;											// distance to the Earth center

	// This calculus comes from matching the unitary vectors i and a of the two velocity frame (see PhD dissertation)
	if (gamma == M_PI/2.0){
		X2[2] = 0;
		X2[3] = -M_PI;
	}else if (gamma == -M_PI/2.0){
		X2[2] = 0;
		X2[3] = 0;
	}else{
		X2[2] = acos(sqrt(sin(gamma)*sin(gamma) + cos(gamma)*cos(gamma)*cos(chi)*cos(chi)));
		if (cos(gamma)*sin(chi) < 0)
			X2[2] = -acos(sqrt(sin(gamma)*sin(gamma) + cos(gamma)*cos(gamma)*cos(chi)*cos(chi))); // sin(theta) = cos(gamma)*sin(chi)  and  0 <= acos(*) <= pi

		real sinPhi = cos(gamma)*cos(chi)/cos(X2[2]);
		if (fabs(sinPhi) < eps && sin(gamma)/cos(X2[2]) < 0) X2[3] = 0; // cos(phi) = - sin(gamma)/cos(theta)
		else if (fabs(sinPhi) < eps && sin(gamma)/cos(X2[2]) > 0) X2[3] = -M_PI;
		else if (sinPhi > 0) X2[3] = acos(-sin(gamma)/cos(X2[2]));
		else X2[3] = -acos(-sin(gamma)/cos(X2[2]));
	}

	real theta = X2[2];
	real phi = X2[3];

	Eigen::Matrix<real,6,6> Jac1, Jac2;
	Eigen::Matrix<real,6,1> p1, p2, pTemp;

	Jac1(0,0) = cos(L)*cos(l); Jac1(1,0) = -r*sin(L)*cos(l); Jac1(2,0) = -r*cos(L)*sin(l); Jac1(3,0) = 0.0; Jac1(4,0) = 0.0; Jac1(5,0) = 0.0;

	Jac1(0,1) = cos(L)*sin(l); Jac1(1,1) = -r*sin(L)*sin(l); Jac1(2,1) = r*cos(L)*cos(l); Jac1(3,1) = 0.0; Jac1(4,1) = 0.0; Jac1(5,1) = 0.0;

	Jac1(0,2) = sin(L); Jac1(1,2) = r*cos(L); Jac1(2,2) = 0.0; Jac1(3,2) = 0.0; Jac1(4,2) = 0.0; Jac1(5,2) = 0.0;

	Jac1(0,3) = 0.0; Jac1(1,3) = (-cos(L)*cos(l)*cos(gamma)*cos(chi) - sin(L)*cos(l)*sin(gamma))*v;
	Jac1(2,3) = (sin(L)*sin(l)*cos(gamma)*cos(chi) - cos(l)*cos(gamma)*sin(chi) - cos(L)*sin(l)*sin(gamma))*v;
	Jac1(3,3) = (sin(L)*cos(l)*sin(gamma)*cos(chi) + sin(l)*sin(gamma)*sin(chi) + cos(L)*cos(l)*cos(gamma))*v;
	Jac1(4,3) = (sin(L)*cos(l)*cos(gamma)*sin(chi) - sin(l)*cos(gamma)*cos(chi))*v;
	Jac1(5,3) = (-sin(L)*cos(l)*cos(gamma)*cos(chi) - sin(l)*cos(gamma)*sin(chi) + cos(L)*cos(l)*sin(gamma))*v;

	Jac1(0,4) = 0.0; Jac1(1,4) = (-cos(L)*sin(l)*cos(gamma)*cos(chi) - sin(L)*sin(l)*sin(gamma))*v;
	Jac1(2,4) = (-sin(L)*cos(l)*cos(gamma)*cos(chi) - sin(l)*cos(gamma)*sin(chi) + cos(L)*cos(l)*sin(gamma))*v;
	Jac1(3,4) = (sin(L)*sin(l)*sin(gamma)*cos(chi) - cos(l)*sin(gamma)*sin(chi) + cos(L)*sin(l)*cos(gamma))*v;
	Jac1(4,4) = (sin(L)*sin(l)*cos(gamma)*sin(chi) + cos(l)*cos(gamma)*cos(chi))*v;
	Jac1(5,4) = (-sin(L)*sin(l)*cos(gamma)*cos(chi) + cos(l)*cos(gamma)*sin(chi) + cos(L)*sin(l)*sin(gamma))*v;

	Jac1(0,5) = 0.0; Jac1(1,5) = (-sin(L)*cos(gamma)*cos(chi) + cos(L)*sin(gamma))*v;
	Jac1(2,5) = 0.0; Jac1(3,5) = (-cos(L)*sin(gamma)*cos(chi) + sin(L)*cos(gamma))*v;
	Jac1(4,5) = -cos(L)*cos(gamma)*sin(chi)*v; Jac1(5,5) = (cos(L)*cos(gamma)*cos(chi) + sin(L)*sin(gamma))*v;


	Jac2(0,0) = cos(L)*cos(l); Jac2(1,0) = -r*sin(L)*cos(l); Jac2(2,0) = -r*cos(L)*sin(l); Jac2(3,0) = 0.0; Jac2(4,0) = 0.0; Jac2(5,0) = 0.0;

	Jac2(0,1) = cos(L)*sin(l); Jac2(1,1) = -r*sin(L)*sin(l); Jac2(2,1) = r*cos(L)*cos(l); Jac2(3,1) = 0.0; Jac2(4,1) = 0.0; Jac2(5,1) = 0.0;

	Jac2(0,2) = sin(L); Jac2(1,2) = r*cos(L); Jac2(2,2) = 0.0; Jac2(3,2) = 0.0; Jac2(4,2) = 0.0; Jac2(5,2) = 0.0;		

	Jac2(0,3) = 0.0; Jac2(1,3) = (-cos(L)*cos(l)*cos(theta)*sin(phi) + sin(L)*cos(l)*cos(theta)*cos(phi))*v;
	Jac2(2,3) = (sin(L)*sin(l)*cos(theta)*sin(phi) - cos(l)*sin(theta) + cos(L)*sin(l)*cos(theta)*cos(phi))*v;
	Jac2(3,3) = (sin(L)*cos(l)*sin(theta)*sin(phi) - sin(l)*cos(theta) + cos(L)*cos(l)*sin(theta)*cos(phi))*v;
	Jac2(4,3) = (-sin(L)*cos(l)*cos(theta)*cos(phi) + cos(L)*cos(l)*cos(theta)*sin(phi))*v;
	Jac2(5,3) = (-sin(L)*cos(l)*cos(theta)*sin(phi) - sin(l)*sin(theta) - cos(L)*cos(l)*cos(theta)*cos(phi))*v;

	Jac2(0,4) = 0.0; Jac2(1,4) = (-cos(L)*sin(l)*cos(theta)*sin(phi) + sin(L)*sin(l)*cos(theta)*cos(phi))*v;
	Jac2(2,4) = (-sin(L)*cos(l)*cos(theta)*sin(phi) - sin(l)*sin(theta) - cos(L)*cos(l)*cos(theta)*cos(phi))*v;
	Jac2(3,4) = (sin(L)*sin(l)*sin(theta)*sin(phi) + cos(l)*cos(theta) + cos(L)*sin(l)*sin(theta)*cos(phi))*v;
	Jac2(4,4) = (-sin(L)*sin(l)*cos(theta)*cos(phi) + cos(L)*sin(l)*cos(theta)*sin(phi))*v;
	Jac2(5,4) = (-sin(L)*sin(l)*cos(theta)*sin(phi) + cos(l)*sin(theta) - cos(L)*sin(l)*cos(theta)*cos(phi))*v;

	Jac2(0,5) = 0.0; Jac2(1,5) = (-sin(L)*cos(theta)*sin(phi) - cos(L)*cos(theta)*cos(phi))*v;
	Jac2(2,5) = 0.0; Jac2(3,5) = (-cos(L)*sin(theta)*sin(phi) + sin(L)*sin(theta)*cos(phi))*v;
	Jac2(4,5) = (cos(L)*cos(theta)*cos(phi) + sin(L)*cos(theta)*sin(phi))*v;
	Jac2(5,5) = (cos(L)*cos(theta)*sin(phi) - sin(L)*cos(theta)*cos(phi))*v;


	p1(0) = X1[6]; 
	p1(1) = X1[10]; 
	p1(2) = X1[11]; 
	p1(3) = X1[8]; 
	p1(4) = X1[9]; 
	p1(5) = X1[7];

	pTemp = Jac1.lu().solve(p1);
	p2 = Jac2*pTemp;

	X2[6] = p2(0);
	X2[7] = p2(5);
	X2[8] = p2(3);
	X2[9] = p2(4);
	X2[10] = p2(1);
	X2[11] = p2(2);

	return X2;

}

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::mstate  interceptor::ConversionState21(mstate const& X2) const{
	mstate X1 = X2;

	const real eps = 1e-18;

	// current state value (state + costate)
	real h = 		X2[0];		//h
	real v = 		X2[1];		//v
	real theta = 	X2[2];		//theta
	real phi = 		X2[3];		//phi
	real L = 		X2[4];		//L
	real l = 		X2[5];		//l

	real r = h + data->R_Earth;											// distance to the Earth center

	// This calculus comes from matching the unitary vectors i and a of the two velocity frame (see PhD dissertation)
	if (theta == M_PI/2.0)
	{
		X1[2] = 0;
		X1[3] = M_PI/2.0;
	}
	else if (theta == -M_PI/2.0)
	{
		X1[2] = 0;
		X1[3] = -M_PI/2.0;
	}
	else
	{
		X1[2] = acos(sqrt(sin(theta)*sin(theta) + cos(theta)*cos(theta)*sin(phi)*sin(phi)));
		if (cos(theta)*cos(phi) > 0)
			X1[2] = -acos(sqrt(sin(theta)*sin(theta) + cos(theta)*cos(theta)*sin(phi)*sin(phi))); // sin(gamma) = -cos(theta)*sin(phi)  and  0 <= acos(*) <= pi

		real sinChi = sin(theta)/cos(X1[2]);
		if (fabs(sinChi) < eps && sin(phi)*cos(theta)/cos(X1[2]) > 0) X1[3] = 0; // cos(chi) = sin(phi)*cos(theta)/cos(gamma)
		else if (fabs(sinChi) < eps && sin(phi)*cos(theta)/cos(X1[2]) < 0) X1[3] = -M_PI;
		else if (sinChi > 0) X1[3] = acos(sin(phi)*cos(theta)/cos(X1[2]));
		else X1[3] = -acos(sin(phi)*cos(theta)/cos(X1[2]));
	}

	real gamma = X1[2];
	real chi = X1[3];

	Eigen::Matrix<real,6,6> Jac1, Jac2;
	Eigen::Matrix<real,6,1> p1, p2, pTemp;


	Jac1(0,0) = cos(L)*cos(l); Jac1(1,0) = -r*sin(L)*cos(l); Jac1(2,0) = -r*cos(L)*sin(l); Jac1(3,0) = 0.0; Jac1(4,0) = 0.0; Jac1(5,0) = 0.0;

	Jac1(0,1) = cos(L)*sin(l); Jac1(1,1) = -r*sin(L)*sin(l); Jac1(2,1) = r*cos(L)*cos(l); Jac1(3,1) = 0.0; Jac1(4,1) = 0.0; Jac1(5,1) = 0.0;

	Jac1(0,2) = sin(L); Jac1(1,2) = r*cos(L); Jac1(2,2) = 0.0; Jac1(3,2) = 0.0; Jac1(4,2) = 0.0; Jac1(5,2) = 0.0;

	Jac1(0,3) = 0.0; Jac1(1,3) = (-cos(L)*cos(l)*cos(gamma)*cos(chi) - sin(L)*cos(l)*sin(gamma))*v;
	Jac1(2,3) = (sin(L)*sin(l)*cos(gamma)*cos(chi) - cos(l)*cos(gamma)*sin(chi) - cos(L)*sin(l)*sin(gamma))*v;
	Jac1(3,3) = (sin(L)*cos(l)*sin(gamma)*cos(chi) + sin(l)*sin(gamma)*sin(chi) + cos(L)*cos(l)*cos(gamma))*v;
	Jac1(4,3) = (sin(L)*cos(l)*cos(gamma)*sin(chi) - sin(l)*cos(gamma)*cos(chi))*v;
	Jac1(5,3) = (-sin(L)*cos(l)*cos(gamma)*cos(chi) - sin(l)*cos(gamma)*sin(chi) + cos(L)*cos(l)*sin(gamma))*v;

	Jac1(0,4) = 0.0; Jac1(1,4) = (-cos(L)*sin(l)*cos(gamma)*cos(chi) - sin(L)*sin(l)*sin(gamma))*v;
	Jac1(2,4) = (-sin(L)*cos(l)*cos(gamma)*cos(chi) - sin(l)*cos(gamma)*sin(chi) + cos(L)*cos(l)*sin(gamma))*v;
	Jac1(3,4) = (sin(L)*sin(l)*sin(gamma)*cos(chi) - cos(l)*sin(gamma)*sin(chi) + cos(L)*sin(l)*cos(gamma))*v;
	Jac1(4,4) = (sin(L)*sin(l)*cos(gamma)*sin(chi) + cos(l)*cos(gamma)*cos(chi))*v;
	Jac1(5,4) = (-sin(L)*sin(l)*cos(gamma)*cos(chi) + cos(l)*cos(gamma)*sin(chi) + cos(L)*sin(l)*sin(gamma))*v;

	Jac1(0,5) = 0.0; Jac1(1,5) = (-sin(L)*cos(gamma)*cos(chi) + cos(L)*sin(gamma))*v;
	Jac1(2,5) = 0.0; Jac1(3,5) = (-cos(L)*sin(gamma)*cos(chi) + sin(L)*cos(gamma))*v;
	Jac1(4,5) = -cos(L)*cos(gamma)*sin(chi)*v; Jac1(5,5) = (cos(L)*cos(gamma)*cos(chi) + sin(L)*sin(gamma))*v;


	Jac2(0,0) = cos(L)*cos(l); Jac2(1,0) = -r*sin(L)*cos(l); Jac2(2,0) = -r*cos(L)*sin(l); Jac2(3,0) = 0.0; Jac2(4,0) = 0.0; Jac2(5,0) = 0.0;

	Jac2(0,1) = cos(L)*sin(l); Jac2(1,1) = -r*sin(L)*sin(l); Jac2(2,1) = r*cos(L)*cos(l); Jac2(3,1) = 0.0; Jac2(4,1) = 0.0; Jac2(5,1) = 0.0;

	Jac2(0,2) = sin(L); Jac2(1,2) = r*cos(L); Jac2(2,2) = 0.0; Jac2(3,2) = 0.0; Jac2(4,2) = 0.0; Jac2(5,2) = 0.0;		

	Jac2(0,3) = 0.0; Jac2(1,3) = (-cos(L)*cos(l)*cos(theta)*sin(phi) + sin(L)*cos(l)*cos(theta)*cos(phi))*v;
	Jac2(2,3) = (sin(L)*sin(l)*cos(theta)*sin(phi) - cos(l)*sin(theta) + cos(L)*sin(l)*cos(theta)*cos(phi))*v;
	Jac2(3,3) = (sin(L)*cos(l)*sin(theta)*sin(phi) - sin(l)*cos(theta) + cos(L)*cos(l)*sin(theta)*cos(phi))*v;
	Jac2(4,3) = (-sin(L)*cos(l)*cos(theta)*cos(phi) + cos(L)*cos(l)*cos(theta)*sin(phi))*v;
	Jac2(5,3) = (-sin(L)*cos(l)*cos(theta)*sin(phi) - sin(l)*sin(theta) - cos(L)*cos(l)*cos(theta)*cos(phi))*v;

	Jac2(0,4) = 0.0; Jac2(1,4) = (-cos(L)*sin(l)*cos(theta)*sin(phi) + sin(L)*sin(l)*cos(theta)*cos(phi))*v;
	Jac2(2,4) = (-sin(L)*cos(l)*cos(theta)*sin(phi) - sin(l)*sin(theta) - cos(L)*cos(l)*cos(theta)*cos(phi))*v;
	Jac2(3,4) = (sin(L)*sin(l)*sin(theta)*sin(phi) + cos(l)*cos(theta) + cos(L)*sin(l)*sin(theta)*cos(phi))*v;
	Jac2(4,4) = (-sin(L)*sin(l)*cos(theta)*cos(phi) + cos(L)*sin(l)*cos(theta)*sin(phi))*v;
	Jac2(5,4) = (-sin(L)*sin(l)*cos(theta)*sin(phi) + cos(l)*sin(theta) - cos(L)*sin(l)*cos(theta)*cos(phi))*v;

	Jac2(0,5) = 0.0; Jac2(1,5) = (-sin(L)*cos(theta)*sin(phi) - cos(L)*cos(theta)*cos(phi))*v;
	Jac2(2,5) = 0.0; Jac2(3,5) = (-cos(L)*sin(theta)*sin(phi) + sin(L)*sin(theta)*cos(phi))*v;
	Jac2(4,5) = (cos(L)*cos(theta)*cos(phi) + sin(L)*cos(theta)*sin(phi))*v;
	Jac2(5,5) = (cos(L)*cos(theta)*sin(phi) - sin(L)*cos(theta)*cos(phi))*v;

	p2(0) = X2[6]; 
	p2(1) = X2[10]; 
	p2(2) = X2[11]; 
	p2(3) = X2[8]; 
	p2(4) = X2[9]; 
	p2(5) = X2[7];

	pTemp = Jac2.lu().solve(p2);
	p1 = Jac1*pTemp;

	X1[6] = p1(0);
	X1[7] = p1(5);
	X1[8] = p1(3);
	X1[9] = p1(4);
	X1[10] = p1(1);
	X1[11] = p1(2);

	return X1;

}

////////////////////////////////////////////////////////////////////////////////////////////
void interceptor::InitAnalytical(real const& ti, mstate & Xi, real const& tf, mstate & Xf) const{
	real h = 		Xi[0];		//h
	real v = 		Xi[1];		//v
	real gamma = 	Xi[2];		//gamma
	real chi = 		Xi[3];		//chi
	real L = 		Xi[4];		//L
	real l = 		Xi[5];		//l

	real hf = 		Xf[0];		//h
	real vf = 		Xf[1];		//v
	real gammaf = 	Xf[2];		//gamma
	real chif = 	Xf[3];		//chi
	real Lf = 		Xf[4];		//L
	real lf = 		Xf[5];		//l

	// temporary variables
	real mass = ComputeMass(ti, Xi);									//mass
	real c_max = data->parameters.c0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;	// max curvature at alatitude h
	real d = data->parameters.d0*exp(-h/data->parameters.hr)*(data->parameters.propellant_mass+data->parameters.empty_mass)/mass;		// max drag at altitude h
	real r = h + data->R_Earth;											// distance to the Earth center
	real rf = hf + data->R_Earth;										// distance to the Earth center
	real g = data->mu0/r/r*data->parameters.mu_gft;						// gravity
	real eta = data->parameters.eta;
	real b = sqrt(c_max*d/(2*eta));
	real hr = data->parameters.hr;
	real alpha_max = data->parameters.alpha_max;
	real R_Earth = data->R_Earth;
	real R = sqrt((rf*cos(Lf)*cos(lf) - r*cos(L)*cos(l))*(rf*cos(Lf)*cos(lf) - r*cos(L)*cos(l)) +
                                    (rf*cos(Lf)*sin(lf) - r*cos(L)*sin(l))*(rf*cos(Lf)*sin(lf) - r*cos(L)*sin(l)) +
																		(rf*sin(Lf) - r*sin(L))*(rf*sin(Lf) - r*sin(L))); // distance to the RDV point
	real Rdot = -( (rf*cos(Lf)*cos(lf) - r*cos(L)*cos(l))*(sin(gamma)*cos(L)*cos(l) - cos(gamma)*cos(chi)*sin(L)*cos(l) - cos(gamma)*sin(chi)*sin(l)) +
                                (rf*cos(Lf)*sin(lf) - r*cos(L)*sin(l))*(sin(gamma)*cos(L)*sin(l) - cos(gamma)*cos(chi)*sin(L)*sin(l) + cos(gamma)*sin(chi)*cos(l)) +
                                (rf*sin(Lf) - r*sin(L))*(sin(gamma)*sin(L) + cos(L)*cos(gamma)*cos(chi)) ) / R;
	real bdot = -c_max*d*sin(gamma)*sqrt(2*eta/(c_max*d))/(2*eta*hr);
	real N1 = exp(b*R) - exp(-b*R) - 2*b*R, D1 = 4 + exp(b*R)*(b*R - 2) - exp(-b*R)*(b*R + 2);
	real dN1 = (bdot*R + b*Rdot)*(exp(b*R) + exp(-b*R) - 2), dD1 = (bdot*R + b*Rdot)*(exp(b*R)*(b*R - 2) + (exp(-b*R)*(b*R + 2)) + exp(b*R) - exp(-b*R));
	real k1 = b*R*(exp(b*R) - exp(-b*R) - 2*b*R)/(4 + exp(b*R)*(b*R - 2) - exp(-b*R)*(b*R + 2));
	real k1dot = (bdot*R + b*Rdot)*N1/D1 + (b*R*(dN1*D1 - dD1*N1)/(D1*D1));
	real N2 = exp(b*R)*(b*R - 1) + exp(-b*R)*(b*R + 1), D2 = 4 + exp(b*R)*(b*R - 2) - exp(-b*R)*(b*R + 2);
	real dN2 = b*R*(bdot*R + b*Rdot)*(exp(b*R) - exp(-b*R)), dD2 = (bdot*R + b*Rdot)*(exp(b*R)*(b*R - 2) + (exp(-b*R)*(b*R + 2)) + exp(b*R) - exp(-b*R));
	real k2 = b*R*(exp(b*R)*(b*R - 1) + exp(-b*R)*(b*R + 1))/(4 + exp(b*R)*(b*R - 2) - exp(-b*R)*(b*R + 2));
	real k2dot = (bdot*R + b*Rdot)*N2/D2 + (b*R*(dN2*D2 - dD2*N2)/(D2*D2));
	real k3 = 2 + k1 - k2;
	real k3dot = k1dot - k2dot;
	// compute lambda_1
	real l1 = 0;
	real dist = fabs(rf*(cos(L)*cos(Lf)*cos(l)*cos(lf) +
                    cos(L)*cos(Lf)*sin(l)*sin(lf) +
                    sin(L)*sin(Lf)) - r);
	real x_E_R = -cos(L)*cos(l)*(rf*cos(Lf)*cos(lf) - r*cos(L)*cos(l))
            - cos(L)*sin(l)*(rf*cos(Lf)*sin(lf) - r*cos(L)*sin(l))
            - sin(L)*(rf*sin(Lf) - r*sin(L));
	if (R == 0) l1 = gammaf;
	else if (dist/R >= 1 && x_E_R > 0) l1 = -M_PI/2.0;
	else if (dist/R >= 1) l1 = M_PI/2.0; 
	else if (x_E_R > 0) l1 = -asin(dist/R);
	else l1 = asin(dist/R);
	// compute lambda_2
	real l2;
	real tlam2 = r - rf*(cos(Lf)*cos(lf)*cos(L)*cos(l) + cos(Lf)*sin(lf)*cos(L)*sin(l) + sin(Lf)*sin(L));
	real normProj = sqrt((rf*cos(Lf)*cos(lf) + (tlam2-r)*cos(L)*cos(l))*(rf*cos(Lf)*cos(lf) + (tlam2-r)*cos(L)*cos(l))
                    + (rf*cos(Lf)*sin(lf) + (tlam2-r)*cos(L)*sin(l))*(rf*cos(Lf)*sin(lf) + (tlam2-r)*cos(L)*sin(l))
                    + (rf*sin(Lf) + (tlam2-r)*sin(L))*(rf*sin(Lf) + (tlam2-r)*sin(L)));
	real prodScal = -(rf*cos(Lf)*cos(lf) + (tlam2-r)*cos(L)*cos(l))*sin(L)*cos(l)
                - (rf*cos(Lf)*sin(lf) + (tlam2-r)*cos(L)*sin(l))*sin(L)*sin(l)
                + (rf*sin(Lf) + (tlam2-r)*sin(L))*cos(L);
	real coordProj_el = -sin(l)*(rf*cos(Lf)*cos(lf) + (tlam2-r)*cos(L)*cos(l)) + cos(l)*(rf*cos(Lf)*sin(lf) + (tlam2-r)*cos(L)*sin(l));
	if (normProj == 0) l2 = 0;
	else if (prodScal/normProj <= -1) l2 = M_PI;
	else if (prodScal/normProj >= 1) l2 = 0;
	else if (coordProj_el >= 0) l2 = acos(prodScal/normProj);
	else l2 = -acos(prodScal/normProj);
	// compute u1
	real u1 = -(k1*(gammaf - l1)/R + k2*sin(gamma - l1)/R + k3*cos(gamma)/(2*hr))/c_max;
	// compute u2
	real u2 = -(k1*(chif - l2)*cos(gamma)/R + k2*sin(chi - l2)*cos(gamma)/R)/c_max;
	// compute du1
	real d1DivC = sin(gamma)/(c_max*hr);
	real du1 = d1DivC*c_max*u1 -
        (k1dot*(gammaf - l1)/R + k1*sin(gamma - l1)/(R*R) - k1*(gammaf - l1)*Rdot/(R*R) + k2dot*sin(gamma - l1)/R +
        k2*cos(gamma - l1)*(c_max*u1 + sin(gamma - l1)/R)/R - k2*sin(gamma - l1)*Rdot/(R*R) + k3dot*cos(gamma)/(2*hr) -
        c_max*u1*k3*sin(gamma)/(2*hr))/c_max;
	// compute du2
	real du2 = d1DivC*c_max*u2 -
        (k1dot*cos(gamma)*(chif - l2)/R - k1*c_max*u1*sin(gamma)*(chif - l2)/R + k1*cos(gamma)*sin(chi - l2)/(R*R) -
        k1*cos(gamma)*(chif - l2)*Rdot/(R*R) + k2dot*cos(gamma)*sin(chi - l2)/R - k2*sin(gamma)*sin(chi - l2)*c_max*u1/R +
        k2*cos(gamma)*cos(chi - l2)*(c_max*u2/cos(gamma) + sin(chi - l2)/R)/R - k2*cos(gamma)*sin(chi - l2)*Rdot/(R*R))/c_max;

	Xi[7] =		-1;							//p_v
	Xi[8] =		2*eta*u1;					//p_gamma
	Xi[9] =	2*eta*u2*cos(gamma);			//p_chi
	Xi[6] =		(-sin(gamma)*c_max*u1*Xi[8]*cos(gamma) - sin(gamma)*eta*c_max*cos(gamma)*u1*u1 - sin(gamma)*d*cos(gamma) -
        sin(gamma)*eta*c_max*cos(gamma)*u2*u2 - sin(gamma)*c_max*u2*Xi[9] - c_max*u2*tan(gamma)*Xi[9]*cos(gamma) +
        2*du1*eta*cos(gamma)*cos(gamma)) / cos(gamma);				//p_h
	Xi[10] =	r*(sin(gamma)*c_max*u2*tan(gamma)*Xi[9]*cos(chi) - 2*sin(gamma)*sin(gamma)*sin(chi)*eta*cos(gamma)*du2 +
        2*sin(gamma)*sin(gamma)*sin(chi)*eta*cos(gamma)*c_max*u1*u2*tan(gamma) - 2*sin(gamma)*du1*eta*cos(gamma)*cos(chi) -
        c_max*u1*Xi[8]*cos(gamma)*cos(gamma)*cos(chi) - eta*c_max*cos(gamma)*cos(gamma)*u1*u1*cos(chi) -
        2*cos(gamma)*cos(gamma)*cos(gamma)*sin(chi)*eta*du2 + 2*cos(gamma)*cos(gamma)*cos(gamma)*sin(chi)*eta*c_max*u1*u2*tan(gamma) -
        d*cos(gamma)*cos(gamma)*cos(chi) - eta*c_max*cos(gamma)*cos(gamma)*u2*u2*cos(chi) - c_max*u2*Xi[9]*cos(chi)*cos(gamma)) / cos(gamma);		//p_L
	Xi[11] =	-r*cos(L)*(-sin(chi)*sin(gamma)*c_max*u2*tan(gamma)*Xi[9] + 2*sin(chi)*sin(gamma)*du1*eta*cos(gamma) +
                     sin(chi)*c_max*u1*Xi[8]*cos(gamma)*cos(gamma) + sin(chi)*eta*c_max*cos(gamma)*cos(gamma)*u1*u1 +
                     sin(chi)*d*cos(gamma)*cos(gamma) + sin(chi)*eta*c_max*cos(gamma)*cos(gamma)*u2*u2 + sin(chi)*c_max*u2*Xi[9]*cos(gamma) -
                     2*eta*cos(gamma)*cos(gamma)*cos(gamma)*du2*cos(chi) - 2*eta*cos(gamma)*du2*cos(chi)*sin(gamma)*sin(gamma) +
                     2*eta*cos(gamma)*cos(gamma)*cos(gamma)*c_max*u1*u2*tan(gamma)*cos(chi) +
                     2*eta*cos(gamma)*c_max*u1*u2*tan(gamma)*cos(chi)*sin(gamma)*sin(gamma)) / cos(gamma);											//p_l
	Xi[8] = v*Xi[8];
	Xi[9] = v*Xi[9];
	Xi[6] = v*Xi[6];
	Xi[10] = v*Xi[10];
	Xi[11] = v*Xi[11];

}

////////////////////////////////////////////////////////////////////////////////////////////
interceptor::mstate interceptor::SetChart(real const& t, mstate const& X){
	mstate Xs = X;

	// determine which chart to be used
	if (data->currentChart==1){
		if (fabs(cos(Xs[2]))>=data->chartLimit){
			// state remains in chart 1
		}else{
			// switch chart
			Xs = ConversionState12(Xs);
			data->currentChart = 2;
		}
	}else{
		if ((fabs(cos(Xs[2]))>=data->chartLimit)){
			// state remains in chart 2
		}else{
			// switch chart
			Xs = ConversionState21(Xs);
			data->currentChart = 1;
		}
	}

	return Xs;
}

////////////////////////////////////////////////////////////////////////////////////////////
real interceptor::ComputeMass(real const& t, mstate const& X) const {
	// if mode = 1 : propulsion, if mode = 0 : no propulsion
	int mode = GetMode(t, X);

	// temporary variables
	real qm = data->parameters.q*data->parameters.mu_gft;		// mass flow rate
	real t1 = data->parameters.propellant_mass / data->parameters.q;

	real mass = data->parameters.empty_mass + data->parameters.propellant_mass;
	if (mode==1)
		mass = data->parameters.empty_mass + data->parameters.propellant_mass - qm*t;
	else
		mass = data->parameters.empty_mass + data->parameters.propellant_mass - qm*t1;

	return mass;
}