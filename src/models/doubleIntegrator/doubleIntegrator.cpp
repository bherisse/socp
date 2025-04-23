/*
 * doubleIntegrator.cpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE
 */
#include <fstream>

#include <math.h>

#include "doubleIntegrator.hpp"

/**
* Vehicle data
*/
struct doubleIntegrator::data_struct{
	doubleIntegrator::parameters_struct parameters;			///< those parameters can be used for continuation
	int  stepNbr;											///< step number for ModelInt
	std::string strFileTrace;								///< trace file
	std::ofstream fileTrace;								///< trace file stream
};

/**
* Constructor
*/
doubleIntegrator::doubleIntegrator(int modelOrder, std::string the_fileTrace) : model(6, modelOrder, 30, the_fileTrace) {

	// vehicle data
	data = new data_struct;
	data->parameters.u_max = 1;				// max normalized control
	data->parameters.a_max = 1;				// max acceleration
	data->parameters.muT = 0.01;			// weight for time cost

};

/**
* Destructor
*/
doubleIntegrator::~doubleIntegrator(){

	// close trace file
  	data->fileTrace.close();
	delete(data);
};

/**
* State model of the vehicle
*/
doubleIntegrator::mstate doubleIntegrator::Model(real const& t, mstate const& X, int isJac) const{
	mstate Xdot(2 * dim);
	
	if (isJac == 0) {
		Xdot = ModelState(t, X);
	}
	else {
		Xdot.resize(4 * dim * dim);
		Xdot = ModelJacobian(t, X);
	}

	return Xdot;

};

/**
* State model of the vehicle
*/
doubleIntegrator::mstate doubleIntegrator::ModelState(real const& t, mstate const& X) const {
	mstate Xdot(X.size());

	// current state value
	real x = X[0];
	real y = X[1];
	real z = X[2];
	real vx = X[3];
	real vy = X[4];
	real vz = X[5];
	real p_x = X[6];
	real p_y = X[7];
	real p_z = X[8];
	real p_vx = X[9];
	real p_vy = X[10];
	real p_vz = X[11];

	real norm_v = sqrt(vx*vx + vy * vy + vz * vz);

	// temporary variables
	real a_max = data->parameters.a_max;

	// control computation
	mcontrol u = Control(t, X);

	// state and costate equations
	Xdot[0] = vx;
	Xdot[1] = vy;
	Xdot[2] = vz;
	Xdot[3] = a_max * u[0];
	Xdot[4] = a_max * u[1];
	Xdot[5] = a_max * u[2];
	Xdot[6] = 0;
	Xdot[7] = 0;
	Xdot[8] = 0;
	Xdot[9] = -p_x;
	Xdot[10] = -p_y;
	Xdot[11] = -p_z;

	return Xdot;

};

/**
* State model of the vehicle
*/
doubleIntegrator::mstate doubleIntegrator::ModelJacobian(real const& t, mstate const& X) const {
	mstate XJdot;
	mstate f(2 * dim, 0);
	mstate dR_dt(4 * dim * dim, 0);
	mstate df_dX;
	mstate dfx_dX(2 * dim, 0);
	mstate dfy_dX(2 * dim, 0);
	mstate dfz_dX(2 * dim, 0);
	mstate dfvx_dX(2 * dim, 0);
	mstate dfvy_dX(2 * dim, 0);
	mstate dfvz_dX(2 * dim, 0);
	mstate dfpx_dX(2 * dim, 0);
	mstate dfpy_dX(2 * dim, 0);
	mstate dfpz_dX(2 * dim, 0);
	mstate dfpvx_dX(2 * dim, 0);
	mstate dfpvy_dX(2 * dim, 0);
	mstate dfpvz_dX(2 * dim, 0);


	// current state value
	real x = X[0];
	real y = X[1];
	real z = X[2];
	real vx = X[3];
	real vy = X[4];
	real vz = X[5];
	real p_x = X[6];
	real p_y = X[7];
	real p_z = X[8];
	real p_vx = X[9];
	real p_vy = X[10];
	real p_vz = X[11];

	real norm_v = sqrt(vx*vx + vy * vy + vz * vz);

	// temporary variables
	real a_max = data->parameters.a_max;

	// control computation
	mcontrol u = Control(t, X);

	// state and costate equations
	f[0] = vx;
	f[1] = vy;
	f[2] = vz;
	f[3] = a_max * u[0];
	f[4] = a_max * u[1];
	f[5] = a_max * u[2];
	f[6] = 0;
	f[7] = 0;
	f[8] = 0;
	f[9] = -p_x;
	f[10] = -p_y;
	f[11] = -p_z;

	dfx_dX.assign({ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0});
	dfy_dX.assign({ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 });
	dfz_dX.assign({ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 });
	dfvx_dX.assign({ 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0 }); // a revoir en cas de saturation sur la commande
	dfvy_dX.assign({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0 }); // a revoir en cas de saturation sur la commande
	dfvz_dX.assign({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1 }); // a revoir en cas de saturation sur la commande
	dfpx_dX.assign({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
	dfpy_dX.assign({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
	dfpz_dX.assign({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
	dfpvx_dX.assign({ 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0 });
	dfpvy_dX.assign({ 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0 });
	dfpvz_dX.assign({ 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0 });
	df_dX.insert(df_dX.end(), dfx_dX.begin(), dfx_dX.end());
	df_dX.insert(df_dX.end(), dfy_dX.begin(), dfy_dX.end());
	df_dX.insert(df_dX.end(), dfz_dX.begin(), dfz_dX.end());
	df_dX.insert(df_dX.end(), dfvx_dX.begin(), dfvx_dX.end());
	df_dX.insert(df_dX.end(), dfvy_dX.begin(), dfvy_dX.end());
	df_dX.insert(df_dX.end(), dfvz_dX.begin(), dfvz_dX.end());
	df_dX.insert(df_dX.end(), dfpx_dX.begin(), dfpx_dX.end());
	df_dX.insert(df_dX.end(), dfpy_dX.begin(), dfpy_dX.end());
	df_dX.insert(df_dX.end(), dfpz_dX.begin(), dfpz_dX.end());
	df_dX.insert(df_dX.end(), dfpvx_dX.begin(), dfpvx_dX.end());
	df_dX.insert(df_dX.end(), dfpvy_dX.begin(), dfpvy_dX.end());
	df_dX.insert(df_dX.end(), dfpvz_dX.begin(), dfpvz_dX.end());

	for (int i = 0; i < 2 * dim; i++) {
		for (int j = 0; j < 2 * dim; j++) {
			for (int k = 0; k < 2 * dim; k++) {
				dR_dt[2 * dim * i + j] += df_dX[2 * dim *i + k] * X[2 * dim * (k + 1) + j];
				
			}
		}
	}

	//std::cout << "Xdot = " << dR_dt.size() << std::endl;
	//for (int i = 0; i < dR_dt.size(); i++) {
	//	std::cout << dR_dt[i] << ' ';
	//}
	//std::cout << std::endl;

	XJdot.insert(XJdot.end(), f.begin(), f.end());
	XJdot.insert(XJdot.end(), dR_dt.begin(), dR_dt.end());

	return XJdot;

};

/**
* Control model of the vehicle
*/
doubleIntegrator::mcontrol doubleIntegrator::Control(real const& t, mstate const& X) const{
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
	mcontrol u(3);
	u[0] = -p_vx/a_max;
	u[1] = -p_vy/a_max; 
	u[2] = -p_vz/a_max;

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

/**
* Hamiltonian of the vehicle
*/
doubleIntegrator::mstate doubleIntegrator::Hamiltonian(real const& t, mstate const& X, int isJac) const{
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
	mcontrol u = Control(t, X);
	real norm_u = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);

	mstate H(1);
	if (isJac == 0) {
		// H is computed here
		H[0] = data->parameters.muT + a_max * a_max*norm_u*norm_u / 2 + p_x * vx + p_y * vy + p_z * vz + a_max * (p_vx*u[0] + p_vy * u[1] + p_vz * u[2]);
	}
	else {
		H.resize(2 * dim + 1);
		H.assign({ 0., 0., 0., p_x, p_y, p_z, vx, vy, vz, -p_vx , -p_vy, -p_vz, 0.});
	}

	return H;
}


/**
* Get parameters pointer
*/
doubleIntegrator::parameters_struct  & doubleIntegrator::GetParameterData(){
	return data->parameters;
}

/**
* Set number of integration steps
*/
void doubleIntegrator::SetStepNumber(int step){
	data->stepNbr = step;	
}