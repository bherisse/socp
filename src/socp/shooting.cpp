/*
 * shooting.cpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE & Riccardo BONALLI (ONERA/DTIS)
 */

#include <iostream>
#include <math.h>
#include <algorithm>
#include <assert.h>

#include "mThread.hpp"      // thread for parallel computing
#include "shooting.hpp"

#include <cminpack.h>

/**
* Shooting data
*/
struct shooting::data_struct{
	std::vector<real> timed;				///< desired vector of times
	std::vector<real> time;					///< current vector of times (for continuation)
	std::vector<real> time_prec;			///< previous vector of times (for continuation)
	std::vector<model::mstate> Xd;			///< desired vector of states
	std::vector<model::mstate> X;			///< current vector of states (for continuation)
	std::vector<model::mstate> X_prec; 		///< previous vector of states (for continuation)
	int dim;								///< state dimension of the model
	int numMulti;							///< number of trajectories for multiple shooting
	std::vector<real> tab_param;			///< parameters to find
	std::vector<real> tab_param_temp;		///< temporary parameters to find
	real continuationStepMin;				///< real parameter for homotopy
	std::vector< std::vector<int> > mode_X;	///< 0 for fixed X, 1 for free X, 2 for continuity
	std::vector<int> mode_t;				///< 0 for fixed time, 1 for free time, 2 for continuity
	int numParam;							///< number of parameters to be found

	int numThread;							///< Number of threads used for multithreading
	std::vector< std::vector<void*> > userThreadData;	///< vector of user data for multithreading
	std::vector< mThread > mulThread;		///< vector of threads
	std::vector< int > threadNum;			///< vector of thread numbers

	int n;									///< number of parameters for the solver
	int maxfev; 							///< max call fcn for the solver
	real xtol;								///< relative tolerance for the solver
	real epsfcn; 							///< for forward-difference approximation for the solver
	int scalingMode; 						///< scaling (1: by the fcn, 2: by xscal) for the solver
	real factor; 							///< initial step bound for the solver
	int nprint; 							///< for the solver
	int info; 								///< returned info (1 if OK) by the solver
	int nfev; 								///< iteration number of the solver
	int stopFlag;							///< flag to stop the solver

};

// Constructor
shooting::shooting(model & model, int numMulti, int numThread): myModel(model){
	// Shooting data
	data = new data_struct;

	// Verify that parameters are correct
	try
	{
		if (numMulti<1) 
			throw std::string("ERROR : numMulti should be superior or equal to 1"); 
		else 
			data->numMulti = numMulti;
		if (numThread<1) 
			throw std::string("ERROR : numThread should be superior or equal to 1"); 
		else 
			data->numThread = numThread;	
	}
	catch(std::string const& chaine)
    {
		std::cerr << std::endl << chaine.c_str() << std::endl;
		exit(1);
    }

	data->time_prec = std::vector<real>(data->numMulti+1);
	data->time = std::vector<real>(data->numMulti+1);
	data->timed = std::vector<real>(data->numMulti+1);
	data->X_prec = std::vector<model::mstate>(data->numMulti+1);
	data->X = std::vector<model::mstate>(data->numMulti+1);
	data->Xd = std::vector<model::mstate>(data->numMulti+1);
	data->dim = myModel.GetDim();
	const int n_param_max = 2*data->dim*data->numMulti + (1+data->numMulti); // initial costate dim + junction states for multiple shooting + time vector
	data->tab_param.reserve(n_param_max);
	data->tab_param_temp.reserve(n_param_max);
	data->continuationStepMin = 1e-12;										// if b is lower than continuationStepMin, homotopy is stoped
	data->mode_X = std::vector< std::vector<int> >(data->numMulti+1);		// 0 for fixed X, 1 for free X
	data->mode_t = std::vector<int>(data->numMulti+1);						// 0 for fixed time, 1 for free time
	data->numParam = n_param_max;											// number of parameters to be found
	
	// Parameters used by the hybrd algorithm
	data->maxfev = 10000; 									// max call fcn
	data->xtol = 1e-8;										// relative tolerance
	myModel.SetODEIntPrecision(data->xtol);			// Set default precision for model integration
	data->epsfcn = 1e-15; 									// for forward-difference approximation
	data->scalingMode = 1; 									// scaling (1: by the fcn, 2: by xscal)
	data->factor = 1; 										// initial step bound
	data->nprint = 0; 										// for iflag
	data->info = 0; 										// returned info (1 if OK)
	data->nfev = 0; 										// iteration number
	data->stopFlag = 0;
	
	// For multithreading
	data->userThreadData = std::vector< std::vector<void*> >(data->numThread);
	data->mulThread = std::vector< mThread >(data->numThread); //the thread
	data->threadNum = std::vector<int>(data->numThread); 
	for (int i=0; i<data->numThread; i++){
		data->userThreadData[i] = std::vector<void*>(6);
	}

};

// Destructor
shooting::~shooting(){
	delete data;
};

// Resize
void shooting::Resize(int numMulti, int numThread) const {
	// Verify that parameters are correct
	try
	{
		if (numMulti<1)
			throw std::string("ERROR : numMulti should be superior or equal to 1");
		else
			data->numMulti = numMulti;
		if (numThread<1)
			throw std::string("ERROR : numThread should be superior or equal to 1");
		else
			data->numThread = numThread;
	}
	catch (std::string const& chaine)
	{
		std::cerr << std::endl << chaine.c_str() << std::endl;
		exit(1);
	}

	data->time_prec.resize(data->numMulti + 1);
	data->time.resize(data->numMulti + 1);
	data->timed.resize(data->numMulti + 1);
	data->X_prec.resize(data->numMulti + 1);
	data->X.resize(data->numMulti + 1);
	data->Xd.resize(data->numMulti + 1);
	const int n_param_max = 2 * data->dim*data->numMulti + (1 + data->numMulti); // initial costate dim + junction states for multiple shooting + time vector
	data->tab_param.reserve(n_param_max);
	data->tab_param_temp.reserve(n_param_max);
	data->mode_X.resize(data->numMulti + 1);		// 0 for fixed X, 1 for free X
	data->mode_t.resize(data->numMulti + 1);		// 0 for fixed time, 1 for free time
	data->numParam = n_param_max;					// number of parameters to be found

													// For multithreading
	data->userThreadData.resize(data->numThread);
	data->mulThread = std::vector< mThread >(data->numThread);
	data->threadNum.resize(data->numThread);
	for (int i = 0; i<data->numThread; i++) {
		data->userThreadData[i] = std::vector<void*>(6);
	}
};

// Set mode
void shooting::SetMode(int const& mode_tf, std::vector<int> const& mode_Xf) const{
	// by default, initial state is fixed with fixed initial time
	data->mode_t[0] = model::FIXED;
	data->mode_X[0] = std::vector<int>(data->dim, model::FIXED);
	// by default, continuity conditions for interior points
	for (int i=1; i<data->numMulti; i++){
		data->mode_X[i] = std::vector<int>(data->dim, model::CONTINUOUS);
		data->mode_t[i] = model::CONTINUOUS;
	}
	// define final mode
	data->mode_X[data->numMulti] = mode_Xf;
	data->mode_t[data->numMulti] = mode_tf;

	// update internal data
	data->numParam = 2*data->dim*data->numMulti + mode_tf;
	data->tab_param.resize(data->numParam);
	data->tab_param_temp.resize(data->numParam);
};

// Set mode
void shooting::SetMode(std::vector<int> const& mode_t, std::vector< std::vector<int> > const& mode_X) const{
	for (int i=0; i<data->numMulti+1; i++){
		data->mode_X[i] = mode_X[i];
	}
	data->mode_t = mode_t;

	// update internal data
	int sum_mode_t = 0;
	for(std::vector<int>::iterator it = data->mode_t.begin(); it != data->mode_t.end(); ++it){
		if(*it==1) sum_mode_t += *it;
	}
	data->numParam = 2*data->dim*data->numMulti + sum_mode_t;
	data->tab_param.resize(data->numParam);
	data->tab_param_temp.resize(data->numParam);
};

// Init Shooting
void shooting::InitShooting(real const& ti, model::mstate const& Xi, real const& tf, model::mstate const& Xf) const
{
	int nMultiple = data->numMulti;
	int numParam = data->numParam;

	// defaut values for initial state, final state (used to initialize the shooting algorithm)
	// time
	for (int i = 0; i<nMultiple + 1; i++) {
		data->time[i] = ti + i*(tf - ti) / (nMultiple);
		data->timed[i] = ti + i*(tf - ti) / (nMultiple);
		data->time_prec[i] = ti + i*(tf - ti) / (nMultiple);
	}
	// state
	data->X[0] = Xi;
	data->Xd[0] = Xi;
	data->X_prec[0] = Xi;
	for (int i = 1; i<nMultiple; i++) {
		data->X[i] = Move(ti, Xi, data->time[i]);
		data->Xd[i] = Move(ti, Xi, data->timed[i]);
		data->X_prec[i] = Move(ti, Xi, data->time_prec[i]);
	}
	data->X[nMultiple] = Xf;
	data->Xd[nMultiple] = Xf;
	data->X_prec[nMultiple] = Xf;

	// guess for optimal parameters
	for (int i = 0; i<2 * data->dim; i++) {
		data->tab_param[i] = data->X[0][i];	
	}
	for (int j = 1; j<nMultiple; j++) {
		for (int i = 0; i<2 * data->dim; i++) {
			data->tab_param[j * 2 * data->dim + i] = data->X[j][i];		// state and costate vector as parameters from t(1) to t(nMultiple-1)
		}
	}

	int nbrParam = 2 * data->dim*data->numMulti;
	for (int j = 0; j<nMultiple + 1; j++) {
		if (data->mode_t[j] == model::FREE) {
			nbrParam += 1;
			data->tab_param[nbrParam - 1] = data->time[j];				// if t[k] (k>0) is free, it becomes a parameter 
		}
	}

};

// Init Shooting
void shooting::InitShooting(std::vector<real> const& vt, std::vector<model::mstate> const& vX) const
{
	int nMultiple = vt.size() - 1;
	int numParam = data->numParam;

	// defaut values for initial state, final state (used to initialize the shooting algorithm)
	// time
	data->time_prec = vt;
	data->time = vt;
	data->timed = vt;
	for (int i = 0; i<nMultiple + 1; i++) {
		data->time_prec[i] = vt[i];
		data->time[i] = vt[i];
		data->timed[i] = vt[i];
	}
	// state
	data->X_prec = vX;
	data->X = vX;
	data->Xd = vX;
	for (int i = 0; i<nMultiple + 1; i++) {
		data->X_prec[i] = vX[i];
		data->X[i] = vX[i];
		data->Xd[i] = vX[i];
	}

	// guess for optimal parameters
	for (int i = 0; i<2 * data->dim; i++) {
		data->tab_param[i] = data->X[0][i];		// initial costate vector is the parameter
	}
	for (int j = 1; j<nMultiple; j++) {
		for (int i = 0; i<2 * data->dim; i++) {
			data->tab_param[j * 2 * data->dim + i] = data->X[j][i];		// state and costate vector as parameters from t(1) to t(nMultiple-1)
		}
	}

	int nbrParam = 2 * data->dim*data->numMulti;
	for (int j = 0; j<nMultiple + 1; j++) {
		if (data->mode_t[j] == model::FREE) {
			nbrParam += 1;
			data->tab_param[nbrParam - 1] = data->time[j];				// if t[k] (k>0) is free, it becomes a parameter 
		}
	}

};

// Set desired state
void shooting::SetDesiredState(real const& ti, model::mstate const& Xi, real const& tf, model::mstate const& Xf) const {
	// initial time
	data->timed[0] = ti;
	// initial state
	data->Xd[0] = Xi;			
	// final time
	data->timed[data->numMulti] = tf;
	// final state
	data->Xd[data->numMulti] = Xf;
};

// Set desired state
void shooting::SetDesiredState(std::vector<real> const& vt, std::vector<model::mstate> const& vX) const {
	int nMultiple = vt.size() - 1;
	for (int i = 0; i<nMultiple + 1; i++) {
		data->timed[i] = vt[i];
		data->Xd[i] = vX[i];
	}
};

// Solve OCP
int shooting::SolveOCP(real const& continuationStep) const {
	if (continuationStep <= 0) {
		// no continuation
		return SolveShooting();
	}
	else {
		// discrete continuation with step "continuationStep"
		return SolveShootingContinuation(continuationStep);
	}

	return 0;
}

// Solve OCP with timeout
int shooting::SolveOCP(real const& continuationStep, double const& timeoutMS) const {
	std::vector<void*> OCPargs(3);
	int info = 0;

	OCPargs[0] = (void*) this;
	OCPargs[1] = (void*)&continuationStep;
	OCPargs[2] = (void*)&info;

	mThread OCPthread;
	OCPthread.create(TimeoutOCPFunction, (void*)OCPargs.data());

	if (OCPthread.join_for(timeoutMS) == -1) {
		data->stopFlag = -1;
		OCPthread.~mThread();
		data->stopFlag = 0;
		return -1;
	}
	else
		return info;
}

// Solve OCP with data continuation
int shooting::SolveOCP(real const& continuationStep, real & Rdata, real const& Rgoal) const {
	if (continuationStep <= 0) {
		// no continuation
		return SolveShootingContinuation(1.0, Rdata, Rgoal);
	}
	else {
		// discrete continuation with step "continuationStep"
		return SolveShootingContinuation(continuationStep, Rdata, Rgoal);
	}

	return 0;
}

// Function to move the model from ti to tf
model::mstate shooting::Move(real const& ti, model::mstate const& Xi, real const& tf) const {

	// Compute trajectory from t0 to tf
	model::mstate X_tf = myModel.ComputeTraj(ti, Xi, tf, 0);

	return X_tf;

}

// Function to move the model from ti to tf
void shooting::Move(real const& ti, model::mstate const& Xi, real const& tf, model::mstate & Xf) const {

	// Compute trajectory from t0 to tf
	Xf = myModel.ComputeTraj(ti, Xi, tf, 0);

}

// Function to move the model to tf
model::mstate shooting::Move(real const& tf) const {
	int numParam = data->numParam;

	// Initial model state
	model::mstate X0 = data->X[0];
	for (int i = 0; i<2 * data->dim; i++)	X0[i] = data->tab_param[i];

	// initial time
	real t0 = data->time[0];
	if (data->mode_t[0] == model::FIXED) {
		t0 = data->time[0];
	}
	else {
		t0 = data->tab_param[2 * data->numMulti*data->dim];
	}
	// final time
	real t_f, tf_traj;
	if (data->mode_t[data->numMulti] == model::FIXED) {
		tf_traj = data->time[data->numMulti];
	}
	else {
		tf_traj = data->tab_param[numParam - 1];
	}
	t_f = tf_traj;
	if (tf >= t0 && tf <= tf_traj) {
		t_f = tf;
	}

	// Get the timeline
	std::vector<real> timeLine(data->numMulti + 1);
	ComputeTimeLine(data->tab_param, timeLine);

	// Find which piece of trajectory include t_f
	real t1, t2;
	int iter = 0;
	t1 = timeLine[0];
	t2 = timeLine[1];
	while (t2<t_f) {
		iter += 1;
		t1 = timeLine[iter];
		t2 = timeLine[iter + 1];
	}
	// Compute trajectory from t1 to tf
	model::mstate X1 = X0;
	int index = 2 * iter*data->dim;
	if (iter > 0 && iter < data->numMulti) {
		// update data
		for (int i = 0; i<2 * data->dim; i++)	X1[i] = data->tab_param[index + i];
	}

	model::mstate X_tf = Move(t1, X1, t_f);

	return X_tf;

}

// Function to move the model to tf
void shooting::Move(real const& tf, model::mstate & Xf) const {

	Xf = Move(tf);

}

// Set relative tolerance for the solver
void shooting::SetPrecision(real const& xtol) const {
	data->xtol = xtol;
	myModel.SetODEIntPrecision(data->xtol);			// Set same precision for model integration
}

// Set minimum continuation step
void shooting::SetContinuationMinStep(real const& step) const{
	data->continuationStepMin = step;
}

// Get current kth parameter
real shooting::GetParameters(int const& k) const {
	return data->tab_param[k];
};

// Get current parameters as a pointer to an array
real *shooting::GetParameters() const {
	real *param = new real[data->numParam];
	for (int k = 0; k<data->numParam; k++) {
		param[k] = data->tab_param[k];
	}

	return param;
};

// Get current parameters as a vector
void shooting::GetParameters(std::vector<real> & paramVector) const {
	paramVector.resize(data->numParam);
	for (int k = 0; k<data->numParam; k++) {
		paramVector[k] = data->tab_param[k];
	}
};

// Get solution
void shooting::GetSolution(std::vector<real>& vt, std::vector<model::mstate>& vX) const {
	// update from computed solution (from data->tab_param to data->time and data->X)
	UpdateSolution();
	vt = data->time;
	for (int i = 0; i<data->numMulti + 1; i++) {
		vX[i] = data->X[i];
	}
}

// Get the number of calls to the function to be solved
double shooting::GetCallNumber() const {
	return data->nfev;
}

// Trace
void shooting::Trace() const
{
	int numParam = data->numParam;

	UpdateSolution();

	// Initial model state
	model::mstate X0 = data->X[0];
	for (int i = 0; i<2 * data->dim; i++)	X0[i] = data->tab_param[i];

	// initial time
	real t0 = data->time[0];
	if (data->mode_t[0] == model::FIXED) {
		t0 = data->time[0];
	}
	else {
		t0 = data->tab_param[2 * data->numMulti*data->dim];
	}

	// final time
	real tf;
	if (data->mode_t[data->numMulti] == model::FIXED) {
		tf = data->time[data->numMulti];
	}
	else {
		tf = data->tab_param[numParam - 1];
	}

	// Get the timeline
	std::vector<real> timeLine(data->numMulti + 1);
	ComputeTimeLine(data->tab_param, timeLine);

	// Compute trajectory from t0 to tf
	model::mstate X1 = X0;
	real t1, t2;
	int index;
	for (int i = 0; i<data->numMulti; i++) {
		t1 = timeLine[i];
		t2 = timeLine[i + 1];
		index = 2 * (i + 1)*data->dim;
		// Compute trajectory from t1 to t2
		myModel.ComputeTraj(t1, X1, t2, 1);
		if (i<(data->numMulti - 1)) {
			// update data
			for (int j = 0; j<2 * data->dim; j++)	X1[j] = data->tab_param[index + j];
		}
	}

}

// Get Model
model & shooting::GetModel() const {
	return myModel;
}

// Static function for solving OCP with timeout
void shooting::TimeoutOCPFunction(void *arg){
	void **userArgs = (void**) arg;
	shooting* This = (shooting*) userArgs[0];
	real* continuationStep = (real*) userArgs[1];
	int* info = (int*) userArgs[2];

	if (*continuationStep<=0){
		// no continuation
		*info = This->SolveShooting();
	}else{
		// discrete continuation with step "continuationStep"
		*info = This->SolveShootingContinuation(*continuationStep);
	}
}

// Solve Shooting
int shooting::SolveShooting() const
{
	int numParam = data->numParam;

	// temporary values are initialized
	for (int k=0;k<numParam;k++){
		data->tab_param_temp[k] = data->tab_param[k];
	}

	for (int i=0; i<data->numMulti+1; i++){
		// time
		data->time[i] = data->timed[i];
		// state
		for(int j=0;j<data->dim;j++) data->X[i][j] = data->Xd[i][j];
	}
	
	// solve the shooting problem
	int ret = SolveShootingFunction(numParam, data->tab_param_temp);

	// update data
	if(ret == 1){
		for (int k=0;k<numParam;k++){
			data->tab_param[k] = data->tab_param_temp[k];
		}
	}

	return ret;
}

// Solve Shooting with continuation step continuationStep
int shooting::SolveShootingContinuation(real const& continuationStep) const
{
	int numParam = data->numParam;

	// parameter for spatial continuation (discrete continuation algortihm)
	real bStep = continuationStep;
	real b = (std::min<real>)(continuationStep, 1.0);  // 1 is the max step
	real b_prec = 0;

	for (int i=0; i<data->numMulti+1; i++){
		// time
		data->time[i] = (1-b) * data->time_prec[i] + b * data->timed[i];
		// state
		for(int j=0;j<data->dim;j++) data->X[i][j] = (1-b)*data->X_prec[i][j] + b*data->Xd[i][j];
	}

	// temporary values are initialized
	for (int k=0;k<numParam;k++){
		data->tab_param_temp[k] = data->tab_param[k];
	}

	bool condition = true;
	int ret;
	while(condition){
		//printf("b = %lf \n",b);

		// solve the shooting problem
		ret = SolveShootingFunction(numParam, data->tab_param_temp);

		if (ret != 1){
			if (fabs(b-b_prec)<data->continuationStepMin){
				// continuation do not give satisfactory results
				//std::cout << "bmin reached (b = " << b << ")" << std::endl;
				condition = false;
			}

			// continuation parameters is decreased
			b = b_prec + (b-b_prec)/2;

			// data->tab_param_temp[k] is re-initialized
			for (int k=0;k<numParam;k++){
				data->tab_param_temp[k] = data->tab_param[k];
			}

			for (int i=0; i<data->numMulti+1; i++){
				// time
				data->time[i] = (1-b) * data->time_prec[i] + b * data->timed[i];
				// state
				for(int j=0;j<data->dim;j++) data->X[i][j] = (1-b)*data->X_prec[i][j] + b*data->Xd[i][j];
			}

		}

		if(ret == 1){
			if (b==1){
				// continuation worked
				condition = false;
			}else{
				// continuation is not finished
				b_prec = b;
				b = (std::min<real>)(b+bStep,1.0);
				for (int k=0;k<numParam;k++){
					data->tab_param[k] = data->tab_param_temp[k];
				}

				for (int i=0; i<data->numMulti+1; i++){
					// time
					data->time[i] = (1-b) * data->time_prec[i] + b * data->timed[i];
					// state
					for(int j=0;j<data->dim;j++) data->X[i][j] = (1-b)*data->X_prec[i][j] + b*data->Xd[i][j];
				}
				
			}
		}

	}

	// update data
	if(ret == 1){
		for (int k=0;k<numParam;k++){
			data->tab_param[k] = data->tab_param_temp[k];
		}

		for (int i=0; i<data->numMulti+1; i++){
			// time
			data->time_prec[i] = data->timed[i];
			// state
			for(int j=0;j<data->dim;j++) data->X_prec[i][j] = data->Xd[i][j];
		}

	}

	// return
	return ret;
}

// Solve Shooting for any other real data Rdata with continuation step continuationStep
int shooting::SolveShootingContinuation(real const& continuationStep, real & Rdata, real const& Rgoal) const
{
	real Rstart = Rdata; // starting value

	int numParam = data->numParam;

	// continuation parameter (discrete continuation algortihm)
	real bStep = continuationStep;
	real b = (std::min<real>)(continuationStep,1.0);  // 1 is the max step
	real b_prec = 0;

	// initial data value
	Rdata = (1-b)*Rstart + b*Rgoal;

	// temporary values are initialized
	for (int k=0;k<numParam;k++){
		data->tab_param_temp[k] = data->tab_param[k];
	}

	for (int i = 0; i<data->numMulti + 1; i++) {
		// time
		data->time[i] = data->timed[i];
		// state
		for (int j = 0; j<data->dim; j++) data->X[i][j] = data->Xd[i][j];
	}

	bool condition = true;
	int ret;
	while(condition){
		//printf("data = %lf \n",Rdata);
		//printf("b = %lf \n",b);

		// solve the shooting problem
		ret = SolveShootingFunction(numParam, data->tab_param_temp);

		if (ret != 1){
			if (fabs(b-b_prec)<data->continuationStepMin){
				// continuation do not give satisfactory results
				//std::cout << "bmin reached (b = " << b << ")" << std::endl;
				condition = false;
			}

			// continuation parameters is decreased
			b = b_prec + (b-b_prec)/2;

			// data->tab_param_temp[k] is re-initialized
			for (int k=0;k<numParam;k++){
				data->tab_param_temp[k] = data->tab_param[k];
			}

			// *Rdata is updated
			Rdata = (1-b)*Rstart + b*Rgoal;

		}

		if(ret == 1){
			if (b==1){
				// continuation worked
				condition = false;
			}else{
				// continuation is not finished
				b_prec = b;
				b = (std::min<real>)(b+bStep,1.0);

				for (int k=0;k<numParam;k++){
					data->tab_param[k] = data->tab_param_temp[k];
				}

				Rdata = (1-b)*Rstart + b*Rgoal;
			}
		}

	}

	// update data
	if(ret == 1){
		for (int k=0;k<numParam;k++){
			data->tab_param[k] = data->tab_param_temp[k];
		}
	}

	// return
	return ret;
}

// Solve root finding
int shooting::SolveShootingFunction(int const & numParam, std::vector<real> & param) const {

	// Parameters used by the hybrd algorithm
	std::vector<real> xscal(numParam);
	for (int i = 0; i<numParam; i++) {
		xscal[i] = 1;
	}
	std::vector<real> fvec(numParam);
	int ml = numParam - 1; 									// n-1 a least if the jacobian is not "lower band"
	int mu = numParam - 1; 									// n-1 a least if the jacobian is not "upper band"
	std::vector<real> fjac(numParam*numParam); 				// jacobian (Q)
	int ldfjac = numParam; 									// leading dimension
	int lr = numParam*(numParam + 1) / 2; 						// integer not less than (n*(n+1))/2
	std::vector<real> r(lr); 								// from qr factorization
	std::vector<real> qtf(numParam);						// (Q transpose)*fvec
	std::vector<real> wa1(numParam);						// work arrays
	std::vector<real> wa2(numParam);
	std::vector<real> wa3(numParam);
	std::vector<real> wa4(numParam);

	// root finding of the function "StaticShootingFunction"
	data->info = __cminpack_func__(hybrd)(StaticShootingFunction,
		(void*) this,
		numParam,
		param.data(),
		fvec.data(),
		data->xtol,
		data->maxfev,
		ml,
		mu,
		data->epsfcn,
		xscal.data(),
		data->scalingMode,
		data->factor,
		data->nprint,
		&data->nfev,
		fjac.data(),
		ldfjac,
		r.data(),
		lr,
		qtf.data(),
		wa1.data(),
		wa2.data(),
		wa3.data(),
		wa4.data());

	return data->info;

}

// Static function for the solver
int shooting::StaticShootingFunction(void *userData, int n, const real *param, real *fvec, int iflag) {		// for cminpack
	shooting* This = (shooting*)userData;
	std::vector<real> param_vec(n), fvec_vec(n);
	param_vec.assign(param, param + n);

	if (This->data->numThread <= 1) {
		This->ShootingFunction(n, param_vec, fvec_vec, iflag);
	}
	else {
		This->ShootingFunctionParallel(n, param_vec, fvec_vec, iflag);
	}

	for (int j = 0; j<n; j++) fvec[j] = fvec_vec[j];

	return This->data->stopFlag;
};

// Function for root finding
void shooting::ShootingFunction(int n, std::vector<real> const& param, std::vector<real> & fvec, int iflag) const
{
	int numParam = data->numParam;

	// Initial model state
	model::mstate X_t0 = data->X[0];
	for (int i = 0; i<2 * data->dim; i++)	X_t0[i] = param[i];

	// Get the timeline
	std::vector<real> timeLine(data->numMulti + 1);
	ComputeTimeLine(param, timeLine);

	// Process multiple shooting
	model::mstate X1 = X_t0;		// current state
	model::mstate X_tf = X_t0;	// next state to be computed
	model::mstate Xp = X1;
	real t1, t2;
	int index;
	int nbrParam = 2 * data->dim*data->numMulti;
	model::mstate funcMS(2 * data->dim);
	for (int i = 0; i<data->numMulti; i++) {
		// current time
		t1 = timeLine[i];
		t2 = timeLine[i + 1];
		// Compute trajectory from t1 to t2
		X_tf = Move(t1, X1, t2);
		index = 2 * (i + 1)*data->dim;
		if (i == 0) {
			// Compute the function to be solved
			std::vector<real> func;
			if (data->mode_t[0] == model::FIXED) {
				func = std::vector<real>(data->dim);
				myModel.InitialFunction(timeLine[0], X1, data->X[0], data->mode_X[0], func);
				for (int k = 0; k<data->dim; k++) fvec[k] = func[k];
			}
			else {
				func = std::vector<real>(data->dim + 1);
				myModel.InitialHFunction(timeLine[0], X1, data->X[0], data->mode_X[0], func);
				for (int k = 0; k<data->dim; k++) fvec[k] = func[k];
				fvec[2 * data->dim*data->numMulti] = func[data->dim];
				nbrParam += 1;
			}
		}
		if (i<(data->numMulti - 1)) {
			for (int k = 0; k<2 * data->dim; k++)	Xp[k] = param[index + k];
			// Switching function for free intermediate times
			if (data->mode_t[i + 1] == model::FREE) {
				real Si = myModel.SwitchingTimesFunction(t2, X_tf);
				fvec[nbrParam] = Si;
				nbrParam += 1;
			}
			// continuity function for multiple shooting
			MultipleShootingFunction(t2, X_tf, Xp, data->X[i + 1], data->mode_X[i + 1], funcMS);
			for (int k = 0; k<2 * data->dim; k++) fvec[index + k] = funcMS[k];
			X1 = Xp;
		}
		if (i == (data->numMulti - 1)) {
			// Compute the function to be solved
			std::vector<real> func;
			if (data->mode_t[data->numMulti] == model::FIXED) {
				func = std::vector<real>(data->dim);
				myModel.FinalFunction(t2, X_tf, data->X[data->numMulti], data->mode_X[data->numMulti], func);
				for (int k = 0; k<data->dim; k++) fvec[k + data->dim] = func[k];
			}
			else {
				func = std::vector<real>(data->dim + 1);
				myModel.FinalHFunction(t2, X_tf, data->X[data->numMulti], data->mode_X[data->numMulti], func);
				for (int k = 0; k<data->dim; k++) fvec[k + data->dim] = func[k];
				fvec[nbrParam] = func[data->dim];
				nbrParam += 1;
				assert(nbrParam == data->numParam);
			}
		}

	}

}

// Function for root finding (Parallelized)
void shooting::ShootingFunctionParallel(int n, std::vector<real> const& param, std::vector<real> & fvec, int iflag) const
{
	int numParam = data->numParam;

	// Get the timeline
	std::vector<real> timeLine(data->numMulti + 1);
	ComputeTimeLine(param, timeLine);

	// Process multiple shooting
	for (int i = 0; i<data->numThread; i++) {
		data->threadNum[i] = i;
		data->userThreadData[i][0] = (void*) this;
		data->userThreadData[i][1] = (void*)&n;
		data->userThreadData[i][2] = (void*)param.data();
		data->userThreadData[i][3] = (void*)fvec.data();
		data->userThreadData[i][4] = (void*)&data->threadNum[i];
		data->userThreadData[i][5] = (void*)&timeLine;

		// Create thread for multi-shooting
		data->mulThread[i].create(StaticShootingParallelFunction, (void*)data->userThreadData[i].data());
	}

	for (int i = 0; i<data->numThread; i++) {
		data->mulThread[i].join();
	}
}

// Static function for parallel shooting
void shooting::StaticShootingParallelFunction(void *arg) {
	void **userArgs = (void**)arg;
	shooting* This = (shooting*)userArgs[0];
	int* n = (int*)userArgs[1];
	real* param = (real*)userArgs[2];
	real* fvec = (real*)userArgs[3];
	int* threadNum = (int*)userArgs[4];
	std::vector<real> *timeLine = (std::vector<real>*) userArgs[5];

	This->ShootingParallelFunctionThread(*n, param, fvec, *threadNum, *timeLine);
}

// Function for parallel shooting
void shooting::ShootingParallelFunctionThread(int n, real *param, real *fvec, int threadNum, std::vector<real> const& timeLine) const {
	// Data
	model::mstate X1(data->X[0].size());
	model::mstate X2(data->X[0].size());
	int index;
	int nMulti;
	int j, start;

	nMulti = floor((real)data->numMulti / data->numThread);
	int remainder = data->numMulti%data->numThread;
	if (threadNum < remainder) {
		nMulti += 1; 	// allocate remaining pieces to threads
		start = threadNum*nMulti;
	}
	else {
		start = threadNum*nMulti + remainder;
	}

	// parameter to get current switching time
	int nbrParam = 2 * data->dim*data->numMulti;
	for (int k = 0; k <= start; k++) {											// get first switching time processed by current thread
		if (data->mode_t[k] == model::FREE) {
			nbrParam += 1;
		}
	}

	model::mstate Xp = X1;
	model::mstate funcMS(2 * data->dim);
	real t1, t2;
	for (int i = 0; i<nMulti; i++) {
		j = start + i; 		// index of the considered traj
							// current time
		t1 = timeLine[j];
		t2 = timeLine[j + 1];

		// Initial state
		if (i == 0) {
			for (int k = 0; k<2 * data->dim; k++)	X1[k] = param[2 * j*data->dim + k];
		}

		// Compute trajectory from t1 to t2
		X2 = Move(t1, X1, t2);

		index = 2 * (j + 1)*data->dim;
		if (j == 0) {
			// Compute the function to be solved
			std::vector<real> func;
			if (data->mode_t[0] == model::FIXED) {
				func = std::vector<real>(data->dim);
				myModel.InitialFunction(timeLine[0], X1, data->X[0], data->mode_X[0], func);
				for (int k = 0; k<data->dim; k++) fvec[k] = func[k];
			}
			else {
				func = std::vector<real>(data->dim + 1);
				myModel.InitialHFunction(timeLine[0], X1, data->X[0], data->mode_X[0], func);
				for (int k = 0; k<data->dim; k++) fvec[k] = func[k];
				fvec[2 * data->dim*data->numMulti] = func[data->dim];
			}
		}
		if (j<(data->numMulti - 1)) {
			for (int k = 0; k<2 * data->dim; k++)	Xp[k] = param[index + k];
			// Switching function for free intermediate times
			if (data->mode_t[j + 1] == model::FREE) {
				real Si = myModel.SwitchingTimesFunction(t2, X2);
				fvec[nbrParam] = Si;
				nbrParam += 1;
			}
			// continuity function for multiple shooting
			MultipleShootingFunction(t2, X2, Xp, data->X[j + 1], data->mode_X[j + 1], funcMS);
			for (int k = 0; k<2 * data->dim; k++) fvec[index + k] = funcMS[k];
			X1 = Xp;
		}
		if (j == (data->numMulti - 1)) {
			// Compute the function to be solved
			std::vector<real> func;
			if (data->mode_t[data->numMulti] == model::FIXED) {
				func = std::vector<real>(data->dim);
				myModel.FinalFunction(t2, X2, data->X[data->numMulti], data->mode_X[data->numMulti], func);
				for (int k = 0; k<data->dim; k++) fvec[k + data->dim] = func[k];
			}
			else {
				func = std::vector<real>(data->dim + 1);
				myModel.FinalHFunction(t2, X2, data->X[data->numMulti], data->mode_X[data->numMulti], func);
				for (int k = 0; k<data->dim; k++) fvec[k + data->dim] = func[k];
				fvec[nbrParam] = func[data->dim];
				nbrParam += 1;
				assert(nbrParam == data->numParam);
			}
		}

	}

};

// Update solution
void shooting::UpdateSolution() const{
	int numParam = data->numParam;

	// Initial model state
	model::mstate X0 = data->X[0]; 
	for(int i=0;i<2*data->dim;i++)	X0[i] = data->tab_param[i];

	// initial time
	real t0 = data->time[0];
	if (data->mode_t[0]==model::FIXED){
		t0 = data->time[0];
	}else{
		t0 = data->tab_param[2*data->numMulti*data->dim];
	}
	// final time
	real t_f, tf_traj;
	if (data->mode_t[data->numMulti]==model::FIXED){
		tf_traj = data->time[data->numMulti];
	}else{
		tf_traj = data->tab_param[numParam-1];
	}
	t_f = tf_traj;

	// Get the timeline
	std::vector<real> timeLine(data->numMulti+1);
	ComputeTimeLine(data->tab_param, timeLine);

	// Compute trajectory from t0 to tf
	model::mstate X1 = X0;
	real t1, t2;
	int index;
	for (int i=0; i<data->numMulti+1; i++){
		// update time
		data->time[i] = timeLine[i];
		// update state
		data->X[i] = X1;

		index = 2*(i+1)*data->dim;
		if (i<(data->numMulti-1)){
			// update data
			for(int j=0;j<2*data->dim;j++)	X1[j] = data->tab_param[index+j];
		}else{
			X1 = Move(t_f);
		}
	}

}

// Multiple Shooting function
void shooting::MultipleShootingFunction(real const& t, model::mstate const& X, model::mstate const& Xp, model::mstate const& Xd, std::vector<int> const& mode_X, model::mstate & fvec) const{
	int dim = data->dim;
	for (int j=0;j<dim;j++){
		switch (mode_X[j]){
			case model::FIXED :		fvec[j] = X[j] - Xd[j];
								fvec[j+dim] = Xp[j] - Xd[j];
								break;
			case model::FREE :			myModel.SwitchingStateFunction(t, j, X, Xp, Xd, fvec);
								break;
			case model::CONTINUOUS :	fvec[j] = X[j] - Xp[j];
								fvec[j+dim] = X[j+dim] - Xp[j+dim];
								break;
			default :			fvec[j] = X[j] - Xp[j];
								fvec[j+dim] = X[j+dim] - Xp[j+dim];
								break;
		}
	}
};

// Compute timeline
void shooting::ComputeTimeLine(std::vector<real> const& param, std::vector<real> & timeLine) const{
	assert(timeLine.size() == data->numMulti+1);

	// Get the timeline
	std::vector<real> SwitchingTimes;
	int nbrParam = 2*data->dim*data->numMulti;
	int currentJunction = 0;
	for(int j=0;j<data->numMulti+1;j++){
		// check mode tj
		if(data->mode_t[j]==model::FIXED){
			// if mode_tj=0 (fixed time), tj is loaded from curent time vector
			timeLine[j] = data->time[j];
			// update timeline between currentJunction and j
			for (int k=currentJunction+1;k<j;k++){
				// in that case, data->mode_t[j] != 0 or 1, tk is defined between currentJunction and j
				timeLine[k] = timeLine[currentJunction] + (k-currentJunction)*(timeLine[j]-timeLine[currentJunction])/(j-currentJunction);
			}
			currentJunction = j;
		}
		if(data->mode_t[j]==model::FREE){
			// if mode_tj=1 (free time), tj is a parameter 
			nbrParam+=1;
			//switchNbr+=1;
			timeLine[j] = param[nbrParam-1];
			// if j < data->numMulti, tj is a switching time (the model or the costate can be discontinuous at this time)
			if(j<data->numMulti) SwitchingTimes.push_back(timeLine[j]);
			// update timeline between currentJunction and j
			for (int k=currentJunction+1;k<j;k++){
				// in that case, data->mode_t[j] != 0 or 1, tk is defined between currentJunction and j
				timeLine[k] = timeLine[currentJunction] + (k-currentJunction)*(timeLine[j]-timeLine[currentJunction])/(j-currentJunction);
			}
			currentJunction = j;
		}
		
	}
	// update switching times in the model (to manage discontinuities)
	myModel.SwitchingTimesUpdate(SwitchingTimes);	

}
