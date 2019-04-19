/*
 * shooting.hpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE & Riccardo BONALLI (ONERA/DTIS)
 */

#include "commonType.hpp"

#include "model.hpp"

#ifndef _SHOOTING_H_
#define _SHOOTING_H_


/**************			shooting class			******************************/
class shooting
{

public:

	/**
	* Constructor
	* @param model a reference to the model used for shooting
	* @param numMulti a number for multiple shooting (from 1 shooting to infinity)
	* @param numThread a number for multithreading (from 1 thread to infinity)
	*/
	shooting(model & model, int numMulti, int numThread);

	/**
	* Destructor
	*/
	~shooting();

	/**
	* Init Shooting
	* @param ti initial time
	* @param Xi initial state
	* @param tf final time
	* @param Xf final state
	*/
	void InitShooting(real const& ti, model::mstate const& Xi, real const& tf, model::mstate const& Xf) const;

	/**
	* Init Shooting using multiple shooting
	* @param vt vector of times
	* @param vX vector of states
	*/
	void InitShooting(std::vector<real> const& vt, std::vector<model::mstate> const& vX) const;

	/**
	* Set desired state
	* @param ti initial time
	* @param Xi initial state
	* @param tf final time
	* @param Xf final state
	*/
	void SetDesiredState(real const& ti, model::mstate const& Xi, real const& tf, model::mstate const& Xf) const;

	/**
	* Set desired state using multiple shooting
	* @param vt vector of times
	* @param vX vector of states
	*/
	void SetDesiredState(std::vector<real> const& vt, std::vector<model::mstate> const& vX) const;

	/**
	* Solve OCP with discrete time and state continuation
	* @param continuationStep step of continuation as a ratio : 0 for no continuation, 0 < continuationStep <= 1 for continuation
	* @return result from the solver : 1 if success
	*/
	int SolveOCP(real const& continuationStep) const;

	/**
	* Solve OCP with discrete continuation on a real parameter from Rdata to Rgoal
	* @param continuationStep step of continuation as a ratio : 0 for no continuation, 0 < continuationStep <= 1 for continuation
	* @param Rdata real data for continuation
	* @param Rgoal real number to achieve
	* @return result from the solver : 1 if success (Rgoal is achieved)
	*/
	int SolveOCP(real const& continuationStep, real & Rdata, real const& Rgoal) const;

	/**
	* Function to move the model from ti to tf
	* @param ti initial time
	* @param Xi initial state
	* @param tf final time
	* @return final state
	*/
	model::mstate Move(real const& ti, model::mstate const& Xi, real const& tf) const;

	/**
	* Function to move the model from ti to tf
	* @param ti initial time
	* @param Xi initial state
	* @param tf final time
	* @param Xf final state
	*/
	void Move(real const& ti, model::mstate const& Xi, real const& tf, model::mstate & Xf) const;

	/**
	* Function to move the model to tf using initial or previously solved solution and handling multiple shooting
	* @return final state
	*/
	model::mstate Move(real const& tf) const;

	/**
	* Function to move the model to tf using initial or previously solved solution and handling multiple shooting
	* @param final state
	*/
	void Move(real const& tf, model::mstate & Xf) const;

	/**
	* Set mode for fixed/free final time and final state
	* @param mode_tf mode for final time : 0 for fixed final time, 1 for free final time
	* @param mode_Xf mode for final state : 0 for fixed final state, 1 for free final state
	*/
	void SetMode(int const& mode_tf, std::vector<int> const& mode_Xf) const;

	/**
	* Set mode for fixed/free time and state using multiple shooting
	* @param mode_t mode for final times : 0 for fixed time, 1 for free time, 2 for dependent time
	* @param mode_X vector of mode for states : 0 for fixed state, 1 for free state, 2 for continuity state
	*/
	void SetMode(std::vector<int> const& mode_t, std::vector< std::vector<int> > const& mode_X) const;

	/**
	* Set relative tolerance for the solver
	* @param xtol realtive tolerance
	*/
	void SetPrecision(real const& xtol);

	/**
	* Set minimum continuation step
	* @param step minimum continuation step
	*/
	void SetContinuationMinStep(real const& step);

	/**
	* Get current kth guess parameter
	* @param k the parameter number to be returned
	* @return the kth guess parameter
	*/
	real GetParameters(int const& k) const;

	/**
	* Get current guess parameters
	* @return a pointer to the array of guess parameters 
	*/
	real *GetParameters() const;

	/**
	* Get current guess parameters
	* @param paramVector a vector of guess parameters
	*/
	void GetParameters(std::vector<real> & paramVector) const;

	/**
	* Get the current solution
	* @param vt vector of times
	* @param vX vector of states
	*/
	void GetSolution(std::vector<real> & vt, std::vector<model::mstate> & vX) const;

	/**
	* Get the number of calls to the function to be solved
	* @return the number of calls to the function to be solved
	*/
	double GetCallNumber() const;

	/**
	* Trace in the specified file (in model) by handling multiple shooting
	*/
	void Trace() const;							

private:

	/**
	* Compute Shooting without continuation
	* @return result from the solver : 1 if success
	*/
	int SolveShooting() const;

	/**
	* Compute Shooting with discrete time and state continuation
	* @param continuationStep step of continuation as a ratio : 0 for no continuation, 0 < continuationStep <= 1 for continuation
	* @return result from the solver : 1 if success
	*/
	int SolveShootingContinuation(real const& continuationStep) const;

	/**
	* Compute Shooting with discrete continuation on parameter Rdata to Rgoal
	* @param continuationStep step of continuation as a ratio : 0 for no continuation, 0 < continuationStep <= 1 for continuation
	* @param Rdata real data for continuation
	* @param Rgoal real number to achieve
	* @return result from the solver : 1 if success (Rgoal is achieved)
	*/
	int SolveShootingContinuation(real const& continuationStep, real & Rdata, real const& Rgoal) const; 

	/**
	* Static function for the solver
	* @param userdata a pointer to user data
	* @param n the number of parameters
	* @param param a pointer to the parameters
	* @param fvec a pointer to the function values
	* @param iflag an optional flag
	*/
	static int StaticShootingFunction(void *userdata, int n, const real *param, real *fvec, int iflag);		// for cminpack

	/**
	* Function to be solved
	* @param n the number of parameters
	* @param param a vector to the parameters
	* @param fvec a pointer to the function values
	* @param iflag a pointer to an optional flag
	*/
	void ShootingFunction(int n, std::vector<real> const& param, std::vector<real> & fvec, int iflag) const;

	/**
	* Function to be solved using multithreading
	* @param n a pointer to the number of parameters
	* @param param a vector to the parameters
	* @param fvec a pointer to the function values
	* @param iflag a pointer to an optional flag
	*/
	void ShootingFunctionPar(int n, std::vector<real> const& param, std::vector<real> & fvec, int iflag) const;

	/**
	* Static function for the thread
	* @param arg a pointer to the arguments passed to the thread
	*/
	static void StaticShootingPar(void *arg);

	/**
	* Compute a part of the shooting function on a single thread
	* @param n the number of parameters
	* @param param a pointer to the parameters
	* @param fvec a pointer to the function values
	* @param iflag a pointer to an optional flag
	* @param threadNum thread number
	* @param timeLine time vector considered by the thread
	*/
	void ShootingFunctionParThread(int n, real *param, real *fvec, int threadNum, std::vector<real> const& timeLine) const;
	
	/**
	* Shooting data structure
	*/
	struct data_struct;

	/**
	* Shooting data pointer
	*/
	data_struct *data;

	/**
	* Solve Function
	* @param numParam number of parameters
	* @param param a vector to the parameters
	* @return result from the solver : 1 if success
	*/
	int SolveShootingFunction(int const & numParam, std::vector<real> & param) const;

	/**
	* Update solution
	*/
	void UpdateSolution() const;

	/**
	* Multiple Shooting function
	* @param t the free time
	* @param X the state at t
	* @param Xp the parameter state at t
	* @param Xd the desired state at t
	* @param mode_X the mode (free/fixed state)
	* @param fvec a pointer to the values of the function
	*/
	void MultipleShootingFunction(real const& t, model::mstate const& X, model::mstate const& Xp, model::mstate const& Xd, std::vector<int> const& mode_X, model::mstate & fvec) const;

	/**
	* Compute timeline
	* @param param a vector to the parameters
	* @param timeLine a vector to the timeline
	*/
	void ComputeTimeLine(std::vector<real> const& param, std::vector<real> & timeLine) const;

};

#endif //_SHOOTING_H_
