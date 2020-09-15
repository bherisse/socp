/*
 * shooting.hpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE & Riccardo BONALLI (ONERA/DTIS)
 */
#include <unordered_map>

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
	* @param numMulti an integer for multiple shooting (from 1 shooting to infinity)
	* @param numThread an integer for multithreading (from 1 thread to infinity)
	*/
	shooting(model & model, int numMulti = int(1), int numThread = int(1));

	/**
	* Destructor
	*/
	~shooting();

	/**
	* Resize the OCP (optional, to be called before SetMode() function)
	* @param numMulti an integer for multiple shooting (from 1 shooting to infinity)
	* @param numThread an integer for multithreading (from 1 thread to infinity)
	*/
	void Resize(int numMulti, int numThread) const;

	/**
	* Set mode for fixed/free final time and final state, default values for other components (to be called before InitShooting() function)
	* @param mode_tf mode for final time : 0 for fixed final time, 1 for free final time
	* @param mode_Xf mode for final state : 0 for fixed final state, 1 for free final state
	*/
	void SetMode(int const& mode_tf, std::vector<int> const& mode_Xf) const;

	/**
	* Set mode for fixed/free time and state (to be called before InitShooting() function)
	* @param mode_t mode for final times : 0 for fixed time, 1 for free time, 2 for continuity
	* @param mode_X vector of mode for states : 0 for fixed state, 1 for free state, 2 for continuity
	*/
	void SetMode(std::vector<int> const& mode_t, std::vector< std::vector<int> > const& mode_X) const;

	/**
	* Init Shooting (to be called before SetDesiredState() function)
	* @param ti initial time
	* @param Xi initial state
	* @param tf final time
	* @param Xf final state
	*/
	void InitShooting(real const& ti, model::mstate const& Xi, real const& tf, model::mstate const& Xf) const;

	/**
	* Init Shooting using multiple shooting (to be called before SetDesiredState() function)
	* @param vt vector of times
	* @param vX vector of states
	*/
	void InitShooting(std::vector<real> const& vt, std::vector<model::mstate> const& vX) const;

	/**
	* Set desired state (optional, to be called before SolveOCP() function)
	* @param ti initial time
	* @param Xi initial state
	* @param tf final time
	* @param Xf final state
	*/
	void SetDesiredState(real const& ti, model::mstate const& Xi, real const& tf, model::mstate const& Xf) const;

	/**
	* Set desired state using multiple shooting (optional, to be called before SolveOCP() function)
	* @param vt vector of times
	* @param vX vector of states
	*/
	void SetDesiredState(std::vector<real> const& vt, std::vector<model::mstate> const& vX) const;

	/**
	* Solve OCP with discrete time and state continuation
	* @param continuationStep step of continuation as a ratio : 0 for no continuation, 0 < continuationStep <= 1 for continuation
	* @return result from the solver: 1 if success
	*/
	int SolveOCP(real const& continuationStep) const;

	/**
	* Solve OCP with discrete time and state continuation and timeout in ms
	* @param continuationStep step of continuation as a ratio : 0 for no continuation, 0 < continuationStep <= 1 for continuation
	* @param timeoutMS timeout in ms
	* @return result from the solver: 1 if success, -1 if timeout exceeded
	*/
	int SolveOCP(real const& continuationStep, double const& timeoutMS) const;

	/**
	* Solve OCP with discrete continuation on a real parameter from Rdata to Rgoal
	* @param continuationStep step of continuation as a ratio : 0 for no continuation, 0 < continuationStep <= 1 for continuation
	* @param Rdata real data for continuation
	* @param Rgoal real number to achieve
	* @return result from the solver : 1 if success (Rgoal is achieved)
	*/
	int SolveOCP(real const& continuationStep, real & Rdata, real const& Rgoal) const;

	/**
	* Solve OCP with discrete continuation on a real parameter from Rdata to Rgoal
	* @param continuationStep step of continuation as a ratio : 0 for no continuation, 0 < continuationStep <= 1 for continuation
	* @param Rdata string data for continuation
	* @param Rgoal real number to achieve
	* @return result from the solver : 1 if success (Rgoal is achieved)
	*/
	int SolveOCP(real const& continuationStep, std::string const& Rdata, real const& Rgoal) const {
		if (myModel.parameters.find(Rdata) != myModel.parameters.end()) {
			return SolveOCP(continuationStep, myModel.parameters[Rdata], Rgoal);
		}
		else {
			std::cout << std::endl << std::endl << "Data " << Rdata << " does not exist!" << std::endl << std::endl;
		}
		
	};

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
	* @param tf the final time
	* @return the final state
	*/
	model::mstate Move(real const& tf) const;

	/**
	* Function to move the model to tf using initial or previously solved solution and handling multiple shooting
	* @param tf the final time
	* @param Xf the final state
	*/
	void Move(real const& tf, model::mstate & Xf) const;

	/**
	* Set relative tolerance for the solver
	* @param xtol realtive tolerance
	*/
	void SetPrecision(real const& xtol) const;

	/**
	* Set minimum continuation step
	* @param step minimum continuation step
	*/
	void SetContinuationMinStep(real const& step) const;

	/**
	* Get current kth parameter
	* @param k the parameter number to be returned
	* @return the kth parameter
	*/
	real GetParameters(int const& k) const;

	/**
	* Get current parameters
	* @return a pointer to the array of parameters 
	*/
	real *GetParameters() const;

	/**
	* Get current parameters
	* @param paramVector a vector of parameters
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

	/**
	* Get Model
	* @return the model
	*/
	model & GetModel() const;

private:

	/**
	* Model of the Optimal Control Problem
	*/
	model & myModel;

	/**
	* Shooting data structure
	*/
	struct data_struct;

	/**
	* Shooting data pointer
	*/
	data_struct *data;

	/**
	* Static function for solving OCP with timeout
	* @param arg a pointer to the arguments passed to the thread
	*/
	static void TimeoutOCPFunction(void *arg);

	/**
	* Solve Shooting without continuation
	* @return result from the solver : 1 if success
	*/
	int SolveShooting() const;

	/**
	* Solve Shooting with discrete time and state continuation
	* @param continuationStep step of continuation as a ratio : 0 for no continuation, 0 < continuationStep <= 1 for continuation
	* @return result from the solver : 1 if success
	*/
	int SolveShootingContinuation(real const& continuationStep) const;

	/**
	* Solve Shooting with discrete continuation on parameter Rdata to Rgoal
	* @param continuationStep step of continuation as a ratio : 0 for no continuation, 0 < continuationStep <= 1 for continuation
	* @param Rdata real data for continuation
	* @param Rgoal real number to achieve
	* @return result from the solver : 1 if success (Rgoal is achieved)
	*/
	int SolveShootingContinuation(real const& continuationStep, real & Rdata, real const& Rgoal) const; 

	/**
	* Solve Function
	* @param numParam number of parameters
	* @param param a vector to the parameters
	* @return result from the solver : 1 if success
	*/
	int SolveShootingFunction(int const & numParam, std::vector<real> & param) const;

	/**
	* Static function for the solver
	* @param userdata a pointer to user data
	* @param n the number of parameters
	* @param param a pointer to the parameters
	* @param fvec a pointer to the function values
	* @param iflag an optional flag
	*/
	static int StaticShootingFunction(void *userdata, int n, const real *param, real *fvec, int iflag);	
	
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
	void ShootingFunctionParallel(int n, std::vector<real> const& param, std::vector<real> & fvec, int iflag) const;
	
	/**
	* Static function for shooting with multithreading
	* @param arg a pointer to the arguments passed to the thread
	*/
	static void StaticShootingParallelFunction(void *arg);
	
	/**
	* Compute a part of the shooting function on a single thread
	* @param n the number of parameters
	* @param param a pointer to the parameters
	* @param fvec a pointer to the function values
	* @param iflag a pointer to an optional flag
	* @param threadNum thread number
	* @param timeLine time vector considered by the thread
	*/
	void ShootingParallelFunctionThread(int n, real *param, real *fvec, int threadNum, std::vector<real> const& timeLine) const;
	
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
