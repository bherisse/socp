/*
 * model.hpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */
#include <fstream>
#include <map>

#include "commonType.hpp"

#include "odeTools.hpp"

#ifndef _MODEL_H_
#define _MODEL_H_

/**************			model class (abstract)		******************************/
class model : public odeTools
{

public:

	/**
	* vector for model state (model state and/or costate)
	*/
	typedef odeTools::odeVector mstate;

	/**
	* vector for model control
	*/
	typedef odeTools::odeVector mcontrol;

	/**
	* Constructor
	* @param stateDim the state dimension
	*/
	model(int const& stateDim, int _stepNbr = 10, std::string _fileTrace = std::string("")) :
		dim(stateDim),
		stepNbr(_stepNbr),
		strFileTrace(_fileTrace)
	{
		parameters.clear();
		// trace file
		std::ofstream fileTrace;
		fileTrace.open(strFileTrace.c_str(), std::ios::trunc);	// erase file
		fileTrace.close();
	};

	/**
	* Destructor
	*/
	virtual ~model(){};

	/**
	* Get state dimension
	* @return the dimension
	*/
	virtual int GetDim() const {return dim;};

	/**
	* Compute trajectory from t0 to tf with initial state X0
	* @param t0 the initial time
	* @param X0 the initial state
	* @param tf the final time
	* @param isTrace trace flag (0 if no trace required)
	* @return the final state and costate
	*/
	virtual mstate ComputeTraj(real const& t0, mstate const& X0, real const& tf, int isTrace){
		return ModelInt(t0, X0, tf, isTrace);
	};

	/**
	* Final function value for the considered Optimal Control Problem
	* @param tf the final time
	* @param X_tf the final state
	* @param Xf the final desired state
	* @param mode_X the mode (free/fixed state)
	* @param fvec a vector of the values of the function
	*/
	virtual void FinalFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const{	
		// Compute the function to solve
		for (int j=0;j<dim;j++){
			if (mode_X[j]==1){
				// continuity => X(k+n) = 0 (transversality condition)
				fvec[j] = X_tf[j+dim];
			}else{
				// continuity => X(k) = Xf(k)
				fvec[j] = X_tf[j] - Xf[j];
			}
		}
	};

	/**
	* Function value for the considered control problem with free final time
	* @param tf the final time
	* @param X_tf the final state
	* @param Xf the final desired state
	* @param mode_X the mode (free/fixed state)
	* @param fvec a vector of the values of the function
	*/
	virtual void FinalHFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const{	
		// Compute the function to solve
		for (int j=0;j<dim;j++){
			if (mode_X[j]==1){
				// continuity => X(k+n) = 0 (transversality condition)
				fvec[j] = X_tf[j+dim];
			}else{
				// continuity => X(k) = Xf(k)
				fvec[j] = X_tf[j] - Xf[j];
			}
		}
		fvec[dim] = Hamiltonian(tf,X_tf);
	};

	/**
	* Final function value for the considered Optimal Control Problem
	* @param t0 the final time
	* @param X_t0 the final state
	* @param X0 the final desired state
	* @param mode_X the mode (free/fixed state)
	* @param fvec a vector of the values of the function
	*/
	virtual void InitialFunction(real const& t0, mstate const& X_t0, mstate const& X0, std::vector<int> const& mode_X, std::vector<real> & fvec) const{	
		// Compute the function to solve
		for (int j=0;j<dim;j++){
			if (mode_X[j]==1){
				// continuity => X(k+n) = 0 (transversality condition)
				fvec[j] = X_t0[j+dim];
			}else{
				// continuity => X(k) = Xf(k)
				fvec[j] = X_t0[j] - X0[j];
			}
		}
	};

	/**
	* Function value for the considered control problem with free final time
	* @param t0 the final time
	* @param X_t0 the final state
	* @param X0 the final desired state
	* @param mode_X the mode (free/fixed state)
	* @param fvec a vector of the values of the function
	*/
	virtual void InitialHFunction(real const& t0, mstate const& X_t0, mstate const& X0, std::vector<int> const& mode_X, std::vector<real> & fvec) const{	
		// Compute the function to solve
		for (int j=0;j<dim;j++){
			if (mode_X[j]==1){
				// continuity => X(k+n) = 0 (transversality condition)
				fvec[j] = X_t0[j+dim];
			}else{
				// continuity => X(k) = Xf(k)
				fvec[j] = X_t0[j] - X0[j];
			}
		}
		fvec[dim] = Hamiltonian(t0,X_t0);
	};

	/**
	* Switching function 
	* @param t the free time
	* @param X the state at t
	* @param fvec value of the function
	*/
	virtual real SwitchingTimesFunction(real const& t, mstate const& X) const{
		// function to implement if mode_t = 1
		return Hamiltonian(t,X);		// by default
	};

	/**
	* Switching function 
	* @param t the free time
	* @param stateID the index of the state component
	* @param X the state at t
	* @param Xp the state at t+
	* @param fvec a vector of the values of the function
	*/
	virtual void SwitchingStateFunction(real const& t, int const& stateID, mstate const& X, mstate const& Xp, mstate const& Xd, mstate & fvec) const{
		// function to implement if mode_X = 1
	};

	/**
	* Update Switching times
	* @param switchingTimes vector of switching times
	*/
	virtual void SwitchingTimesUpdate(std::vector<real> const& switchingTimes){};

	/**
	* Set relative tolerance if dopri5 is used
	* @param xtol value of the precision
	*/
	virtual void SetODEIntPrecision(real const& xtol) {
		odeIntTol = xtol;
	};
	
//protected:

	int dim;									///< state dimension

	std::map<std::string,real> parameters;		///< parameters for continuation

	std::string strFileTrace;					///< trace file

	int  stepNbr;								///< step number for ModelInt

	/**
	* Control model of the vehicle : U = Control(t, X)
	* @param t the time
	* @param X the state
	* @return the control U as a control state
	*/
	virtual mcontrol Control(real const& t, mstate const& X) const = 0;

	/**
	* Hamiltonian of the vehicle
	* @param t the time
	* @param X the state
	* @return the Hamiltonian value
	*/
	virtual real Hamiltonian(real const& t, mstate const& X) const = 0;

	/**
	* Integrate state equations
	* @param X the state
	* @param t0 the initial time
	* @param tf the final time
	* @param isTrace trace flag (0 if no trace required)
	* @return the state at tf
	*/
	virtual mstate ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace) {
		std::stringstream ss;

		real dt = (tf - t0) / stepNbr;		// time step
		mstate Xs = X;

		if (isTrace) {
			odeTools::integrate(modelStruct(this), Xs, t0, tf, dt, observerStruct(this,ss));
			// write in trace file
			std::ofstream fileTrace;
			fileTrace.open(strFileTrace.c_str(), std::ios::app);
			fileTrace << ss.str();
			fileTrace.close();
		}
		else {
			odeTools::integrate(modelStruct(this), Xs, t0, tf, dt);
		}

		return Xs;
	}

	/**
	* Trace state
	* @param t the time
	* @param X the state
	* @param file the file (ofstream)
	*/
	virtual void Trace(real const& t, mstate const& X, std::ofstream & file) const{
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
		file << H << std::endl;
	}

	/**
	* Trace state
	* @param t the time
	* @param X the state
	* @param file the file (sstream)
	*/
	virtual void Trace(real const& t, mstate const& X, std::stringstream & file) const{
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
		file << H << std::endl;
	}

	/**
	* Get Vehicle mode (for hybrid systems)
	* @param t the time
	* @param X the state
	* @return the mode
	*/
	virtual int GetMode(real const& t, mstate const& X) const {return 0;};

private:

	/**
	* Default constructor
	*/
	model() {};

};

#endif //_MODEL_H_

