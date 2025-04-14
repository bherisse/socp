/*
 * model.hpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */
#include <fstream>
#include <map>

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
	* Enumerate optimization mode for times and components of the state vector
	*/
	enum {
		FIXED,		///< fixed time or state 
		FREE,		///< free time or state (constraints and transversality conditions are defined in user model functions)
		CONTINUOUS	///< use standard continuity conditions when doing multiple shooting (for interior points only)
	};

	/**
	* Constructor
	* @param stateDim the state dimension
	*/
	model(int const& _stateDim, int _modelOrder = 0, int _stepNbr = 10, std::string _fileTrace = std::string("")) :
		dim(_stateDim),
		modelOrder(_modelOrder),
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
	* @param isJac a flag (0 for state only, 1 for state+jacobian)
	* @return the final state and costate
	*/
	virtual mstate ComputeTraj(real const& t0, mstate const& X0, real const& tf, int isTrace, int isJac){
		return ModelInt(t0, X0, tf, isTrace, isJac);
	};

	/**
	* Final function value for the considered Optimal Control Problem
	* @param tf the final time
	* @param X_tf the final state
	* @param Xf the final desired state
	* @param mode_X the mode (free/fixed state)
	* @param fvec a vector of the values of the function
	* @param isJac a flag (0 for state only, 1 for state+jacobian)
	*/
	virtual void FinalFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec, int isJac) const{
		if (isJac == 0) {
			// Compute the function to solve
			for (int j = 0;j < dim;j++) {
				if (mode_X[j] == FREE) {
					// continuity => X(j+n) = 0 (transversality condition)
					fvec[j] = X_tf[j + dim];
				}
				else {
					// continuity => X(j) = Xf(j)
					fvec[j] = X_tf[j] - Xf[j];
				}
			}
		}
		else {
			// Compute the jacobian of the function to solve
			for (int j = 0;j < dim;j++) {
				if (mode_X[j] == FREE) {
					// continuity => X(j+n) = 0 (transversality condition)
					for (int i = 0; i < 2 * dim; i++) {
						fvec[2 * dim * j + i] = X_tf[2 * dim * (j + dim + 1) + i];
					}
				}
				else {
					// continuity => X(j) = Xf(j)
					for (int i = 0; i < 2 * dim; i++) {
						fvec[2 * dim * j + i] = X_tf[2 * dim * (j + 1) + i];
					}
				}
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
	* @param isJac a flag (0 for state only, 1 for state+jacobian)
	*/
	virtual void FinalHFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec, int isJac) const{
		if (isJac == 0) {
			// Compute the function to solve
			for (int j = 0;j < dim;j++) {
				if (mode_X[j] == FREE) {
					// continuity => X(j+n) = 0 (transversality condition)
					fvec[j] = X_tf[j + dim];
				}
				else {
					// continuity => X(j) = Xf(j)
					fvec[j] = X_tf[j] - Xf[j];
				}
			}

			fvec[dim] = Hamiltonian(tf, X_tf, 0)[0];
		}
		else {
			// Compute the jacobian of the function to solve
			mstate State_tf = { X_tf.begin(), X_tf.begin() + 2 * dim };
			mstate fxtf = Model(tf, State_tf, 0);
			for (int j = 0;j < dim;j++) {
				if (mode_X[j] == FREE) {
					// continuity => X(j+n) = 0 (transversality condition)
					for (int i = 0; i < 2 * dim; i++) {
						fvec[(2 * dim + 1) * j + i] = X_tf[2 * dim * (j + dim + 1) + i];
					}
					fvec[(2 * dim + 1) * j + 2 * dim] = fxtf[j + dim];
				}
				else {
					// continuity => X(j) = Xf(j)
					for (int i = 0; i < 2 * dim; i++) {
						fvec[(2 * dim + 1) * j + i] = X_tf[2 * dim * (j + 1) + i];
					}
					fvec[(2 * dim + 1) * j + 2 * dim] = fxtf[j];
				}
			}

			mstate dH_dX = Hamiltonian(tf, State_tf, 1);
			for (int i = 0; i < 2 * dim; i++) {
				fvec[(2 * dim + 1) * dim + i] = 0;
				for (int k = 0; k < 2 * dim; k++) {
					fvec[(2 * dim + 1) * dim + i] += dH_dX[k] * X_tf[2 * dim * (k + 1) + i];
				}
			}
			fvec[(2 * dim + 1) * dim + 2 * dim] = 0;
			for (int k = 0; k < 2 * dim; k++) {
				fvec[(2 * dim + 1) * dim + 2 * dim] += dH_dX[k] * fxtf[k];
			}
			fvec[(2 * dim + 1) * dim + 2 * dim] += dH_dX[2 * dim];
			
		}
		
	};

	/**
	* Final function value for the considered Optimal Control Problem
	* @param t0 the final time
	* @param X_t0 the final state
	* @param X0 the final desired state
	* @param mode_X the mode (free/fixed state)
	* @param fvec a vector of the values of the function
	* @param isJac a flag (0 for state only, 1 for state+jacobian)
	*/
	virtual void InitialFunction(real const& t0, mstate const& X_t0, mstate const& X0, std::vector<int> const& mode_X, std::vector<real> & fvec, int isJac) const{
		if (isJac == 0) {
			// Compute the function to solve
			for (int j = 0;j < dim;j++) {
				if (mode_X[j] == FREE) {
					// continuity => X(j+n) = 0 (transversality condition)
					fvec[j] = X_t0[j + dim];
				}
				else {
					// continuity => X(j) = Xf(j)
					fvec[j] = X_t0[j] - X0[j];
				}
			}
		}
		else {
			// Compute the jacobian of the function to solve
			for (int j = 0;j < dim;j++) {
				if (mode_X[j] == FREE) {
					// continuity => X(j+n) = 0 (transversality condition)
					for (int i = 0; i < 2 * dim; i++) {
						fvec[2 * dim  * j + i] = X_t0[2 * dim * (j + dim + 1) + i];
					}
				}
				else {
					// continuity => X(j) = Xf(j)
					for (int i = 0; i < 2 * dim; i++) {
						fvec[2 * dim * j + i] = X_t0[2 * dim * (j + 1) + i];
					}
				}
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
	* @param isJac a flag (0 for state only, 1 for state+jacobian)
	*/
	virtual void InitialHFunction(real const& t0, mstate const& X_t0, mstate const& X0, std::vector<int> const& mode_X, std::vector<real> & fvec, int isJac) const{
		if (isJac == 0) {
			// Compute the function to solve
			for (int j = 0;j < dim;j++) {
				if (mode_X[j] == FREE) {
					// continuity => X(j+n) = 0 (transversality condition)
					fvec[j] = X_t0[j + dim];
				}
				else {
					// continuity => X(j) = Xf(j)
					fvec[j] = X_t0[j] - X0[j];
				}
			}

			fvec[dim] = Hamiltonian(t0, X_t0, 0)[0];
		}
		else {
			// Compute the jacobian of the function to solve
			mstate State_t0 = { X_t0.begin(), X_t0.begin() + 2 * dim };
			mstate fxt0 = Model(t0, State_t0, 0);
			for (int j = 0;j < dim;j++) {
				if (mode_X[j] == FREE) {
					// continuity => X(j+n) = 0 (transversality condition)
					for (int i = 0; i < 2 * dim; i++) {
						fvec[(2 * dim + 1) * j + i] = X_t0[2 * dim * (j + dim + 1) + i];
					}
					fvec[(2 * dim + 1) * j + 2 * dim] = fxt0[j + dim];
				}
				else {
					// continuity => X(j) = Xf(j)
					for (int i = 0; i < 2 * dim; i++) {
						fvec[(2 * dim + 1) * j + i] = X_t0[2 * dim * (j + 1) + i];
					}
					fvec[(2 * dim + 1) * j + 2 * dim] = fxt0[j];
				}
			}

			mstate dH_dX = Hamiltonian(t0, State_t0, 1);
			for (int i = 0; i < 2 * dim; i++) {
				fvec[(2 * dim + 1) * dim + i] = 0;
				for (int k = 0; k < 2 * dim; k++) {
					fvec[(2 * dim + 1) * dim + i] += dH_dX[k] * X_t0[2 * dim * (k + 1) + i];
				}
			}
			fvec[(2 * dim + 1) * dim + 2 * dim] = 0;
			for (int k = 0; k < 2 * dim; k++) {
				fvec[(2 * dim + 1) * dim + 2 * dim] += dH_dX[k] * fxt0[k];
			}
			fvec[(2 * dim + 1) * dim + 2 * dim] += dH_dX[2 * dim];

		}
	};

	/**
	* Switching function 
	* @param t the free time
	* @param X the state at t
	* @param fvec value of the function
	* @param isJac a flag (0 for state only, 1 for state+jacobian)
	*/
	virtual mstate SwitchingTimesFunction(real const& t, mstate const& X, mstate const& Xp, int isJac) const{
		// function to implement if mode_t = 1
		if (isJac == 0) {
			mstate fvec(1, 0);
			fvec[0] = Hamiltonian(t, X, 0)[0] - Hamiltonian(t, Xp, 0)[0];
			return fvec;		// by default
		}
		else {
			// Compute the jacobian of the function to solve
			mstate fvec(4 * dim + 1, 0);
			mstate State_t = { X.begin(), X.begin() + 2 * dim };
			mstate fxt = Model(t, State_t, 0);
			mstate dH_dX = Hamiltonian(t, State_t, 1);
			mstate State_p = { Xp.begin(), Xp.begin() + 2 * dim };
			mstate fxp = Model(t, State_p, 0);
			mstate dHp_dX = Hamiltonian(t, State_p, 1);
			for (int i = 0; i < 2 * dim; i++) {
				for (int k = 0; k < 2 * dim; k++) {
					fvec[i] += dH_dX[k] * X[2 * dim * (k + 1) + i];
					fvec[2 * dim + i] -= dHp_dX[k] * Xp[2 * dim * (k + 1) + i];
				}
			}
			for (int k = 0; k < 2 * dim; k++) {
				fvec[4 * dim] += dH_dX[k] * fxt[k] - dHp_dX[k] * fxp[k];		
			}
			fvec[4 * dim] += dH_dX[2 * dim] - dHp_dX[2 * dim];
			return fvec;
		}
		
	};

	/**
	* Switching function 
	* @param t the free time
	* @param stateID the index of the state component
	* @param X the state at t
	* @param Xp the state at t+
	* @param fvec a vector of the values of the function
	* @param isJac a flag (0 for state only, 1 for state+jacobian)
	*/
	virtual void SwitchingStateFunction(real const& t, int const& stateID, mstate const& X, mstate const& Xp, mstate const& Xd, mstate & fvec, int isJac) const{
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

	int modelOrder;								///< 0 if jacobian is not provided, 1 otherwise

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
	* @param isJac a flag (0 for state only, 1 for state+jacobian)
	* @return the Hamiltonian value
	*/
	virtual mstate Hamiltonian(real const& t, mstate const& X, int isJac) const = 0;

	/**
	* Integrate state equations
	* @param X the state
	* @param t0 the initial time
	* @param tf the final time
	* @param isTrace trace flag (0 if no trace required)
	* @param isJac a flag (0 for state only, 1 for state+jacobian)
	* @return the state at tf
	*/
	virtual mstate ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace, int isJac) {
		std::stringstream ss;

		real dt = (tf - t0) / stepNbr;		// time step
		mstate Xs = X;

		if (isTrace) {
			integrate(modelStruct(this, isJac), Xs, t0, tf, dt, observerStruct(this,ss));
			// write in trace file
			std::ofstream fileTrace;
			fileTrace.open(strFileTrace.c_str(), std::ios::app);
			fileTrace << ss.str();
			fileTrace.close();
		}
		else {
			integrate(modelStruct(this, isJac), Xs, t0, tf, dt);
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
		real H = Hamiltonian(t, X, 0)[0];

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
		real H = Hamiltonian(t, X, 0)[0];

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

