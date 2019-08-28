/*
 * model.hpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */
#include <fstream>
#include <sstream>

#include "commonType.hpp"

#include "odeTools.hpp"

// if boost installed
#ifdef _USE_BOOST
	#include "boost/numeric/odeint.hpp"
	//using namespace boost::numeric::odeint;
	typedef boost::numeric::odeint::runge_kutta_dopri5< odeTools::odeVector > dopri_stepper_type;
#endif

#ifndef _MODEL_H_
#define _MODEL_H_

/**************			model class (abstract)		******************************/
class model
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
	*/
	model(){};

	/**
	* Constructor
	* @param stateDim the state dimension
	*/
	model(int const& stateDim) {dim = stateDim; odeIntTol = 1e-8;};

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
	virtual void SwitchingTimesFunction(real const& t, mstate const& X, real & fvec) const{
		// function to implement if mode_t = 1
		fvec = Hamiltonian(t,X);		// by default
	};

	/**
	* Switching function 
	* @param t the free time
	* @param stateID the index of the state component
	* @param X the state at t
	* @param Xp the state at t+
	* @param fvec a vector of the values of the function
	*/
	virtual void SwitchingStateFunction(real const& t, int const& stateID, mstate const& X, mstate const& Xp, mstate & fvec) const{
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
	
protected:

	int dim;							///< state dimension

	real odeIntTol;						///< precision if dopri5 is used

	/**
	* Model structure
	*/
	struct modelStruct : public odeTools::odeStruct
	{
		model* m_model;

		modelStruct( model* model ) : m_model( model ) { }

		virtual void operator()( model::mstate const& X , model::mstate& dXdt , real const& t ) const
		{
			dXdt = m_model->Model(t, X);
		}
	};

	/**
	* Observer structure
	*/
	struct observerStruct
	{
		model* m_model;
		std::stringstream & m_file;

		observerStruct(model* const model, std::stringstream & file) : m_model( model ), m_file( file ) { }

		virtual void operator()(model::mstate const& X, double const t) const
		{
			 m_model->Trace(t, X, m_file);
		}
	};

	/**
	* State model of the vehicle : dX/dt = Model(t, X)
	* @param t the time
	* @param X the state
	* @return dX/dt as a model state
	*/
	virtual mstate Model(real const& t, mstate const& X) const = 0;

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
	virtual mstate ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace) = 0;

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
	
	/**
	* Integrate model from t0 to tf with an observer
	* @param model structure of the model
	* @param X the state
	* @param t0 initial time
	* @param tf final time
	* @param dt the step of integration
	* @param observer structure of the observer
	*/
	static void integrate(modelStruct const& model, mstate & X, double const& t0, double const& tf, double const& dt, observerStruct const& observer){
		#ifdef _USE_BOOST
			// if boost used, use dopri 5
			//boost::numeric::odeint::integrate(model, X, t0, tf, dt, observer);
			boost::numeric::odeint::integrate_const( boost::numeric::odeint::make_dense_output< dopri_stepper_type >(model.m_model->odeIntTol, model.m_model->odeIntTol), model, X, t0, tf, dt, observer);
		#else
			real t = t0;							// time
			observer(X,t);
			while (t < (tf-dt/2)) {
				if (t + dt > tf) {
					odeTools::RK4(t, X, tf-t, model);
				}
				else {
					odeTools::RK4(t, X, dt, model);
				}
				t += dt;
				observer(X,t);	
			}
		#endif
	};

	/**
	* Integrate model from t0 to tf
	* @param model structure of the model
	* @param X the state
	* @param t0 initial time
	* @param tf final time
	* @param dt the step of integration
	*/
	static void integrate(modelStruct const& model, mstate & X, double const& t0, double const& tf, double const& dt){
		#ifdef _USE_BOOST
			// if boost used, use dopri 5
			//boost::numeric::odeint::integrate(model, X, t0, tf, dt);
			boost::numeric::odeint::integrate_adaptive( boost::numeric::odeint::make_dense_output< dopri_stepper_type >(model.m_model->odeIntTol, model.m_model->odeIntTol), model, X, t0, tf, dt);
		#else
			real t = t0;			// time
			while (t < (tf-dt/2)) {
				if (t+dt > tf){
					odeTools::RK4(t, X, tf-t, model);
				}else{
					odeTools::RK4(t, X, dt, model);
				}				
				t += dt;
			}
		#endif
	};

};

#endif //_MODEL_H_

