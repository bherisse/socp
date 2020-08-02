/*
 * odeTools.hpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */

#include "commonType.hpp"

#include <sstream>
#include <vector>

#ifndef _ODETOOLS_H_
#define _ODETOOLS_H_

/**************			odeTools class	******************************/
class odeTools
{

public:

	/**
	* odeVector
	*/
	typedef std::vector<real> odeVector;

	/**
	* Model structure
	*/
	struct modelStruct
	{
		odeTools* m_ode;

		modelStruct(odeTools* ode) : m_ode(ode) { }

		virtual void operator()(odeVector const& X, odeVector& dXdt, real const& t) const
		{
			dXdt = m_ode->Model(t, X);
		}
	};

	/**
	* Observer structure
	*/
	struct observerStruct
	{
		odeTools* m_ode;
		std::stringstream & m_file;

		observerStruct(odeTools* ode, std::stringstream & file) : m_ode(ode), m_file(file) { }

		virtual void operator()(odeVector const& X, double const t) const
		{
			m_ode->Trace(t, X, m_file);
		}
	};

	/**
	* Constructor
	*/
	odeTools() : odeIntTol(1e-8) {};


	real odeIntTol;								///< precision if dopri5 is used

	/**
	* Set relative tolerance if dopri5 is used
	* @param xtol value of the precision
	*/
	virtual void SetODEIntPrecision(real const& xtol) {
		odeIntTol = xtol;
	};

	/**
	* State model of the vehicle : dX/dt = Model(t, X)
	* @param t the time
	* @param X the state
	* @return dX/dt as a model state
	*/
	virtual odeVector Model(real const& t, odeVector const& X) const = 0;

	/**
	* Trace state
	* @param t the time
	* @param X the state
	* @param file the file (sstream)
	*/
	virtual void Trace(real const& t, odeVector const& X, std::stringstream & file) const = 0;

	/**
	* Multiply an odeVector X by a real
	* @param a a real
	* @param X a vector
	* @return a vector
	*/
	static odeVector MultState(real a, odeVector const& X);

	/**
	* Add two vectors
	* @param X a vector
	* @param Y a vector
	* @return a vector
	*/
	static odeVector AddState(odeVector const& X, odeVector const& Y);

	/**
	* compute one step RK1
	* @param t time
	* @param X an input vector
	* @param step a time step
	* @param function the function to integrate
	* @param context a user argument
	* @return a vector
	*/
	static odeVector RK1(real const& t, odeVector const& X, real const& step, odeVector (*function)(real const&,odeVector const&, void*), void* context);

	/**
	* compute one step RK1
	* @param t time
	* @param X an input vector
	* @param step a time step
	* @param ode structure of the ode
	*/
	static void RK1(real const& t, odeVector & X, real const& step, modelStruct const& ode);

	/**
	* compute one step RK2
	* @param t time
	* @param X an input vector
	* @param step a time step
	* @param function the function to integrate
	* @param context a user argument
	* @return a vector
	*/
	static odeVector RK2(real const& t, odeVector const& X, real const& step, odeVector (*function)(real const&,odeVector const&, void*), void* context);

	/**
	* compute one step RK2
	* @param t time
	* @param X an input vector
	* @param step a time step
	* @param ode structure of the ode
	*/
	static void RK2(real const& t, odeVector & X, real const& step, modelStruct const& ode);

	/**
	* compute one step RK4
	* @param t time
	* @param X an input vector
	* @param step a time step
	* @param function the function to integrate
	* @param context a user argument
	* @return a vector
	*/
	static odeVector RK4(real const& t, odeVector const& X, real const& step, odeVector (*function)(real const&,odeVector const&, void*), void* context);

	/**
	* compute one step RK4
	* @param t time
	* @param X an input vector
	* @param step a time step
	* @param ode structure of the ode
	*/
	static void RK4(real const& t, odeVector & X, real const& step, modelStruct const& ode);

	/**
	* Integrate model from t0 to tf with an observer
	* @param model structure of the model
	* @param X the state
	* @param t0 initial time
	* @param tf final time
	* @param dt the step of integration
	* @param observer structure of the observer
	*/
	static void integrate(modelStruct const& _model, odeVector & X, double const& t0, double const& tf, double const& dt, observerStruct const& _observer);

	/**
	* Integrate model from t0 to tf
	* @param model structure of the model
	* @param X the state
	* @param t0 initial time
	* @param tf final time
	* @param dt the step of integration
	*/
	static void integrate(modelStruct const& _model, odeVector & X, double const& t0, double const& tf, double const& dt);

private:

};

#endif //_ODETOOLS_H_

