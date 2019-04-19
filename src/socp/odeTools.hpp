/*
 * odeTools.hpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */

#include "commonType.hpp"

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
	* ODE structure 
	*/
	struct odeStruct
	{
		odeStruct( ) { }

		virtual void operator()( odeVector const& X , odeVector& dXdt , real const& t ) const {	}
	};

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
	static void RK1(real const& t, odeVector & X, real const& step, odeStruct const& ode);

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
	static void RK2(real const& t, odeVector & X, real const& step, odeStruct const& ode);

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
	static void RK4(real const& t, odeVector & X, real const& step, odeStruct const& ode);

private:

};

#endif //_ODETOOLS_H_

