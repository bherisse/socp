/*
 * goddard.hpp
 *
 *  Created on: February 02, 2017
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */
/*
 * The algorithm is presented in the paper "Singular Arcs in the Generalized Goddard's Problem", F. Bonnans, P. Martinon, E. Trélat (J Optim Theory Appl (2008) 139: 439-461)
 */

#include "../../socp/model.hpp"
#include "../../socp/map.hpp"

#include <iostream>

#ifndef _GODDARD_H_
#define _GODDARD_H_

/**************			goddard class			******************************/
class goddard:public model
{

public:

	/**
	* Vehicle parameters
	*/
	struct parameters_struct{
		real C;					///< coefficient for thrust
		real b;					///< coefficient for mass flow rate
		real KD;				///< coefficient for drag
		real kr;				///< coefficient for density of air
		real u_max;				///< max normalized control
		real mu1;				///< weight for the cost on the norm of the control in [0,1]
		real mu2;				///< weight for the quadratic cost on the control in [0,1]
		real singularControl;	///< singular control value
	};

	/**
	* Constructor
	*/
	goddard(std::string the_fileTrace = std::string(""));

	/**
	* Destructor
	*/
	virtual ~goddard();

	/**
	* Get data
	* @return parameters structure
	*/
	parameters_struct & GetParameterData();

	/**
	* Switching function 
	*/
	virtual void SwitchingTimesFunction(real const& t, mstate const& X, real & fvec) const;

	/**
	* Update Switching times
	*/
	virtual void SwitchingTimesUpdate(std::vector<real> const& switchingTimes);

	/**
	* set parameter data by name
	* @param name the name of the parameter
	* @param value the value of the parameter
	*/
	void SetParameterDataName(std::string name, real value);

	/**
	* Get parameter data by name
	* @param name the name of the parameter
	* @return the value of the parameter
	*/
	real & GetParameterDataName(std::string name);

private:

	/**
	* Vehicle data
	*/
	struct data_struct;
	data_struct *data;

	/**
	* State model of the vehicle
	*/
	virtual mstate Model(real const& t, mstate const& X) const;

	/**
	* Control model of the vehicle
	*/
	virtual mcontrol Control(real const& t, mstate const& X) const;

	/**
	* Hamiltonian of the vehicle
	*/
	virtual real Hamiltonian(real const& t, mstate const& X) const;

	/**
	* Integrate state equations with a RK4
	*/
	virtual mstate ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace);

	/**
	* Trace state
	*/
	virtual void Trace(real const& t, mstate const& X, std::stringstream & file) const;

	/**
	* Get Singular control
	*/
	real GetSingularControl(real t, mstate const& X) const;

};

#endif //_GODDARD_H_

