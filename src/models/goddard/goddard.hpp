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
		real C = 3.5;					///< coefficient for thrust
		real b = 7.0;					///< coefficient for mass flow rate
		real KD = 310.0;				///< coefficient for drag
		real kr = 500.0;				///< coefficient for density of air
		real u_max = 1.0;				///< max normalized control
		real mu1 = 1.0;					///< weight for the cost on the norm of the control in [0,1]
		real mu2 = 0.0;					///< weight for the quadratic cost on the control in [0,1]
		real singularControl = -1;		///< singular control value
	};

	/**
	* Constructor
	*/
	goddard(std::string the_fileTrace = std::string(""), int stepNbr = 10);

	/**
	* Destructor
	*/
	virtual ~goddard();

	/**
	* Switching function 
	*/
	virtual mstate SwitchingTimesFunction(real const& t, mstate const& X, mstate const& Xp, int isJac) const;

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
	virtual mstate Model(real const& t, mstate const& X, int isJac) const;

	/**
	* Control model of the vehicle
	*/
	virtual mcontrol Control(real const& t, mstate const& X) const;

	/**
	* Hamiltonian of the vehicle
	*/
	virtual mstate Hamiltonian(real const& t, mstate const& X, int isJac) const;

	/**
	* Integrate state equations with a RK4
	*/
	virtual mstate ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace, int isJac);

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

