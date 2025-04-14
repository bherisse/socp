/*
* covid19.hpp
*
*  Created on: April 19, 2020
*      Author: Bruno HERISSE (ONERA/DTIS)
*/
/*
* The model used here is a SEIR model with estimated covid19 parameters
*/

#include "../../socp/model.hpp"
#include "../../socp/map.hpp"

#include <iostream>

#ifndef _COVID19_H_
#define _COVID19_H_

/**************			covid19 class			******************************/
class covid19:public model
{

public:

	/**
	* Vehicle parameters
	*/
	struct parameters_struct{
		real R0;					///< the number of secondary infections each infected individual produces
		real Tinf;					///< duration patient is infectious
		real Tinc;					///< incubation period
		real N;						///< size of population
		real Imax;					///< maximum Infectious permitted
		real muI;					///< penalization coefficient to constrain Infectious
		real umin;					///< minimum control
		real umax;					///< maximum control
	};

	/**
	* Constructor
	*/
	covid19(std::string the_fileTrace = std::string(""));

	/**
	* Destructor
	*/
	virtual ~covid19();

	/**
	* Get data
	* @return parameters structure
	*/
	parameters_struct & GetParameterData();

private:

	/**
	* Vehicle data
	*/
	struct data_struct;
	data_struct *data;

	/**
	* State model of the vehicle
	*/
	virtual mstate Model(real const& t, mstate const& X, int isJac = 0) const;

	/**
	* Control model of the vehicle
	*/
	virtual mcontrol Control(real const& t, mstate const& X) const;

	/**
	* Hamiltonian of the vehicle
	*/
	virtual mstate Hamiltonian(real const& t, mstate const& X, int isJac = 0) const;

	/**
	* Integrate state equations with a RK4
	*/
	virtual mstate ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace, int isJac = 0);

};

#endif //_COVID19_H_

