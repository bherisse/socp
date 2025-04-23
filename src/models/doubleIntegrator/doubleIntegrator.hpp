/*
 * doubleIntegrator.hpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE
 */

#include "../../socp/model.hpp"

#include <iostream>

#ifndef _DOUBLEINTEGRATOR_H_
#define _DOUBLEINTEGRATOR_H_

/**************			doubleIntegrator class			******************************/
class doubleIntegrator : public model
{

public:

	/**
	* Vehicle parameters
	*/
	struct parameters_struct{
		real u_max;				// max normalized control
		real a_max;				// max acceleration
		real muT;				// weight for time cost
	};

	/**
	* Constructor
	*/
	doubleIntegrator(int modelOrder, std::string the_fileTrace);

	struct odeStruct;
	odeStruct *my_odeStruct;

	/**
	* Destructor
	*/
	virtual ~doubleIntegrator();

	/**
	* Get data pointer
	*/
	parameters_struct & GetParameterData();

	/**
	* Set number of integration steps
	*/
	void SetStepNumber(int step);

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
	mstate ModelState(real const& t, mstate const& X) const;
	mstate ModelJacobian(real const& t, mstate const& X) const;

	/**
	* Control model of the vehicle
	*/
	virtual mcontrol Control(real const& t, mstate const& X) const;

	/**
	* Hamiltonian of the vehicle
	*/
	virtual mstate Hamiltonian(real const& t, mstate const& X, int isJac) const;

};

#endif //_DOUBLEINTEGRATOR_H_

