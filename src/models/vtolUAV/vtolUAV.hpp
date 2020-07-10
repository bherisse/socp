/*
 * vtolUAV.hpp 
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */

#include "../../socp/model.hpp"
#include "../../socp/map.hpp"

#include <iostream>

#ifndef _VTOLUAV_H_
#define _VTOLUAV_H_

/**************			vtolUAV class			******************************/
class vtolUAV:public model
{

public:

	/**
	* Vehicle parameters
	*/
	struct parameters_struct{
		real u_max;				///< max normalized control
		real a_max;				///< max acceleration
		real alphaT;			///< weight for time cost
		real alphaV;			///< weight for Vd
		real invSigmaXwp;		///< weight for cost at intermediate points 
		real Vd;				///< desired velocity
		real ca;				///< drag coefficient
	};

	/**
	* Constructor
	*/
	vtolUAV(map & the_map, std::string the_fileTrace = std::string(""));

	/**
	* Destructor
	*/
	virtual ~vtolUAV();

	/**
	* Get data
	* @return parameters structure
	*/
	parameters_struct & GetParameterData();

	/**
	* Final function value for the considered Optimal Control Problem
	*/
	virtual void FinalFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const;

	/**
	* Function value for the considered control problem with free final time
	*/
	virtual void FinalHFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const;

	/**
	* Update Switching times
	*/
	virtual void SwitchingTimesUpdate(std::vector<real> const& switchingTimes);

	/**
	* Switching function
	*/
	virtual void SwitchingStateFunction(real const& t, int const& stateID, mstate const& X, mstate const& Xp, mstate const& Xd, mstate & fvec) const;

	/**
	* Get Model
	* @return the map
	*/
	map & GetMap() const;

private:

	/**
	* Map
	*/
	map & myMap;

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

};

#endif //_VTOLUAV_H_

