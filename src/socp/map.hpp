/*
 * map.hpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */
#include "commonType.hpp"

#ifndef _MAP_H_
#define _MAP_H_

/**************			map class			******************************/
class map
{

public:

	/**
	* Constructor
	*/
	map() {};

	/**
	* Destructor
	*/
	virtual ~map() {};

	/**
	* Obstacle function
	*/
	virtual void ComputeObstacles(real *position_) const = 0;

	/**
	* Obstacle function
	*/
	virtual void ComputeFuncObstacles(real *position_, real *funcObs) const = 0;

	/**
	* Obstacle gradient
	*/
	virtual void ComputeGradObstacles(real *position_, real *gradObs) const = 0;

	/**
	* Obstacle function
	*/
	virtual real GetObsFunc() const = 0;

	/**
	* Obstacle gradient
	*/
	virtual real *GetObsGrad() const = 0;

private:


};

#endif //_MAP_H_

