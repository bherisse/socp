/*
 * map.hpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */
#include <vector>

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
	* Penalization function
	* @param position the position to be considered
	* @param penFunc the value of the penalization function
	*/
	virtual void PenalizationFunction(std::vector<real> const& position, real & penFunc) const = 0;

	/**
	* Penalization gradient
	* @param position the position to be considered
	* @param penGrad the value of the penalization gradient
	*/
	virtual void PenalizationGradient(std::vector<real> const& position, std::vector<real> & penGrad) const = 0;

private:


};

#endif //_MAP_H_

