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
	* @param state the state to be considered
	* @param func the value of the function
	*/
	virtual void Function(std::vector<real> const& state, real & func) const = 0;

	/**
	* Penalization gradient
	* @param state the state to be considered
	* @param grad the value of the gradient
	*/
	virtual void Gradient(std::vector<real> const& state, std::vector<real> & grad) const = 0;

private:


};

#endif //_MAP_H_

