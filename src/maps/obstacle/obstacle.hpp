/*
 * obstacle.hpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE
 */
#include "../../socp/map.hpp"

#ifndef _OBSTACLE_H_
#define _OBSTACLE_H_

/**************			map class			******************************/
class obstacle:public map
{

public:

	/**
	* Map parameters
	*/
	struct parameters_struct{
		real phiObs;          	///< weight of obstacle function
		real psiWP;          	///< weight of waypoint function
		real muObs;  		///< parameter for penalization dispersion
		real sigmaWP;			///< parameter for penalization dispersion
	};

	/**
	* Constructor
	*/
	obstacle(std::string the_fileObstacles = std::string(""), std::string the_fileWP = std::string(""));

	/**
	* Destructor
	*/
	virtual ~obstacle();

	/**
	* Total function
	*/
	virtual void Function(std::vector<real> const& position, real & funcTot) const;

	/**
	* Total gradient
	*/
	virtual void Gradient(std::vector<real> const& position, std::vector<real> & gradTot) const;

	/**
	* Obstacle function
	*/
	void ObstaclePenalizationFunction(std::vector<real> const& position, real & funcObs) const;

	/**
	* Obstacle gradient
	*/
	void ObstaclePenalizationGradient(std::vector<real> const& position, std::vector<real> & gradObs) const;

	/**
	* WayPoints function
	*/
	void WPPenalizationFunction(std::vector<real> const& position, real & funcWP) const;

	/**
	* WayPoints gradient
	*/
	void WPPenalizationGradient(std::vector<real> const& position, std::vector<real> & gradWP) const;

	/**
	* Get data pointer
	*/
	parameters_struct & GetParameterData();

	/**
	* Get path as a vector
	*/
	const std::vector<std::vector<real>> & GetPath();

private:

	/**
	* Map data
	*/
	struct data_struct;
	data_struct *data;

	/**
	* Read parameters in obstacle file
	*/
	void ReadObstacleInput();

	/**
	* Read parameters in waypoint file
	*/
	void ReadWPInput();

};

#endif //_OBSTACLE_H_

