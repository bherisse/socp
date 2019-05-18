/*
 * interceptor.hpp
 *
 *  Created on: June 30, 2016
 *      Author: Riccardo BONALLI & Bruno HERISSE (ONERA/DTIS)
 */
/*
 * The algorithm is presented in the paper "Analytical Initialization of a Continuation-Based Indirect
Method for Optimal Control of Endo-Atmospheric Launch Vehicle Systems", R. Bonalli, B. Hérissé, E. Trélat (IFAC WC 2017)
 */

#include "../../socp/model.hpp"
#include "../../socp/map.hpp"

#include <iostream>

#ifndef _INTERCEPTOR_H_
#define _INTERCEPTOR_H_

/**************			interceptor class			******************************/
class interceptor:public model
{

public:

	/**
	* Vehicle parameters
	*/
	struct parameters_struct{
		real c0;				///< max curvature at ground level (1/m)
		real hr;				///< reference altitude (m)
		real d0;				///< drag at ground level (1/m)	
		real eta;				///< coeff of efficiency
		real propellant_mass;	///< propellant mass (kg)
		real empty_mass; 		///< empty mass (kg)
		real q; 				///< mass flow rate (kg/s)
		real ve; 				///< gas speed (m/s)
		real alpha_max; 		///< max angle of attack (rd)
		real u_max;				///< max of the normalized control (for saturation)
		real a_max;				///< max acceleration allowed (for saturation)
		real r_2p;				///< ratio of used gas for the second propulsion phase
		real t_2p;				///< start time for the second propulsion phase
		real mu_gft;			///< parameter for considering gravity and propulsion
		real muT;				///< weight for time cost
		real muV;				///< weight for velocity cost
		real muC;				///< weight for quadratic control cost
	};

	/**
	* Constructor
	*/
	interceptor(std::string the_fileTrace = std::string(""));

	/**
	* Destructor
	*/
	virtual ~interceptor();

	/**
	* Get data
	* @return parameters structure
	*/
	parameters_struct & GetParameterData();

	/**
	* Compute trajectory from t0 to tf with initial state X0
	*/
	virtual mstate ComputeTraj(real const& t0, mstate const& X0, real const& tf, int isTrace);

	/**
	* Function value for the considered control problem
	*/
	virtual void FinalFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const;
	virtual void FinalHFunction(real const& tf, mstate const& X_tf, mstate const& Xf, std::vector<int> const& mode_X, std::vector<real> & fvec) const;

	/**
	* Guess from analytical solution
	* @param ti the initial time
	* @param Xi the initial state
	* @param tf the final time
	* @param Xf the final state
	*/
	void InitAnalytical(real const& ti, mstate & Xi, real const& tf, mstate & Xf) const;

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
	static mstate static_Model(real const& t, interceptor::mstate const& X, void *model);
	mstate Model_1(real const& t, mstate const& X) const;
	mstate Model_2(real const& t, mstate const& X) const;

	/**
	* Control model of the vehicle
	*/
	virtual mcontrol Control(real const& t, mstate const& X) const;
	mcontrol Control_1(real const& t, mstate const& X) const;
	mcontrol Control_2(real const& t, mstate const& X) const;

	/**
	* Hamiltonian of the vehicle
	*/
	virtual real Hamiltonian(real const& t, mstate const& X) const;
	real Hamiltonian_1(real const& t, mstate const& X) const;
	real Hamiltonian_2(real const& t, mstate const& X) const;

	/**
	* Integrate state equations with a RK4
	*/
	virtual mstate ModelInt(real const& t0, mstate const& X, real const& tf, int isTrace);

	/**
	* Trace state
	*/
	virtual void Trace(real const& t, mstate const& X, std::stringstream & file) const;

	/**
	* Get Vehicle mode (for hybrid systems)
	*/
	virtual int GetMode(real const& t, mstate const& X) const;

	/**
	* Switch charts
	*/
	mstate SetChart(real const& t, mstate const& X);
	mstate  ConversionState12(mstate const& X1) const;
	mstate  ConversionState21(mstate const& X2) const;

	/**
	* Compute Mass
	*/
	real ComputeMass(real const& t, mstate const& X) const;

};

#endif //_INTERCEPTOR_H_

