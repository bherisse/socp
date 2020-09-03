/*
 * odeTools.cpp
 *
 *  Created on: June 30, 2016
 *      Author: Bruno HERISSE (ONERA/DTIS)
 */

#include "odeTools.hpp"

// if boost installed
#ifdef _USE_BOOST
#include "boost/numeric/odeint.hpp"
 //using namespace boost::numeric::odeint;
typedef boost::numeric::odeint::runge_kutta_dopri5< odeTools::odeVector > dopri_stepper_type;
#endif

/**
* Multiply a odeVector X by a real
*/
odeTools::odeVector odeTools::MultState(real a, odeVector const& X){
	odeVector Y(X.size());

	for (int i = 0; i < X.size(); i++) {
		Y[i] = a*X[i];
	}

	return Y;
};

/**
* Add two vectors
*/
odeTools::odeVector odeTools::AddState(odeVector const& X, odeVector const& Y){
	odeVector Z(X.size());

	for (int i = 0; i < X.size(); i++) {
		Z[i] = X[i] + Y[i];
	}

    return Z;
};

/**
* compute one step RK1
*/
odeTools::odeVector odeTools::RK1(real const& t, odeVector const& X, real const& step, odeVector (*function)(real const&,odeVector const&, void*), void* context){
	return AddState(X,MultState(step,function(t, X, context)));
}

void odeTools::RK1(real const& t, odeVector & X, real const& step, modelStruct const& ode) {
	odeVector F1 = X;

	ode(X, F1, t);

	X = AddState(X, MultState(step, F1));
}

/**
* compute one step RK2
*/
odeTools::odeVector odeTools::RK2(real const& t, odeVector const& X, real const& step, odeVector (*function)(real const&,odeVector const&, void*), void* context){
	odeVector F1 = function(t, X, context);
	odeVector F2 = function(t + step/2.0, AddState(X,MultState(step/2.0,F1)), context);

	return AddState(X,MultState(step,F2));
}

void odeTools::RK2(real const& t, odeVector & X, real const& step, modelStruct const& ode) {
	odeVector F1 = X, F2 = X;

	ode(X, F1, t);
	ode(AddState(X, MultState(step / 2.0, F1)), F2, t + step / 2.0);

	X = AddState(X, MultState(step, F2));
}

/**
* compute one step RK4
*/
odeTools::odeVector odeTools::RK4(real const& t, odeVector const& X, real const& step, odeVector (*function)(real const&,odeVector const&, void*), void* context){
	odeVector F1 = function(t, X, context);
	odeVector F2 = function(t + step/2.0, AddState(X,MultState(step/2.0,F1)), context);
    odeVector F3 = function(t + step/2.0, AddState(X,MultState(step/2.0,F2)), context);
    odeVector F4 = function(t + step, AddState(X,MultState(step,F3)), context);	

	return AddState(X,MultState(step/6.0,AddState(F1,AddState(F4,MultState(2.0, AddState(F2,F3))))));
}

void odeTools::RK4(real const& t, odeVector & X, real const& step, modelStruct const& ode){
	odeVector F1=X, F2=X, F3=X, F4=X;
	
	ode(X, F1, t);
	ode(AddState(X,MultState(step/2.0,F1)), F2, t + step/2.0);
	ode(AddState(X,MultState(step/2.0,F2)), F3, t + step/2.0);
	ode(AddState(X,MultState(step,F3)), F4, t + step);		

	X = AddState(X,MultState(step/6.0,AddState(F1,AddState(F4,MultState(2.0, AddState(F2,F3))))));
}

/**
* Integrate model from t0 to tf with an observer
*/
void odeTools::integrate(modelStruct const& _model, odeVector & X, double const& t0, double const& tf, double const& dt, observerStruct const& _observer) {
#ifdef _USE_BOOST
	// if boost used, use dopri 5
	//boost::numeric::odeint::integrate(model, X, t0, tf, dt, observer);
	//boost::numeric::odeint::integrate_const(boost::numeric::odeint::make_dense_output< dopri_stepper_type >(_model.m_model->odeIntTol, _model.m_model->odeIntTol), _model, X, t0, tf, dt, _observer);
	boost::numeric::odeint::integrate_const(boost::numeric::odeint::make_dense_output< dopri_stepper_type >(odeIntTol, odeIntTol), _model, X, t0, tf, dt, _observer);
#else
	real t = t0;							// time
	_observer(X, t);
	while (t < (tf - dt / 2)) {
		if (t + dt > tf) {
			odeTools::RK4(t, X, tf - t, _model);
		}
		else {
			odeTools::RK4(t, X, dt, _model);
		}
		t += dt;
		_observer(X, t);
	}
#endif
};

/**
* Integrate model from t0 to tf
*/
void odeTools::integrate(modelStruct const& _model, odeVector & X, double const& t0, double const& tf, double const& dt) {
#ifdef _USE_BOOST
	// if boost used, use dopri 5
	//boost::numeric::odeint::integrate(model, X, t0, tf, dt);
	//boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_dense_output< dopri_stepper_type >(_model.m_model->odeIntTol, _model.m_model->odeIntTol), _model, X, t0, tf, dt);
	boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_dense_output< dopri_stepper_type >(odeIntTol, odeIntTol), _model, X, t0, tf, dt);
#else
	real t = t0;			// time
	while (t < (tf - dt / 2)) {
		if (t + dt > tf) {
			odeTools::RK4(t, X, tf - t, _model);
		}
		else {
			odeTools::RK4(t, X, dt, _model);
		}
		t += dt;
	}
#endif
};
