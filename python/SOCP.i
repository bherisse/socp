%module(directors="1", allprotected="1") pySOCP
%include "std_string.i"
%include "std_vector.i"
%include "carrays.i"
%{

#include "socp/commonType.hpp"
#include "socp/shooting.hpp"
#include "socp/odeTools.hpp"
#include "socp/map.hpp"
#include "socp/model.hpp"

%}

%include "socp/commonType.hpp"
%include "socp/shooting.hpp"
%include "socp/odeTools.hpp"
%include "socp/map.hpp"
%include "socp/model.hpp"

%array_class(double, doubleArray);

namespace std {
   %template(iVector) vector<int>;
   %template(dVector) vector<double>;
   %template(iMatrix) vector<vector<int>>;
   %template(dMatrix) vector<vector<double>>;
};