%module(directors="1", allprotected="1") pyGoddard
%include "std_string.i"
%include "std_vector.i"
%{

#include "socp/commonType.hpp"
#include "socp/map.hpp"
#include "socp/model.hpp"

#include "models/goddard/goddard.hpp"

%}

%include "socp/commonType.hpp"
%include "socp/map.hpp"
%include "socp/model.hpp"

%include "models/goddard/goddard.hpp"
