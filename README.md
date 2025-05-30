[![Build Status](https://travis-ci.com/bherisse/socp.svg?branch=master)](https://travis-ci.com/bherisse/socp)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/bherisse/socp)](https://ci.appveyor.com/project/bherisse/socp)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

SOCP (Shooting for Optimal Control Problems)
===========
	
This software is a C++ library for optimal control using shooting methods.
	
Installing SOCP
----------

To install the software use cmake tool and build CMinPack first. In the SOCP folder, do 
- mkdir build
- cd build
- cmake ..
- make CMinPack
- cmake ..
- make

To use Boost library (ode integration with adaptive stepsize), do not forget to set Boost_USE to ON (adding cmake option -DBoost_USE:BOOL=ON). 

To build Python modules (Release mode only), add cmake option -DBUILD_PYTHON:BOOL=ON. 
Note that you need SWIG (minimum version: 3) and PYTHON (minimum version: 3) installed on your computer. 

Using SOCP
---------

The best way to start using SOCP is to have a look at the example provided in the "tests" folder. 

To build Goddard model, add cmake option -DBUILD_GODDARD:BOOL=ON. To build Goddard test, add cmake option -DBUILD_TESTING:BOOL=ON.
If Python modules are built, tests/testGoddard.py can be used.

License
--------

SOCP is released under the 3-clause BSD license. It is free software, and is released as an open-source package. 
Please however be aware that a solver from the CMinpack package within SOCP is licensed using a different 
BSD-like agreement to SOCP, thus ensure you are familiar with it before using SOCP in any commercial work, 
or generating products based on SOCP.
CMinpack License file is copied in "documentation" folder.
