'''
Created on Apr 9, 2019

@author: Bruno Herisse (ONERA/DTIS)

The algorithm is presented in the paper "Singular Arcs in the Generalized Goddard’s Problem", F. Bonnans, P. Martinon, E. Trélat (J Optim Theory Appl (2008) 139: 439–461)
'''
import math
import sys
import time

sys.path.insert(0, '../binaries/python/')
import pySOCP

sys.path.insert(0, '../binaries/python/goddard')
import pyGoddard

#################################################################################################################################
print("Initialize the OCP: ")

# new model object
my_fileTrace = "../trace/goddard/trace.dat"					# file full path for trace
my_goddard = pyGoddard.goddard(my_fileTrace)				# The OCP is defined in the goddard class
goddardDim = pyGoddard.goddard.GetDim(my_goddard)			# GetDim return the dimension of the dynamic model

# new shooting object
nMulti = 6;													# multiple shooting parameter (>0)
my_shooting = pySOCP.shooting(my_goddard, nMulti, 1);		# model / number for multiple shooting (1 if simple shooting) / number of threads (1 for one thread)
pySOCP.shooting.SetPrecision(my_shooting, 1e-6);			# relative tolerance

# vectors to read the trajectory
vt = pySOCP.dVector(nMulti+1);								# time line (t0 .. t(nMulti+1))
vX = pySOCP.dMatrix(nMulti+1);								# state and costate at times vt
							
# set mode for final time and state
mode_tf = 1;												# 0 for fixed time, 1 for free optimized time, 2 otherwise
mode_Xf = pySOCP.iVector(goddardDim,0);						# 0 for fixed state, 1 for free optimized state, 2 otherwise
mode_Xf[3] = 1;	# vxf is free
mode_Xf[4] = 1;	# vyf is free
mode_Xf[5] = 1;	# vzf is free
mode_Xf[6] = 1;	# mf is free
pySOCP.shooting.SetMode(my_shooting, mode_tf, mode_Xf);		# setting mode to the solver

# set initial guess for the shooting (trivial guess)
ti = 0;											# initial time
Xi = pySOCP.dVector(2*goddardDim);				# with the adjoint vector, the full initial state dimension is twice the model dimension
Xi[0] =		0.999949994;						# x
Xi[1] =		0.0001;								# y
Xi[2] =		0.01;								# z
Xi[3] =		1e-10;								# vx			// !=0 to avoid dividing by 0
Xi[4] =		1e-10;								# vy
Xi[5] =		1e-10;								# vz	
Xi[6] =		1.0;								# mass	
Xi[7] =		0.1;								# p_x			// the costate vector is the guess
Xi[8] =		0.1;								# p_y				
Xi[9] =		0.1;								# p_z
Xi[10] =	0.1;								# p_vx
Xi[11] =	0.1;								# p_vy
Xi[12] =	0.1;								# p_vz
Xi[13] =	0.1;								# p_mass
tf = 0.1;										# initial guess for tf
Xf = pySOCP.dVector(2*goddardDim);				# no need to specify the costate vector for Xf
Xf[0] =		1.01;								# x
Xf[1] =		0.0;								# y	
Xf[2] =		0.0;								# z
Xf[3] =		0.0;								# vx (free)
Xf[4] =		0.0;								# vy (free)
Xf[5] =		0.0;								# vz (free)
Xf[6] =		0.0;								# mass (free)
pySOCP.shooting.InitShooting(my_shooting, ti, Xi, tf, Xf);	# init shooting with this guess

print("	done.")

#################################################################################################################################
print("Solve the OCP: ")

# to return error flag
info = 1;

# to return computing time
time1 = time.time();

# start with no drag to compute a better guess from the trivial guess
pyGoddard.goddard.SetParameterDataName(my_goddard, "KD", 0.0);

# solve OCP
print("	Solve OCP without drag (KD = 0): ", end="")
if info==1:
	info = pySOCP.shooting.SolveOCP(my_shooting, 0.0);		# solve OCP without continuation (continuationstep = 0)
print("OK = " + str(info))

# then, start a continuation on KD
print("	Continuation on parameter KD (KD = 0 to 310): ", end="")
if info==1:
	info = pySOCP.shooting.SolveOCP(my_shooting, 1.0, pyGoddard.goddard.GetParameterDataName(my_goddard, "KD"), 310);		# solve OCP with continuation step 1
print("OK = " + str(info))

# then, start a continuation on muL2
print("	Continuation on parameter muL2 (muL2 = 1 to 0.2): ", end="")
if info==1:
	info = pySOCP.shooting.SolveOCP(my_shooting, 1.0, pyGoddard.goddard.GetParameterDataName(my_goddard, "muL2"), 0.2);		# solve OCP with continuation step 1
print("OK = " + str(info))

# to set muL2 = 0, we need to change the structure of the solution (Bang-Singular-Off)
pySOCP.shooting.GetSolution(my_shooting, vt, vX)		# guess from previous solution
tf = vt[nMulti]											# tf from last computed solution
SwitchingTimes_1 = 0.0227; 								# Define switching times for singular arcs (estimated from observation of the previous solution)
SwitchingTimes_2 = 0.08;	
vt[0] = ti; 											# Define a new timeline according to switching times
vt[1] = SwitchingTimes_1/2; 
vt[2] = SwitchingTimes_1; 
vt[3] = (SwitchingTimes_2+SwitchingTimes_1)/2; 
vt[4] = SwitchingTimes_2; 
vt[5] = (SwitchingTimes_2+tf)/2;
vt[6] = tf;		
for i in range(0, nMulti):
	vX[i] = pySOCP.shooting.Move(my_shooting, vt[i]);	# Compute the guess on this new timeline

# update mode (define the new structure of the trajectory)
mode_t = pySOCP.iVector(nMulti+1, 2);					# define time modes on the timeline				
mode_t[0] = 0;	# ti (fixed)
mode_t[2] = 1;	# SwitchingTimes_1 (free)
mode_t[4] = 1;	# SwitchingTimes_2 (free)
mode_t[nMulti] = mode_tf;	# tf
mode_X = pySOCP.iMatrix(nMulti+1); 						# define state modes on the timeline
mode_X[0] = pySOCP.iVector(goddardDim,0);				# Xi is fixed
mode_X[nMulti] = mode_Xf;								
for i in range(1, nMulti-1):
	mode_X[i] = pySOCP.iVector(goddardDim,2);
pySOCP.shooting.SetMode(my_shooting, mode_t, mode_X);	# setting mode to the solver

# init
pySOCP.shooting.InitShooting(my_shooting, vt, vX);		# init shooting with the guess

# change cost of the model (L1 norm only to maximize the mass)
pyGoddard.goddard.SetParameterDataName(my_goddard, "muL2", 0);
# set singular control
#pyGoddard.goddard.SetParameterDataName(my_goddard, "singularControl", 0.6);	# constant approximation for the singular control
pyGoddard.goddard.SetParameterDataName(my_goddard, "singularControl", -1);		# if singularControl < 0, then true explicit singular control

# solve OCP
print("	Solve OCP with singular control (muL2 = 0): ", end="");
if info==1:
	info = pySOCP.shooting.SolveOCP(my_shooting, 0.0);							# solve OCP with continuation step 0 (no continuation possible here)
print("OK = " + str(info));

# print computing time
time2 = time.time();
time = (time2 - time1)
print("	Computing time: " + str(time) + " secondes.") 

#################################################################################################################################
print("Trace the solution: ")
print("	Writting trace/goddard/trace.dat")
pySOCP.shooting.Trace(my_shooting)

#################################################################################################################################
input("\nPress enter to exit !")