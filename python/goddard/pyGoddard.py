'''
Created on Jul 31, 2020

@author: Bruno Herisse (ONERA/DTIS)

The algorithm is presented in the paper "Singular Arcs in the Generalized Goddard’s Problem", F. Bonnans, P. Martinon, E. Trélat (J Optim Theory Appl (2008) 139: 439–461)
'''
import math
import sys
import time

import os.path
from os import path

import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

sys.path.insert(0, '../../binaries/python/')
import pySOCP

class parameters_struct(object):
    def __init__(self, C, b, KD, kr, u_max, mu1, mu2, singularControl):
        self.C = C;                             # coefficient for thrust
        self.b = b;                             # coefficient for mass flow rate
        self.KD = KD;                           # coefficient for drag
        self.kr = kr;                           # coefficient for density of air
        self.u_max = u_max;                     # max normalized control
        self.mu1 = mu1;                         # weight for the cost on the norm of the control in [0,1]
        self.mu2 = mu2;                         # weight for the quadratic cost on the control in [0,1]
        self.singularControl = singularControl; # singular control value

class data_struct(object):
    def __init__(self, switchingTimes):
        self.switchingTimes = switchingTimes    # switching times for singular control

class goddard(pySOCP.model):
    def __init__(self, the_fileTrace, stepNbr):
        super().__init__(7, stepNbr, the_fileTrace)         # call the constructor of model class
        paramStruct = parameters_struct(3.5, 7, 310, 500, 1, 1, 0, -1)
        self.parameters["C"] = paramStruct.C;								# C: coefficient for thrust
        self.parameters["b"] = paramStruct.b;								# b: coefficient for mass flow rate
        self.parameters["KD"] = paramStruct.KD;							# KD: coefficient for drag
        self.parameters["kr"] = paramStruct.kr;							# kr: coefficient for density of air
        self.parameters["u_max"] = paramStruct.u_max;						# u_max: max normalized control
        self.parameters["mu1"] = paramStruct.mu1;							# mu1: weight for the cost on the norm of the control in [0,1]
        self.parameters["mu2"] = paramStruct.mu2;							# mu2: weight for the quadratic cost on the control in [0,1]
        self.parameters["singularControl"] = paramStruct.singularControl;	# singularControl: a constant approximation of the singular arc
        switchingTimes = [0.0227, 0.08]
        self.data = data_struct(switchingTimes)
        
    def Model(self, t, X):
        # current state and costate value
        x = X[0]
        y = X[1]
        z = X[2]
        vx = X[3]
        vy = X[4]
        vz = X[5]
        mass = X[6]
        p_x = X[7]
        p_y = X[8]
        p_z = X[9]
        p_vx = X[10]
        p_vy = X[11]
        p_vz = X[12]
        p_mass = X[13]
        
        r = math.sqrt(x*x + y*y + z*z)
        v = math.sqrt(vx*vx + vy*vy + vz*vz)
        pvdotv = p_vx*vx + p_vy*vy + p_vz*vz
        b = self.parameters["b"]
        C = self.parameters["C"]
        KD = self.parameters["KD"]
        kr = self.parameters["kr"]
        g = 1 / r / r                                           # normalized gravity (g=1 for r=1 that corresponds to Earth radius)
        norm_pv = math.sqrt(p_vx*p_vx + p_vy*p_vy + p_vz*p_vz)
        
        # control computation
        u = self.Control(t, X)
        norm_u = math.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2])
        pvdotu = p_vx*u[0] + p_vy*u[1] + p_vz*u[2]
        
        # state and costate equations
        # Xdot = list(range(len(X)))
        Xdot = pySOCP.dVector(2*self.dim);
        Xdot[0] = vx;
        Xdot[1] = vy;
        Xdot[2] = vz;
        Xdot[3] = -KD*v*vx*math.exp(-kr*(r - 1)) / mass - g*x / r + C*u[0] / mass;
        Xdot[4] = -KD*v*vy*math.exp(-kr*(r - 1)) / mass - g*y / r + C*u[1] / mass;
        Xdot[5] = -KD*v*vz*math.exp(-kr*(r - 1)) / mass - g*z / r + C*u[2] / mass;
        Xdot[6] = -b*norm_u;
        Xdot[7] = -kr*KD / mass*v*math.exp(-kr*(r - 1))*x / r*pvdotv + g*(p_vx*(1 - 3 * x*x / r / r) / r - p_vy * 3 * x*y / r / r / r - p_vz * 3 * x*z / r / r / r);
        Xdot[8] = -kr*KD / mass*v*math.exp(-kr*(r - 1))*y / r*pvdotv + g*(-p_vx * 3 * y*x / r / r / r + p_vy*(1 - 3 * y*y / r / r) / r - p_vz * 3 * y*z / r / r / r);
        Xdot[9] = -kr*KD / mass*v*math.exp(-kr*(r - 1))*z / r*pvdotv + g*(-p_vx * 3 * z*x / r / r / r - p_vy * 3 * z*y / r / r / r + p_vz*(1 - 3 * z*z / r / r) / r);
        Xdot[10] = -p_x + KD / mass*math.exp(-kr*(r - 1))*(pvdotv*vx / v + p_vx*v);
        Xdot[11] = -p_y + KD / mass*math.exp(-kr*(r - 1))*(pvdotv*vy / v + p_vy*v);
        Xdot[12] = -p_z + KD / mass*math.exp(-kr*(r - 1))*(pvdotv*vz / v + p_vz*v);
        Xdot[13] = -KD*math.exp(-kr*(r - 1)) / mass / mass*v*pvdotv + C / mass / mass*pvdotu;

        return Xdot
    
    @staticmethod
    def staticModel(X, t, goddardModel):
        vecX = goddardModel.Model(t,X)
        listX = list(range(vecX.size()))
        for i in range(len(listX)):
            listX[i] = vecX[i]
        return listX
    
    def Control(self, t, X):
        # current state and costate value
        x = X[0]
        y = X[1]
        z = X[2]
        vx = X[3]
        vy = X[4]
        vz = X[5]
        mass = X[6]
        p_x = X[7]
        p_y = X[8]
        p_z = X[9]
        p_vx = X[10]
        p_vy = X[11]
        p_vz = X[12]
        p_mass = X[13]
        
        r = math.sqrt(x*x + y*y + z*z)
        v = math.sqrt(vx*vx + vy*vy + vz*vz)
        pvdotv = p_vx*vx + p_vy*vy + p_vz*vz
        b = self.parameters["b"]
        C = self.parameters["C"]
        KD = self.parameters["KD"]
        kr = self.parameters["kr"]
        g = 1 / r / r                                           # normalized gravity (g=1 for r=1 that corresponds to Earth radius)
        norm_pv = math.sqrt(p_vx*p_vx + p_vy*p_vy + p_vz*p_vz)
        
        # control computation (minimisation of the Hamiltoninan)
        alpha_u = 0;
        u = list(range(3));
        Switch = self.parameters["mu1"] - b*p_mass - C / mass*norm_pv;        # switching function
        
        if (self.parameters["mu2"] > 0):
            # if a quadratic cost is used 
            if (Switch<0):
                alpha_u = -Switch / 2 / self.parameters["mu2"];
            else:       #(Switch>=0)
                alpha_u = 0;
        else:
            # the strcture of the control is imposed : Bang-Singular-Off
            if (t <= (self.data.switchingTimes[0])):
                alpha_u = 1.0;
            elif (t > self.data.switchingTimes[0] and t <= self.data.switchingTimes[1]):
                if (self.parameters["singularControl"] < 0):
                    alpha_u = self.GetSingularControl(t, X)     # true singular control 
                else:
                    alpha_u = self.parameters["singularControl"];     # approximation of the singular control by a constant parameter
            else:
                alpha_u = 0;
        
        u[0] = -p_vx*alpha_u / norm_pv;
        u[1] = -p_vy*alpha_u / norm_pv;
        u[2] = -p_vz*alpha_u / norm_pv;

        # saturation
        norm_u = math.fabs(alpha_u);
        u_max = self.parameters["u_max"]
        if (norm_u > u_max):
            u[0] = u[0] / norm_u*u_max;
            u[1] = u[1] / norm_u*u_max;
            u[2] = u[2] / norm_u*u_max;
            norm_u = u_max;

        # control = list(range(3));
        control = pySOCP.dVector(3);
        control[0] = u[0];
        control[1] = u[1];
        control[2] = u[2];

        return control;
    
    def GetSingularControl(self, t, X):
        # current state and costate value
        x = X[0]
        y = X[1]
        z = X[2]
        vx = X[3]
        vy = X[4]
        vz = X[5]
        mass = X[6]
        p_x = X[7]
        p_y = X[8]
        p_z = X[9]
        p_vx = X[10]
        p_vy = X[11]
        p_vz = X[12]
        p_mass = X[13]
        
        r = math.sqrt(x*x + y*y + z*z)
        v = math.sqrt(vx*vx + vy*vy + vz*vz)
        rdotv = x*vx + y*vy + z*vz
        pvdotv = p_vx*vx + p_vy*vy + p_vz*vz
        b = self.parameters["b"]
        C = self.parameters["C"]
        KD = self.parameters["KD"]
        kr = self.parameters["kr"]
        g = 1 / r / r                                           # normalized gravity (g=1 for r=1 that corresponds to Earth radius)
        norm_pv = math.sqrt(p_vx*p_vx + p_vy*p_vy + p_vz*p_vz)
        D = KD*math.exp(-kr*(r - 1))
        
        p_xdot = -kr*KD / mass*v*math.exp(-kr*(r - 1))*x / r*pvdotv + g*(p_vx*(1 - 3 * x*x / r / r) / r - p_vy * 3 * x*y / r / r / r - p_vz * 3 * x*z / r / r / r)
        p_ydot = -kr*KD / mass*v*math.exp(-kr*(r - 1))*y / r*pvdotv + g*(-p_vx * 3 * y*x / r / r / r + p_vy*(1 - 3 * y*y / r / r) / r - p_vz * 3 * y*z / r / r / r)
        p_zdot = -kr*KD / mass*v*math.exp(-kr*(r - 1))*z / r*pvdotv + g*(-p_vx * 3 * z*x / r / r / r - p_vy * 3 * z*y / r / r / r + p_vz*(1 - 3 * z*z / r / r) / r)
        p_vxdot = -p_x + KD / mass*math.exp(-kr*(r - 1))*(pvdotv*vx / v + p_vx*v)
        p_vydot = -p_y + KD / mass*math.exp(-kr*(r - 1))*(pvdotv*vy / v + p_vy*v)
        p_vzdot = -p_z + KD / mass*math.exp(-kr*(r - 1))*(pvdotv*vz / v + p_vz*v)
        
        prdotdotv = p_xdot*vx + p_ydot*vy + p_zdot*vz
        prdotdotpv = p_xdot*p_vx + p_ydot*p_vy + p_zdot*p_vz
        prdotpvdot = p_x*p_vxdot + p_y*p_vydot + p_z*p_vzdot
        prdotpv = p_x*p_vx + p_y*p_vy + p_z*p_vz
        pvdotdotg = p_vxdot*g*x / r + p_vydot*g*y / r + p_vzdot*g*z / r
        pvdotdotv = p_vxdot*vx + p_vydot*vy + p_vzdot*vz
        pvdotdotpv = p_vxdot*p_vx + p_vydot*p_vy + p_vzdot*p_vz
        vdotg = vx*g*x / r + vy*g*y / r + vz*g*z / r
        pvdotg = p_vx*g*x / r + p_vy*g*y / r + p_vz*g*z / r
        prdotv = p_x*vx + p_y*vy + p_z*vz
        prdotg = p_x*g*x / r + p_y*g*y / r + p_z*g*z / r
        
        # singular control computation (minimisation of the Hamiltoninan)
        alpha_u = 0
        au = 2 * norm_pv*C / mass*pvdotv\
            + 2 * pvdotv*C / mass*norm_pv\
            - b / mass*(2 * pvdotv*pvdotv + norm_pv*norm_pv*v*v)\
            - b / D*v*prdotpv - C / D*prdotpv / v*pvdotv / norm_pv
        bu = -2 * norm_pv*norm_pv*(vdotg + D / mass*v*v*v) + 2 * v*v*pvdotdotpv\
            - 2 * pvdotv*(pvdotg + D / mass*v*pvdotv - pvdotdotv)\
            + b / C*(2 * norm_pv*pvdotv*(vdotg + D / mass*v*v*v) + norm_pv*v*v*(pvdotg + D / mass*v*pvdotv - pvdotdotv) - v*v*pvdotv / norm_pv*pvdotdotpv)\
            - mass / D*kr*rdotv / r*v*prdotpv + mass / D*prdotpv / v*(vdotg + D / mass*v*v*v) - mass / D*v*(prdotdotpv + prdotpvdot)
        alpha_u = bu / au
        
        return alpha_u;
    
    def Hamiltonian(self, t, X):
        # current state and costate value
        x = X[0]
        y = X[1]
        z = X[2]
        vx = X[3]
        vy = X[4]
        vz = X[5]
        mass = X[6]
        p_x = X[7]
        p_y = X[8]
        p_z = X[9]
        p_vx = X[10]
        p_vy = X[11]
        p_vz = X[12]
        p_mass = X[13]
        
        r = math.sqrt(x*x + y*y + z*z)
        v = math.sqrt(vx*vx + vy*vy + vz*vz)
        pvdotv = p_vx*vx + p_vy*vy + p_vz*vz
        b = self.parameters["b"]
        C = self.parameters["C"]
        KD = self.parameters["KD"]
        kr = self.parameters["kr"]
        g = 1 / r / r                                           # normalized gravity (g=1 for r=1 that corresponds to Earth radius)
        norm_pv = math.sqrt(p_vx*p_vx + p_vy*p_vy + p_vz*p_vz)
        
        u = self.Control(t, X)
        norm_u = math.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2])

        # H is computed
        H = self.parameters["mu1"]*norm_u + self.parameters["mu2"]*norm_u*norm_u\
            + p_x*vx + p_y*vy + p_z*vz\
            + p_vx*(-KD*v*vx*math.exp(-kr*(r - 1)) / mass - g*x / r + C*u[0] / mass)\
            + p_vy*(-KD*v*vy*math.exp(-kr*(r - 1)) / mass - g*y / r + C*u[1] / mass)\
            + p_vz*(-KD*v*vz*math.exp(-kr*(r - 1)) / mass - g*z / r + C*u[2] / mass)\
            - p_mass*b*norm_u;

        return H;
    
    # Optional redefinition of ModelInt
    # def ModelInt(self, t0, X0, tf, isTrace):
    #     dt = (tf-t0) / self.stepNbr;     # time step
    #     Xs = list(range(X0.size()))
    #     for i in range(X0.size()):
    #         Xs[i] = X0[i]
        
    #     t = np.linspace(t0, tf, self.stepNbr+1)      # the points of evaluation of the solution
    #     X = odeint(goddard.staticModel, Xs, t, args=(self,),rtol=self.odeIntTol, atol=self.odeIntTol)
                   
    #     Xs = X[-1,:]
        
    #     if (isTrace):
    #         fileTrace  = open(self.strFileTrace, "a+")
    #         for i in range(t.size):
    #             fileTrace.write("%s\t" % t[i])
    #             for item in X[i,:]:
    #                 fileTrace.write("%s\t" % item)
    #             control = self.Control(t[i], X[i,:]);
    #             for item in control:
    #                 fileTrace.write("%s\t" % item)
    #             fileTrace.write("%s\t" % self.Hamiltonian(t[i], X[i,:]))
    #             fileTrace.write("\n")
    #         fileTrace.close();
        
    #     Xr = pySOCP.dVector(2*self.dim);
    #     for i in range(Xr.size()):
    #         Xr[i] = Xs[i]
            
    #     return Xr;
    
    def SwitchingTimesFunction(self, t, X):
        # current state and costate value
        x = X[0]
        y = X[1]
        z = X[2]
        vx = X[3]
        vy = X[4]
        vz = X[5]
        mass = X[6]
        p_x = X[7]
        p_y = X[8]
        p_z = X[9]
        p_vx = X[10]
        p_vy = X[11]
        p_vz = X[12]
        p_mass = X[13]

        b = self.parameters["b"]
        C = self.parameters["C"]
        norm_pv = math.sqrt(p_vx*p_vx + p_vy*p_vy + p_vz*p_vz)

        # Switching function
        #Switch = data->parameters.mu1 - b*p_mass - C/mass*norm_pv;

        # Using Switch or Hamiltonian is equivalent in this case (free final time)
        #fvec = Switch;
        fvec = self.Hamiltonian(t, X);
        
        return fvec;
        
    def SwitchingTimesUpdate(self, switchingTimes):
        # Switching times
        self.data.switchingTimes = list(range(switchingTimes.size()));
        for i in range(switchingTimes.size()):
            self.data.switchingTimes[i] = switchingTimes[i];
            
    def SetParameterDataName(self, name, value):
        self.parameters[name] = value;
            
    def GetParameterDataName(self, name):
        return self.parameters[name]
        return self.parameters[name];

# TESTS
# my_goddard = goddard("../../trace/goddard/trace.dat")
# my_goddard.SetParameterDataName("KD",0)
# print(my_goddard.GetParameterDataName("KD"))
# print(my_goddard.odeIntTol)
# print(my_goddard.GetDim())
# print(my_goddard.GetParameterDataName("C"))
# switchPoints = pySOCP.dVector(3,0)
# switchPoints[0] = 1
# switchPoints[1] = 2
# switchPoints[2] = 4
# my_goddard.SwitchingTimesUpdate(switchPoints)
# print(my_goddard.data.switchingTimes)
# ti = 0;                                         # initial time
# Xi = pySOCP.dVector(2*7);              # with the adjoint vector, the full initial state dimension is twice the model dimension
# Xi[0] =     0.999949994;                        # x
# Xi[1] =     0.0001;                             # y
# Xi[2] =     0.01;                               # z
# Xi[3] =     1e-10;                              # vx            // !=0 to avoid dividing by 0
# Xi[4] =     1e-10;                              # vy
# Xi[5] =     1e-10;                              # vz    
# Xi[6] =     1.0;                                # mass  
# Xi[7] =     0.1;                                # p_x           // the costate vector is the guess
# Xi[8] =     0.1;                                # p_y               
# Xi[9] =     0.1;                                # p_z
# Xi[10] =    0.1;                                # p_vx
# Xi[11] =    0.1;                                # p_vy
# Xi[12] =    0.1;                                # p_vz
# Xi[13] =    0.1;                                # p_mass
# tf = 0.1;                                       # initial guess for tf
# print(my_goddard.Control(ti,Xi)[0])
# print(my_goddard.Model(ti,Xi)[0])
# print(my_goddard.Hamiltonian(ti,Xi))
# print(my_goddard.SwitchingTimesFunction(ti,Xi))
# Xf = my_goddard.ModelInt(ti,Xi,tf,1)
# print(Xf[0])
# Xf = my_goddard.ComputeTraj(ti,Xi,tf,0)
# print(Xf[0])
# print(my_goddard.parameters["C"])

