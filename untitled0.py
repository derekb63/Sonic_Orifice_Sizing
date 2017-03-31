#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sympy solver for calculating the unknown variables in the mass flow rate
through a sonic orifice

@author: beande
"""

from sympy.solvers import solve
from sympy import Symbol
import cantera as ct
import numpy as np

m_dot = Symbol('m')
A = Symbol('A')
P_o = Symbol('P')
k = Symbol('k')
C_d = Symbol('Cd')
[rho, k, MW] = Calc_Props(Gas, T, P_u)
R = 8314
T = 298

def m_dot(Orifice, or_unit, P_u, p_unit, T, Gas, C_d=.99):
    Orifice = conv_in_m(Orifice, or_unit, 'm')
    P_u = conv_Pa_psi(P_u, p_unit, 'Pa')
    A = A_orf(Orifice)
    [rho, k, MW] = Calc_Props(Gas, T, P_u)
    R = ct.gas_constant
    m_dot = A*P_u*k*C_d*np.sqrt((2/(k+1))**((k+1)/(k-1)))/np.sqrt((R*T)/(k*MW))
    return m_dot


# Calc_Props takes the input termperature and pressure of the specified gas and
# outputs the properties required for the calculation of the mass flow rate
# through the sonic orifice using Cantera and GriMech 3.0
def Calc_Props(Gas, T, P):
    gas = ct.Solution('gri30.cti')
    gas.TPX = T, P, '{0}:1'.format(Gas)
    rho = gas.density
    k = gas.cp_mass/gas.cv_mass
    MW = gas.mean_molecular_weight
    return rho, k, MW
