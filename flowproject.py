#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 13:33:10 2016

@author: derek
"""
import time
import numpy as np
import cantera as ct

from flowprojectfunc import Fuel_Oxidizer_Ratio
from flowprojectfunc import Mix_rho
from flowprojectfunc import A_orf
from flowprojectfunc import conv_in_m
from flowprojectfunc import pressure_orifice_finder

t0 = time.time()
# Input parameters to define the PDE/System variables
fuel = 'C3H8'
ox = 'N2O'
phi = 1.0
T = 298
P = 101325
P_avg = 700000
L = 2
D_tube = 0.0762
Op_freq = 1
p_max_ox = 1.379E6
p_max_fuel = 689467

# Possible orifice sizes with converted to m for inputting into the
# find_closest function
Orifices = np.array(conv_in_m([0.040, 0.047, 0.063, 0.142], 'in'))

# Use the Fuel_Oxidizer_Ratio function to calculate the F/O for the mixture
# at the specified eqivalence ratio
F_O = Fuel_Oxidizer_Ratio(phi, fuel, ox)

# Calculate the necessary operating flow rates of the PDE based on the
# previously defined operating paramters
rho_mix = Mix_rho(fuel, ox, F_O, T, ct.one_atm)
m_dot_tube = A_orf(D_tube)*L*rho_mix*Op_freq

# Calculate the mass flow rate of the oxidizer for the defined conditions
m_dot_ox = m_dot_tube / (1 + F_O)
[Pressure_ox, Orifice_ox] = pressure_orifice_finder(ox, m_dot_ox, T, P_avg,
                                                    Orifices, p_max_ox)

# Calculate the mass flow rate of the oxidizer for the defined conditions
m_dot_fuel = F_O * m_dot_ox
[Pressure_f, Orifice_f] = pressure_orifice_finder(fuel, m_dot_fuel, T, P_avg,
                                                  Orifices, p_max_fuel)


print(round(Pressure_f, 2), round(Orifice_f, 3))
print(round(Pressure_ox, 2), round(Orifice_ox, 3))

t1 = time.time()
total = t1-t0
print(total)
