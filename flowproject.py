#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 13:33:10 2016

@author: derek
"""
import time
import numpy as np
import cantera as ct
from tabulate import tabulate
from sympy.solvers import solve
from sympy import Symbol
import matplotlib.pyplot as plt

from flowprojectfunc import Fuel_Oxidizer_Ratio
from flowprojectfunc import Mix_rho
from flowprojectfunc import A_orf
from flowprojectfunc import conv_in_m
from flowprojectfunc import pressure_orifice_finder
from flowprojectfunc import m_dot
from flowprojectfunc import conv_Pa_psi

t0 = time.time()

## Input Parameters ##
'Note: all units have to be in mksK (meters,kilogram,seconds,Kelvin) units, or their derived units(ie Pa,N, etc)'
# Fuel/Ox molecular formulas
fuel = 'C3H8'
ox = 'N2O'

# Fuel Oxidizer Ratio
phi = 1.0
# ambient conditions
T = 298
P = 101325
# ??
P_avg = 700000 #7 atm
# tube physical dimensions
L = conv_in_m(12,'in','m')
D_tube = conv_in_m(0.5,'in','m')

# the amount of time you want the valves open, in Hz
Op_freq = 20

' Input variables to iterate across to find optimal orifices'
# min and max pressures allowed for the fuel and oxidizer
p_max_ox = conv_Pa_psi(130, 'psi','Pa')
p_max_fuel = conv_Pa_psi(80, 'psi','Pa')
p_min_gas = conv_Pa_psi(20, 'psi','Pa') 

# Possible orifice sizes
Orifices = np.array(conv_in_m(np.arange(0.0125, 0.250, 0.001), 'in', 'm'))

## End Input Parameters ##



fuel_error = []
ox_error = []

# Use the Fuel_Oxidizer_Ratio function to calculate the F/O for the mixture
# at the specified eqivalence ratio
F_O = Fuel_Oxidizer_Ratio(phi, fuel, ox)

# Calculate the necessary operating flow rates of the PDE based on the
# previously defined operating paramters
rho_mix = Mix_rho(fuel, ox, F_O, T, ct.one_atm)
m_dot_tube = A_orf(D_tube)*L*rho_mix*Op_freq

# Calculate the mass flow rate of the oxidizer for the defined conditions
m_dot_ox = m_dot_tube / (1 + F_O)
# Use the matching function to find the nearest orifice and pressure
# combination for the oxidizer
# Orifices = np.array(conv_in_m([.040,0.142], 'in', 'm'))
[Pressure_ox, Orifice_ox] = pressure_orifice_finder(ox, m_dot_ox, T, P_avg,
                                                    Orifices, p_max_ox,
                                                    p_min_gas)

# Calculate the mass flow rate of the oxidizer for the defined conditions
m_dot_fuel = m_dot_tube - m_dot_ox
# Use the matching function to find the nearest orifice and pressure
# combination for the fuel
# Orifices = np.array(conv_in_m([0.047], 'in', 'm'))
[Pressure_f, Orifice_f] = pressure_orifice_finder(fuel, m_dot_fuel, T, P_avg,
                                                  Orifices, p_max_fuel,
                                                  p_min_gas)

# Calculate the error of the output variables to double check the values
m_dot_ox_check = m_dot(Orifice_ox, 'in', Pressure_ox, 'psi', T, ox)
error_ox = np.divide(m_dot_ox-m_dot_ox_check, m_dot_ox) * 100

m_dot_fuel_check = m_dot(Orifice_f, 'in', Pressure_f, 'psi', T, fuel)
error_fuel = np.divide(m_dot_fuel-m_dot_fuel_check, m_dot_fuel) * 100

print(tabulate([[fuel, round(Pressure_f, 2), Orifice_f,
               '{0:.2E}'.format(error_fuel)],
                [ox, round(Pressure_ox, 2), Orifice_ox,
                '{0:.2E}'.format(error_ox)]],
               headers=['Molecule', 'Pressure (psi)', 'Diameter (in)', 'Error (%)']))
# print('FO Error: ', ((m_dot_fuel_check/m_dot_ox_check) - F_O)/F_O * 100, '%')
# print('Fuel: ')
# print('Pressure (psi):', 'Orifice Diameter (in):', 'Error (%)')
# print(round(Pressure_f, 2), round(Orifice_f, 3), round(error_fuel, 10))
# print('Ox: ', round(Pressure_ox, 2), round(Orifice_ox, 3), round(error_ox, 10))

t1 = time.time()
total = t1-t0
#print('Time: ', total)

m_dot_diluent = Symbol('m_dot_diluent')
a =[]
for species_dilution in np.arange(0.05, 0.4, 0.005):

    a.append((solve((m_dot_diluent/(m_dot_diluent+m_dot_fuel+m_dot_ox))-species_dilution,
                   m_dot_diluent), species_dilution))
c = [(m_dot(0.052, 'in', x, 'psi', 298, 'CO2'), x) for x in np.linspace(60, 200, num=55)]

plt.plot([i[1] for i in c], [i[0] for i in c])
plt.xlabel('Upstream Pressure (psi)')
plt.ylabel('m_dot diluent')

#plt.show()







