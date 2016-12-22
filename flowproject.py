#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 13:33:10 2016

@author: derek
"""
import time
t0 = time.time()
import numpy as np
import cantera as ct


from flowprojectfunc import Fuel_Oxidizer_Cals
from flowprojectfunc import Calc_Props
from flowprojectfunc import Mix_rho
from flowprojectfunc import A_orf
from flowprojectfunc import APu_prod
from flowprojectfunc import conv_in_m
from flowprojectfunc import find_closest


fuel = 'C3H8'
ox = 'N2O'
phi = 1.0
F_O = Fuel_Oxidizer_Cals(phi, fuel, ox)
T = 298
P = 101325
P_guess = 700000
L = 2
D_tube = 0.0762
Op_freq = 1
p_max_ox = 1.379E6
p_max_fuel = 689467
Orifices = np.array(conv_in_m([0.040, 0.047, 0.063, 0.142], 'in'))


[rho_fuel, k_fuel, MW_fuel] = Calc_Props(fuel, T, P)
[rho_ox, k_ox, MW_ox] = Calc_Props(ox, T, P)

rho_mix = Mix_rho(fuel, ox, F_O, T, ct.one_atm)
m_dot_tube = A_orf(D_tube)*L*rho_mix*Op_freq

m_dot_ox = m_dot_tube / (1 + F_O)
m_dot_fuel = F_O * m_dot_ox

APu_fuel = APu_prod(m_dot_fuel, T, fuel, P_guess)
APu_ox = APu_prod(m_dot_ox, T, ox, P_guess)

[Pressure_f, Orifice_f] = find_closest(p_max_fuel, APu_fuel, Orifices)
[Pressure_ox, Orifice_ox] = find_closest(p_max_ox, APu_ox, Orifices)

print(round(Pressure_f, 2), round(Orifice_f, 3))
print(round(Pressure_ox, 2), round(Orifice_ox, 3))

t1 = time.time()
total = t1-t0
print(total)