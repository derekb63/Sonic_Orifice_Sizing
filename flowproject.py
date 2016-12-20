#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 13:33:10 2016

@author: derek
"""

import numpy as np
import cantera as ct

from flowprojectfunc import Fuel_Oxidizer_Cals
from flowprojectfunc import Calc_Props
from flowprojectfunc import Mix_rho
from flowprojectfunc import m_dot
from flowprojectfunc import A_orf
from flowprojectfunc import APu_prod

fuel = 'C3H8'
ox = 'N2O'
phi = 1.0
F_O = Fuel_Oxidizer_Cals(phi, fuel, ox)
T = 298
P = 101325
L = 2
D_tube = 0.0762
Op_freq = 1

[rho_fuel, k_fuel, MW_fuel] = Calc_Props(fuel, T, P)
[rho_ox, k_ox, MW_ox] = Calc_Props(ox, T, P)

rho_mix = Mix_rho(fuel, ox, F_O, T, ct.one_atm)
m_dot_tube = A_orf(D_tube)*L*rho_mix*Op_freq
print(str(m_dot_tube) + ' (kg/s)')

m_dot_ox = m_dot_tube / (1 + F_O)
m_dot_fuel = F_O * m_dot_ox
