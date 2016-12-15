#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 13:33:10 2016

@author: derek
"""

import numpy as np
import pandas as pd
import cantera as ct


# Fuel_Oxidizer_Cals takes the imput strings of the fuel and oxidizer as well
# as the desired equivalence ratio and calulates the fuel/oxidizer mass ratio
def Fuel_Oxidizer_Cals(phi, fuel, ox):
    Elements = ['C', 'H', 'O', 'N']
    ox_val = [0, 0, 0, 0]
    if ox == 'Air':
        ox = ['O2', 'N2']
        ox_val = [0, 0, 2, 7.52]
        MW_ox = 28.8
    elif ox == 'N2O':
        ox_val = [0, 0, 1, 2]
        MW_ox = 44
    elif ox == 'O2':
        ox_val = [0, 0, 2, 0]
        MW_ox = 32
    else:
        print('Your Oxidizer is not reconized')
    if fuel == 'H2':
        fuel_val = [0, 2, 0, 0]
        MW_fuel = 2
    elif fuel == 'CH4':
        fuel_val = [1, 4, 0, 0]
        MW_fuel = 16.04
    elif fuel == 'C3H8':
        fuel_val = [3, 8, 0, 0]
        MW_fuel = 44
    else:
        print('Your Fuel is not reconized')
    react_names = [fuel]
    react_names += [ox]
    product_vals = [(1, 0, 2, 0), (0, 2, 1, 0), (0, 0, 0, 2)]
    product_names = ['CO2', 'H2O', 'N2']
    names = [ox] + product_names
    A = pd.DataFrame(np.transpose(np.vstack([ox_val, product_vals])),
                     index=Elements, columns=names)
    coeffs = np.abs(np.linalg.solve(A[:][:], [-x for x in fuel_val]))
    F_O_s = (1*MW_fuel)/(coeffs[0]*MW_ox)
    F_O = phi*F_O_s
    return F_O


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


def Mix_rho(fuel, ox, F_O, T, P):
    gas = ct.Solution('gri30.cti')
    X = {'{0}'.format(fuel): F_O, '{1}'.format(ox): (1-F_O)}
    gas.TPX = T, P, X
    rho_mix = gas.density
    return rho_mix

fuel = 'H2'
ox = 'O2'
phi = 0.5
F_O = Fuel_Oxidizer_Cals(phi, fuel, ox)
T = 298
P = 101325
L = 2
D_tube = 0.0762
Op_freq = 1

[rho_fuel, k_fuel, MW_fuel] = Calc_Props(fuel, T, P)
[rho_ox, k_ox, MW_ox] = Calc_Props(ox, T, P)

rho_mix = Mix_rho(fuel, ox, F_O, T, ct.one_atm)
m_dot_tube = np.pi*D_tube**2*L*rho_mix*Op_freq
