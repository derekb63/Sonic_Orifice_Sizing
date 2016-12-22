# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 17:27:48 2016

@author: beande
"""

# Function definitions for the flowproject file
# Import the required packages
import numpy as np
import pandas as pd
import cantera as ct


# Fuel_Oxidizer_Cals takes the imput strings of the fuel and oxidizer as well
# as the desired equivalence ratio and calulates the fuel/oxidizer mass ratio
def Fuel_Oxidizer_Ratio(phi, fuel, ox):
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


# Mix_rho caclulates the density of the completely burned products to get a
# more accurate extiamte of the flow rates that the PDE will need rather than
# just using air
def Mix_rho(fuel, ox, F_O, T, P):
    gas = ct.Solution('gri30.cti')
    X = {'{0}'.format(fuel): F_O, '{0}'.format(ox): (1-F_O)}
    gas.TPX = T, P, X
    rho_mix = gas.density
    return rho_mix


# Return the mass flow rate through the sonic orifice assuming the ideal gas
def m_dot(A, P_u, T, Gas):
    [rho, k, MW] = Calc_Props(Gas, T, P_u)
    R = ct.gas_constant
    m_dot = A * P_u * k * np.sqrt((2/(k+1))**((k+1)/(k-1)))/np.sqrt((k*R*T)/MW)
    return m_dot


# Caclualte the are of the orifice based on the orifice diameter
def A_orf(D):
    A_orf = np.pi / 4 * D**2
    return A_orf


# Calculate the product of the area and upstream pressure since both of these
# variables are varaible in the experiment and constrined betwwe nsome values
# the product of the two will be used to find the optimal orifice and pressure
def APu_prod(m_dot, T, Gas, P_guess):
    [rho, k, MW] = Calc_Props(Gas, T, P_guess)
    R = ct.gas_constant
    APu = (m_dot*np.sqrt((k*R*T)/MW))/np.sqrt((2/(k+1))**((k+1)/(k-1)))
    return APu


# Convert from in to m
def conv_in_m(measurement_to_convert, starting_unit):
    if starting_unit=='in':
        output = np.multiply(measurement_to_convert, 0.0254)
    elif starting_unit=='m':
        output = np.divide(measurement_to_convert, 0.0254)
    else:
        print('Unit is not recognized')
    return output


# Find the index of the closest value in an array to the input variable
def find_closest(p_max, value, Orifices):
   possible = np.linspace(350000, p_max, num=200)
   Areas = A_orf(Orifices)
   APu_poss = pd.DataFrame(np.einsum('i,j-> ji', Areas, possible),
                    index=possible, columns=Orifices)
   idx = APu_poss.sub(value).abs().min().idxmin()
   val = APu_poss.sub(value).abs().min(axis=1).idxmin()
   Pressure = conv_Pa_psi(val,'Pa')
   Orifice_size = conv_in_m(idx, 'm')
   return (Pressure, Orifice_size)

# Convert from Pa to Psi and from psi to Pa
def conv_Pa_psi(value,starting_unit):
    if starting_unit=='psi':
        output = np.multiply(value, 6894.75728)
    elif starting_unit=='Pa':
        output = np.multiply(value, 0.000145037738007)
    else:
        print('Unit is not recognized')
    return output
    
