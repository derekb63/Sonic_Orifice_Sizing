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

if __name__ == '__main__':
    F_O=Fuel_Oxidizer_Ratio(1,'CH4','O2')


    
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


# Mix_rho caclulates the density of the fuel and oxidizer to get a
# more accurate estimate of the flow rates that the PDE will need rather than
# just using air
def Mix_rho(fuel, ox, F_O, T, P):
    gas = ct.Solution('gri30.cti')
    X = {'{0}'.format(fuel): F_O, '{0}'.format(ox): (1-F_O)}
    gas.TPX = T, P, X
    rho_mix = gas.density
    return rho_mix


# Return the mass flow rate through the sonic orifice assuming the ideal gas
def m_dot(Orifice, or_unit, P_u, p_unit, T, Gas):
    Orifice = conv_in_m(Orifice, or_unit, 'm')
    P_u = conv_Pa_psi(P_u, p_unit, 'Pa')
    A = A_orf(Orifice)
    [rho, k, MW] = Calc_Props(Gas, T, P_u)
    R = ct.gas_constant
    m_dot = A * P_u * k * np.sqrt((2/(k+1))**((k+1)/(k-1)))/np.sqrt((R*T)/(k*MW))
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
    APu = (m_dot*np.sqrt((k*R*T)/MW))/np.sqrt((2/(k+1))**((k+1)/(k-1)))/k
    return APu


# Convert from in to m
def conv_in_m(measurement_to_convert, starting_unit, ending_unit):
    if starting_unit == ending_unit:
        output = measurement_to_convert
    elif starting_unit == 'in' and ending_unit == 'm':
        output = np.multiply(measurement_to_convert, 0.0254)
    elif starting_unit == 'm' and ending_unit == 'in':
        output = np.divide(measurement_to_convert, 0.0254)
    else:
        print('Unit combination is not recognized')
    return output


# Convert from Pa to Psi and from psi to Pa
def conv_Pa_psi(value, starting_unit, ending_unit):
    if starting_unit == ending_unit:
        output = value
    elif starting_unit == 'psi' and ending_unit == 'Pa':
        output = np.multiply(value, 6894.75728)
    elif starting_unit == 'Pa' and ending_unit == 'psi':
        output = np.multiply(value, 0.000145037738007)
    else:
        print('Unit combination is not recognized')
    return output


# Calculate the required pressure and orifice size for the prescribed
# conditions
def pressure_orifice_finder(gas, m_dot_gas, T, P_avg, Orifices, p_max_gas,
                            p_min_gas):
    # APu is the product of the orifice area and upstream pressure calculated
    # by rearranging the mass flow rate equation for a sonic orifice
    APu_gas = APu_prod(m_dot_gas, T, gas, P_avg)
    # Find the index of the closest value in an array to the input variable
    possible_pressure = np.linspace(p_min_gas, p_max_gas, num=10000)
    Areas = A_orf(Orifices)
    APu_poss = pd.DataFrame(np.einsum('i,j-> ji', Areas, possible_pressure),
                            index=possible_pressure, columns=Orifices)
    idx = APu_poss.sub(APu_gas).abs().min().idxmin()
    val = APu_poss.sub(APu_gas).abs().min(axis=1).idxmin()
    print(val)
    Pressure = conv_Pa_psi(val, 'Pa', 'psi')
    print(Pressure)
    Orifice = conv_in_m(idx, 'm', 'in')
    
    print(APu_poss)
    '''
    print()
    print(idx)
    x=APu_poss.sub(APu_gas).abs().min()
    print(x)
    print()
    '''
    return (Pressure, Orifice)






