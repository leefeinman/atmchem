#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 16:18:11 2023

@author: leefeinman

A module containing various back-of-envelope-type calculations for use in an
atmospheric chemistry lab setting. 
"""
#%% Common imports 
import numpy as np


#%% chi_water (mol fraction of water)

def chi_water(T_celc,RH_percent, P_Torr, chi_h20_dry_air = False):
    '''
    chi_water takes air temperature (degrees C), relative humidity (%), and P (Torr)
    and calculates the mole fraction of water in the atmosphere...
    
    
    Created on Fri Apr  7 09:51:46 2023

    Author: leefeinman
    from Andrew Lindsay, translated from .m to .py

    Parameters
    ----------
    T_celc : numpy array, list, tuple, or float 
        Temperature data in centigrade. If held constant, use float if desired. 
    RH_percent : numpy array, list, tuple, or float 
        % relative humidity data. If held constant, use float if desired.
    P_Torr : numpy array, list, tuple, or float 
        Total pressure data. If held constant, use float if desired.
    chi_h20_dry_air : boolean, optional
        If looking for mole fraction of water expressed vs. the amount in dry air. The default is False.

    Returns
    -------
    x_h2o : numpy array or float
        returns the mole fraction of water.

    '''
    
      # changing the arguments to numpy arrays if they aren't already unless
      # they are float/int
    if not isinstance(T_celc, (float, int, np.ndarray)):
        T_celc = np.array(T_celc)
    if not isinstance(RH_percent, (float, int, np.ndarray)):
        RH_percent = np.array(RH_percent)        
    if not isinstance(P_Torr, (float, int, np.ndarray)):
        P_Torr = np.array(P_Torr)
    
    
      # conversion formulas are from Vaisala
      # https://www.vaisala.com/sites/default/files/documents/Humidity_Conversion_Formulas_B210973EN-F.pdf
    A = 6.116441
    m = 7.591386
    Tn = 240.7263 # all constants for T between -20C and +50C (0.08% error at most)
    
      
      # H2O saturation vapor pressure in hectopascal (eq6)
    H2O_sat_P = A * 10**((m * T_celc) / (T_celc + Tn))  # Yes, T is in celcius... this is not an error.
      
    # water pressure in hPa
    Pw = 0.01 * RH_percent * H2O_sat_P
    
      # converting the total pressure (one of the function's arguments)
      # to hPa from torr
    P_hpa = P_Torr * (10**3 / 760)
    
      # note x_h2o is in other contexts expressed vs. the amount in dry air,
      # in that case the denominator should be changed to (P - Pw).
    if chi_h20_dry_air:
        x_h2o = Pw / (P_hpa - Pw)
    else:
        x_h2o = Pw / P_hpa
        
    
    return x_h2o
#%% p_to_c (partial pressure to concentration in ppb)

def p_to_c(P, Punit = "mbar", Cunit = 'ppb', Psystem = 1013.25):
    """
    Calculates the concentration of a substance in a plumbing system given the
    vapor pressure (P) of a pure liquid or partial pressure (also P) of a
    substance from an initial solution headspace.
    
    
    Created on Wed May  3 15:13:46 2023

    Author: leefeinman

    Parameters
    ----------
    P : int, float, numpy array
        Partial pressure or vapor pressure of pure liquid (see function description).
    Punit : str, optional
        Accepts "mbar" (default), "torr", "atm", and "Pa".
    Cunit : str, optional
        Unit of concentration desired. Accepts "ppm", "ppb" (default), and "ppt".
    Psystem : int, float, numpy array, optional
        Pressure of the tubing where the substance will be introduced. Use the 
        same units as Punit. Default pressure of tubing is atmospheric pressure
        (1013.25 mbar).

    Raises
    ------
    ValueError
        Error raised when units are not properly provided.

    Returns
    -------
    TYPE
        int, float, numpy array.
    
    """
    
    
    if Punit == 'Pa':
        P = P * 0.01
        Psystem = Psystem * 0.01
        
    elif Punit == 'atm':
        P = P * 1013.25
        Psystem = Psystem * 1013.25
        
    elif Punit == 'torr':
        P = P * (1013.25/760)
        Psystem = Psystem * (1013.25/760)
        
    else:
        None
        
        
    if Cunit == 'ppm':
        factor = 10 **6
        
    elif Cunit == 'ppb':
        factor = 10 **9
        
    elif Cunit == 'ppt':
        factor = 10 **12

    else:
        raise ValueError("Provided units or unit formats are not supported. " +
                         "See documentation for correct inputs."
                         )
        
        
        
        
    return (P / Psystem) * factor

#%% raoults (finds the partial pressure of substance above its liquid mixture)

def raoults(pvapsolute, conc, mw_solute, mw_solvent = 18.02, density = 1):
    """
    Using the vapor pressure of a pure substance, this function uses Raoult's
    law to calculate the the partial pressure in headspace of a solute in a
    bicomponent mixture given: the solute's concentration in mg/L, the
    molecular weights of solvent & solute (g/mol), and the solution's density. 
    
    Default mw of solvent is water's. Default relative density (to water)
    of solvent is 1 (g/mL).
    
    Created on Fri May  5 15:19:57 2023

    Author: leefeinman
    
    
    Parameters
    ----------
    pvapsolute : int, float, numpy array
        Vapor pressure of the solute.
    conc : int, float, numpy array
        Concentration of solute in mg/L
        (miligrams solute per liters of solution).
    mw_solute : int, float, numpy array
        Molecular weight of solute.
    mw_solvent : int, float, numpy array | optional
        Molecular weight of solvent. The default is water: 18.02.
    density : int, float, numpy array | optional
        Relative density of the solution compared to water. This is uniteless.
        The default is water: 1.

    Returns
    -------
    int, float, numpy array
        Returns the partial pressure.

    """
 
    
    
    mol_solute = (conc * 10**-3) / mw_solute
    mol_solvent = 1000 * density / mw_solvent
    
    chi_solute = mol_solute / (mol_solute + mol_solvent)
    
    
    
    return pvapsolute * chi_solute

