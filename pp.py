#!/usr/bin/env python

"""This script sets thermophysical properties of phases and species.

The definition of phases must be done here. This file is necessary to 
other scripts.

@author = Rodolfo Rodrigues
@contact = rodolfo.rodrigues@ufsm.br
@data = April, 2012, rev.: June, 2013 (adapted to use cython Cantera)
"""
#==============================================================================
# import libraries
#==============================================================================
import cantera as ct
import numpy as np

#==============================================================================
# set parameters
#==============================================================================
To = 298.195 # K (reference temperature)
R = ct.gas_constant # 8314.47215 Pa*m^3/K/kmol

#==============================================================================
# load components to phases
#==============================================================================
s = ct.Solution('data/min-species.xml','solid')
g = ct.Solution('data/min-species.xml','gas')
#print s.species_names
air = ct.Solution('data/air.cti','air')
f = ct.Mixture([(s,1),(g,1)])
# preallocate variable
nsp = f.n_species

#==============================================================================
# get standard state properties
#==============================================================================
# molecular weight
Mw_s = s.molecular_weights
Mw_g = g.molecular_weights
Mw = np.concatenate((Mw_s, Mw_g))
# standard enthalpy of formation
Hfo_s = s.standard_enthalpies_RT*To*R
Hfo_g = g.standard_enthalpies_RT*To*R
Hfo = np.concatenate((Hfo_s, Hfo_g))
## find important indices
# fuel species
i_C = f.species_index('solid','C(gr)')
i_C_ = f.species_index('gas','C')
i_H = f.species_index('gas','H')
i_O = f.species_index('gas','O')
i_N = f.species_index('gas','N')
i_S = f.species_index('gas','S')
i_Cl = f.species_index('gas','CL')
# product species
i_Ar = f.species_index('gas','Ar')
i_H2 = f.species_index('gas','H2')
i_CH4 = f.species_index('gas','CH4')
#i_C2H6 = f.species_index('gas','C2H6')
i_H2O = f.species_index('gas','H2O')
i_CO = f.species_index('gas','CO')
i_CO2 = f.species_index('gas','CO2')
i_O2 = f.species_index('gas','O2')
i_N2 = f.species_index('gas','N2')
i_SO2 = f.species_index('gas','SO2')
i_ClO = f.species_index('gas','CLO')
#i_H2Ol = f.species_index('solid','H2O(L)')
# ash species
i_SiO2 = f.species_index('solid','SiO2(hqz)')
i_CaO = f.species_index('solid','CaO(s)')
i_Al2O3 = f.species_index('solid','AL2O3(a)')
i_Fe2O3 = f.species_index('solid','Fe2O3(s)')
i_Na2O = f.species_index('solid','Na2O(c)')
i_K2O = f.species_index('solid','K2O(s)')
i_MgO = f.species_index('solid','MgO(s)')
i_TiO2 = f.species_index('solid','TiO2(ru)')
i_Cr2O3 = f.species_index('solid','Cr2O3(s)')
i_P2O5 = f.species_index('gas','P2O5')
i_SO3 = f.species_index('gas','SO3')

## get specific properties
# standard enthalpy of formation (J/kmol)
Hfo_C     = Hfo[i_C]
Hfo_H2    = Hfo[i_H2]
Hfo_CH4   = Hfo[i_CH4]
Hfo_CO    = Hfo[i_CO]
Hfo_CO2   = Hfo[i_CO2]
Hfo_O2    = Hfo[i_O2]
Hfo_N2    = Hfo[i_N2]
Hfo_Ar     = Hfo[i_Ar]
Hfo_H2Og  = Hfo[i_H2O]
Hfo_SO2   = Hfo[i_SO2]
Hfo_ClO   = Hfo[i_ClO]
#Hfo_H2Ol  = Hfo[i_H2Ol]
Hfo_H2Ol  = -283970115.359
Hfo_SiO2  = Hfo[i_SiO2]
Hfo_CaO   = Hfo[i_CaO]
Hfo_Al2O3 = Hfo[i_Al2O3]
Hfo_Fe2O3 = Hfo[i_Fe2O3]
Hfo_Na2O  = Hfo[i_Na2O]
Hfo_K2O   = Hfo[i_K2O]
Hfo_MgO   = Hfo[i_MgO]
Hfo_P2O5  = Hfo[i_P2O5]
Hfo_TiO2  = Hfo[i_TiO2]
Hfo_SO3   = Hfo[i_SO3]
Hfo_Cr2O3 = Hfo[i_Cr2O3]
# water enthalpy of vaporization
H_vap     = 4.0667e7 # J/kmol
#H_vap_mass= H_vap/Mw[i_H2Ol] # J/kg
H_vap_mass= H_vap/Mw[i_H2O] # J/kg
Hfo_air   = 0.23211606*Hfo[i_O2] + 0.75507754*Hfo[i_N2] + 0.0128064*Hfo[i_Ar]
# molecular weight (kg/kmol)
Mw_f = [Mw[i_C], Mw[i_H], Mw[i_O], Mw[i_N], Mw[i_S], Mw[i_Cl],
        Mw[i_SiO2], Mw[i_CaO], Mw[i_Al2O3], Mw[i_Fe2O3],
        Mw[i_Na2O], Mw[i_K2O], Mw[i_MgO], Mw[i_P2O5],
        Mw[i_TiO2], Mw[i_SO3], Mw[i_Cr2O3],
#        Mw[i_H2Ol]]
        Mw[i_H2O]]
# air molecular weight (kg/kmol)
Mw_air = air.mean_molecular_weight
