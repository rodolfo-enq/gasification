#!/usr/bin/python
import feedstock
import gasifier
import numpy as np

#===============================================================================
# predefine variables
#===============================================================================
fuel = ['HighAshCoal','Average']
frac = feedstock.fraction(fuel)['mass']
my_species = ['C(gr)','H2','CO','CH4','CO2','N2']#,'H2O']
one_bar = 0.986923267 # atm

#==============================================================================
# Fig. 4.2, Chap. 4, RODRIGUES (2015), URL: http://hdl.handle.net/10183/140478
#==============================================================================
gasifier.coprocessing(frac,
                      fuel_id=fuel, moisture=np.array([0.15]),
                      blend=np.array([1.0]), # 100% biomassa (average)
                      T=np.array([1000.0, 1100.0]) - 273.15, # degC
                      P=one_bar, # atm
                      ER=np.linspace(0.0, 1.0, 21),
                      db='y', normalized='y',
                      species=my_species)
