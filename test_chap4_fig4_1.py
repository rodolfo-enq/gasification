#!/usr/bin/python
#==============================================================================
# import libraries/files
#==============================================================================
import feedstock
import pp
import gasifier
import numpy as np

#===============================================================================
# predefine variables
#===============================================================================
fuel1 = ['SMC','PS2']
fuel2 = ['SMC','RS']
frac1 = feedstock.fraction(fuel1)['mass']
frac2 = feedstock.fraction(fuel2)['mass']
my_species = ['CO2','H2','N2','CO','CH4']
ns = len(my_species)

#==============================================================================
# Tab. 2, LI; ZHANG; BI (2010), doi: 10.1016/j.ijhydene.2009.04.046
# Fig. 4.1, Chap. 4, RODRIGUES (2015), URL: http://hdl.handle.net/10183/140478
#==============================================================================
# SMC-PS blending
# Run: 1, 3, 9, 10
gasifier.coprocessing1(frac1,
                       fuel_id=fuel1, moisture=np.array([0, 0, 0, 0]),
                       blend=np.array([1.0/3.0, 1.0/3.0, 0.6/3.0, 1.0/3.0]),
                       T=np.array([990.0, 1021.0, 1045.0, 953.0]), # C
                       P=1.0, # Pa
                       ER=np.array([0.30, 0.42, 0.37, 0.34]),
                       SR=np.array([0.41, 0.41, 0.48, 0.51])/ \
                                   (pp.Mw[pp.i_H2O]/pp.Mw[pp.i_C]),
                       species=my_species)
# SMC-RS blending
# Run: 12
gasifier.coprocessing1(frac2,
                       fuel_id=fuel2, moisture=0,
                       blend=np.array([0.6/3.0]),
                       T=np.array([940.0]), # C
                       P=1.0, # Pa
                       ER=np.array([0.32]),
                       SR=np.array([0.49])/(pp.Mw[pp.i_H2O]/pp.Mw[pp.i_C]),
                       species=my_species)