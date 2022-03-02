#!/usr/bin/env python

"""This script defines functions to get values from a fuel database.

The fuel database is a predefined CSV file ('fuels.csv'). The file can be 
edited using any spreadsheet editor to add new fuel compounds.

@author = Rodolfo Rodrigues
@contact = rodolfo.rodrigues@ufsm.br
@data = April, 2012, rev.: June, 2013 (adapted to use cython Cantera)
"""

# FIXME: 2013/05/25 Not possible to use a list of single value as parameter to 
# functions.

#==============================================================================
# import libraries
#==============================================================================
import numpy as np
import csv
import pp
fuels = 'fuels.csv'

with open(fuels,'r') as f:
    myID = []
    rownum = 0
    for row in csv.reader(f):
        if rownum == 0:
            header = row
            typeIndex = header.index('Type')
            categoryIndex = header.index('Category')
            lmoistIndex = header.index('Lower moisture')
            hmoistIndex = header.index('Higher moisture')
            fcarbIndex = header.index('Fixed carbon')
            volmatIndex = header.index('Volatile matter')
            ashIndex = header.index('Ash')
            hhvIndex = header.index('HHV')
            lhvIndex = header.index('LHV')
            cellIndex = header.index('Cellulose')
            hcelIndex = header.index('Hemicellulose')
            ligIndex = header.index('Lignin')
            cIndex = header.index('C')
            hIndex = header.index('H')
            oIndex = header.index('O')
            nIndex = header.index('N')
            sIndex = header.index('S')
            clIndex = header.index('Cl')
            sio2Index = header.index('SiO2')
            caoIndex = header.index('CaO')
            al2o3Index = header.index('Al2O3')
            fe2o3Index = header.index('Fe2O3')
            na2oIndex = header.index('Na2O')
            k2oIndex = header.index('K2O')
            mgoIndex = header.index('MgO')
            p2o5Index = header.index('P2O5')
            tio2Index = header.index('TiO2')
            so3Index = header.index('SO3')
            cr2o3Index = header.index('Cr2O3')
        else:
            myID = np.append(myID,row[0])
        rownum += 1
    f.close()

def moisture(self):
    """
    Get the moisture fraction value for a list of fuels. 
    These fuels must be available in the database (file: 'fuels.csv').

    Parameters
    ----------
    self : list of string
        List of fuels

    Returns
    -------
    moist : float | array of float
        Moisture fraction [kg/kg]
    """
    nf = len(self) # number of fuels at the list (self)
    with open(fuels,'r') as f:
        idIndex = np.zeros((nf), dtype=object) # start value
        moist = np.zeros(nf) # start value
        for i in range(nf):
            # add 1 to take into account the header
            idIndex[i] = np.where(myID == self[i])[0] + 1
        rownum = 0
        for row in csv.reader(f):
            for i in range(nf):
                if rownum == idIndex[i]:
                    if row[cIndex] != '': # verify if values is empty
                        moist[i] = row[lmoistIndex]
                        break
            rownum += 1
        # define type of variable as scalar or array
        if nf == 1:
            return float(moist)/100
        else:
            return moist.astype('float64')/100
    f.close()   
    
    
def ash(self):
    """
    Get the ash fraction value for a list of fuels. 
    These fuels must be available in the database (file: 'fuels.csv').
    
    Default values are used if there are not available. Those values are
    chosen from VASSILEV et al. (2013) according to type and category of fuel.

    Parameters
    ----------
    self : string|list
        List of fuels

    Returns
    -------
    ash : float
        Ash fraction [kg/kg]
    composition : ndarray
        Composition of ash fraction [kg/kg]
    
    References
    ----------
    VASSILEV, S. V.; BAXTER, D.; ANDERSEN, L. K.; VASSILEVA, C. G. An overview 
    of the composition and application of biomass ash. Part 1. Phase-mineral 
    and chemical composition and classification. Fuel. v. 105, p. 40-76, 2013.
    """
    nf = len(self) # number of fuels at the list (self)
    with open(fuels,'r') as f:
        idIndex = np.zeros((nf), dtype=object) # start value
        frac = np.zeros(nf) # start value
        comp = np.zeros((nf,11))
        for i in range(nf):
            # add 1 due to header
            idIndex[i] = np.where(myID == self[i])[0] + 1
        rownum = 0
        for row in csv.reader(f):
            for i in range(nf):
                if rownum == idIndex[i]:
                    frac[i] = row[ashIndex]
                    # SiO2 CaO Al2O3 Fe2O3 Na2O K2O MgO P2O5 TiO2 SO3 Cr2O3
                    if row[sio2Index] == '' and row[caoIndex] == '':
                        if row[typeIndex] == 'Coal':
                            if row[categoryIndex] == 'Subbituminous':
                                ## mean compostion [Vassilev2013]
                                comp[i,:] = np.array([54.74, 7.05, 22.86, 5.30,
                                                       1.09, 1.67,  2.14, 0.08,
                                                       1.00, 4.07,     0])
                            else:
                                ## mean compostion [Vassilev2013]
                                comp[i,:] = np.array([54.06, 6.57, 23.18, 6.85,
                                                       0.82, 1.60,  1.83, 0.50,
                                                       1.05, 3.54,     0])
                        if row[typeIndex] == 'Biomass':
                            if row[categoryIndex] == 'Wood' or \
                               row[categoryIndex] == 'Charcoal':
                                ## mean compostion [Vassilev2013]
                                comp[i,:] = np.array([22.22, 43.03, 5.09, 3.44,
                                                       2.85, 10.75, 6.07, 3.48,
                                                       0.29,  2.78, 0])
                            if row[categoryIndex] == 'Straw':
                                ## mean compostion [Vassilev2013]
                                comp[i,:] = np.array([43.94, 14.13, 2.71, 1.42,
                                                       1.35, 24.49, 4.66, 4.13,
                                                       0.16,  3.01, 0])
                            if row[categoryIndex] == 'Grass':
                                ## mean compostion [Vassilev2013]
                                comp[i,:] = np.array([46.18, 11.23, 1.39, 0.98,
                                                       1.25, 24.59, 4.02, 6.62,
                                                       0.08,  3.66, 0])
                            if row[categoryIndex] == 'Agricultural':
                                ## mean compostion [Vassilev2013]
                                comp[i,:] = np.array([33.39, 14.86, 3.66, 3.26,
                                                       2.29, 26.65, 5.62, 6.48,
                                                       0.18,  3.61, 0])
                            else:
                                ## mean compostion [Vassilev2013]
                                comp[i,:] = np.array([29.76, 25.27, 5.51, 4.00,
                                                       2.48, 17.91, 5.42, 5.71,
                                                       0.66,  3.28, 0])
                    else:
                        comp[i,0] = row[sio2Index]
                        comp[i,1] = row[caoIndex]
                    if row[al2o3Index] != '':
                        comp[i,2] = row[al2o3Index]
                    if row[fe2o3Index] != '':
                        comp[i,3] = row[fe2o3Index]
                    if row[na2oIndex] != '':
                        comp[i,4] = row[na2oIndex]
                    if row[k2oIndex] != '':
                        comp[i,5] = row[k2oIndex]
                    if row[mgoIndex] != '':
                        comp[i,6] = row[mgoIndex]
                    if row[p2o5Index] != '':
                        comp[i,7] = row[p2o5Index]
                    if row[tio2Index] != '':
                        comp[i,8] = row[tio2Index]
                    if row[so3Index] != '':
                        comp[i,9] = row[so3Index]
                    if row[cr2o3Index] != '':
                        comp[i,10] = row[cr2o3Index]
                    break
            rownum += 1
        # define type of variable as scalar or array
        if nf == 1:
            frac = float(frac)/100
            comp = comp/comp.sum
        else:
            frac = frac.astype('float64')/100
            comp /= comp.sum(axis=1)[:, np.newaxis]
    f.close()
    return {'fraction':frac, 'composition':comp}
    
def fixed_carbon(self):
    """
    Get the fixed carbon fraction value for a list of fuels. 
    These fuels must be available in the database (file: 'fuels.csv').

    Parameters
    ----------
    self : string|list
        List of fuels

    Returns
    -------
    fixed_carb : float
        Fixed carbon fraction [kg/kg]
    """
    nf = len(self) # number of fuels at the list (self)
    with open(fuels,'r') as f:
        idIndex = np.zeros(nf) # start value
        fixedcarb = np.zeros(nf) # start value
        for i in range(nf):
            # add 1 to take into account the header
            idIndex[i] = np.where(myID==self[i])[0]+1
        rownum = 0
        for row in csv.reader(f):
            for i in range(nf):
                if rownum == idIndex[i]:
                    if row[cIndex] != '': # verify if values is empty
                        fixedcarb[i] = row[fcarbIndex]
                    break
            rownum += 1
        # define type of variable as scalar or array
        if nf == 1:
            return float(fixedcarb)/100
        else:
            return fixedcarb.astype('float64')/100
    f.close()

def volatile(self):
    nf = len(self) # number of fuels at the list (self)
    with open(fuels,'r') as f:
        idIndex = np.zeros(nf) # start value
        volat = np.zeros(nf) # start value
        for i in range(nf):
            # add 1 to take into account the header
            idIndex[i] = np.where(myID==self[i])[0]+1 
        rownum = 0
        for row in csv.reader(f):
            for i in range(nf):
                if rownum == idIndex[i]:
                    volat[i] = row[volmatIndex]
                    break
            rownum += 1
        # define type of variable as scalar or array
        if nf == 1:
            return float(volat)/100
        else:
            return volat.astype('float64')/100
    f.close()

def HHV(self):
    """
    Get the higher heating value for a list of fuels. 
    
    These fuels must be available in the database (file: 'fuels.csv').

    Parameters
    ----------
    self : string|list
        List of fuels

    Returns
    -------
    hhv : float
        Higher heating value, HHV [MJ/kg, d.b.]
    """
    nf = len(self) # number of fuels at the list (self)
    with open(fuels,'r') as f:
        idIndex = np.zeros(nf) # start value
        hhv = np.zeros(nf) # start value
        for i in range(nf):
            # add 1 to take into account the header
            idIndex[i] = np.where(myID==self[i])[0]+1 
        rownum = 0
        for row in csv.reader(f):
            for i in range(nf):
                if rownum == idIndex[i]:
                    if row[hhvIndex] == '':
                        hhv[i] = 0
                    else:
                        hhv[i] = row[hhvIndex]
                    break
            rownum += 1
        # define type of variable as scalar or array
        if nf == 1:
            hhv = float(hhv)
        else:
            hhv = hhv.astype('float64')
    f.close()
    return hhv

def LHV(self):
    """
    Get the lower heating value for a list of fuels.
    
    These fuels must be available in the database (file: 'fuels.csv').

    Parameters
    ----------
    self : list of string
        List of fuels

    Returns
    -------
    lhv : float
        Lower heating value, d.b., LHV [MJ/kg, d.b.]
    """
    nf = len(self) # number of fuels at the list (self)
    with open(fuels,'r') as f:
        idIndex = np.zeros(nf) # start value
        lhv = np.zeros(nf) # start value
        for i in range(nf):
            # add 1 to take into account the header
            idIndex[i] = np.where(myID==self[i])[0]+1 
        rownum = 0
        for row in csv.reader(f):
            for i in range(nf):
                if rownum == idIndex[i]:
                    if row[lhvIndex] != '': # verify if values is empty
                        lhv[i] = row[lhvIndex]
                        break
            rownum += 1
        # define type of variable as scalar or array
        if nf == 1:
            lhv = float(lhv)
        else:
            lhv = lhv.astype('float64')
    f.close()
    return lhv

def heating_values(self, moisture=0):
    '''
    d.b. if moisture=0
    '''
    n = np.size(self) # number of fuels
    lhv = LHV(self) # LHV of fuels
    hhv = HHV(self) # HHV of fuels
    h = H(self) # hydrogen mass fractions of fuels, d.b.
    if np.size(moisture) == 1:
        m = moisture*np.ones(n) # moisture fraction of fuels
    else:
        m = moisture # moisture fraction of fuels
    for i in range(n):
        if lhv[i] == 0 and hhv[i] == 0:
            lhv[i] = 20 # FIXME: Add LHV estimation calculation
            hhv[i] = 20 # FIXME: Add HHV estimation calculation
        elif lhv[i] != 0 and hhv[i] == 0:
            hhv[i] = (lhv[i] + 2.258*m[i])/(1 - m[i]) + 20.1790296248*h[i]
        elif lhv[i] == 0 and hhv[i] != 0:
            lhv[i] = (hhv[i] - 20.1790296248*h[i])*(1 - m[i]) - 2.258*m[i]
        lhv[i] *= 1 - m[i]
        hhv[i] *= 1 - m[i]
    return {'LHV': lhv, 'HHV': hhv}

def heat_of_formation(self, basis='kmol'):
    """
    Return the enthalpy (heat) of formation based on AspenTech's heat of 
    combustion-based correlation.

    Parameters
    ----------
    self : list of string
        List of fuels from CSV file
    basis : string
        Mole ('kmol', default) or mass ('kg') basis

    Returns
    -------
    ho : float
        Enthalpy of formation [J/kmol, J/kg]
    
    Reference
    ---------
    ASPEN TECHNOLOGY. Aspen physical property system: Physical property models. 
    version 7.3.2. Burlington, USA, 2012. (Chap. 4, p. 310)
    ZHU, Y. Evaluation of gas turbine and gasifier-based power generation 
    system. PhD-thesis. Depart. of Civil, Construction, and Environmental 
    Engineering. NC State University. Raleigh, USA, 2004. (Appendix A, p. 229)
    """
    ho = heating_values(self)['HHV'] - (1.418e6*H(self) + 3.278e6*C(self) 
        + 9.264e4*S(self) - 2.418e6*N(self) - 1.426e4*Cl(self))*1e4 # [kJ/kg]
    if basis == 'kg':
        return ho/1e3 # [J/kg]
    else:
        return ho*molecular_weight(self)/1e3 # [J/kmol]

def C(self):
    """
    Get the carbon fraction value for a list of fuels. 
    These fuels must be available in the database (file: 'fuels.csv').

    Parameters
    ----------
    self : list of string
        List of fuels

    Returns
    -------
    c : float
        Carbon fraction [kg/kg]
    """
    nf = len(self) # number of fuels at the list (self)
    with open(fuels,'r') as f:
        idx = np.zeros((nf), dtype=object) # start value
        c = np.zeros(nf) # start value
        for i in range(nf):
            # add 1 to take into account the header
            idx[i] = np.where(myID == self[i])[0] + 1 
        rownum = 0
        for row in csv.reader(f):
            for i in range(nf):
                if rownum == idx[i]:
                    if row[cIndex] != '': # verify if values is empty
                        c[i] = row[cIndex]
                        break
            rownum += 1
        # define type of variable as scalar or array
        if nf == 1:
            return float(c)/100
        else:
            return c.astype('float64')/100
    f.close()

def H(self):
    nf = len(self) # number of fuels at the list (self)
    with open(fuels,'r') as f:
        idx = np.zeros((nf), dtype=object) # start value
        h = np.zeros(nf) # start value
        for i in range(nf):
            # add 1 to take into account the header
            idx[i] = np.where(myID == self[i])[0] + 1
        rownum = 0
        for row in csv.reader(f):
            for i in range(nf):
                if rownum == idx[i]:
                    if row[hIndex] != '': # verify if values is empty
                        h[i] = row[hIndex]
                        break
            rownum += 1
        # define type of variable as scalar or array
        if nf == 1:
            return float(h)/100
        else:
            return h.astype('float64')/100
    f.close()

def O(self):
    nf = len(self) # number of fuels at the list (self)
    with open(fuels,'r') as f:
        idx = np.zeros((nf), dtype=object) # start value
        o = np.zeros(nf) # start value
        for i in range(nf):
            # add 1 to take into account the header
            idx[i] = np.where(myID == self[i])[0] + 1
        rownum = 0
        for row in csv.reader(f):
            for i in range(nf):
                if rownum == idx[i]:
                    if row[oIndex] != '': # verify if values is empty
                        o[i] = row[oIndex]
                        break
            rownum += 1
        # define type of variable as scalar or array
        if nf == 1:
            return float(o)/100
        else:
            return o.astype('float64')/100
    f.close()

def N(self):
    nf = len(self) # number of fuels at the list (self)
    with open(fuels,'r') as f:
        idx = np.zeros((nf), dtype=object) # start value
        n = np.zeros(nf) # start value
        for i in range(nf):
            # add 1 to take into account the header
            idx[i] = np.where(myID == self[i])[0] + 1 
        rownum = 0
        for row in csv.reader(f):
            for i in range(nf):
                if rownum == idx[i]:
                    if row[nIndex] != '': # verify if values is empty
                        n[i] = row[nIndex]
                        break
            rownum += 1
        # define type of variable as scalar or array
        if nf == 1:
            return float(n)/100
        else:
            return n.astype('float64')/100
    f.close()

def S(self):
    nf = len(self) # number of fuels at the list (self)
    with open(fuels,'r') as f:
        idx = np.zeros((nf), dtype=object) # start value
        s = np.zeros(nf) # start value
        for i in range(nf):
            # add 1 to take into account the header
            idx[i] = np.where(myID == self[i])[0] + 1
        rownum = 0
        for row in csv.reader(f):
            for i in range(nf):
                if rownum == idx[i]:
                    if row[sIndex] != '': # verify if values is empty
                        s[i] = row[sIndex]
                        break
            rownum += 1
        # define type of variable as scalar or array
        if nf == 1:
            return float(s)/100
        else:
            return s.astype('float64')/100
    f.close()

def Cl(self):
    nf = len(self) # number of fuels at the list (self)
    with open(fuels,'r') as f:
        idx = np.zeros((nf), dtype=object) # start value
        cl = np.zeros(nf) # start value
        for i in range(nf):
            # add 1 to take into account the header
            idx[i] = np.where(myID == self[i])[0] + 1
        rownum = 0
        for row in csv.reader(f):
            for i in range(nf):
                if rownum == idx[i]:
                    if row[clIndex] != '': # verify if values is empty
                        cl[i] = row[clIndex]
                        break
            rownum += 1
        # define type of variable as scalar or array
        if nf == 1:
            return float(cl)/100
        else:
            return cl.astype('float64')/100
    f.close()
    
def fraction(self):
    '''
    Get the mass fraction in d.b.
    '''
    nf = len(self) # number of fuels at the list (self)
    ## get values
    ash_ = np.array(ash(self)['fraction'])
    ash_i = np.array(ash(self)['composition'])
    ## preallocate variable
    massf = np.zeros((nf, 6))
    xmass = np.zeros((nf, 17))
    for i in range(nf):
        massf[i,:] = np.array(( C(self)[i], H(self)[i], 
                                O(self)[i], N(self)[i],
                                S(self)[i], Cl(self)[i] ))
    for i in range(nf):
        for j in range(6):
            massf[i,j] = (1.0 - ash_[i])*massf[i,j]
    ## normalised mass fraction, w.b. [kg/kg]
    for i in range(nf):
        xmass[i,:] = np.hstack((massf[i,:], ash_[i]*ash_i[i,:]))
    xmass /= xmass.sum(axis=1)[:, np.newaxis]
#    for i in range(nf):
#        xmass[i,:] = xmass[i,:]/np.sum(xmass[i,:]) # normalise
    ### TODO: Mole fraction calculation
    xmol = 0
    return {'mass':xmass, 'mole':xmol}
    
def molecular_weight(self):
    '''
    Get the mean molecular weight of fuel, in d.b.
    '''
    mw = pp.Mw_f[:-1]
    return 1/np.sum(fraction(self)['mass']/mw, axis=1)

def biochemical_composition(self, dist=[0.6, 0.8, 0.8]):
    """
    Get cellulose, hemicellulose, and lignin fraction values for a list of 
    fuels. These fuels must be available in the database (file: 'fuels.csv').
    
    If there are not those values they are estimated following the procedure 
    presented in CUOCI et al. (2007).
    
    TODO: The function must get the values from database if they are available. 

    Parameters
    ----------
    self : string
        ID of fuel according to 'fuels.csv' file.
    dist : ndarray
        Mass distribution of surrogate species in terms of cellulose, 
        hemicellulose and lignins [kg/kg]
        Default values are from Ranzi et al. (2008): 60%, 80% and 80%. 
        See supplemental material of Ranzi et al. (2008)

    Returns
    -------
    CELL : float
        Cellulose mass fraction [kg/kg]
    HCE : float
        Hemicellulose mass fraction [kg/kg]
    LIG-C : float
        Lignin-C mass fraction [kg/kg]
    LIC-H : float
        Lignin-H mass fraction [kg/kg]
    LIG-O : float
        Lignin-O mass fraction [kg/kg]
    
    References
    ----------
    Cuoci, A.; Faravelli, T.; Frassoldati, A.; Granata, S.; Migliavacca, G.; 
    Pierucci, S.; Ranzi, E. & Sommariva, S. A general mathematical model of 
    biomass devolatilization. Note 2. Detailed kinetics of volatile species.
    30th Meeting of the Italian Section of The Combustion Institute, 2007.
    Ranzi, E.; Cuoci, A.; Faravelli, T.; Frassoldati, A.; Migliavacca, G.; 
    Pierucci, S. & Sommariva, S. Chemical kinetics of biomass pyrolysis. Energ.
    Fuel., 2008, 22, 4292-4300
    """
    if dist == [0.6,0.8,0.8]:
        # (Cuoci et al., 2007)
        a = np.array([[0.44850160, 0.58942,    0.61653427],
                      [0.06171176, 0.05517644, 0.06825135],
                      [0.48978665, 0.35540356, 0.31521439]]) 
    else:
        # mass fraction
        cellu = np.array([0.44446117, 0.06216388, 0.49337496])
        hemic = np.array([0.45456224, 0.06103358, 0.48440417])
        ## (Cuoci et al., 2007)
        lig_c = np.array([0.677644,   0.05686658, 0.26548942]) 
        lig_h = np.array([0.60125683, 0.07109754, 0.32764563]) 
        lig_o = np.array([0.567364,   0.05475391, 0.37788209])
        # definition of surrogate species
        s1 = dist[0]*cellu + (1-dist[0])*hemic
        s2 = dist[1]*lig_o + (1-dist[1])*lig_c
        s3 = dist[2]*lig_h + (1-dist[2])*lig_c
        # matrix of CHO fractions in terms of s1,s2,s3 surrogate species    
        a = np.array([[s1[0], s2[0], s3[0]],
                      [s1[1], s2[1], s3[1]],
                      [s1[2], s2[2], s3[2]]])
    # get values of fuels
    c = np.array([C(self)])
    h = np.array([H(self)])
    o = np.array([O(self)])
    # CHO normalized mass fraction of fuel
    b = np.array([c,h,o])/sum(np.array([c,h,o]))
    # solve the problem
    x = np.linalg.solve(a,b)
    cell = dist[0]*x[0]
    hcel = (1-dist[0])*x[0]
    ligo = dist[1]*x[1]
    ligh = dist[2]*x[2]
    ligc = (1-dist[1])*x[1] + (1-dist[2])*x[2]
    return 'CELL:%7.5f, HCE:%7.5f, LIGC:%7.5f, LIGH:%7.5f, LIGO:%7.5f'\
            %(cell, hcel, ligc, ligh, ligo)
#    return {'CELL':float(cell), 'HCE':float(hcel), 'LIGC':float(ligc), 
#            'LIGH':float(ligh), 'LIGO':float(ligo)}
#    return np.array([float(cell),hcel,ligc,ligh,ligo])
