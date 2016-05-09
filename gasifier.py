#!/usr/bin/env python

"""This script defines functions to equilibrium simulation of gasification 
processes. It uses some predefined functions from Cantera package.

@author = Rodolfo Rodrigues
@contact = rodolfo@unipampa.edu.br
@data = April, 2012, rev.: June, 2013 (adapted to use cython Cantera)
"""
#==============================================================================
# import libraries/files
#==============================================================================
import pp
import feedstock
import cantera as ct
import numpy as np
import scipy.optimize as opt
import csv

#==============================================================================
# predefine parameters
#==============================================================================
R = ct.gas_constant  # 8314.4621 Pa*m^3/K/kmol
Tn = 273.15  # K
Pn = ct.one_atm  # 101315 Pa
zero = np.zeros(1)
one = np.ones(1)

#==============================================================================
# special functions
#==============================================================================
def get_feed(self, moist=0, fuel=1.0, air=0, o2=0, stm=0):
    """
    This function creates a mixture of phases to denote the fuel.
    The fuel is composed as a mixture of char, gas, ash, and moisture phases.

    Parameters
    ----------
    self : ndarray
        Mass fraction of fuel compounds, d.b. [kg/kg]
    moist : float
        Mass fraction of moisture fuel [kg/kg]
    fuel : float
        Mass amount of fuel, d.b. [kg]
    air : float
        Mass amount of air [kg]
    o2 : float
        Mass amount of pure O2 [kg] (default value is zero)
    stm : float
        Mass amount of steam [kg] (default value is zero)

    Returns
    -------
    feed : object
        Feedstock object [mixture of phases]
    """
    # convert everything to array
    moist *= one
    fuel *= one
    air *= one
    o2 *= one
    stm *= one
    # preallocate variables
    no = np.zeros(pp.nsp)
    # mass amount of fuel, w.b.
    mass = fuel*(1 + moist)*np.append(self*(1 - moist), moist)
    ## NOTE: It's not possible to estimate the molecular weight of a fuel
    ## starting from its mass fraction composition. This parameter is taken
    ## from the whole-number multiple and empirical formula.
    ## attribute values for species
    #
    # mole amount of fuel, w.b.
    mol = mass/pp.Mw_f
    ## attribute values for species
    # mole amount of CHONSCl content
    no[pp.i_C] = mol[0]
#    no[pp.i_C_] = 0.3*mol[0]
    no[pp.i_H] = mol[1]
    no[pp.i_O] = mol[2]
    no[pp.i_N] = mol[3]
    no[pp.i_S] = mol[4]
    no[pp.i_Cl] = mol[5]
    # mole amount of ash content 
    no[pp.i_SiO2] = mol[6]
    no[pp.i_CaO] = mol[7]
    no[pp.i_Al2O3] = mol[8]
    no[pp.i_Fe2O3] = mol[9]
    no[pp.i_Na2O] = mol[10]
    no[pp.i_K2O] = mol[11]
    no[pp.i_MgO] = mol[12]
    no[pp.i_P2O5] = mol[13]
    no[pp.i_TiO2] = mol[14]
    no[pp.i_SO3] = mol[15]
    no[pp.i_Cr2O3] = mol[16]
    # mole amount of moisture content
    no[pp.i_H2O] = mol[17]
    # mole amount of air content 
    # air composition: 23.2%wt O2, 75.47%wt N2, 1.2%wt Ar
    if (o2.all() == 0 and air.any() != 0):
        no[pp.i_O2] = 0.23211606*air/pp.Mw[pp.i_O2]
        no[pp.i_N2] = 0.75507754*air/pp.Mw[pp.i_N2]
        no[pp.i_Ar] = 0.01280640*air/pp.Mw[pp.i_Ar]
    elif (o2.any() != 0 and air.all() == 0):
        no[pp.i_O2] = o2/pp.Mw[pp.i_O2]
    # mole amount of steam content
    no[pp.i_H2O] += stm/pp.Mw[pp.i_H2O]
    ## attribute values for phase
    # mole amount to each phase
    no_s = np.sum(no[:pp.s.n_species]) # solid phase
    no_g = np.sum(no[pp.s.n_species:]) # gas phase
    # set mole amount to each phase
    pp.f.set_phase_moles(pp.f.phase_index('solid'), no_s)
    pp.f.set_phase_moles(pp.f.phase_index('gas'), no_g)
    # set mole amount to each species
    pp.f.species_moles = no
    return pp.f

def get_water(self, moisture=0, fuel=1.0, steam=0):
    """
    This function get the mole amount of water as moisture and steam.

    Parameters
    ----------
    self : ndarray
        Mass fraction of fuel compounds, d.b. [kg/kg]
    moisture : float
        Mass fraction of moisture fuel [kg/kg]
    fuel : float
        Mass amount of fuel, d.b. [kg]
    steam : float
        Mass amount of steam [kg] (default value is zero)

    Returns
    -------
    mole_moisture : float
        Mole amount of moisture [kmol]
    mole_steam : float
        Mole amount of steam [kmol]
    """
    # convert everything to array
    moisture *= one
    fuel *= one
    steam *= one
    # mass amount of fuel, w.b.
    mass = fuel*(1 + moisture)*np.append(self*(1 - moisture), moisture)
    # mole amount of fuel, w.b.
    mol = mass/pp.Mw_f
    # mole amount of moisture content
    mole_moisture = mol[17]
    # mole amount of steam content
    mole_steam = steam/pp.Mw[pp.i_H2O]
    return mole_moisture, mole_steam

def get_enthalpy(self, value='h', duty=0):
    '''
    Return enthalpy (h) and specific heat capacity (cp) of a mixture of phases.
    
    TODO: Add duty term to enthalpy calculation
    '''
    # enthalpy [J] per 1 kg of fuel
    h = (self.phase_moles(self.phase_index('solid')) \
        * self.phase(self.phase_index('solid')).enthalpy_mole \
        + self.phase_moles(self.phase_index('gas')) \
        * self.phase(self.phase_index('gas')).enthalpy_mole \
        )/sum(self.species_moles)
    if value == 'h': return h
    # specific heat capacity [J/kmol/K]
    cp = (self.phase_moles(self.phase_index('solid')) \
        * self.phase(self.phase_index('solid')).cp_mole \
        + self.phase_moles(self.phase_index('gas')) \
        * self.phase(self.phase_index('gas')).cp_mole \
        )/sum(self.species_moles)
    if value == 'cp': return cp
    return  h, cp
    
def equilibrate_tp(self, moisture, fuel, air, o2=zero, steam=zero,
                T=1273, P=ct.one_atm):
    """
    Isothermic multi-phase equilibrium calculation holding temperature and 
    pressure fixed.
    
    The enthalpy of feedstocks/reagents (fuel, air, and steam) does not matter 
    for calculation of products in this approach.

    Parameters
    ----------
    self : ndarray
        Mass fraction of fuel compounds in d.b. [kg/kg]
    moisture : float
        Mass fraction of moisture fuel [kg/kg]
    fuel : float
        Mass amount of fuel in d.b. [kg]
    air : float
        Mass amount of air [kg]
    o2 : float
        Mass amount of oxygen, O2 [kg] (default value is zero)
    steam : float
        Mass amount of steam [kg]
    T : float
        Temperature [K]
    P : float
        Pressure [Pa] (default = 1 atm)

    Returns
    -------
    inlet : float
        Mole amount of inlet species [kmol]
    outlet : float
        Mole amount of outlet species [kmol]
    """
    f = get_feed(self, moisture, fuel, air, o2, steam)
    ## save initial composition
    inlet = f.species_moles
    # set desired condition
    f.T = T
    f.P = P
    # calculate equilibrium
    f.equilibrate('TP')#, solver='vcs')#, estimate_equil=1)
    ## save final composition
    # mole amount
    outlet = f.species_moles
    # FIXME: That is not possible to use labels at phaseMoles function
    return outlet, inlet#, GasMoles, GasComposition

def simple_equilibrate_hp(self, moisture, fuel, air=zero, steam=zero, 
                          P=ct.one_atm, duty=0):
    """
    Adiabatic multi-phase equilibrium calculation holding enthalpy and 
    pressure fixed.
    
    Use `equilibrate_hp' function for nonconventional fuels.

    Parameters
    ----------
    self : ndarray
        Mass fraction of fuel compounds in d.b. [kg/kg]
    moisture : float
        Mass fraction of moisture fuel [kg/kg]
    fuel : float
        Mass amount of fuel in d.b. [kg]
    air : float
        Mass amount of air [kg]
    steam : float
        Mass amount of steam [kg]
    P : float
        Pressure [Pa] (default = 1 atm)
    duty : float
        Duty fraction of outlet energy (default = 0)
        Positive value means lost heat.

    Returns
    -------
    content : object
        Reactor state
    inlet : float
        Mole amount of inlet species [kmol]
    outlet : float
        Mole amount of outlet species [kmol]
    T : float
        Equilibrium temperature [K]
    """
    f = get_feed(self, moisture, fuel, air, steam)
    # save initial composition
    inlet = f.species_moles
    # get enthalpy
    H = f.H
    # set desired condition
    f.P = P
    if duty != 0: f.H = (1-duty)*H
    # calculate equilibrium
    f.equilibrate('HP') #, solver='vcs', max_iter=200, estimate_equil=-1)
    # save final composition
    outlet = f.species_moles
    T = f.T
    return {'content':f, 'outlet':outlet, 'T':T, 'inlet':inlet}

def equilibrate_hp(self, hfo, fuel, mw, moisture=zero, air=zero, steam=zero, 
                   P=ct.one_atm, duty=0, guess=None, solver=0, disp=0):
    '''
    Non-isothermic multi-phase equilibrium calculation holding enthalpy and 
    pressure fixed.
    
    Use `simple_equilibrate_hp' function for conventional fuels.

    Parameters
    ----------
    self : ndarray
        Mass fraction of fuel compounds in d.b. [kg/kg]
    moisture : float
        Mass fraction of moisture fuel [kg/kg]
    fuel : float
        Mass amount of fuel in d.b. [kg]
    mw : float
        Molecular weight of fuel in d.b. [kg/kmol]
    air : float
        Mass amount of air [kg]
    steam : float
        Mass amount of steam [kg]
    P : float
        Pressure [Pa] (default = 1 atm)
    duty : float
        Duty fraction of outlet energy (default = 0)
        Positive value means lost heat.
    guess : float
        Guess value of temperature for equilibrium calculations [K]
    solver : integer
        solver = 0, default calculation
        solver = 1, scipy calculation
    disp : integer
        Display status notification of calculation.
        Default = 0, no notification.

    Returns
    -------
    content : objet
        Reactor state    
    inlet : float
        Mole amount of inlet species [kmol]
    outlet : float
        Mole amount of outlet species [kmol]
    T : float
        Equilibrium temperature [K]
    '''
    f = get_feed(self, moisture, fuel, air, steam)
    mole_moisture, mole_steam = get_water(self, moisture, fuel, steam)
    # save initial composition
    inlet = f.species_moles
    # get moles of fuel
    mole_fuel = fuel/mw
    # get moles of air species
    mole_O2 = inlet[pp.i_O2]
    mole_N2 = inlet[pp.i_N2]
    mole_Ar = inlet[pp.i_Ar]    
    # inlet enthalpy [J/kmol]
    inlet_h = (mole_fuel*hfo + mole_moisture*(pp.Hfo_H2Ol + pp.H_vap) \
                + mole_O2*pp.Hfo_O2 + mole_N2*pp.Hfo_N2 + mole_Ar*pp.Hfo_Ar \
                + mole_steam*pp.H_vap)/(mole_fuel + mole_moisture + mole_O2 \
                + mole_N2 + mole_Ar + mole_steam)
    # use default guess value
    if guess == None: guess = pp.To
    # equilibrium calculation at T and P constant
    def equilibrate_tp(self, T, P):
        self.T = T
        self.P = P
        self.equilibrate('TP')
        return self
    # set phases
    f = equilibrate_tp(f, guess, P)
    # choose solver
    # 0: own solver (default) (adapted from CATON et al., 2009)
    # 1: scipy solver (scipy.optimize.minimize_scalar)
    if solver == 0:
        # default solver (adapted from CATON et al., 2009)
        # set parameters to iterative calculation
        dT = 50 # temperature increment
        tol = 0.01 # tolerance
        iters = 0 # initial iteration
        # first state
        # enthalpy and specific heat of outlet species
        outlet_h, outlet_cp  = get_enthalpy(f,'h,cp')
        # duty
        outlet_h = (1-duty)*outlet_h
        outlet_cp = (1-duty)*outlet_cp
        # define the error
        T_err0 = (outlet_h - inlet_h)/outlet_cp
        # iterative calculation
        # estimate equilibrium temperature and product composition
        while (abs(T_err0) > tol):
            guess += dT
            f = equilibrate_tp(f, guess, P)
            outlet_h, outlet_cp  = get_enthalpy(f,'h,cp')
            # duty
            outlet_h = (1-duty)*outlet_h
            outlet_cp = (1-duty)*outlet_cp
            T_err = (outlet_h - inlet_h)/outlet_cp
            if (cmp(T_err, 0) != cmp(T_err0, 0)): # verify change of sign
                guess -= dT # go back to previous temperature
                dT *= 0.5 # decrease increment
            else:
                # verify change of curve inclination after highest temperature
                if (abs(T_err) > abs(T_err0)):
                    dT *= -1 # change of increment sign
                T_err0 = T_err # update value!
            iters += 1 # counter
            if iters == 200:
                print 'maximum number of iterations reached'
                break
            if disp == 2: 
                print 'T = %4.2f, T_err = %0.4g, iters = %2.0f' %(guess,
                                                                  T_err,iters)
        if disp == 1:
            print 'T = %4.2f, T_err = %0.4g, iters = %2.0f' %(guess,
                                                              T_err,iters)
        T = f.T
        outlet = f.species_moles
    else:
        # alternative solver (it uses minimize_scalar method)
        def residual(x):
            # set phases
            f.T = x
            f.P = P
            f.equilibrate('TP')
            # outlet enthalpy [J/kmol] with duty source
            outlet_h  = (1-duty)*get_enthalpy(f,'h')
            return (outlet_h - inlet_h)**2
        # estimate equilibrium temperature
        res = opt.minimize_scalar(residual,method='bounded',bounds=(200,6000),
                                  bracket=(residual(1200),residual(3000)))
        # estimate equilibrium product composition
        T = res.x[0]
        f = equilibrate_tp(f, T, P)
        outlet = f.species_moles
    return {'content':f, 'outlet':outlet, 'T':T, 'inlet':inlet}

def get_fuel_db(self):
#    fuel = get_feed(self, zero, one, zero) # 1 kg of fuel in d.b.
    fuel = get_feed(self) # 1 kg of fuel in d.b.
    nsp = fuel.n_species
    sp = fuel.species_moles
    # initiate variables
    mol_of_C = 0
    mol_of_H = 0
    mol_of_O = 0
    mol_of_S = 0
    mol_of_Cl = 0
    mol_of_Si = 0
    mol_of_Ca = 0
    mol_of_Al = 0
    mol_of_Fe = 0
    mol_of_Na = 0
    mol_of_K = 0
    mol_of_Mg = 0
    mol_of_P = 0
    mol_of_Ti = 0
    mol_of_Cr = 0
    mol_of_Ar = 0
    mol = 0
    # count moles of C,H,O in fuel species
    # IMPORTANT: I have to count S, Cl and ash species for precise estimation 
    # of stoichiometric oxygen amount. This is important mainly for high ash
    # fuels
    for i in range(nsp):
        if sp[i] != 0:
#            if i != fuel.species_index('gas', 'H2O'):
#                if i != fuel.species_index('gas', 'CO2'):
#                    mol_of_C += sp[i] * fuel.n_atoms(i, 'C')
#                    mol_of_H += sp[i] * fuel.n_atoms(i, 'H')
#                    mol_of_O += sp[i] * fuel.n_atoms(i, 'O')
            mol_of_C += sp[i] * fuel.n_atoms(i, 'C')
            mol_of_H += sp[i] * fuel.n_atoms(i, 'H')
            mol_of_O += sp[i] * fuel.n_atoms(i, 'O')
            mol_of_S += sp[i] * fuel.n_atoms(i, 'S')
            mol_of_Cl += sp[i] * fuel.n_atoms(i, 'Cl')
            mol_of_Si += sp[i] * fuel.n_atoms(i, 'Si')
            mol_of_Ca += sp[i] * fuel.n_atoms(i, 'Ca')
            mol_of_Al += sp[i] * fuel.n_atoms(i, 'Al')
            mol_of_Fe += sp[i] * fuel.n_atoms(i, 'Fe')
            mol_of_Na += sp[i] * fuel.n_atoms(i, 'Na')
            mol_of_K += sp[i] * fuel.n_atoms(i, 'K')
            mol_of_Mg += sp[i] * fuel.n_atoms(i, 'Mg')
            mol_of_P += sp[i] * fuel.n_atoms(i, 'P')
            mol_of_Ti += sp[i] * fuel.n_atoms(i, 'Ti')
            mol_of_Cr += sp[i] * fuel.n_atoms(i, 'Cr')
            mol_of_Ar += sp[i] * fuel.n_atoms(i, 'Ar')
            mol += sp[i]
    # normalise per mole of fuel
    mol_of_C /= mol
    mol_of_H /= mol
    mol_of_O /= mol
    mol_of_S /= mol
    mol_of_Cl /= mol
    mol_of_Si /= mol
    mol_of_Ca /= mol
    mol_of_Al /= mol
    mol_of_Fe /= mol
    mol_of_Na /= mol
    mol_of_K /= mol
    mol_of_Mg /= mol
    mol_of_P /= mol
    mol_of_Ti /= mol
    mol_of_Cr /= mol
    mol_of_Ar /= mol
    # stoichiometric moles of oxygen per mole of fuel
    stoic = mol_of_C + 0.25*mol_of_H - 0.5*mol_of_O + mol_of_S \
            - 0.5*mol_of_Cl + mol_of_Si + 0.5*mol_of_Ca + 3/2*mol_of_Al \
            + 3/2*mol_of_Fe + 0.25*mol_of_Na + 0.25*mol_of_K + 0.5*mol_of_Mg \
            + 2.5*mol_of_P + mol_of_Ti + 3/2*mol_of_Cr
    if stoic < 0:   # FIXME: Figure out the issue of a negative stoic
                    # oxygen. This happens when there is a fuel with high
                    # oxygen content, that is, 
                    # 0.5*mol_of_O > mol_of_C + 0.25*mol_of_H
        stoic += 0.5*mol_of_O
    return fuel, stoic

def enthalpy_of_formation(self, hhv):
    '''
    Estimate the standard enthalpy of formation of fuel [J/kg] from higher 
    heating value and species composition.
    
    Parameters
    ----------
    self : ndarray

    Returns
    -------
    hfo : ndarray
        standard enthalpy of formation of fuel [J/kg]
    '''
    f, stoic = get_fuel_db(self)
    mol = f.species_moles # kmol
    Mw = sum(mol*pp.Mw)
    # standard enthalpy of formation [J/kg]
    return (mol[pp.i_C]*pp.Hfo_CO2 + mol[pp.i_H]/2*pp.Hfo_H2Ol \
            + mol[pp.i_N]*pp.Hfo_N2 + mol[pp.i_S]*pp.Hfo_SO2 \
            + mol[pp.i_Cl]*pp.Hfo_ClO + mol[pp.i_SiO2]*pp.Hfo_SiO2 \
            + mol[pp.i_CaO]*pp.Hfo_CaO + mol[pp.i_Al2O3]*pp.Hfo_Al2O3 \
            + mol[pp.i_Fe2O3]*pp.Hfo_Fe2O3 + mol[pp.i_Na2O]*pp.Hfo_Na2O \
            + mol[pp.i_K2O]*pp.Hfo_K2O + mol[pp.i_MgO]*pp.Hfo_MgO \
            + mol[pp.i_P2O5]*pp.Hfo_P2O5 + mol[pp.i_TiO2]*pp.Hfo_TiO2 \
            + mol[pp.i_SO3]*pp.Hfo_SO3 + mol[pp.i_Cr2O3]*pp.Hfo_Cr2O3 \
            - stoic*pp.Hfo_O2 + hhv*1e6*Mw)/mol[pp.i_C]
            
def mass_of_air(self, fuel, ER=1.0):
    fuel_db, stoic = get_fuel_db(self)
    mol_of_fuel = fuel * np.sum(self/pp.Mw_f[:-1])
    # mole amount of gasifying agent
    mol_of_air = ER * stoic * mol_of_fuel/0.21
    # mass amount of gasifying agent
    return mol_of_air * pp.Mw_air

def equivalence_ratio(self, fuel, air, o2=0):
    fuel_db, stoic = get_fuel_db(self)
    mol_of_fuel = fuel * np.sum(self/pp.Mw_f[:-1])
    if air!=0 and o2==0:
        mol_of_O2 = 0.21 * (air/pp.Mw_air)
    elif air==0 and o2!=0:
        mol_of_O2 = o2/pp.Mw[pp.i_O2]
    else:
        mol_of_O2 = 0.21 * (air/pp.Mw_air) + o2/pp.Mw[pp.i_O2]
    return mol_of_O2/(stoic * mol_of_fuel)

def steam_to_carbon_ratio(self, fuel, steam):
    mol = chon_moles(self, 0, fuel, 0, 0, 0)
    mol_of_C = mol[0]
    mol_of_steam = steam / pp.Mw[pp.i_H2O]
    return mol_of_steam / mol_of_C
    
def mass_of_steam(self, fuel, SR=0):
    mol = chon_moles(self, 0, fuel, 0, 0, 0)
    mol_of_C = mol[0]
    mol_of_steam = SR * mol_of_C
    return mol_of_steam * pp.Mw[pp.i_H2O]

def chon_moles(self, moist, fuel, air, o2, stm):
    f = get_feed(self, moist, fuel, air, o2, stm)
    nsp = f.n_species
    sp = f.species_moles
    # initiate variables
    mol_of_C = 0
    mol_of_H = 0
    mol_of_O = 0
    mol_of_N = 0
    # count moles of C,H,O in fuel species
    for i in range(nsp):
        if sp[i] != 0:
            mol_of_C += sp[i] * f.n_atoms(i, 'C')
            mol_of_H += sp[i] * f.n_atoms(i, 'H')
            mol_of_O += sp[i] * f.n_atoms(i, 'O')
            mol_of_N += sp[i] * f.n_atoms(i, 'N')
    return mol_of_C, mol_of_H, mol_of_O, mol_of_N
    
def ohc_ratio(self, moist, fuel, air, o2, stm):
    C, H, O, N = chon_moles(self, moist, fuel, air, o2, stm)
    return H/C, O/C

def gas_yield(self, basis='vol', db='y'):
    """
    Gas yield of reactor outlet.

    Parameters
    ----------
    self : ndarray
        Mole of products [kmol]
    basis : string
        Mole amount ('kmol')
        Mass amount ('kg')
        Normal volume amount ('Nm3')
        Normal condition at 273.15K and 1 atm.
    db : string
        Dry basis ('y', default) or wet basis ('n')
    

    Returns
    -------
    yield : float
        Syngas yield [kmol] [kg] [Nm3]
    """
    # mole of gas species
    mol = self[pp.s.n_species:]
    # wet basis
    if (db == 'n'):        
        if (basis == 'mole'):
            return np.sum(mol) - self[pp.i_N2]
        if (basis == 'mass'):
            return np.sum(mol*pp.Mw_g) - self[pp.i_N2]*pp.Mw[pp.i_N2]
        if (basis == 'vol'):
            return ((np.sum(mol) - self[pp.i_N2])*R*Tn)/Pn
    # dry basis
    if (db == 'y'):
        if (basis == 'mole'):
            return np.sum(mol) - self[pp.i_H2O] - self[pp.i_N2]
        if (basis == 'mass'):
            return np.sum(mol*pp.Mw_g) - self[pp.i_H2O]*pp.Mw[pp.i_H2O] \
                - self[pp.i_N2]*pp.Mw[pp.i_N2]
        if (basis == 'vol'):
            return ((np.sum(mol) - self[pp.i_H2O] - self[pp.i_N2])*R*Tn)/Pn

def get_species(self, species=[], eps=1e-6):
    '''
    Get a list of species which mole fractions in 'self' are higher than 'eps'.
    this function is useful to find a minimum number of species to handle out a
    chemical equilibrium problem.
    '''
    i = 1
    while i < pp.nsp:
        if self[i] > eps:
            species_name = pp.f.species_name(i)
            try:
                species.index(species_name)
            except:
                # exclude liquid species
                if 'L)' not in species_name:
                    species.append(species_name)
        i += 1
    return species
    
def get_fraction(self, species, normalized='n', db='n', eps=None):
    '''
    db : string
        Dry basis ('y') or wet basis ('n', default)
    '''
    ## TODO: Make available for mass fraction calculation
    idx = len(species)
    mole = np.zeros(idx, 'd')
    i = 0
    while i < idx:
        # get values
        try:
            mole[i] = self[pp.f.species_index('solid', species[i])]#/mole_solid
        except:
            mole[i] = self[pp.f.species_index('gas', species[i])]#/mole_gas
        if eps != None:
            # make small values as zero
            if mole[i] < eps:
                mole[i] = 0
        i += 1
    # convert mole amount to mole fraction
    mole /= sum(self)
    if db == 'y':
        mole *= (1 - self[pp.i_H2O])
    if normalized == 'y':
        # normalize values
        mole /= np.sum(mole)
    return mole
    
def h2co_ratio(self):
    h2 = self[pp.f.species_index('gas', 'H2')]
    co = self[pp.f.species_index('gas', 'CO')]
    return h2/co
    
def carbon_conversion(products, reagents):
    return (reagents[pp.i_C] - products[pp.i_C]) / reagents[pp.i_C]

def syngas_hhv(self, fuel_mass=1.0, basis='vol'):
    """
    Higher heating value of gas-phase products (syngas).

    Parameters
    ----------
    self : ndarray
        Mole of products [kmol]
    fuel : float
        Mass of fuel, w.b.
    basis : string
        HHV in mass fraction = 'w', mole fraction = 'm', 
        volume fraction = 'v' (default)

    Returns
    -------
    HHV : float
        Higher heating value in the respective basis (mass, mole, or volume), 
        d.b. [MJ/kg] [MJ/kmol] [MJ/Nm3]
    """
    ns = pp.nsp
    # preallocate variables
    sp = []
    hhv_i = np.zeros(ns) # will be nonzero to 'heating' species
    # find key species
    for i in range(ns):
        if (i == pp.f.species_index('gas','H2') or \
            i == pp.f.species_index('gas','CH4') or \
            i == pp.f.species_index('gas','CO') #or \
#            i == pp.f.species_index('gas','C2H6')
            ):
            sp = np.append(sp, pp.f.species_name(i))
            hhv_i[i] = pp.Hfo[i] + (pp.f.n_atoms(i,'C') \
            + 0.25*pp.f.n_atoms(i,'H'))*pp.Hfo[pp.i_O2] \
            - (pp.f.n_atoms(i,'C'))*pp.Hfo[pp.i_CO2] \
            # FIXME: liquid or gas water?
            - (0.5*pp.f.n_atoms(i,'H'))*pp.Hfo[pp.i_H2O] # [J/kmol]
    # higher heating value
    hhv = np.sum(self*hhv_i)*1e-6 # [MJ]
    if (basis == 'syngas mole'):
        return hhv/gas_yield(self, db='y', basis='mole') # d.b. [MJ/kmol]
    if (basis == 'syngas mass'):
        return hhv/gas_yield(self, db='y', basis='mass') # d.b. [MJ/kg]
    if (basis == 'fuel mass'):
        return hhv/fuel_mass # [MJ/kg]
    if (basis == 'syngas vol'):
        return hhv/gas_yield(self, db='y', basis='vol') # d.b. [MJ/Nm3]

def syngas_lhv(self, fuel_mass=1.0):
    """
    Lower heating value (LHV) of gas-phase products (syngas).

    Parameters
    ----------
    self : ndarray
        Mole of products [kmol]
    fuel : float
        Mass of fuel, w.b.
    basis : string
        LHV in mass fraction = 'w', mole fraction = 'm', 
        volume fraction = 'v' (default)

    Returns
    -------
    lhv : float
        Lower heating value [MJ/kg]
    """
    lhv_CO = 10.160*pp.Mw[pp.i_CO] # MJ/kmol
    lhv_CH4 = 49.855*pp.Mw[pp.i_CH4] # MJ/kmol
#    lhv_C2H6 = 47.208*pp.Mw[pp.i_C2H6] # MJ/kmol
    lhv_H2 = 120.092*pp.Mw[pp.i_H2] # MJ/kmol
    return (lhv_CO*self[pp.i_CO] + lhv_CH4*self[pp.i_CH4] \
#            + lhv_C2H6*self[pp.i_C2H6] 
            + lhv_H2*self[pp.i_H2])*(1 \
            - self[pp.i_H2O]/gas_yield(self, db='n', basis='mole'))
        
def gas_hhv(self, basis='vol'):
    """
    Higher heating value of gas-phase products (fuel gas).

    Parameters
    ----------
    self : ndarray
        Mole of products [kmol]
    basis : string
        HHV in mass fraction = 'w', mole fraction = 'm', 
        volume fraction = 'v' (default)

    Returns
    -------
    HHV : float
        Higher heating value in the respective basis (mass, mole, or volume), 
        d.b. [MJ/kg] [MJ/kmol] [MJ/Nm3]
    """
    ns = pp.nsp
    # preallocate variables
    sp = []
    hhv_i = np.zeros(ns) # will be nonzero to 'heating' species
    # find 'heating' species
    for i in range(ns):
        if (i == pp.f.species_index('gas','H2') or \
            i == pp.f.species_index('gas','CO')):
            # Combustion of hydrogen
            # H2 + 0.5O2 --> H2O + <<HHV>>
            # Combustion of carbon monoxide
            # CO + 0.5O2 --> CO2 + <<HHV>>
            sp = np.append(sp, pp.f.species_name(i))
            hhv_i[i] = pp.Hfo[i] + (pp.f.n_atoms(i,'C') \
            + 0.25*pp.f.n_atoms(i,'H'))*pp.Hfo[pp.i_O2] \
            - (pp.f.n_atoms(i,'C'))*pp.Hfo[pp.i_CO2] \
            # FIXME: liquid or gas water?
            - (0.5*pp.f.n_atoms(i,'H'))*pp.Hfo[pp.i_H2O] # [J/kmol]
        if (pp.f.n_atoms(i,'C') >= 1 and pp.f.n_atoms(i,'H') >= 1):
            if (pp.f.n_atoms(i,'N') == 0 and pp.f.n_atoms(i,'O') == 0 and \
                pp.f.n_atoms(i,'S') == 0):
                # Combustion of hydrocarbons
                # CxHy + (x+0.25y)O2 --> xCO2 + 0.5yH2O + <<HHV>>
                sp = np.append(sp, pp.f.species_name(i))
                hhv_i[i] = pp.Hfo[i] + (pp.f.n_atoms(i,'C') \
                + 0.25*pp.f.n_atoms(i,'H'))*pp.Hfo[pp.i_O2] \
                - (pp.f.n_atoms(i,'C'))*pp.Hfo[pp.i_CO2] \
                # FIXME: liquid or gas water?
                - (0.5*pp.f.n_atoms(i,'H'))*pp.Hfo[pp.i_H2O] # [J/kmol]
    ## N2 H2 CO CH4 CO2 C2H6
    # higher heating value
    hhv = np.sum(self*hhv_i)*1e-6 # [MJ]
    if (basis == 'mole'):
        return hhv/gas_yield(self, db='y', basis='mole') # d.b. [MJ/kmol]
    if (basis == 'mass'):
        return hhv/gas_yield(self, db='y', basis='mass') # d.b. [MJ/kg]
    if (basis == 'vol'):
        return hhv/gas_yield(self, db='y', basis='vol') # d.b. [MJ/Nm3]

def mass_to_mole_fraction(self, Mw1, Mw2):
    """
    Convert mass fraction to mole fraction for dual-fuel blends.
    
    Parameters
    ----------
    self : ndarray
        Mass fraction of fuel #1 [kg/kg]
    Mw1 : float
        Molecular weight of fuel #1 [kg/kmol]
    Mw2 : float
        Molecular weight of fuel #2 [kg/kmol]
    
    Returns
    -------
    mole_fraction : ndarray
        Mole fraction of fuel #1 [kmol/kmol]
    """
    idx = len(self)
    if (self.ndim == 1):
        mole_fraction = self/Mw1/(self/Mw1 + (1.0 - self)/Mw2)
    else:
        mole_fraction = np.zeros(idx,'d')
        for i in range(idx):
            mole_fraction[i] = self[i]/Mw1/(self[i]/Mw1 + (1 - self[i])/Mw2)
    return mole_fraction

def mole_to_mass_fraction(self, Mw1, Mw2):
    """
    Convert mole fraction to mass fraction for dual-fuel blends.
    
    Parameters
    ----------
    self : ndarray
        Mole fraction of fuel #1 [kmol/kmol]
    Mw1 : float
        Molecular weight of fuel #1 [kg/kmol]
    Mw2 : float
        Molecular weight of fuel #2 [kg/kmol]
    
    Returns
    -------
    mass_fraction : ndarray
        Mass fraction of fuel #1 [kg/kg]
    """
    idx = len(self)
    if (self.ndim == 1):
        mass_fraction = self*Mw1/(Mw2 - self(Mw1 - Mw2))
    else:
        mass_fraction = np.zeros(idx,'d')
        for i in range(idx):
            mass_fraction[i] = self[i]*Mw1/(Mw2 - self[i]*(Mw1 - Mw2))
    return mass_fraction

def mixture(f, prop1, prop2):
    n1 = np.size(f)
    if (prop1.ndim <= 0):
        prop3 = np.zeros((n1))
        for i in range(n1):
            prop3[i] = f[i]*prop1 + (1.0 - f[i])*prop2
    else:
        n2 = len(prop1)
        prop3 = np.zeros((n1,n2))
        for i in range(n1):
            for j in range(n2):
                prop3[i,j] = f[i]*prop1[j] + (1.0 - f[i])*prop2[j]
    return prop3

def blending(f, coal, biomass):
    """
    f : float
        %wt biomass in coal-biomass blend
    """
    return (1.0 - f)*coal + (f)*biomass

def avg_error(mes, sim):
    """
    Return average error
    sim : ndarray
        simulated values
    mes: ndarray
        mesuared values
    """
    return np.sum(np.abs(sim-mes)/mes)/len(mes)

def cold_gas_efficiency(self, fuel_lhv, moisture_fuel):
    """
    Return cold gas efficiency of gasification.
    fuel_lhv : ndarray
        Fuel LHV
    moisture_fuel : ndarray
        Fuel moisture
    """
    return (syngas_lhv(self, 1 + moisture_fuel)/fuel_lhv)    

def coprocessing(self, fuel_id, blend, moisture, T, P=1.0,
                 air=0, O2=0, ER=0.4, steam=0, SR=0,
                 small=None, db='n', normalized='n', format_='%',
                 species=['C(gr)','N2','O2','H2','CO','CH4','CO2','H2O']):
    """
    Cogasification calculations for binary blends of fuels.
    
    Parameters
    ----------
    self : ndarray
        Mass fraction of fuel #1 [kg/kg]
    fuel_id : list of strings
        List of ID fuel
    blend : float|ndarray
        Fuel #1 to fuel #2 ratio [kg/kg]
    moisture : float|ndarray
        Moisture mass fraction [kg/kg]
    T : float|ndarray
        Temperature [degC]
    P : float|ndarray
        Pressure [atm] (default is 1.0)
    air : float|ndarray
        Air amount [kg] (default is zero)
    O2 : float|ndarray
        O2 amount [kg] (default is zero)
    ER : float|ndarray
        Equivalence ratio [kmol/kmol]
    steam : float|ndarray
        Steam amount [kg] (default is zero)
    SR : float|ndarray
        Steam to carbon ratio [kmol/kmol] (default is zero)
        basis: 1 kg coal-biomass blend, d.b.
    small : float
        Smallest number to report as a fraction value (default is None)
    db : string
        Dry basis composition ('y') or web basis composition ('n') (default
        is 'n')
    normalized : string
        Normalized compostion ('y') or overall composition ('n') (default
        is 'n')
    format_ : string
        Percentual ('%') or 'ppm' compostion (default is '%')
    species : list of strings
        List of chemical species.
        Default is C(gr), N2, O2, H2, CO, CH4, CO2, H2O
    
    Returns
    -------
    file : csv
        Function return a CSV file as following data: %wt biomass ratio 
        (assuming 1st fuel as coal), %wt moisture, T (degC), P (atm), 
        equivalence ratio, steam-to-carbon ratio, O-to-C ratio, H-to-C ratio, 
        species mole fractions, H2-to-CO ratio, % carbon conversion, 
        gas yield (Nm3/kg), HHV (MJ/kg), % cold gas efficiency        
    """
    # convert all values to array
    blend *= one
    moisture *= one
    T *= one
    P *= one
    air *= one
    O2 *= one
    ER *= one
    steam *= one
    SR *= one
    # default values
    steam_ = 0
    SR_ = 0
    air_ = 0
    o2_ = 0
    ER_ = 0
    # get number of points
    n_0 = np.size(fuel_id)
    n_1 = np.size(blend)
    n_2 = np.size(moisture)
    n_3 = np.size(T)
    n_4 = np.size(P)
    
    if np.size(air) > 1:
        n_5 = np.size(air)
    elif np.size(O2) > 1:
        n_5 = np.size(O2)
    elif np.size(ER) > 1:
        n_5 = np.size(ER)
    else:
        n_5 = 1
        
    if np.size(steam) > 1:
        n_6 = np.size(steam)
    elif np.size(SR) > 1:
        n_6 = np.size(SR)
    else:
        n_6 = 1
        
    if format_ == 'ppm':
        ft = 1e6
    else:
        ft = 1e2
#    # start count minimum number of species
#    minimum_species = []
    # start calculations
    for i in range(n_0-1): # asssumed 1st fuel as coal
        csvfile = open(str(fuel_id[0]) + '-' + str(fuel_id[i+1]) + '.csv','w')
        f = csv.writer(csvfile)
        f.writerow(['% BR','% MC','T (C)','P (atm)','ER','SR','O/C',
                    'H/C'] + species + ['H2/CO','% CC','Y (Nm3/kg)',
                    'HHV (MJ/kg)','% CGE'])
        for j in range(n_1): # %coal-biomass blend
            frac = blending(blend[j], self[0,:], self[i+1,:])
            for k in range(n_2): # moisture
                # get lhv to each moisture content of fuels
                fuel_lhv = feedstock.heating_values(fuel_id,moisture[k])['LHV']
                for l in range(n_3): # temperature
                    for m in range(n_4): # pressure
                        for o in range(n_5): # equivalence ratio
                            if air.any() != 0:
                                air_ = air[o]
                                o2_ = 0
                                ER_ = equivalence_ratio(frac, 1.0, air[o])
                            elif O2.any() != 0:
                                air_ = 0
                                o2_ = O2[o]
                                ER_ = equivalence_ratio(frac, 1.0, 0, O2[o])
                            elif ER.any() != 0:
                                air_ = mass_of_air(frac, 1.0, ER[o])
                                o2_ = 0
                                ER_ = ER[o]
                            for q in range(n_6): # steam-to-carbon ratio
                                if SR.any() != 0:
                                    steam_ = mass_of_steam(frac, 1.0, SR[q])
                                    SR_ = SR[q]
                                elif steam.any() != 0:
                                    steam_ = steam[q]
                                    SR_ = steam_to_carbon_ratio(frac, 1.0, 
                                                                steam[q])
                                hc,oc = ohc_ratio(frac, moisture[k], 1.0, 
                                                  air_, o2_, steam_)
                                p,r = equilibrate_tp(frac, moisture[k], 1.0, 
                                                     air_, o2_, steam_, 
                                                     T[l]+273.15,
                                                     ct.one_atm*P[m])
                                fuel_lhv_ = blending(blend[j], fuel_lhv[0], 
                                                     fuel_lhv[i+1])
                                syngas_lhv_ = syngas_lhv(p, 1 + moisture[k])
                                eff = syngas_lhv_/fuel_lhv_
                                hhv = syngas_hhv(p, basis='fuel mass', 
                                                 fuel_mass=1+moisture[k])
                                h2co = h2co_ratio(p)
                                cc = carbon_conversion(p,r)
                                y = gas_yield(p, basis='vol', db='y') # per kg
                                syngas = get_fraction(p, species, eps=small,
                                                      db=db, 
                                                      normalized=normalized)
                                f.writerow([100*blend[j], 100*moisture[k],
                                            T[l], P[m], ER_,
                                            SR_, oc, hc] + list(ft*syngas) 
                                            + [h2co, 100*cc, y, hhv, 100*eff])
#                                minimum_species = get_species(p, 
#                                                              minimum_species, 
#                                                              eps=1e-6)
        csvfile.close()
#        print minimum_species
        print 'Blend #' + str(i+1) + ' (' + str(fuel_id[0]) + '-' \
              + str(fuel_id[i+1]) + '): DONE'                    
    return None

def coprocessing1(self, fuel_id, blend, moisture, T, P=1.0,
                  air=0, O2=0, ER=0.4,
                  steam=0, SR=0,
                  small=None, db='n', normalized='n',
                  species=['C(gr)','N2','O2','H2','CO','CH4','CO2','H2O']):
    """
    Cogasification calculations for binary blends of fuels.
    
    Parameters
    ----------
    self : ndarray
        Mass fraction of fuel #1 [kg/kg]
    fuel_id : list of strings
        List of ID fuel
    blend : float|ndarray
        Fuel #1 to fuel #2 ratio [kg/kg]
    moisture : float|ndarray
        Moisture mass fraction [kg/kg]
    T : float|ndarray
        Temperature [degC]
    P : float|ndarray
        Pressure [atm] (default is 1.0)
    air : float|ndarray
        Air amount [kg] (default is zero)
    O2 : float|ndarray
        O2 amount [kg] (default is zero)
    ER : float|ndarray
        Equivalence ratio [kmol/kmol]
    steam : float|ndarray
        Steam amount [kg] (default is zero)
    SR : float|ndarray
        Steam to carbon ratio [kmol/kmol] (default is zero)
        basis: 1 kg coal-biomass blend, d.b.
    small : float
        Smallest number to report as a fraction value (default is None)
    db : string
        Get dry basis composition ('y') or web basis composition ('n') (default
        is 'y')
    normalized : string
        Get normalized compostion ('y') or overall composition ('n') (default
        is 'y')
    species : list of strings
        List of chemical species.
        Default is N2, O2, H2, CO, CH4, CO2, H2O
    
    Returns
    -------
    file : csv
        Function return a CSV file as following data: %wt biomass ratio 
        (assuming 1st fuel as coal), %wt moisture, T (degC), P (atm), 
        equivalence ratio, steam-to-carbon ratio, O-to-C ratio, H-to-C ratio, 
        species mole fractions, H2-to-CO ratio, % carbon conversion, 
        gas yield (Nm3/kg), HHV (MJ/kg), % cold gas efficiency        
    """
    # convert all values to array
    blend *= one
    moisture *= one
    T *= one
    air *= one
    O2 *= one
    ER *= one
    steam *= one
    SR *= one
    # get number of points
    n_0 = np.size(fuel_id)
    n_1 = np.size(blend)    
    # start calculations
    for i in range(n_0-1): # asssumed 1st fuel as coal
        csvfile = open(str(fuel_id[0]) + '-' + str(fuel_id[i+1]) + '.csv','w')
        f = csv.writer(csvfile)
        f.writerow(['% BR','% MC','T (C)','P (atm)','ER','SR','O/C',
                    'H/C'] + species + ['H2/CO','% CC','Y (Nm3/kg)',
                    'HHV (MJ/kg)','% CGE'])
        for j in range(n_1): # %coal-biomass blend
            frac = blending(blend[j], self[0,:], self[i+1,:])
            # get lhv to each moisture content of fuels
            fuel_lhv = feedstock.heating_values(fuel_id, moisture[j])['LHV']
            if air.any() != 0:
                air_ = air[j]
                o2_ = 0
                ER_ = equivalence_ratio(frac, 1.0, air[j])
            elif O2.any() != 0:
                air_ = 0
                o2_ = O2[j]
                ER_ = equivalence_ratio(frac, 1.0, 0, O2[j])
            elif ER.any() != 0:
                air_ = mass_of_air(frac, 1.0, ER[j])
                o2_ = 0
                ER_ = ER[j]
            else:
                air_ = 0
                o2_ = 0
                ER_ = 0
            if SR.all() != 0:
                steam_ = mass_of_steam(frac, 1.0, SR[j])
                SR_ = SR[j]
            elif steam.all() != 0:
                steam_ = steam[j]
                SR_ = steam_to_carbon_ratio(frac, 1.0, steam[j])
            else:
                steam_ = 0
                SR_ = 0
            hc, oc = ohc_ratio(frac, moisture[j], 1.0, air_, o2_, steam_)
            p, r = equilibrate_tp(frac, moisture[j], 1.0, air_, o2_, steam_, 
                                  T[j]+273.15, ct.one_atm*P)
            fuel_lhv_ = blending(blend[j], fuel_lhv[0], fuel_lhv[i+1])
            syngas_lhv_ = syngas_lhv(p, 1 + moisture[j])
            eff = 100*syngas_lhv_/fuel_lhv_
            hhv = syngas_hhv(p, basis='fuel mass', fuel_mass=1+moisture[j])
            h2co = h2co_ratio(p)
            cc = 100*carbon_conversion(p, r)
            y = gas_yield(p, basis='vol', db='y') # per kg
            syngas = get_fraction(p, species, eps=small, 
                                  db=db, normalized=normalized)
            f.writerow([100*blend[j], 100*moisture[j], T[j], P, ER_, SR_, 
                        oc, hc] + list(100*syngas) + [h2co, cc, y, hhv, eff])
        csvfile.close()
        print 'Blend #' + str(i+1) + ' (' + str(fuel_id[0]) + '-' \
              + str(fuel_id[i+1]) + '): DONE'                    
    return None
