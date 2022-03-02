#!/usr/bin/env python

"""This script defines functions to chemical reactors. It uses source terms 
provided by Cantera package.

@author = Rodolfo Rodrigues
@contact = rodolfo.rodrigues@ufsm.br
@data = September, 2012, rev.: June, 2013 (adapted to use cython Cantera)
"""
#==============================================================================
# import libraries/files
#==============================================================================
import cantera as ct
import numpy as np
from scipy.integrate import ode
from time import clock
import csv

#==============================================================================
# predefine parameters
#==============================================================================
R = ct.gas_constant # 8314.47215 Pa*m^3/K/kmol
Tn = 273.15 # K
Pn = ct.one_atm # 101315 Pa
air = ct.Solution('data/air.cti')

#==============================================================================
# solver settings
#==============================================================================
## Tolerance for CHEMKIN EQUIL method
#abs_tol=1e-9 # absolute tolerance
#rel_tol=1e-6 # relative tolerance

## Tolerance for CHEMKIN PSR model
#abs_tol=1e-9 # absolute tolerance
#rel_tol=1e-4 # relative tolerance

# Set tolerance
abs_tol=1e-12 # absolute tolerance
rel_tol=1e-8 # relative tolerance

nstep = 300 # maximum number of steps
maxiter = 200 # maximum number of iterations
test_fails = 20 # maximium number of error test fails

def mixer(self, F):
    """
    This functions return a mixture of any number of streams 
    """
    ## get number of streams
    n_streams = np.size(F)
    ## verify if there are more than one stream
    if n_streams > 1: 
        ## get number of species of key stream
        n_species = self[0].n_species
        ## preallocate variable
        P = np.zeros(n_streams)
        h = 0.0
        ## set initial vector of mass (kg/s)
        mass = F[0]*self[0].Y # [kg/s][kg/kg]=[kg/s]        
        for j in range(n_streams):
            ## get pressure of streams
            P[j] = self[j].P
            ## get total energy of streams
            h += F[j]*self[j].enthalpy_mass # [kg/s][J/kg]=[J/s]
            if j > 0:
                for k in range(self[j].n_species):
                    for i in range(n_species):
                        if self[0].species_name(i) == self[j].species_name(k):
                            ## get total mass of each species
                            mass[i] += F[j]*self[j].Y[k] # [kg/s]        
        ## mix stream
        ## mass fraction
        mass = mass/sum(mass) # [kg/kg]
        ## specific enthalpy
        h = h/np.sum(F) # [J/s][s/kg]=[J/kg]
        ## pressure
        P = np.amin(P) # [Pa]
        ## mass flow
        F = np.sum(F)        
        ## redefine stream #1 as mix stream
        self = self[0]        
        ## re-calculate mix stream state
        self.Y = mass
        self.HP = h, P
    return self, F

def mole_of_atoms(self, F, mole0=None):
    mole = np.zeros(self.n_elements)
    # count mols of each elements
    for i in range(self.n_species):
        for j in range(self.n_elements):
            mole[j] += self.n_atoms(self.species_name(i), \
                                    self.element_name(j)) \
                        *(F/self.mean_molecular_weight*self.X[i])
    # verify element balance
    if mole0 != None and np.amax(np.abs((mole-mole0)/mole0)) > 0.01:
        print('Warning: Element balance is inconsistent! (dev>1.0%)')
    return mole

def equilibrium_reactor(self, F, P, T=None, energy='on', 
                        solver='auto', disp=0):
    self, F = mixer(self, F)
    if T == None and energy=='on':
        self.TP = None, P
        self.equilibrate('HP', maxsteps=nstep, rtol=rel_tol, solver=solver, 
                         loglevel=disp)
    else:
        self.TP = T, P
        self.equilibrate('TP', maxsteps=nstep, rtol=rel_tol, solver=solver, 
                         loglevel=disp)
    return self

def slope(y1, y2, dx):
    """
    Get the larger absolute slope value from 2 consecutive points
    """
    return np.max(np.abs((y2 - y1)/dx))

def species_conservation(self, Y0, F, V):
    Y = self.Y # mass fraction
    M = self.molecular_weights # kg/kmol
    rate = self.net_production_rates # kmol/m^3/s
    return F*(Y0 - Y) + rate*M*V
        
def energy_conservation(self, X0, H0):
    X = self.X # mole fraction
    H = self.partial_molar_enthalpies # J/kmol
    return X*H - X0*H0

def cstr0(self, F, P, V, T=None, dt=2.0, t=None, disp=0):
    # mix streams if there is more than one
    self, F = mixer(self, F)
    # set initial conditions
    self.TP = T, P
    # create a new reactor
    if T == None or T == 0:
        r = ct.IdealGasReactor(self, energy='on')
    else:
        r = ct.IdealGasReactor(self, energy='off')
    r.volume = V
    # create up and downstream reservoirs
    upstream = ct.Reservoir(self)
    downstream = ct.Reservoir(self)
    # set mass flow into reactor
    m = ct.MassFlowController(upstream, r, mdot=F)
    # set valve to hold pressure constant
    ct.PressureController(r, downstream, master=m, K=1e-5)
    # create reactor network
    s = ct.ReactorNet([r])
    s.rtol, s.atol, s.max_err_test_fails = rel_tol, abs_tol, test_fails
    # set parameters
    time = 0
    residual = 1
    all_done = False
    # forces reinitialization
    s.set_initial_time(0)
    while not all_done:
        Yo = r.thermo.Y
        To = r.thermo.T
        try:
            time += dt
            s.advance(time)
        except:
            dt *= 0.1
            time += dt
            s.advance(time)
        if t != None and time >= t: all_done = False
        if time > 10 * dt:
            all_done = True
            residual = slope(Yo, r.thermo.Y, dt)
            if T == None:
                residual = np.max(np.append(residual,slope(To,r.thermo.T,dt)))
            if residual > rel_tol:
                all_done = False
                break
        if disp == True: print('Residual: %1.3e' %(residual))
    return self, time, residual
    
def cstr(self, F, P, V, T=None, dt=10.0, t=None, disp=0):
    # mix streams if there is more than one
    self, F = mixer(self, F)
    # get initial conditions
    self.TP = T, P
    # create constant pressure reactor
    if T == None:
#        r = ct.IdealGasConstPressureReactor(contents=self, volume=V)
        r = ct.IdealGasReactor(contents=self, volume=V)
    else:
#        r = ct.IdealGasConstPressureReactor(contents=self, volume=V, 
#                                            energy='off')
        r = ct.IdealGasReactor(contents=self, volume=V, energy='off')
    # set simulation parameters
    s = ct.ReactorNet([r])
    s.rtol, s.atol, s.max_err_test_fails = rel_tol, abs_tol, test_fails
    time = 0
    done = 1
    residual = 1
    while done:
        Y0 = self.Y # mass fraction
        T0 = self.T # K
        # start simulation
        try:
            time += dt
            s.advance(time)
        except:
            dt *= 0.1
            time += dt
            s.advance(time)
        residual = slope(Y0, self.Y, dt)
        if T == None: 
            residual = np.max(np.append(residual, slope(T0, self.T, dt)))
        # stop simulation 
        # at final time
        if t != None and time >= t: done = 0
        # at steady-state (accumulation = zero)
        # or at very large time
        if residual <= abs_tol*1e3 or time >= 1e8: done = 0
        if disp == True: print('Residual: %1.3e' %(residual))
#    print(self.report())
    return self, time, residual

def cstr1(self, F, P, V, T=None, dt=10.0, t=1e20, 
          species = ['H2','H2O','CO','CO2','CH4','O2'], disp=0):
    # mix streams if there is more than one
    self, F = mixer(self, F)
    # get initial conditions
    self.TP = T, P
    # create constant pressure reactor
    if T == None:
#        r = ct.IdealGasConstPressureReactor(contents=self, volume=V)
        r = ct.IdealGasReactor(contents=self, volume=V)
    else:
#        r = ct.IdealGasConstPressureReactor(contents=self, volume=V, 
#                                            energy='off')
        r = ct.IdealGasReactor(contents=self, volume=V, energy='off')
    # set simulation parameters
    s = ct.ReactorNet([r])
    s.rtol, s.atol, s.max_err_test_fails = rel_tol, abs_tol, 20
    time = 0
    done = 1
    data = open('cstr1_dyn'+str(t)+'.csv','w')
    f = csv.writer(data)
    f.writerow(['T (C)','P (atm)','V (m3)','t (s)','residual'] + species)
    f.writerow([self.T-273.15, self.P/ct.one_atm, V, time, 0.0]
                    + list(self[species].X))
    while done:
        Y0 = self.Y # mass fraction
        T0 = self.T # K
        # start simulation
        try:
            time += dt
            s.advance(time)
        except:
            dt *= 0.1
            time += dt
            s.advance(time)
        residual = slope(Y0, self.Y, dt)
        if T == None: 
            residual = np.max(np.append(residual, slope(T0, self.T, dt)))
        f.writerow([self.T-273.15, self.P/ct.one_atm, V, time, residual] 
                    + list(self[species].X))
        # at final time
        if t != None and time >= t: done = 0
        if residual <= abs_tol: done = 0
        if disp == True: print('Residual: %1.3e' %(residual))
    data.close()
    return None

def get_contents(state):
    return np.append(np.array([state.TPY[0], state.TPY[1]]), state.TPY[2])

def set_contents(contents, i):
    return contents[i,0], contents[i,1], contents[i,2:]

def ss_isot_pfr(self, F, P, V, T, L, N):
    # mix streams if there is more than one
    self, F = mixer(self, F)
#    print (self.report())
    # preallocate variables
    if np.size(P) == 1: P *= np.ones(N)
    if np.size(T) == 1: T *= np.ones(N)
    contents = np.zeros((N, self.n_species+2))
    # define equations to solve
    def f(z, y, G, T, P):
        # z : independent variable (reactor length)
        # y : dependent variable (mass fraction)
        self.set_unnormalized_mass_fractions(y) # mass fraction
#        self.TPY = T, P, Y
        self.TP = T, P
        rate = self.net_production_rates # kmol/m^3/s
        M = self.molecular_weights # kg/kmol
        return rate*M/G
    # initial values
    y0 = self.Y # initial y
    z0 = 0.0 # initial z
    # get state at z0
    contents[0,:] = get_contents(self)
    print('N=%3i P=%4.2f T=%7.2f r=%1.3e' %(1, P[0]/Pn, T[0]-273.15, 
                                            np.sum(self.net_production_rates)))
    # choose integrator(solver)
    r = ode(f).set_integrator('vode', method='bdf', 
        atol=abs_tol, rtol=rel_tol, nsteps=nstep, with_jacobian=True)
    # set integrator parameters
    r.set_initial_value(y0, z0).set_f_params(F*L/V, T[0], P[0])
    # set aditional parameters
    length, dz = np.linspace(z0, L, N, retstep=True)
    i = 1 # start counter    
    # start integration
    while r.successful() and i < N:
        r.integrate(r.t+dz)
        # set new state for phase
        self.TPY = T[i], P[i], r.y
        # save values
        contents[i,:] = get_contents(self)
        print('N=%3i P=%4.2f T=%7.2f r=%1.3e' %(i+1, P[i]/Pn, T[i]-273.15, 
                                                np.sum(
                                                self.net_production_rates)))
        i += 1
    return self, length, contents

def ss_nonisot_pfr(self, F, P, V, L, N):
    # mix streams if there is more than one
    self, F = mixer(self, F) 
    # preallocate variables
    if np.size(P) == 1: P *= np.ones(N)
    contents = np.zeros((N, self.n_species+2))
#    contents = N*[0, 0, np.zeros(self.n_species)]
    # define equations to solve
    def f(z, y, G, P):
        # z : independent variable (reactor length)
        # y[0] : dependent variable (temperature)
        # y[1:] : dependent variable (mass fraction)
        self.set_unnormalized_mass_fractions(y[1:]) # mass fraction
        self.TP = y[0], self.P
        rate = self.net_production_rates # kmol/m^3/s
        M = self.molecular_weights # kg/kmol
        cp = self.cp_mass # J/kg/K
        h = self.partial_molar_enthalpies # J/kmol
        return np.append(-np.sum(rate*h)/(G*cp), rate*M/G)
    # initial values
    y0 = np.append(self.T, self.Y) # initial y
    z0 = 0.0 # initial t
    # get state at z=0   
    contents[0,:] = get_contents(self)
#    contents[0,:] = self.TPY
    print('N=%3i P=%4.2f T=%7.2f r=%1.3e' %(1, P[0]/Pn, self.T-273.15, 
                                            np.sum(self.net_production_rates)))
    # choose integrator(solver)
    r = ode(f).set_integrator('vode', method='bdf', 
        atol=abs_tol, rtol=rel_tol, nsteps=nstep, 
        with_jacobian=True)
    # set integrator parameters
    r.set_initial_value(y0, z0).set_f_params(F*L/V, P[0])
    # set aditional parameters
    length, dz = np.linspace(z0, L, N, retstep=True)
    i = 1 # start counter
    # start integration
    while r.successful() and i < N:
        r.integrate(r.t+dz)
        # set new state for phase
        self.TPY = r.y[0], P[i], r.y[1:]
        # save values
        contents[i,:] = get_contents(self)
#        contents[i,:] = self.TPY
        print('N=%3i P=%4.2f T=%7.2f r=%1.3e' %(i+1, P[i]/Pn, self.T-273.15,
                                                np.sum(
                                                self.net_production_rates)))
        i += 1
    return self, length, contents
    
def pfr(self, F, P, V, T=None, L=1.0, N=100):
    if T == None or T == 0:
        return ss_nonisot_pfr(self, F, P, V, L, N)
    else:
        return ss_isot_pfr(self, F, P, V, T, L, N)
        
def series_reactors(self, F, P, V, N, T=None, 
                    species=['CH4','CO','H2','CO2','H2O','O2'],
                    kinetic=1, reactor=0, L=1.0, ND=10, disp=0):
    '''
    T
    None -> Nonisothermal reactor
    kinetic
    0 -> equilibrium
    1 -> kinetic
    reactor
    0 -> cstr
    1 -> pfr
    '''
    # mix streams if there is more than one
    self, F = mixer(self, F)
    # convert single values in vector
    if np.size(kinetic) == 1: kinetic *= np.ones(N)
    if np.size(reactor) == 1: reactor *= np.ones(N)
    if np.size(P) == 1: P *= np.ones(N+1)
    if np.size(V) == 1: V *= np.ones(N+1)
    if np.size(T) == 1: 
        if T == None: T = [(None)]*(N+1) # isothermal reactor
        else: T *= np.ones(N+1) # non-isothermal reactor
    if np.size(L) == 1: L *= np.ones(N)
    if np.size(ND) == 1: ND *= np.ones(N)
    # preallocate variable (N reactors + 1 inlet)
    nsp = len(species)
#    contents = np.zeros((N+1, self.n_species+2))
    X = np.zeros([N+1,nsp])
    Y = np.zeros([N+1,nsp]) 
    t_chem = np.zeros([N+1,nsp])
    tau_ = np.zeros([N+1,nsp])
    tau = np.zeros(N+1)
    Da = np.zeros([N+1,nsp])
    # check up
    for i in range(N):
        if kinetic[i] == 0 and reactor[i] == 1:
            print('N='+str(i+1)+', You cannot have an equilibrium PFR')
            break
#        if reactor[i] == 0:
#            if L[i] != 0:
#                print 'N='+str(i+1)+', Change CSTR L to zero'
#                L[i] = 0
#            if ND[i] != 0:
#                print 'N='+str(i+1)+', Change CSTR ND to zero'
#                ND[i] = 0
    # create upstream and downstream reservoirs
    reservoir = ct.Reservoir(self)
#    contents[0,:] = get_contents(self)
    T[0] = self.T
    X[0,:] = self[species].X
    Y[0,:] = self[species].Y
#    downstream = ct.Reservoir(self)
#    # create mass flow controller from upstream reservoir to reactor
#    m = ct.MassFlowController(upstream, r, mdot=F)
#    # create pressure controller from reactor to downstream reservoir
#    v = ct.PressureController(r, downstream, master=m, K=1e-5)
    # start to measure simulation time spent
    t0_sim = clock()
    # start simulation
    for i in range(N):
        if i > 0:
            self.TDY = reservoir.thermo.TDY
        print('Inlet reactor #'+str(i+1))
#        print(self.report())
        if kinetic[i] == 0:
            self = equilibrium_reactor(self, F, P[i+1], T[i+1])
        else:
            if reactor[i] == 0: # cstr
                self, t, zero = cstr(self, F, P[i+1], V[i], T[i+1])
            else: # pfr
                self, l, content = pfr(self, F, P[i+1], V[i], T[i+1], 
                                       L[i], ND[i])
        print('Outlet reactor #'+str(i+1))
#        print(self.report())
        # save values
#        contents[i+1,:] = get_contents(self)
        T[i+1] = self.T
        X[i+1,:] = self[species].X
        Y[i+1,:] = self[species].Y
        t_chem[i+1,:] = species_time(self, species, is_abs='no')
        tau[i+1] = residence_time(self, F, V[i])
        tau_[i+1,:] = residence_time(self, F, V[i], species)
        Da[i+1,:] = damkohler(tau_[i+1,:], t_chem[i+1,:])
        # save state to reservoir
        reservoir.thermo.TDY = self.TDY
        # report
        if disp == True:
            print('N=%3i,  P=%4.2f atm,  V=%5.3f m^3,  T=%4.2f C' %(i+1, 
                                   self.P/ct.one_atm, V[i], self.T-273.15))
    # stop to measure simulation time spent
    return {'state': self, 'F': F, 'P': P, 'V': V, 'T': T, 'X': X, 'Y': Y,
            'time spent': clock() - t0_sim, 'residence time': tau, 
            'chemical time': t_chem, 'damkohler': Da}
    
def fw_of_gasifying_agent(phase_of_fuel, phase_of_gasifying_agent=air, 
                          fw_of_fuel=None, fw_of_all=1.0, 
                          equivalence_ratio=1.0):
    if fw_of_fuel != None:
        fm_of_fuel = fw_of_fuel/phase_of_fuel.mean_molecular_weight # kmol/s
        # mass flowrate of gasifying agent
        fm_of_gasifying_agent = equivalence_ratio * stoich_oxygen(
        phase_of_fuel) * fm_of_fuel * (1 - phase_of_fuel['H2O'].X) \
        /phase_of_gasifying_agent['O2'].X # kmol/s
        return fm_of_gasifying_agent \
                * phase_of_gasifying_agent.mean_molecular_weight
    else:
        A = np.array([[1.0, 1.0],
                      [equivalence_ratio * stoich_oxygen(phase_of_fuel) \
                      * (1 - phase_of_fuel['H2O'].X) \
                      / phase_of_gasifying_agent['O2'].X \
                      / phase_of_fuel.mean_molecular_weight,
                      - 1.0 / phase_of_gasifying_agent.mean_molecular_weight]])
        b = np.array([fw_of_all, 0])
        x = np.linalg.solve(A, b)        
        return x[1]

def flow(self, t, V):
    return self.density*V/t

def species_time(self, species, is_abs='no'):
    time = np.zeros(len(species))
#    print self[species].net_production_rates
    for i in range(len(species)):
        # chemical time is zero if rate is zero too
        if self[i].net_production_rates != 0:
            time[i] = self[i].concentrations/self[i].net_production_rates
    if is_abs == 'yes':
        time = np.abs(time)
    return time
                
def residence_time(self, F, V, species=None):
#    self, F = mixer(self, F)
    if species==None:
        return self.density_mass*V/F
    else:
        return self[species].concentrations*self[species].molecular_weights \
                *V/(self[species].Y*F)

def damkohler(residence_time, chemical_time):
    nsp = np.size(chemical_time)    
    ratio = np.zeros(nsp)
    for i in range(nsp):
        # damkohler is zero if chemical time is zero too
        if chemical_time[i] != 0:
            ratio[i] = residence_time[i]/chemical_time[i]
    return ratio

def damkohler_number(content, species, rates, flow, volume):
    return np.min(
            np.abs(
                (np.sum(content[species].concentrations)/rates) /\
                (np.sum(content[species].density_mass)*volume /\
                (np.sum(content[species].Y)*flow))
            )
        )

def chemical_time(content, species, rates):
#    return np.max(np.abs(np.sum(content[species].concentrations)/rates))
    return np.sum(content[species].concentrations)/\
            np.min(np.abs(rates[np.nonzero(rates)]))
    
def equivalence_ratio(fw_of_fuel, phase_of_fuel, fw_of_gasifying_agent,
                      phase_of_gasifying_agent=air):
    fm_of_fuel = fw_of_fuel/phase_of_fuel.mean_molecular_weight # kmol/s
    # mass flowrate of gasifying agent
    fm_of_gasifying_agent = \
    fw_of_gasifying_agent/phase_of_gasifying_agent.mean_molecular_weight
    return fm_of_gasifying_agent * phase_of_gasifying_agent['O2'].X/(
    stoich_oxygen(phase_of_fuel) * fm_of_fuel * \
    (1 - phase_of_fuel['H2O'].X))

def fw_of_air(fw_of_fuel, phase_of_fuel, phase_of_gasifying_agent=air, ER=1.0):
    fm_of_fuel = fw_of_fuel/phase_of_fuel.mean_molecular_weight # kmol/s
    # mole amount of gasifying agent
    fm_of_air = ER * stoich_oxygen(phase_of_fuel) * fm_of_fuel * (
                1 - phase_of_fuel['H2O'].X) / phase_of_gasifying_agent['O2'].X
    # mass amount of gasifying agent
    return fm_of_air * phase_of_gasifying_agent.mean_molecular_weight
    
def inverse_equivalence_ratio(fw_of_fuel, phase_of_fuel, fw_of_gasifying_agent,
                              phase_of_gasifying_agent=air):
    return 1.0/equivalence_ratio(fw_of_fuel, phase_of_fuel, 
                                 fw_of_gasifying_agent, 
                                 phase_of_gasifying_agent)

def stoich_oxygen(self):
    nsp = self.n_species
    # initiate variable
    mol = np.zeros(4, dtype=float)
    # compute number of C,H,O moles of fuel species
    for i in range(nsp):
        if self.X[i] != 0:
            if i != self.species_index('H2O'):
                if i != self.species_index('CO2'):
                    mol[0] += self.X[i] \
                    * self.n_atoms(self.species_name(i), 'C')
                    mol[1] += self.X[i] \
                    * self.n_atoms(self.species_name(i), 'H')
                    mol[2] += self.X[i] \
                    * self.n_atoms(self.species_name(i), 'O')
                    mol[3] += self.X[i]
    # normalise for 1 mol of fuel
    mol = mol/mol[3]
    # stoichiometric moles of oxygen per mol of fuel
    ## CxHyOz + (x+0.25y-0.5z)O2 --> xCO2 + (2x+0.5y-0.5z)H2O 
    return mol[0] + 0.25*mol[1] - 0.5*mol[2]
