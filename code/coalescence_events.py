#!/usr/bin/env python3

__appname__ = '[coalescence_event.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp, odeint
from coalescence import uniform_pack
from functions import model
import matplotlib.pylab as plt

## CONSTANTS ##

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; DeltaGATP = 75e3 # J/mol
global T; T = 298 # K 
global E; E = 3e6 # J (Total energy contained in metabolites)

## FUNCTIONS ##

def main(argv):
    '''Main function'''

    #Load communities data 
    composition = pd.read_csv('../data/coal_composition.csv', index_col = 0)
    networks = pd.read_csv('../data/coal_network.csv', index_col = 0)
    time_series = pd.read_csv('../data/coal_time_series.csv')

    #Transform data to long_format
    id_vars = ['n_simulation']
    comp_long = pd.melt(composition, id_vars, 
                        var_name = 'str_met',
                        value_name = 'pop_conc') 
    #Get simulation number for which integration diverges
    ind = comp_long.index[comp_long.pop_conc>1e4].tolist()
    bad_simulations = pd.unique(comp_long.n_simulation.iloc[ind])
    #Get Trues where n_simulation is one of the bad_simulations
    trues = comp_long.n_simulation.isin(bad_simulations)
    #Get index of these occurences
    ind = np.where(trues == True)[0]
    #Get rid of divergent communities
    comp_long.drop(comp_long.index[list(ind)], inplace = True)
    #How many communities do we have?
    num_com = len(np.unique(comp_long.n_simulation))
    #Perform the coalescence of two communities
    #Pick communities 117 and 139
    comp_117 = comp_long[comp_long.n_simulation == 117].reset_index(drop = True)
    comp_139 = comp_long[comp_long.n_simulation == 139].reset_index(drop = True)
    #Get networks of strains in that community
    networks_117 = networks[networks.n_simulation == 117].reset_index(drop=True)
    networks_117 = networks_117.drop(['strain'], axis = 1)
    networks_139 = networks[networks.n_simulation == 139].reset_index(drop=True)
    networks_139 = networks_139.drop(['strain'], axis = 1)
    #Drop metabolites rows
    str_117 = comp_117[comp_117.str_met.str.startswith('s')]
    str_139 = comp_139[comp_139.str_met.str.startswith('s')]
    #Create dataframe of community with surviving strains and networks
    community_117 = pd.merge(str_117, networks_117, on = 'n_simulation', 
                             left_index = True, right_index = True)
    community_139 = pd.merge(str_139, networks_139, on = 'n_simulation', 
                             left_index = True, right_index = True)
    #Get rid of extinct species
    stable_com_117 = community_117[community_117.pop_conc>1]
    stable_com_139 = community_139[community_139.pop_conc>1]
    #Get number of strains in each community
    num_117 = len(stable_com_117)
    num_139 = len(stable_com_139)
    #Put together both communities
    mixed = pd.concat([stable_com_117, stable_com_139]).reset_index(drop = True)
    #Transform from string to numeric vector all elements in substrate and 
    #product colums
    mixed['substrate'] = mixed['substrate'].apply(lambda x: 
                               np.fromstring(x[1:-1], sep=' ', dtype = int))
    mixed['product'] = mixed['product'].apply(lambda x: 
                             np.fromstring(x[1:-1], sep=' ', dtype = int))

    #Get parameters for integration
    #Obtain number of metabolites and strains
    m = len(comp_117) - len(str_117)
    s = len(str_117)
    #Obtain reaction networks as a list of tuples
    nets = [tuple([mixed.loc[i, ['substrate']][0],
                   mixed.loc[i, ['product']][0]]) for i in range(len(mixed))]
    #Choose Eta matrix for each strain
    Eta = 0.5*np.ones(shape = (s,m,m))
    diag_ind = np.diag_indices(m)
    for i in range(s):
        Eta[i][diag_ind] = 0
    q_m = [np.ones(1) for i in range(s)]
    ks = [0.1*q_m[i] for i in range(s)]
    kr = [10*q_m[i] for i in range(s)]

    #Initial conditions                
    z0 = list(mixed.pop_conc) + m*[2]

    #Timespan
    tspan = tuple([1, 1e5])

    #Integrate model
    sol = solve_ivp(lambda t,z: model(t, z, s = len(mixed), m = m, 
                                      kappa = 2*np.ones(m), 
                                      gamma = 1*np.ones(m), 
                                      networks = nets, 
                                      mu = uniform_pack(E, m), 
                                      Eta = Eta, 
                                      q_max = q_m, 
                                      ks = ks, 
                                      kr = kr,
                                      g = np.ones(s),
                                      maint = np.array(mixed['maintenance'])),
                    tspan, z0, method = 'BDF', atol = 1e-3)

    t = sol.t
    z = sol.y
    N = z[0:s]
    C = z[s:s+m]
    
    colors = num_117*['red'] + num_139*['blue']
    for i in range(m):
        plt.plot(t, N[i], label= 'Strain'+str(i), color = colors[i])
        
    plt.xscale('log')
    plt.show()
                                





    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
     
