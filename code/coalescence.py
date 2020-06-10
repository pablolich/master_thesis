#!/usr/bin/env python3

__appname__ = '[coalescence.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
from functions import *
from integration import uniform_pack
from scipy.special import comb
from progressbar import ProgressBar
import pandas as pd

## CONSTANTS ##

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; DeltaGATP = 75e3 # J/mol
global T; T = 298 # K 
global E; E = 3e6 # J (Total energy contained in metabolites)

## FUNCTIONS ##

def names(s, m):
    '''Create a vector of column names for results dataframe'''
    #Strains
    s_seq = [str(i) for i in range(s)]
    names_s = np.repeat('s', s)
    label_s = np.char.add(names_s, s_seq)
    m_seq = [str(i) for i in range(m)]
    names_m = np.repeat('m', m)
    label_m = np.char.add(names_m, m_seq)
    return np.concatenate([label_s, label_m])
    

def main(argv):
    '''Main function'''
    
    #Set fix parameters#
    ####################

    #Number of metabolites
    m = 10
    #Number of strains
    s = 50
    #Generate random uniformly distributed chemical potentials
    mu = uniform_pack(E, m)
    #Set growth proportionality constant
    g = np.ones(s) 
    #Set supply and dilution rate of metabolites
    kappa = 2*np.ones(m) 
    gamma = 1*np.ones(m)
    #Average cost per reaction
    xi0 = 0.005
    #Choose Eta matrix for each strain
    Eta = 0.5*np.ones(shape = (s,m,m))
    diag_ind = np.diag_indices(m)
    for i in range(s):
        Eta[i][diag_ind] = 0
    #Choose q_max, ks and krn
    q_max = [np.ones(1) for i in range(s)]
    ks = [0.1*q_max[i] for i in range(s)]
    kr = [10*q_max[i] for i in range(s)]
    #Total number of posible reactions (only forward net reactions are taking 
    #place)
    N_reac = comb(m, 2, exact = True)
    #Perform n_simul simulations
    n_simul = 1000
    #Generate reaction network for each strain
    tot = s*n_simul
    n_reac_s = np.zeros(tot, dtype = 'int')
    maintenance = np.zeros(tot)
    #Prealocate container of all networks 
    tot_reac_network = tot*[0]
    pbar = ProgressBar()
    print('Creating reaction networks')
    for i in pbar(range(tot)):
        #Randomly set number of reactions in the network
        num_reac = np.random.randint(1, N_reac+1)
        #Get reaction network
        reac_network = network(m, num_reac)
        #Store reaction network 
        tot_reac_network[i] = reac_network
        #Store actual number of reactions
        n_reac_s[i] = len(tot_reac_network[i][0])
        #Set cost based on the number of reactions
        maintenance[i] = xi0*n_reac_s[i]

    #Integrate model#
    #################

    #Time span
    tspan = tuple([1, 1e5])
    #Initial conditions 
    N0 = s*[2]
    C0 = m*[2] 
    z0 = N0 + C0


    #Prealocate storing objects
    empty_z = [[] for i in range(m+s)]
    zn = np.array(empty_z)
    zn_stable = np.zeros(shape = (n_simul, m+s))
    tn = np.array([])
    n_label = np.array([])

    pbar = ProgressBar()
    for  n in pbar(range(n_simul)):

        #Network
        network_n = tot_reac_network[s*n:s*(n+1)]
        #Maintenance vector
        maintenance_n = maintenance[s*n:s*(n+1)]
        #Integrate the whole community
        sol = solve_ivp(lambda t,z: model(t,z,s,m,kappa,gamma,network_n,mu,Eta,
                        q_max,ks,kr,g,maintenance_n), tspan, z0, 
                        method = 'BDF', atol = 1e-3)
        #Extract solution and time evaluations
        z = sol.y
        t = sol.t
        #Store
        zn = np.concatenate([zn, z], axis = 1)
        zn_stable[n,:] = zn[:,-1]
        tn = np.concatenate([tn, t])
        n_label = np.concatenate([n_label, n*np.ones(len(t), dtype = int)])

    #Create dataframe of networks
    df_network = pd.DataFrame(tot_reac_network, 
                              columns = ['substrate', 'product'])
    df_network.insert(0, 'strain', np.tile(np.arange(s), n_simul))
    df_network.insert(0, 'n_simulation', np.repeat(np.arange(n_simul), s))
    #Create dataframe of time series solutions
    cols = names(s, m)
    df_sol = pd.DataFrame(np.transpose(zn), columns = cols)
    df_sol.insert(0, 't', tn)
    df_sol.insert(0, 'n_simulation', n_label)
    #Create dataframe of stable community compositions
    df_stable = pd.DataFrame(zn_stable, columns = cols)
    df_stable.insert(0, 'n_simulation', np.arange(n_simul))
    #Save
    df_sol.to_csv('../data/coal_time_series.csv')
    df_network.to_csv('../data/coal_network.csv')
    df_stable.to_csv('../data/coal_composition.csv')

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

