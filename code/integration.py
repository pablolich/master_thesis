#!/usr/bin/env python3

__appname__ = '[integration.py]'
__author__ = 'Pablo Lechon (pl1619@ic.ac.uk)'
__version__ = '0.0.1'


## IMPORTS ##

import sys
from functions import *
from scipy.special import comb
from matplotlib import gridspec
import pandas as pd
from progressbar import ProgressBar

## CONSTANTS ##

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; DeltaGATP = 75e3 # J/mol
global T; T = 298 # K 
global E; E = 3e6 # J (Total energy contained in metabolites)


## FUNCTIONS ##

def labels(i, j):
    '''Create labels for reaction network'''
    ind = np.transpose(np.indices([i,j]).reshape(2, i*j))
    #label = [str(ind[i][0])+str(ind[i][1]) for i in range(i*j)]
    label = pd.DataFrame()
    label['simulation'] =  ind[:,0]
    label['strain'] = ind[:,1]
    return label

def random_pack(E, m):
    '''Generate random uniformly distributed chemical potentials'''
    #Pack total energy into m intervals (chemical potentials of metabolites)
    pack = np.random.uniform(E, size = m)
    #Sort decreasingly 
    mu = np.sort(pack)[::-1]
    #Fix first and last chemical potentials
    mu[0] = E 
    mu[-1] = 0
    return mu

def uniform_pack(E, m):
    '''Generate equidistant chemical potentials'''
    return np.linspace(E,0,m)

def main(argv):
    '''Integrate model for a set of arbitrary parameters'''

    #############################
    #1. Environmental conditions#
    #############################

    #Number of metabolites
    m = 10
    #Number of strains
    s = 200
    #Generate random uniformly distributed chemical potentials
    mu = uniform_pack(E, m)
    #Set growth and maintenance parameters
    g = np.ones(s) 
    #Set supply and dilution rate of metabolites
    #Diversification
    kappa = np.concatenate([np.array([2]), np.zeros(m-1)])
    #Provide all resources
    #kappa = 2*np.ones(m) 
    gamma = 1*np.ones(m)
    #Average cost per reaction
    xi0 = 0.005

    ###########################################
    #2. Biochemical parameters for each strain#
    ###########################################

    #Choose Eta matrix for each strain
    Eta = 0.5*np.ones(shape = (s,m,m))
    diag_ind = np.diag_indices(m)
    for i in range(s):
        Eta[i][diag_ind] = 0
    #Choose q_max, ks and krn
    q_max = [np.ones(1) for i in range(s)]
    ks = [0.1*q_max[i] for i in range(s)]
    kr = [10*q_max[i] for i in range(s)]
    #Prealocate number of reactions possesed by each strain 
    #Total number of posible reactions (only forward net reactions are taking 
    #place)
    N_reac = comb(m, 2, exact = True)
    #Perform n_simul simulations
    n_simul = 2
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
        #Number of reactions is fixed
        #num_reac = 45
        #Cost is constant
        #maintenance[i] = xi0*N_reac
        #Get reaction network
        reac_network = network(m, num_reac, s = 0.1)
        #Store reaction network 
        tot_reac_network[i] = reac_network
        #Store actual number of reactions
        n_reac_s[i] = len(tot_reac_network[i][0])
        #Set cost based on the number of reactions
        maintenance[i] = xi0*n_reac_s[i]

    n = 0

    ####################
    #3. Integrate model#
    ####################

    #Time span
    tspan = tuple([1, 1e5])
    #Initial conditions 
    N0 = s*[2]
    C0 = m*[0] 
    z0 = N0 + C0


    #Prealocate storing objects
    empty_z = [[] for i in range(m+s)]
    zn = np.array(empty_z)
    zn_stable = np.zeros(shape = (n_simul, m+s))
    tn = np.array([])
    n_label = np.array([])
    #K = np.zeros(shape = (n_simul, s))
    while n < n_simul:

        #Network
        network_n = tot_reac_network[s*n:s*(n+1)]
        #Maintenance vector
        maintenance_n = maintenance[s*n:s*(n+1)]
        ##Get individual performance of each species
        ##Initialize vector of individual performances
        ##Integrate each species separately
        #pbar = ProgressBar()
        #print('Calculating individual performance of each strain')
        #for i in pbar(range(s)):
        #    #Eliminate all the reaction networks except the ith one. 
        #    tspani = tuple([1, 100])
        #    network_ni = [(network_n[i])]
        #    Etai = Eta[i].reshape(1,m,m)
        #    q_m = [np.array(q_max[i])]
        #    kri = [np.array(kr[i])]
        #    ksi = [np.array(ks[i])]
        #    gi = np.array([g[i]])
        #    mi = np.array([maintenance_n[i]])
        #    N0i = 1*[2]
        #    z0i = N0i + C0
        #    sol = solve_ivp(lambda t,z: model(t,z,1,m,kappa,gamma,network_ni,
        #                    mu,Etai,q_m,ksi,kri,gi,mi), tspani, z0i, 
        #                    method = 'BDF', atol = 1e-3)
        #    z = sol.y
        #    K[n,i] = z[0][-1]

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

        print(n)
        n += 1

    #Create dataframe of networks
    df_network = pd.DataFrame(tot_reac_network, 
                              columns = ['substrate', 'product'])
    label = labels(n_simul, s)
    df_network = pd.concat([label, df_network], axis = 1)
    #Create dataframe of time series solutions
    df_sol = pd.DataFrame(np.transpose(zn))
    df_sol.insert(0, 't', tn)
    df_sol.index = n_label
    #Create dataframe of stable community compositions
    df_stable = pd.DataFrame(zn_stable)
    #Create a dataframe of individual performances
    #df_performance = pd.DataFrame(K)
    #Save
    df_sol.to_csv('../data/time_series.csv')
    df_network.to_csv('../data/network.csv')
    df_stable.to_csv('../data/composition.csv')
    df_performance.to_csv('../data/performance.csv')

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

