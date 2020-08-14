#!/usr/bin/env python3

__appname__ = '[fitness_decoupling.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
from functions import *
from scipy.special import comb
from progressbar import ProgressBar
import pandas as pd

## CONSTANTS ##

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; DeltaGATP = 75e3 # J/mol
global T; T = 298 # K 
global E; E = 3e6 # J (Total energy contained in metabolites)
global seed; seed = 69

## FUNCTIONS ##

def main(argv):
    '''Main function'''

    #Set fix parameters#
    ####################

    #Number of metabolites
    m = 15
    #Number of strains
    s = 10 
    #Generate random uniformly distributed chemical potentials
    np.random.seed(seed)
    #mu = random_pack(E, m)
    mu = decreasing_pack(E, m)
    #mu = uniform_pack(E, m)
    #Set growth proportionality constant
    g = np.ones(s) 
    #Set supply and dilution rate of metabolites
    kappa = 2*np.ones(m) 
    gamma = 1*np.ones(m)
    #Average cost per reaction
    xi0 = 0.005
    epsilon = 1e-1
    #Choose Eta matrix for each strain
    Eta = np.zeros(shape = (s,m,m))
    triu_ind = np.triu_indices(m,1)
    tril_ind = np.tril_indices(m,-1)
    eta_values_u = (triu_ind[1]-triu_ind[0])/(m-1)
    eta_values_l = (tril_ind[0]-tril_ind[1])/(m-1)
    for i in range(s):
        Eta[i][triu_ind] = eta_values_u
        Eta[i][tril_ind] = eta_values_l
    #Choose q_max, ks and krn
    q_max = [np.ones(1) for i in range(s)]
    ks = [0.1*q_max[i] for i in range(s)]
    kr = [10*q_max[i] for i in range(s)]
    #Total number of posible reactions (only forward net reactions are taking 
    #place)
    N_reac = comb(m, 2, exact = True)
    #Perform n_simul simulations
    n_simul = 50
    #Generate reaction network for each strain
    tot = s*n_simul
    n_reac_s = np.zeros(tot, dtype = 'int')
    maintenance = np.zeros(tot)
    #Initialize redundance
    redundance = np.zeros(tot)
    #Prealocate vector for storing initial resource surpluses
    surplus = np.zeros(tot)

    #Prealocate container of all networks 
    tot_reac_network = tot*[0]
    pbar = ProgressBar()
    print('Creating reaction networks')
    for i in pbar(range(tot)):
        #Randomly set number of reactions in the network
        num_reac = 0
        np.random.seed(i+seed)
        #while num_reac == 0: 
        #    num_reac = round(np.random.beta(1,5)*N_reac)
        num_reac = np.random.randint(1, m + 1)#N_reac+1)
        #Get reaction network
        reac_network = network(m, num_reac)
        #Store reaction network 
        tot_reac_network[i] = reac_network
        #Store actual number of reactions
        n_reac_s[i] = len(tot_reac_network[i][0])
        #Set cost based on the number of reactions
        rand_cost = np.random.normal(0,1)
        while rand_cost < -1/epsilon:
            rand_cost = np.random.normal(0,1)
        #Free energy gap of reactions
        mult = sum(tot_reac_network[i][1]-tot_reac_network[i][0])
        maintenance[i] = xi0*mult#*(1+epsilon*rand_cost)
        redundance[i] = 2*len(np.unique(tot_reac_network[i][0])) - \
                len(tot_reac_network[i][0])

    ####################
    #3. Integrate model#
    ####################

    #Time span
    tspan = tuple([1, 1e5])
    #Initial conditions 
    N0 = s*[1]
    C0 = m*[2] 
    z0 = N0 + C0


    #Prealocate storing objects
    empty_z = [[] for i in range(m+s)]
    zn = np.array(empty_z)
    zn_stable = np.zeros(shape = (n_simul, m+s))
    tn = np.array([])
    n_label = np.array([])
    #Initialize vector of individual performances (carrying capacities)
    K = np.zeros(shape = (n_simul, s))
    pbar = ProgressBar()
    for  n in pbar(range(n_simul)):

        #Network
        network_n = tot_reac_network[s*n:s*(n+1)]
        #Maintenance vector
        maintenance_n = maintenance[s*n:s*(n+1)]
        #Get s matrices of zeros: reaction presence/absence for each strain
        rates = np.zeros(shape = (s,m,m))
        #Calculate initial surplus for each strain given environment C0
        for j in range(s):
            network_t = network_n[j];
            eta = Eta[j]
            q_m = q_max[j]
            k_s = ks[j]
            k_r = kr[j]
            rates[j] = rate_matrix(network_t, eta, q_m, k_s, k_r, np.array(C0),
                                   mu, m)
            surplus[s*n + j] = jgrow(rates[j], eta) - maintenance_n[j]
        #Get individual performance of each species
        #Integrate each species separately
        for i in range(s):
            #Eliminate all the reaction networks except the ith one. 
            network_ni = [(network_n[i])]
            Etai = Eta[i].reshape(1,m,m)
            q_m = [np.array(q_max[i])]
            kri = [np.array(kr[i])]
            ksi = [np.array(ks[i])]
            gi = np.array([g[i]])
            mi = np.array([maintenance_n[i]])
            N0i = 1*[2]
            z0i = N0i + C0
            sol = solve_ivp(lambda t,z: model(t,z,1,m,kappa,gamma,network_ni,
                            mu,Etai,q_m,ksi,kri,gi,mi), tspan, z0i, 
                            method = 'BDF', atol = 1e-3)
            z = sol.y
            K[n,i] = z[0][-1]

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
    #Create dataframe with k, and my proxy of fitness
    df = pd.DataFrame({'K':K.flatten(), 'fitness':surplus+redundance})
    df.to_csv('../data/abundance_fitness.csv')

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

