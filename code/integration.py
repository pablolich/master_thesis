#!/usr/bin/env python3

__appname__ = '[analysis.py]'
__author__ = 'Pablo Lechon (pl1619@ic.ac.uk)'
__version__ = '0.0.1'


## IMPORTS ##

import sys
from implementation import *
from scipy.special import comb
from matplotlib import gridspec


## CONSTANTS ##

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; DeltaGATP = 75e3 # J/mol
global T; T = 298 # K 
global E = 3e6 # J (Total energy contained in metabolites)


## FUNCTIONS ##

def main(argv):
    '''Integrate model for a set of arbitrary parameters'''

    #############################
    #1. Environmental conditions#
    #############################

    #Number of metabolites
    m = 5
    #Number of strains
    s = 3
    #Pack total energy into m intervals (chemical potentials of metabolites)
    pack = np.random.uniform(E, size = m)
    #Sort decreasingly 
    mu = np.sort(pack)[::-1]
    #Fix first and last chemical potentials
    mu[0] = E 
    mu[-1] = 0
    #Set growth and maintenance parameters
    g = np.ones(s) 
    maintenance = 0.2*np.ones(s)
    #Set supply and dilution rate of metabolites
    kappa = 8*np.ones(m) 
    gamma = 1*np.ones(m)


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
    #Initialize reaction rates matrix
    reaction_rates = np.zeros(shape = (s,m,m))
    #Prealocate the reaction tensor
    reactions = np.zeros(shape = (s,m,m), dtype = 'int')
    #Initialize reaction network container for one strain
    reac_network = np.array([[], []], dtype = 'int32')
    #Prealocate container of all networks 
    tot_reac_network = s*[reac_network]
    #Prealocate number of reactions possesed by each strain 
    n_reac_s = np.zeros(s, dtype = 'int')
    #Total number of posible reactions (only forward net reactions are taking 
    #place)
    N_reac = comb(m, 2, exact = True)
    #Generate reaction network fro each strain
    for i in range(s):
        #Set number of reactions in the network
        num_reac = np.random.randint(1, N_reac+1)
        #Get reaction network
        reac_network = network(m, num_reac)
        #Store reaction network information for later use
        tot_reac_network[i] = reac_network
        reactions[i][reac_network] = 1
        non_zero = reactions[i][reactions[i].nonzero()]
        n_reac_s[i] = len(non_zero)

    #Initial conditions 
    N0 = s*[2]
    C0 = m*[0] 
    z0 = N0 + C0
    #Time span
    tspan = tuple([1, 500])
    #Integration
    sol = solve_ivp(lambda t,z: model(t,z,s,m,kappa,gamma,tot_reac_network, mu,
                    Eta, q_max, ks, kr, g, maintenance), tspan, z0, 
                    method = 'BDF')
    t = sol.t
    z = sol.y
    #Store solutions 
    N = z[0:s]
    C = z[s:s+m]

    #Save  reactions, N, C for later analysis
    _type = ['classic', 'new']
    np.save('../data/reactions_' + _type[0] + '.npy', reactions)
    np.save('../data/population_' + _type[0] + '.npy', N)
    np.save('../data/concentration_' + _type[0] + '.npy', C)

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

