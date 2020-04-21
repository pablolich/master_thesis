#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
from model_functions import *
from scipy.special import comb

## CONSTANTS ##

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; DeltaGATP = 75e3 # J/mol
global T; T = 298 # K 

## FUNCTIONS ##

def main(argv):
    '''Main function'''

    #Number of timepoints
    n = 500
    #Timepoints
    t = np.linspace(1, 10, n)
    #Number of metabolites
    m = 4
    #Chemical potentials the interval is 10 because I like it
    interval = 10
    mu = np.sort(np.random.uniform(interval, size = m))
    mu[-1] = interval
    #Number of posible reactions
    n_reac = comb(m, 2, exact = True)
    #Number of strains
    s = 1
    #Initial conditions
    Cinit = [0]*m
    Cinit[0] = 1
    Ninit = [1]
    #Prealocate matrices for solutions
    C = np.zeros(shape = (m, len(t)))
    N = np.zeros(shape = (s, len(t)))
    #Add initial conditions to the first element of the solution
    C[:,0] = Cinit
    N[:,0] = Ninit
    
    
    ################
    #SET PARAMETERS#
    ################

    #1. eta is chosen randomly for each reaction
    ############################################
    #For now we set all of them to 0.5
    Eta = np.ones(shape = (m,m))
    #Get indices of upper triangular elements with an offset of 1 to exclude
    #diagonal elements.
    tri_ind = np.triu_indices(m, 1)
    #Randomly (not for now) choose etas
    eta_set = [0.5]*n_reac
    Eta[tri_ind] = eta_set

    #Numerically integrate the model
    for i in range(1, n):
        #Declare timespan for this iteration of numerical integration
        tspan = [t[i-1], t[i]]

        #2. Rates are calculated according to enzyme kinetics
        #####################################################
        #Note that rates change after every iteration of the integration
        #because the concentrations change.
        #2a. Set reaction network
        reactions = np.zeros(shape = (m,m))
        tri_ind = np.triu_indices(m, 1)
        reaction_set = np.random.randint(2, size = n_reac)
        reactions[tri_ind] = reaction_set
        #Avoid that all elements are 0 (ie, the microbe has at least one 
        #reaction)
        if not np.any(reactions):
            #Choose one element of the matrix
            ind = np.random.randint(n_reac)
            element = (tri_ind[0][ind], tri_ind[1][ind]) 
            #Assign the value 1 to that matrix element
            reactions[element] = 1
        #2b. Calculate rates
        #Get non-zero element indices (ie, labels of metabolites present in
        #reaction network)
        reac_network = np.nonzero(reactions)
        #Get rates of each reaction in the network
        #Substrate and product  according to C and reac_network concentrations
        S = C[reac_network[0], i-1]
        P = C[reac_network[1], i-1]
        #The following calculations are performed for all reactions at once.
        #Calculate reaction quotients
        Q = r_quotient(S, P)
        #Calculate Gibbs free energy changes
        mu_S = mu[reac_network[0]]
        mu_P = mu[reac_network[1]]
        DeltaG = mu_P*P - mu_S*S
        #Calculate equilibrium constant
        import ipdb; ipdb.set_trace(context = 20)
        #Get etas for reaction network
        eta = Eta[reac_networks]
        Keq = K_equilibrium(DeltaG, eta, DeltaGATP, R, T)
        #Calculate equilibrium constants
        R = np.matrix([[0, 0.5], [0, 0]])

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

