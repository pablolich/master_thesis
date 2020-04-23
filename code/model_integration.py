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
    m = 100
    #Chemical potentials 
    interval = 3e6 
    mu_rev = np.sort(np.random.uniform(interval, size = m))
    mu = mu_rev[::-1]
    #Fix first and last chemical potentials
    mu[0] = interval 
    mu[-1] = 0
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
        #Select only reactions that go forward
        tri_ind = np.triu_indices(m, 1)
        #Initialite reaction network
        reac_network = np.array([[], []], dtype = 'int32')
        #List of present metabolite (specified as input argument)
        list_m_i = [0] 
        #Initialize reaction sample
        keep_sampling = True
        while keep_sampling: 
            #Sample one reaction
            #Choose a substrate from list of present substrates
            m_i = int(np.random.choice(list_m_i))
            #Choose a product from list of metabolites with a lower chemical
            #potential, so that the reaction takes place.
            m_j = int(np.random.choice(np.arange(m_i+1, m)))
            #Note that doing m_i + 1 avoids having reactions i-->i
            #Create the tuple representing the reaction
            r_sample = np.array([[m_i], [m_j]])
            #Is it an energetically valid reaction?
            #Check if Gibbs energy change is not too big in order to assure 
            #that the end metabolite is reached through multiple steps.
            if mu[r_sample[1]] - mu[r_sample[0]] > -4*DeltaGATP:
                #Add the reaction to reac_network
                reac_network = np.concatenate((reac_network, r_sample),axis = 1)
                #Eliminate repeated reactions
                reac_network = np.unique(reac_network, axis = 1)
                #Add product to list of substrates
                list_m_i.append(m_j)
                if m-1 in reac_network[1]:
                    #When the last metabolite is reached, stop sampling
                    keep_sampling = False
        
        #Transform to indices in reactions matrix
        reac_network = tuple(reac_network)
        #Switch on the matrix elements in reactions contained in reac_network
        reactions[reac_network] = 1
        #import matplotlib.pylab as plt
        #plt.imshow(reactions)
        #plt.show()

        #2b. Calculate rates
        #Get rates of each reaction in the network
        #Substrate and product concentrations from vector C for metabolites
        #taking part in reac_network 
        S = C[reac_network[0], i-1]
        P = C[reac_network[1], i-1]
        #The following calculations are performed for all reactions at once.
        #Calculate reaction quotients
        Q = r_quotient(S, P)
        #Calculate Gibbs free energy changes
        mu_S = mu[reac_network[0]]
        mu_P = mu[reac_network[1]]
        #Gibs energy change is the difference between chemical potentials 
        #(stechiometric coefficients are all 1)
        DeltaG = mu_P - mu_S
        #Calculate equilibrium constant
        #Get etas for reaction network
        eta = Eta[reac_network]
        Keq = K_equilibrium(DeltaG, eta, DeltaGATP, R, T)
        #Calculate thetas
        theta = Theta(Q, Keq)
        #Calculate rates
        q_max = np.ones(len(reac_network[0]))
        ks = 0.1*q_max
        kr = 10*q_max
        #Include reaction rates in reaction network matrix
        q_reac = rate(q_max, theta, ks, kr, S)
        reactions[reac_network] = q_reac
        import ipdb; ipdb.set_trace(context = 20)


    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

