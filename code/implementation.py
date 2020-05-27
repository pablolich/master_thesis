#!/usr/bin/env python3

__appname__ = '[implementation]'
__author__ = 'Pablo Lechon (pl1619@ic.ac.uk)'
__version__ = '0.0.1'

## CONSTANTS ##

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; DeltaGATP = 75e3 # J/mol
global T; T = 298 # K 

## IMPORTS ##

import numpy as np
import matplotlib.pylab as plt
from functions import *
import scipy.stats as stats

## FUNCTIONS ##

def kth_diag_indices(a, k):
    '''Returns indices of the kth diagonal of matrix a'''
    rows, cols = np.diag_indices_from(a)
    return rows[:-k], cols[k:]

def choose_reaction(m):
    '''
    Chooses one reaction with a probability p that decreases the further 
    away we are from the main diagonal
    '''
    #Probabilities follow a truncated normal distribution N(1, sqrt(m))
    sigma = np.sqrt((m))
    k = round(abs(np.random.normal(1, sigma)))
    #Avoid k being in the main diagonal
    while (k < 1):
        k = round(abs(np.random.normal(1, sigma)))
    #Get indices of that diagonal
    ind = kth_diag_indices(np.ones(shape = (m,m)), k)
    #Choose uniformly at random one of the elements in this diagonal
    try: 
        el = np.random.randint(m-k)
        #Substrate
        sub = ind[0][el]
        #Product
        prod = ind[1][el]
    except: 
        sub = 0
        prod = m-1
    #Create reaction
    reaction = np.array([[sub], [prod]])

    return reaction


def generate_network(s, m, nATP, mu, num_reac):

    '''
    Create reaction network of feasible reactions (those which associated
    Gibbs energy change is lower than nATP*D(G_ATP)) based on the
    number of metabolites,  microbial strains and chemical potentials.

    Parameters:
        
        s (int): number of strains
        m (int): number of metabolites
        nATP (float): moles of ATP generated per reaction
        mu (array): vector of chemical potentials for each metabolite

    Returns:

        tot_reac_network (list): reaction network for each strain
        reactions (array): s matrices of 0 and 1, where 1 represents that
        reaction i-->j is present in strain s.

    '''

    #Initialize output
    reac_network = np.array([[], []], dtype = 'int')
    #list of present substrates (at the start only substrate 0 is present, so 
    #the network will always start at 0)
    list_m_i = [0] 
    #Initialize reaction counter
    num = 0
    while num < num_reac: 
        #Choose a substrate from list of present substrates (at the start
        #only substrate 0 is present, so the network will always start at
        #0)
        m_i = int(np.random.choice(list_m_i))
        #Avoid choosing the last metabolite
        while m_i+1 == m:
            m_i = int(np.random.choice(list_m_i))
        #Create list of possible products of the reaction compatible where
        #Delta G < 0
        products = np.arange(m_i+1, m, dtype = int)
        #Note that we start from m_i to satisfy the latter condition. 
        #doing m_i + 1 avoids having reactions i-->i
        #Choose a product from list of metabolites
        m_j = int(np.random.choice(products))
        #Create the tuple representing the sampled reaction
        r_sample = choose_reaction(m)
        #r_sample = np.array([[m_i], [m_j]])
        #Is it an energetically valid reaction?
        #Check if Gibbs energy change is not too big in order to assure 
        #that the end metabolite is reached through multiple steps.
        if mu[r_sample[1]] - mu[r_sample[0]] > -nATP*DeltaGATP:
            #Add the reaction to reac_network
            reac_network = np.concatenate((reac_network, r_sample),axis = 1)
            #Eliminate repeated reactions
            reac_network = np.unique(reac_network, axis = 1)
            #Add product to list of substrates
            list_m_i.append(m_j)
            #When the last metabolite is reached, stop sampling. The 
            #reaction network has reached the end product (that with the 
            #lowest chemical potential)
        #Increase reaction counter
        num += 1
    reac_network= tuple(reac_network)
    
    return reac_network


def model_integration(t, s, m, tot_reac_network, mu, Eta, q_max, ks, kr, 
                      reaction_rates, g, N, C, maintenance, kappa, gamma):

    '''
    Integrate jacob's model.

    Parameters:
        t (1xn array): time vector
        s (int): number of strains
        m (int): number of metabolites
        tot_reac_network (s tuples): reaction network of each strain. Each 
                                     tuple stores two arrays for substrates and
                                     products.
        mu (1xm array): chemical potential of each metabolite
        Eta (mxmxs array): energy associated to each of the possible single 
                           reactions for each strain.
        q_max (sx1 array): maximal reaction rate. Each element of the array is
                           an array with the maximal reaction rates for each
                           reaction in the network of strain i.
        ks (sx1 array): disassociation constant. Each element of the array is
                        an array with the disassociation constants for each
                        reaction in the network of strain i.
        kr (sx1 array): reversibility factor. Each element of the array is
                        an array with the reversibility factors for each
                        reaction in the network of strain i.
        reaction_rates (mxmxs array): reaction rates of reactions possesed by 
                                      each strain.
        g (1xs array): proportionality constant relating energy to population 
                       of each strain.
        maintenance (1xs array): maintenance energy requirement of each strain.
        kappa (1xm array): supply rate of each metabolite.
        gamma (1xm array): dilution rate of each metabolite.
        N (sxn array): population vector time series for each strain.
        C (mxn array): metabolite concentration time series vector for each 
                       metabolite.

    '''

    #Get number of timepoints
    n = len(t)
    #Loop through time
    for i in range(1, n):
        #Declare timespan for this iteration of numerical integration
        tspan = [t[i-1], t[i]]

        #1. Calculate  rates and asociated quatnitiesfor time step i
        ############################################################

        #Note that rates change after every iteration of the integration
        #because the concentrations change, that is why they need to be 
        #recalculated after each integration step.
        #Prealocate elements that will be used later on
        G = np.zeros(s)
        M = np.zeros(s)
        flux_in_out = np.zeros([s, m])
        #Loop through all strains
        for j in range(s):
            #Get concentrations from those metabolites taking part in reaction 
            #network 
            S = C[tot_reac_network[j][0], i-1]
            P = C[tot_reac_network[j][1], i-1]
            #Get chemical potentials of those metabolites
            mu_S = mu[tot_reac_network[j][0]]
            mu_P = mu[tot_reac_network[j][1]]

            #Calculate quantities necesary for the integration
            ###################################################################
            #The following calculations are performed for all reactions at once
            ###################################################################
            #Calculate reaction quotients
            Q = r_quotient(S, P)
            #Gibs free energy change is the difference between chemical 
            #potentials (stechiometric coefficients are all 1)
            DeltaG = mu_P - mu_S
            #Get etas for reaction network
            eta = Eta[j][tot_reac_network[j]]
            #Calculate equilibrium constant
            Keq = K_equilibrium(DeltaG, eta, DeltaGATP, R, T)
            #Calculate thetas
            theta = Theta(Q, Keq)
            #Calculate rates
            q_reac = rate(q_max[j], theta, ks[j], kr[j], S)
            #Turn nans to 0
            nans = np.isnan(q_reac)
            q_reac[nans] = 0
            #Include reaction rates in reaction network matrix
            reaction_rates[j][tot_reac_network[j]] = q_reac
            #Calculate growth
            Jgrow = jgrow(reaction_rates[j], Eta[j])
            #Calculate flux in - flux out
            flux_in_out[j,:] = vin_out(reaction_rates[j])
            #Calculate Growth and Maintenance vectors
            G[j] = Grow(g[j], N[j, i-1], Jgrow)
            M[j] = Maintenance(g[j], N[j, i-1], maintenance[j])
            ###################################################################

        #2. Integrate model
        ###################

        #Initial conditions for poplation and concentrations
        z0 = list(N[:, i-1])+list(C[:, i-1])
        #Integration
        z = odeint(model, z0, tspan, args=(s,m,G,M,kappa,gamma,flux_in_out))
        #Store solutions 
        N[:,i] = z[1][0:s]
        C[:,i] = z[1][s:s+m]
        #Transform negative elements to 0
        C[:,i] = C[:,i].clip(min=0)
        #Next initial condition
        z0 = z[1]

    return N, C
