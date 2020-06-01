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

def trunc_normal(mu, sigma, a, b):
    '''Sample from a [a, b]-truncated normal distribution N(mu, sigma)'''
    q = np.random.normal(mu, sigma)
    while (q < a) | (q > b):
        q = (np.random.normal(1, sigma))
    return q

def diag_element(m, s):
    '''

    Chooses a matrix element of the form (i, i+k), where k represents the
    kth upper diagonal, and its sampled from a truncated normal distribution,
    and i one of the kth diagonal elements, and its sampled from a uniform 
    distribution.

    Parameters:
        m (int): number of metabolites.
        s (float): shrinkage factor of normal distribution's variance. 

    Returns:
        (2x1 array): indices of the ith element from the kth diagonal.

    '''

    #Pick k from a normal truncated distribution between [1, m-1], N(1, sigma)
    mu = 1
    sigma = np.sqrt(m-1)/s
    a = 1
    b = m-1
    sample = trunc_normal(mu, sigma, a - 0.5, b + 0.5)
    k = round(sample)
    #Pick i from a uniform distribution U(0, m-k-1)
    i = np.random.randint(m-k)
    #Return matrix element (i, i+k)
    return np.array([[i],[i+k]])

##Test of the above function##
#reac_network = np.array([[], []], dtype = 'int')
#m = 50
#n_reac = 5000
#for i in range(n_reac):
#    reac = diag_element(m, 1)
#    reac_network = np.concatenate((reac_network, reac),axis = 1)
#reac_network = tuple(reac_network)
#reactions = np.zeros(shape = (m,m))
#reactions[reac_network] = 1
#import matplotlib.pylab as plt
#plt.imshow(reactions)
#plt.show()

def network(m, n_reac, s = 1):
    '''
    Creates reaction network by sampling n_reac. Note that the number of 
    reactions in the network, N_reac <= n_reac, due to repeated reactions. 

    Parameters:
        m (int): number of metabolites.
        n_reac (int): number of reactions to be sampled.
        s (float): scale of normal distribution of the reaction network

    Returns: 
        network (tuple): reaction network.
    '''

    #Counter for number of reactions
    n = 0
    #Network container
    network = np.array([[],[]], dtype = int)
    while n < n_reac:
        #Add 1 sampled reaction to the network
        network = np.concatenate((network, diag_element(m, s)), axis = 1)
        n += 1
    #Eliminate repeated reactions
    network = np.unique(network, axis = 1)
    return tuple(network)

##Test network
#network = network(10, 10)
#reactions = np.zeros(shape = (10, 10))
#reactions[network] = 1
#
#import matplotlib.pylab as plt
#plt.imshow(reactions)
#plt.show()

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
    Integrate Jacob's model.

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

    return N, C
