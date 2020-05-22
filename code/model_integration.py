#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import matplotlib.pylab as plt
from matplotlib import gridspec
from model_functions import *
from scipy.special import comb

## CONSTANTS ##

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; DeltaGATP = 75e3 # J/mol
global T; T = 298 # K 

## FUNCTIONS ##
def generate_network(s, m, nATP, mu):

    '''
    Create reaction network of feasible reactions (those which associated
    Gibbs energy change is lower than nATP*D(G_ATP)) based on the
    number of metabolites,  microbial strains and chemical potentials

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
    keep_sampling = True
    while keep_sampling: 
        #Choose a substrate from list of present substrates (at the start
        #only substrate 0 is present, so the network will always start at
        #0)
        m_i = int(np.random.choice(list_m_i))
        #Create list of possible products of the reaction compatible where
        #Delta G < 0
        products = np.arange(m_i+1, m, dtype = int)
        #Note that we start from m_i to satisfy the latter condition. 
        #doing m_i + 1 avoids having reactions i-->i
        #Choose a product from list of metabolites
        m_j = int(np.random.choice(products))
        #Create the tuple representing the sampled reaction
        r_sample = np.array([[m_i], [m_j]])
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
            if m-1 in reac_network[1]:
                keep_sampling = False
    reac_network= tuple(reac_network)
    
    return reac_network


def main(argv):
    '''Main function'''


    ################
    #SET PARAMETERS#
    ################

    #1. Initial conditions
    ######################

    #Number of timepoints
    n = 1000
    #Timepoints
    t = np.linspace(1, 80, n)
    #Number of metabolites
    m = 3
    #Chemical potentials 
    #Total energy contained in metabolitess
    interval = 1e5
    #Pack total energy into m intervals (chemical potentials of metabolites)
    pack = np.random.uniform(interval, size = m)
    #Sort decreasingly 
    mu = np.sort(pack)[::-1]
    #Fix first and last chemical potentials
    mu[0] = interval 
    mu[-1] = 0
    nATP = 10
    #Total number of posible reactions if only forward net reactions taking 
    #place
    n_reac = comb(m, 2, exact = True)
    #Number of strains
    s = 2
    #Initial conditions
    #All metabolite concentrations equal 0, and all strains have 1 individual
    Cinit = [0]*m
    Ninit = [1]*s
    #Prealocate matrices for storing solution's time evolution
    C = np.zeros(shape = (m, len(t)))
    N = np.zeros(shape = (s, len(t)))
    #Add initial conditions to the first column of the solution's time series
    C[:,0] = Cinit
    N[:,0] = Ninit

    #2. Choose Eta matrix
    #####################

    #For now we set all of them to 0.5
    Eta = np.ones(shape = (s,m,m))
    #Get indices of upper triangular elements with an offset of 1 to exclude
    #diagonal elements. This will constrain reactions going forward
    tri_ind = np.triu_indices(m, 1)
    #Randomly (not for now) choose etas
    eta_set = [0.5]*n_reac
    #Set upper diagonal values of Eta matrixes to chosen etas for the s o
    #strains
    for i in range(s):
        Eta[i][tri_ind] = eta_set

    #3. Choose reactions matrix for all strains
    ###########################################
    
    #Prealocate the reaction tensor
    reactions = np.zeros(shape = (s,m,m), dtype = 'int')
    #Initialize reaction network container for one strain
    reac_network = np.array([[], []], dtype = 'int32')
    #Prealocate reaction network list where reaction networks for all strains
    #will be stored
    tot_reac_network = s*[reac_network]
    #Create one reaction network for each microbial strain
    for i in range(s):
        import ipdb; ipdb.set_trace(context = 20)
        reac_network = generate_network(s, m, nATP, mu)
        tot_reac_network[i] = reac_network
        reactions[i][reac_network] = 1

    #Initialize reaction rates matrix
    reaction_rates = np.zeros(shape = (s,m,m))

    #Set parameter values that depend on each strain
    q_max = np.ones(len(tot_reac_network[j][0]))
    ks = 0.1*q_max
    kr = 10*q_max
    #Set growth and maintenance parameters
    g = 1 
    maintenance = 0.2
    #Set injection and dilution rate of metabolites
    kappa = 0.05*np.ones(m) 
    gamma = 0.5

    ###################
    #MODEL INTEGRATION#
    ###################

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
            #The following calculations are performed for all reactions at once.
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
            q_reac = rate(q_max, theta, ks, kr, S)
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
            G[j] = Grow(g, N[j, i-1], Jgrow)
            M[j] = Maintenance(g, N[j, i-1], maintenance)
            
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
    
    #Some plotting for meeting
    ##########################
    
    #Create panel layout
    G = gridspec.GridSpec(2,3)
    
    plt.figure(figsize = (7,3.5))
    #These 2 lines will be useful to present on monday
    for i in range(s):
        axes = plt.subplot(G[i,0])
        plt.imshow(reactions[i])

    axes2 = plt.subplot(G[:,1:])
    for i in range(len(N)):
        plt.plot(t, N[i], label= 'Strain'+str(i))

    for i in range(len(C)):
        plt.plot(t, C[i], linestyle = '--', label = 'Metabolite'+str(i))
        
    plt.legend()
    plt.savefig('../results/reac_dynamics_'+option+'.pdf', bbox_inches='tight')

    return 0


## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

