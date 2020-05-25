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
def generate_network(s, m, nATP, mu, num_reac):

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


def main(argv):
    '''Main function'''


    ################
    #SET PARAMETERS#
    ################

    #1. Initial conditions
    ######################

    #Number of timepoints
    #Timepoints
    t = np.linspace(1, 100, 1000)
    #Number of metabolites
    m = 5
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
    nATP = 1 
    #Total number of posible reactions if only forward net reactions taking 
    #place
    n_reac = comb(m, 2, exact = True)
    #Number of strains
    s = 10
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
    #Prealocate number of reactions possesed by each strain 
    n_reac_s = np.zeros(s, dtype = 'int')
    #Create one reaction network for each microbial strain
    for i in range(s):
        #Set number of reactions in the network
        num_reac = np.random.randint(1, n_reac+1)
        #Get reaction network
        reac_network = generate_network(s, m, nATP, mu, num_reac)
        #Store reaction network information for later usk
        tot_reac_network[i] = reac_network
        reactions[i][reac_network] = 1
        non_zero = reactions[i][reactions[i].nonzero()]
        n_reac_s[i] = len(non_zero)

    #Initialize reaction rates matrix
    reaction_rates = np.zeros(shape = (s,m,m))

    #Set parameter values that depend on each strain
    q_max = [np.ones(n_reac_s[i]) for i in range(s)]
    ks = [0.1*q_max[i] for i in range(s)]
    kr = [10*q_max[i] for i in range(s)]
    #Set growth and maintenance parameters
    g = np.ones(s) 
    maintenance = 0.2*np.ones(s)
    #Set injection and dilution rate of metabolites
    kappa = 0.1*np.ones(m) 
    gamma = 0.1*np.ones(m)

    ###################
    #MODEL INTEGRATION#
    ###################
    N, C = model_integration(t, s, m, tot_reac_network, mu, Eta, q_max, ks, kr, 
                             reaction_rates, g, N, C, maintenance, kappa, 
                             gamma)


    #Some plotting for meeting
    ##########################
    
    #Create panel layout
    G = gridspec.GridSpec(2,4)
    
    plt.figure(figsize = (14,7))
    #These 2 lines will be useful to present on monday
    mean_reactions = sum(reactions)/len(reactions)
    axes = plt.subplot(G[:,0:2])
    plt.imshow(mean_reactions, extent = [-0.5,m-0.5,m-0.5,-0.5], 
               origin = 'upper')
    plt.colorbar()

    axes2 = plt.subplot(G[:,2:])
    for i in range(len(N)):
        plt.plot(t, N[i], label= 'Strain'+str(i))

    for i in range(len(C)):
        plt.plot(t, C[i], linestyle = '--', label = 'Metabolite'+str(i))
        
    plt.legend()
    plt.show()
    #plt.savefig('../results/reac_dynamics_'+option+'.pdf', bbox_inches='tight')

    return 0


## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

