#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
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

## FUNCTIONS ##

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


    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

