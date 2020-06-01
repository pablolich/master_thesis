#!/usr/bin/env python3

__appname__ = '[analysis]'
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
    #Total energy contained in metabolitess
    interval = 3e6
    #Pack total energy into m intervals (chemical potentials of metabolites)
    pack = np.random.uniform(interval, size = m)
    #Sort decreasingly 
    mu = np.sort(pack)[::-1]
    #Fix first and last chemical potentials
    mu[0] = interval 
    mu[-1] = 0
    #Limit energy jump in a single reaction to n times of the energy of 1 mol
    #of ATP
    nATP = 1
    #Total number of posible reactions (only forward net reactions are taking 
    #place)
    n_reac = comb(m, 2, exact = True)
    #Set growth and maintenance parameters
    g = np.ones(s) 
    maintenance = 0.2*np.ones(s)
    #maintenance = np.array([0.1, 0.2, 0.3])
    #Set supply and dilution rate of metabolites
    kappa = 8*np.ones(m) 
    gamma = 1*np.ones(m)


    ###########################
    #3. Biochemical parameters#
    ###########################

    #Choose Eta matrix
    #For now we set all of them to 0.5
    Eta = 0.5*np.ones(shape = (s,m,m))
    diag_ind = np.diag_indices(m)
    for i in range(s):
        Eta[i][diag_ind] = 0
    #Initialize reaction rates matrix
    reaction_rates = np.zeros(shape = (s,m,m))
    #Prealocate the reaction tensor
    reactions = np.zeros(shape = (s,m,m), dtype = 'int')
    #Initialize reaction network container for one strain
    reac_network = np.array([[], []], dtype = 'int32')
    #Prealocate reaction network list where reaction networks for all strains
    #will be stored
    tot_reac_network = s*[reac_network]
    #Prealocate number of reactions possesed by each strain 
    n_reac_s = np.zeros(s, dtype = 'int')
    #Set parameter values that depend on each strain
    q_max = [np.ones(1) for i in range(s)]
    ks = [0.1*q_max[i] for i in range(s)]
    kr = [10*q_max[i] for i in range(s)]
    for i in range(s):
        #Set number of reactions in the network
        #num_reac = np.random.randint(1, n_reac+1)
        num_reac = 10
        #Get reaction network
        reac_network = network(m, 5)
        #Store reaction network information for later use
        tot_reac_network[i] = reac_network
        reactions[i][reac_network] = 1
        non_zero = reactions[i][reactions[i].nonzero()]
        n_reac_s[i] = len(non_zero)


    #########################
    #Integrate model k times#
    #########################
    
    #k_max = 1
    #n = len(t)
    ##Prealocate results storing elements
    #Nn = np.zeros(shape = (k_max, s, n))
    #Cn = np.zeros(shape = (k_max, m, n))
    #Rn = k_max*[0]
    #k = 0
    #Iterate k_max times generating a different reaction network each time
    #while k < k_max:
    #Choose reactions matrix for all strains
    #Create one reaction network for each microbial strain
    #Integrate model
    ################

    #Initial conditions for poplation and concentrations
    N0 = s*[2]
    C0 = m*[0] 
    z0 = N0 + C0
    #Time span
    tspan = tuple([1, 500])
    #Integration
    #Using old API. To use this, switch order of t, z in model function.
    #z = odeint(model, z0, tspan, args=(s,m,G,M,kappa,gamma,flux_in_out))
    #Using new API
    sol = solve_ivp(lambda t,z: model(t,z,s,m,kappa,gamma,tot_reac_network, mu,
                    Eta, q_max, ks, kr, g, maintenance), tspan, z0, 
                    method = 'BDF')
    t = sol.t
    z = sol.y
    #Store solutions 
    N = z[0:s]
    C = z[s:s+m]
    #Find negative elements
    neg = np.where(C[:,i-1] < 0)[0]
    #Set those indices to 0
    #print(C[0, i-5:i+2])
    #print(N[0, i-5:i+2])        
    #C[neg, i-1] = 0
    #Set populations at those indices to the previous one
    #N[neg, i-1] = N[neg, i-2]
    #Transform negative elements to 0
    C[:,i] = C[:,i].clip(min=0)
    #If species population are below 1, declare extinct
    N[:, i][N[:, i] < 1] = 0

    #Save  Rn, Nn, Cn for later analysis
    _type = ['classic', 'new']
    np.save('../data/reactions_' + _type[0] + '.npy', reactions)
    np.save('../data/population_' + _type[0] + '.npy', N)
    np.save('../data/concentration_' + _type[0] + '.npy', C)

    #Some plotting for meeting
    ##########################
    
    #Create panel layout
    G = gridspec.GridSpec(2,4)
    
    plt.figure(figsize = (14,7))
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
     

