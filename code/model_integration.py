#!/usr/bin/env python3

__appname__ = '[model_integration]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import matplotlib.pylab as plt
from model_functions import *
from scipy.special import comb
from matplotlib.lines import Line2D

## CONSTANTS ##

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; DeltaGATP = 75e3 # J/mol
global T; T = 298 # K 

## FUNCTIONS ##
def eliminate_posibilities(substrate, products, network, tot_metabolites):
    #Get indices of reactions in the network already with that substrate
    ind = np.where(network[0,:] == substrate)[0]
    #Get products of those reactions
    products_ind = network[1,ind]
    #Get indices of original product vector where those products are
    remove_ind = np.where(products == products_ind)[0]
    #Remove those indices
    products = np.delete(products, remove_ind)
    if not len(products):
        import ipdb; ipdb.set_trace(context = 20)
    return products


def get_reac_network(s, m,  nATP, mu):
    '''Create reaction network of feasible reactions (those which associated
    Gibbs energy change is lower than nATP*D(G_ATP)) based on the
    number of metabolites and microbial strains'''

    #Maximum number of reactions in the network
    n_max = comb(m, 2, exact = True)
    #Prealocate the reaction matrix
    reactions = np.zeros(shape = (s,m,m))
    #Choose number of reactions in the network
    reac_network = np.array([[], []], dtype = 'int32')
    num_reac = np.random.randint(1, n_max + 1)
    #Initialize reaction network container 
    tot_reac_network = np.empty(shape = (2,num_reac), dtype = int)
    #Create reaction network for microbial strain i
    for i in range(s):
        print(i)
        #The number of reactions in the network fixes a subset of metabolites
        #from which the first substrate is chosen. The upper limit of
        #the set of possible first metabolites is calculated according to the
        #notes on this matter. 
        m_up = np.floor(m + 0.5*(1 - np.sqrt(1+8*num_reac)))
        #Create array of first metabolite candidates
        list_first = np.arange(m_up, dtype = int)
        #Choose the first metabolite to start reactions with and include it in
        #the list of present substrates
        list_m_i = [np.random.choice(list_first)] 
        #Initialize counter for number of reactions
        n = 0
        while n < num_reac: 
            print(n)
            #Sample one reaction from all the posibilities
            #Choose a substrate from list of present substrates
            m_i = np.random.choice(list_m_i)
            #Create list of possible products compatible where Delta G < 0 and 
            #the energetic distance between product and substrate is less than
            #nATP*DeltaGATP
            #Vector of energy distances between selected substrate and all
            #metabolites
            d_energy = mu - mu[m_i]
            #Apply first condition, only Delta G < 0 
            m_1 = np.where(d_energy < 0)[0]
            #Apply second condition distance < nATP*DeltaGATP
            products = m_1[list(abs(d_energy[m_1]) < nATP*DeltaGATP)]
            #Eliminate products from list of products that would create a 
            #repeated reaction
            products = eliminate_posibilities(m_i, products, tot_reac_network, 
                                              m)
            print(products)
            m_j = np.random.choice(products)
            #Pick one product from list
            #Add reaction to tot_reac_network
            tot_reac_network[:,n] = [m_i, m_j]
            #Check if the added column is repeated by getting indices of unique
            #columns
            el, ind = np.unique(tot_reac_network, return_index = True,
                                axis = 1)
            #If it isn't repeated, keep going onto the next column
            if n in ind: n += 1
            #Otherwise, set it to 0 and sample a reaction again for that column
            else: tot_reac_network[:,n] = [0, 0]
            #Add product to list of substrates only if the added product is 
            #not the last metabolite, because it can't be a substrate
            if m_j != m-1:
                list_m_i.append(m_j)
        
        #Transform to indices in reactions matrix
        reac_network = tuple(reac_network)
        #Switch on the matrix elements in reactions contained in reac_network
        tot_reac_network[i] = reac_network
        reactions[i][reac_network] = 1
        #Reset reaction network for next iteration.
        reac_network = np.array([[], []], dtype = 'int32')

    return(reactions, tot_reac_network)

def main(argv):
    '''Main function'''

    #Number of timepoints
    n = 100
    #Timepoints
    t = np.linspace(1, 30, n)
    #Number of metabolites
    m = 100
    #Chemical potentials 
    #Total energy contained in metabolitess
    interval = 3e6
    #Sort decreasingly 
    mu_rev = np.sort(np.random.uniform(interval, size = m))
    mu = mu_rev[::-1]
    #Fix first and last chemical potentials
    mu[0] = interval 
    mu[-1] = 0
    nATP = 4
    #Number of posible reactions given thermodynamic constraint of only
    #forward net reactions taking place
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
    
    
    ################
    #SET PARAMETERS#
    ################

    #1. Choose Eta matrix
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

    #2. Choose reactions matrix for each strain
    ###########################################
    reactions, tot_reac_network = get_reac_network(s, m, nATP, mu)
    #plt.imshow(sum(reactions)/s)
    #plt.colorbar()
    #plt.show()

    #Numerically integrate the model
    for i in range(1, n):
        #Declare timespan for this iteration of numerical integration
        tspan = [t[i-1], t[i]]

        #2. Rates are calculated according to enzyme kinetics
        #####################################################
        #Note that rates change after every iteration of the integration
        #because the concentrations change, that is why they need to be 
        #recalculated after each integration step.

        #These 2 lines will be useful to present on monday
        #plt.imshow(reactions)
        #plt.show()

        #2b. Calculate rates
        #Get rates of each reaction in the network
        #Get concentrations metabolites from vector C taking part in reaction 
        #networks as substrate or products
        #Prealocate elements that will be used later on
        G = np.zeros(s)
        M = np.zeros(s)
        flux_in_out = np.zeros([s, m])
        for j in range(s):
            S = C[tot_reac_network[j][0], i-1]
            P = C[tot_reac_network[j][1], i-1]
            #The following calculations are performed for all reactions at once.
            #Calculate reaction quotients
            Q = r_quotient(S, P)
            #Calculate Gibbs free energy changes
            mu_S = mu[tot_reac_network[j][0]]
            mu_P = mu[tot_reac_network[j][1]]
            #Gibs energy change is the difference between chemical potentials 
            #(stechiometric coefficients are all 1)
            DeltaG = mu_P - mu_S
            #Calculate equilibrium constant
            #Get etas for reaction network
            eta = Eta[j][tot_reac_network[j]]
            Keq = K_equilibrium(DeltaG, eta, DeltaGATP, R, T)
            #Calculate thetas
            theta = Theta(Q, Keq)
            #Calculate rates
            q_max = np.ones(len(tot_reac_network[j][0]))
            ks = 0.1*q_max
            kr = 10*q_max
            #Include reaction rates in reaction network matrix
            q_reac = rate(q_max, theta, ks, kr, S)
            #Turn nans to 0
            nans = np.isnan(q_reac)
            q_reac[nans] = 0
            reactions[j][tot_reac_network[j]] = q_reac
            #Calculate growth
            Jgrow = jgrow(reactions[j], Eta[j])
            #Calculate flux in - flux out
            flux_in_out[j,:] = vin_out(reactions[j])
            #Set growth and maintenance parameters
            g = 1 
            maintenance = 0.2*Jgrow
            G[j] = Grow(g, N[j, i-1], Jgrow)
            M[j] = Maintenance(g, N[j, i-1], maintenance)
            
        #Integrate model
        #Initial conditions for poplation and concentrations
        z0 = list(N[:, i-1])+list(C[:, i-1])
        kappa = np.ones(m) 
        gamma = 0.5
        z = odeint(model, z0, tspan, args=(s,m,G,M,kappa,gamma,flux_in_out))
        #Store solutions 
        N[:,i] = z[1][0:s]
        C[:,i] = z[1][s:s+m]
        #Transform negative elements to 0
        C[:,i] = C[:,i].clip(min=0)
        #Next initial condition
        z0 = z[1]
    
    for i in range(len(N)):
        plt.plot(t, N[i])

    for i in range(len(C)):
        plt.plot(t, C[i], linestyle = '--')
        
    custom_lines = [Line2D([0], [0], color = 'black'),
                    Line2D([0], [0], color = 'black', linestyle = '--')]
    plt.legend(custom_lines, ['Population', 'Concentration'])
    plt.xlabel('t')
    plt.show()

    return 0


## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

