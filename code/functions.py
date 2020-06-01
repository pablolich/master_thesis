import numpy as np
from scipy.integrate import solve_ivp, odeint
import matplotlib.pyplot as plt

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; DeltaGATP = 75e3 # J/mol
global T; T = 298 # K 
global E = 3e6 # J (Total energy contained in metabolites)

#Functions to calculate quantities in the model

def r_quotient(S, P, Keq):
    '''Computes reaction quotient from concentration of (S)ubstrates and
    (P)roducts'''
    Q = P/S
    #Turn nans to 0
    nans = np.isnan(Q)
    Knans = Keq[nans]
    Q[nans] = Knans
    return (Q)

def K_equilibrium(DeltaG, eta, DeltaGATP, R, T):
    '''Computes equilibrium constant from gibbs energies, R and T'''
    return (np.exp((-DeltaG - eta*DeltaGATP)/(R*T)))

def Theta(Q, Keq):
    '''Computes theta: how far is my reaction from equilibrium'''
    return (Q/Keq)

def rate(q_max, theta, ks, kr, S):
    '''Computes rate of reaction'''
    return ((q_max*S*(1-theta))/(ks + S*(1+kr*theta)))

def energies(eta, q):
    '''Computes matrix of energies available for each one-step biochemical 
    reaction that strain s can possess'''
    return (np.multiply(eta, q))

def jgrow(R, Eta):
    '''Calculates the energy available for growth'''
    return (np.matmul(R.transpose(),Eta).trace())

def vin_out(R):
    '''Calculates the difference between rates of production and consumption 
    of metabolite beta'''
    return ((R.transpose()-R).sum(axis=1))

def Grow(g, N, Jgrow):
    '''Calculates total energy available'''
    return (g*N*Jgrow)

def Maintenance(g, N, m):
    '''Calculates energy required for maintenance'''
    return (g*N*m)

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

def model(t, z, s, m, kappa, gamma, networks, mu, Eta, q_max, ks, kr, g, maint):
    '''
    Implementation of Jacob Cook's model.

    Parameters:
        t (tuple): time span 
        s (int): number of strains
        m (int): number of metabolites
        networks (s tuples): reaction network of each strain. Each 
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
        maint (1xs array): maintenance energy requirement of each strain.
        kappa (1xm array): supply rate of each metabolite.
        gamma (1xm array): dilution rate of each metabolite.
        N (sxn array): population vector time series for each strain.
        C (mxn array): metabolite concentration time series vector for each 
                       metabolite.

    '''

    #Separate variables
    N = z[0:s]
    C = z[s:m+s]

    #1. Calculate  rates and asociated quatnities
    #############################################

    #Prealocate elements that will be used later on
    G = np.zeros(s)
    M = np.zeros(s)
    flux_in_out = np.zeros([s, m])
    rates = np.zeros(shape = (s,m,m))
    #Loop through all strains
    for j in range(s):
        #Get concentrations from those metabolites taking part in reaction 
        #network 
        S = C[networks[j][0]]
        P = C[networks[j][1]]
        #The following calculations are performed for all reactions at once
        ###################################################################
        #Get chemical potentials of those metabolites
        mu_S = mu[networks[j][0]]
        mu_P = mu[networks[j][1]]
        #Gibs free energy change is the difference between chemical 
        #potentials (stechiometric coefficients are all 1)
        DeltaG = mu_P - mu_S
        #Get etas for reaction network
        eta = Eta[j][networks[j]]
        #Calculate equilibrium constant
        Keq = K_equilibrium(DeltaG, eta, DeltaGATP, R, T)
        #Calculate reaction quotients
        Q = r_quotient(S, P, Keq)
        #Calculate thetas
        theta = Theta(Q, Keq)
        #Calculate rates
        q_reac = rate(q_max[j], theta, ks[j], kr[j], S)
        #Turn nans to 0
        nans = np.isnan(q_reac)
        q_reac[nans] = 0
        #Include reaction rates in reaction network matrix
        rates[j][networks[j]] = q_reac
        #Calculate growth
        Jgrow = jgrow(rates[j], Eta[j])
        #Calculate flux in - flux out
        flux_in_out[j,:] = vin_out(rates[j])
        #Calculate Growth and Maintenance vectors
        G[j] = Grow(g[j], N[j], Jgrow)
        M[j] = Maintenance(g[j], N[j], maintenance[j])


    #Vectorized model equations
    dNdt = list(G - M)
    #Reshape N in order to multiply the flux vectors
    N_mul = np.repeat(N, m).reshape(s, m)
    dCdt = list(kappa - gamma*C + sum(flux_in_out*N_mul))
    dzdt = dNdt + dCdt

    return dzdt 
