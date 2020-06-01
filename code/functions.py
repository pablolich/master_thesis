import numpy as np
from scipy.integrate import solve_ivp, odeint
import matplotlib.pyplot as plt

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; DeltaGATP = 75e3 # J/mol
global T; T = 298 # K 

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

def model(t, z, s, m, kappa, gamma, networks, mu, Eta, q_max, ks, kr, g, 
          maintenance):
    '''Implementation of the model'''

    #Variables
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


    #Model equations
    dNdt = list(G - M)
    #Reshape N in order to multiply the flux vectors
    N_mul = np.repeat(N, m).reshape(s, m)
    dCdt = list(kappa - gamma*C + sum(flux_in_out*N_mul))
    dzdt = dNdt + dCdt

    return dzdt 
