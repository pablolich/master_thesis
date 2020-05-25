import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#Functions to calculate quantities in the model

def r_quotient(S, P):
    '''Computes reaction quotient from concentration of (S)ubstrates and
    (P)roducts'''
    Q = P/S
    #Turn nans to 0
    nans = np.isnan(Q)
    Q[nans] = 0
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

def model(z, t, s, m, G, M, kappa, gamma, vin_out):
    '''Implementation of the model'''

    #Variables
    N = z[0:s]
    C = z[s:m+s]

    #Model equations
    dNdt = list(G - M)
    #Reshape N in order to multiply the flux vectors
    N_mul = np.repeat(N, m).reshape(s, m)
    dCdt = list(kappa - gamma*C + sum(vin_out*N_mul))
    dzdt = dNdt + dCdt

    return dzdt 
