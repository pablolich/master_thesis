import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#Functions to calculate quantities in the model

def r_cocient(S, P):
    '''Computes reaction cocient from concentration of (S)ubstrates and
    (P)roducts'''
    return (P/S)

def K_equilibrium(DeltaG, eta, DeltaGATP, R, T):
    '''Computes equilibrium constant from gibbs energies, R and T'''
    return (np.exp((-DeltaG - eta*DeltaGATP)/(R*T)))

def theta(Q, Keq):
    '''Computes theta: how far is my reaction from equilibrium'''
    return (Q/Keq)


def rate(q_max, theta, ks, kr, C):
    '''Computes rate of reaction'''
    return ((q_max*C*(1-theta))/(ks + C*(1+kr*theta)))

def energies(eta, q):
    '''Computes matrix of energies available for each one-step biochemical 
    reaction that strain s can possess'''
    return (np.multiply(eta, q))

def jgrow(R, Eta):
    '''Calculates the energy available for growth'''
    return ((R.transpose()*Eta).trace())

def vin_out(R, Eta, beta):
    '''Calculates the difference between rates of production and consumption 
    of metabolite beta'''
    return ((R.transpose()-R).sum(axis=1))

def Grow(g, N, Jgrow):
    '''Calculates total energy available'''
    return (g*N*Jgrow)

def Maintenance(g, N, m):
    '''Calculates energy required for maintenance'''
    return (g*N*m)

def model(z, t, par):
    '''Implementation of the model'''

    #Variables
    N = z[0]
    C = z[1]

    #Parameters
    G = par[0]
    M = par[2]
    kappa = par[3]
    gamma = par[4]
    vin_out = par[5]

    #Model equations
    dNdt = G - M 
    dCdt = kappa - gamma*C + sum(vin_out)*N
    dzdt = [dNdt, dCdt]

    return dzdt 
