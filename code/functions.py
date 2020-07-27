import numpy as np
from scipy.integrate import solve_ivp, odeint
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from itertools import combinations, product, permutations
from scipy.special import comb
from progressbar import ProgressBar

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; DeltaGATP = 75e3 # J/mol
global T; T = 298 # K 
global E; E = 3e6 # J (Total energy contained in metabolites)

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

def tuple2matrix(tuples, abundances, s, m):
    '''Transform a list of tuples into adjacency matrices''' 
    #Initialize s mxm matrices to store weighted adjacency matrix of each 
    #reaction network
    network_adj = np.zeros(shape = (s, m, m))
    #Transform tuples to matrices
    for j in range(len(tuples)):
        network_adj[j][tuples[j]] = abundances[j]
    return network_adj



def networks2community(networks, abundances, s, m):
    '''
    Transform a list of tuples into adjacency matrices and sum them weighted
    by abundances

    Parameters: 
        networks (list of tuples): List of s tuples storing reaction network of
                                   each strain.
        abundances (1xs array): Array with abundance of each strain.

    '''
    network_adj = tuple2matrix(networks, abundances, s, m)
    #Add up all networks to form the community adjacency matrix
    community_adj = np.sum(network_adj, axis = 0)
    community_nx = nx.from_numpy_matrix(community_adj, 
                                        create_using = nx.DiGraph)
    return community_nx



def model(t, z, s, m, kappa, gamma, networks, mu, Eta, q_max, ks, kr, g, maint):
    '''
    Implementation of Jacob Cook's model.

    Parameters:
        t (tuple): time span 
        z (s+m list): concatenated N [(sxn array): population for each strain.]
                      and C [(mxn array): concentration for each metabolite.]
        s (int): number of strains.
        m (int): number of metabolites.
        kappa (1xm array): supply rate of each metabolite.
        gamma (1xm array): dilution rate of each metabolite.
        networks (s tuples): reaction network of each strain. Each 
                             tuple stores two arrays for substrates and
                             products.
        mu (1xm array): chemical potential of each metabolite.
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
        M[j] = Maintenance(g[j], N[j], maint[j])


    #Vectorized model equations
    dNdt = list(G - M)
    #Reshape N in order to multiply the flux vectors
    N_mul = np.repeat(N, m).reshape(s, m)
    dCdt = list(kappa - gamma*C + sum(flux_in_out*N_mul))
    dzdt = dNdt + dCdt

    return dzdt 

def F_calculator(z, t, m, s, net):
    '''Calculate F = -log(R/T)'''
    t_points = len(t)
    #Create community network and evaluate centrality measures for each
    #timepoint
    F = np.zeros(t_points)
    for i in range(t_points):
        #Don't worry about divergent cases
        if any(z[:, -1]>1e4):
            pass
        else:
            #Abundance of each strain at timepoint i
            abundances = z[0:s,i]
            #Metabolite concentrations at timepoint i
            concentrations = z[s:s+m, i]
            #Number of individuals demanding each resource at timepoint i
            sub = [int(round(abundances[p]))*list(net[p][0]) for 
                    p in range(s)]
            flat_sub = [item for sublist in sub for item in sublist]
            #Get frequency of consumption of each metabolite
            T = np.bincount(flat_sub)
            #Remove 0 from T
            T = T[T != 0]
            #Identify metabolites that are not being consumed at all
            missing = np.setxor1d(np.array(flat_sub), np.arange(m))
            #Remove unused metabolites from vector of concentrations
            R = np.delete(concentrations, missing)
            #Calculate efficiency in simultaneous depletion
            F[i] = sum(np.log(abs(R)/T))
    return F

def coalescence_event(C1, C2, m, s):
    '''Perform the coalescence of two communities'''

    #Put together both communities
    mixed = pd.concat([C1, C2]).reset_index(drop = True)

    #Get parameters for integration
    #Obtain reaction networks as a list of tuples
    substrates = mixed['substrate']
    products = mixed['product']
    nets = vector2tuple(substrates, products)
    #Choose Eta matrix for each strain
    Eta = 0.5*np.ones(shape = (s,m,m))
    diag_ind = np.diag_indices(m)
    for i in range(s):
        Eta[i][diag_ind] = 0
    q_m = [np.ones(1) for i in range(s)]
    ks = [0.1*q_m[i] for i in range(s)]
    kr = [10*q_m[i] for i in range(s)]

    #Initial conditions                
    z0 = list(mixed['stable.state']) + m*[0]

    #Timespan
    tspan = tuple([1, 1e5])

    #Integrate model
    sol = solve_ivp(lambda t,z: model(t, z, s = len(mixed), m = m, 
                                      kappa = 2*np.ones(m), 
                                      gamma = 1*np.ones(m), 
                                      networks = nets, 
                                      #mu = uniform_pack(E, m),
                                      mu = decreasing_pack(E, m), 
                                      Eta = Eta, 
                                      q_max = q_m, 
                                      ks = ks, 
                                      kr = kr,
                                      g = np.ones(s),
                                      maint = np.array(mixed['maintenance'])),
                    tspan, z0, method = 'BDF', atol = 1e-3)

    #Unpack results
    t = sol.t
    z = sol.y

    return t, z, nets

def random_pack(E, m):
    '''Generate random uniformly distributed chemical potentials'''
    #Pack total energy into m intervals (chemical potentials of metabolites)
    pack = np.random.uniform(E, size = m)
    #Sort decreasingly 
    mu = np.sort(pack)[::-1]
    #Fix first and last chemical potentials
    mu[0] = E 
    mu[-1] = 0
    return mu

def uniform_pack(E, m):
    '''Generate equidistant chemical potentials'''
    return np.linspace(E,0,m)

def second_order_polynomial(x):
        return(-x**0.5)

def decreasing_pack(E, m):
    '''Generate quadraticaly decreasing chemical potentials'''
    #x = np.linspace(0,1,m)
    #import ipdb; ipdb.set_trace(context = 20)
    #distances = second_order_polynomial(x)+1
    #d_norm = E*distances/sum(distances)
    #mu = [E - sum(d_norm[0:i]) for i in range(len(d_norm))]
    x = np.arange(m)
    mu = E*(1 - np.sqrt(x/(m)))
    return mu

def rate_matrix(network, Eta, q_max, ks, kr, C, mu, m):
    '''Calculate matrix of reaction rates'''

    #Initialize matrix of rates
    r = np.zeros(shape=(m,m))
    #Get concentrations from those metabolites taking part in reaction 
    #network 
    S = C[network[0]]
    P = C[network[1]]
    #Get chemical potentials of those metabolites
    mu_S = mu[network[0]]
    mu_P = mu[network[1]]
    #Gibs free energy change is the difference between chemical 
    #potentials (stechiometric coefficients are all 1)
    DeltaG = mu_P - mu_S
    #Get etas for reaction network
    eta = Eta[network]
    #Calculate equilibrium constant
    Keq = K_equilibrium(DeltaG, eta, DeltaGATP, R, T)
    #Calculate reaction quotients
    Q = r_quotient(S, P, Keq)
    #Calculate thetas
    theta = Theta(Q, Keq)
    #Calculate rates
    q_reac = rate(q_max, theta, ks, kr, S)
    #Turn nans to 0
    nans = np.isnan(q_reac)
    q_reac[nans] = 0
    #Include reaction rates in reaction network matrix
    r[network] = q_reac
    return r

#def cooperation_index(adj_matrix):
#    '''
#    Calculate the cooperation index from the adjacency matrix
#
#    Parameters: 
#        
#        adj_matrix (mxm array): adjacency matrix
#
#    Returns:
#         
#        coop_idx (float): Cooperation index. A measure of how similar is the 
#                          adjacency matrix to a matrix where all the rows have 
#                          none or one non-zero element.
#    '''
#    #Get dimension of matrix
#    dimension = len(adj_matrix)
#    #Prealocate cooperation index vector for adjacency matrix rows
#    coop_ind_vec = np.zeros(dimension)
#    #Loop through rows of matrix
#    for i in range(dimension):
#        #Number of non_zero elements
#        non_zero = sum(np.nonzero(adj_matrix[i,:]))
#        #Set cooperation row index to 0 if all elements are 0
#        if  non_zero == 0:
#            coop_ind_vec[i] = 0
#        else:
#            #Calculate an index of cooperation
#            coop_ind_vec[i] = 1 - non_zero/dimension
#    coop_idx = mean(coop_ind_vec)
#    return coop_idx

def facilitation_matrix(networks, m):
    '''Create facilitation matrix'''
    
    #Get number of strains
    s = len(networks)
    #Get matrix indices (species pairs)
    pairs = permutations(range(s),2)
    #Preallocate facilitation matrix
    fac_mat = np.zeros(shape=(s,s))
    for i, j in pairs:
        #Get products of species j and substrates by species i
        sub_i = networks[i][0]
        prod_j = networks[j][1]
        #Get shared elements
        shared = np.intersect1d(prod_j, sub_i)
        #Calculate reliability of strain i in each of its substrates
        rel_subs = reliability_dict(sub_i, m, facilitation = False)
        #Get reliability in facilitated substrates for species i
        rel_fac = [rel_subs[str(i)] for i in shared] 
        #Get how much of the metabolites of the other are being facilitated
        rel_prods = reliability_dict(prod_j, m, facilitation = False)
        #Get reliablility in facilitated products for species j
        rel_fac_prod = [rel_prods[str(i)] for i in shared]
        #Sum to get interaction level
        tot_facilitation = 0.5*(np.array(rel_fac) + np.array(rel_fac_prod))
        #Facilitation degree between species i and j
        fac_mat[i,j] = sum(tot_facilitation) 

    return fac_mat

def competition_matrix(networks, m):
    '''Create competition matrix'''
    
    #Get number of strains
    s = len(networks)
    #Get matrix indices (species pairs)
    pairs = permutations(range(s), 2)
    #Preallocate facilitation matrix
    comp_mat = np.zeros(shape=(s,s))
    for i, j in pairs:
        #get used substrates by species i and products used by species j
        sub_i = networks[i][0]
        sub_j = networks[j][0]
        #get shared elements
        shared = np.intersect1d(sub_i, sub_j)
        #calculate reliability of strain i in each of its metaboltes
        rel = reliability_dict(sub_i, m)
        #get reliability in competing metabolites with species j
        rel_comp = [rel[str(i)] for i in shared] 
        #calculate reliablitiy of strain j in each of its metabllites
        rel_j = reliability_dict(sub_j, m)
        #get reliability in competing metabolites with species i
        rel_comp_j = [rel_j[str(i)] for i in shared] 
        #sum to get total competition level, and normalize
        rel_comp_tot = 0.5*(np.array(rel_comp) + np.array(rel_comp_j))
        #Facilitation degree between species i and j
        comp_mat[i,j] = sum(rel_comp_tot)

    return comp_mat

def reliability_dict(vec, m, facilitation = False):
    '''Calculates the normalized index of repetition of each element in vec'''
    if facilitation:
        #Eliminate elements of vector that cannot ever be facilitated
        #ie, first and last metabolite.
        vec = np.delete(vec, np.where((vec == 0) | (vec == m-1))[0])
    #Initialize keys of dictionary
    keys = [str(i) for i in np.unique(vec)]
    #Calculate number of repetitions of each element
    values = np.bincount(vec)
    #Eliminate 0s
    values = values[values != 0]/len(vec)
    #Create dictionary
    reliability = dict(zip(keys, values))
    return reliability 

def individual_fitness(_in, out_, competition):
    return (_in - out_)/competition
    
def cooperation_index(cooperation_matrix):
    '''
    Calculate the cooperation index from the cooperation matrix. It is a 
    measure of how overlapping the products of species are with the substrates
    of others. 

    Parameters:
        
        cooperation_matrix (sxs array): Matrix of cooperation between all
                                        possible pair of species

    Returns:
         
       cooperation_idx (float)
     '''

    #Get number of strains
    s = len(cooperation_matrix)
    cooperation_idx = np.zeros(s)
    #Get autocooperation of each species
    autocoop = np.diag(cooperation_matrix)
    for i in range(s):
        cooperation_idx[i] = sum(cooperation_matrix[i,:]) + \
                             sum(cooperation_matrix[:,i]) + \
                             2 * cooperation_matrix[i,i] 
    return cooperation_idx

def competition_index(comp_matrix):
    return sum(comp_matrix)

def autofacilitation_index(networks, m):
    '''
    Calcualte autofacilitation index, a measure of how many elements are both 
    substrates and products in the reaction network of a given strain
    '''

    #Get number of strains
    s = len(networks)
    #Preallocate vector of autofacilitation indexes
    auto_idx = np.zeros(s)
    #Loop through strains
    for i in range(s):
        #Facilitation degree
        substrates = np.unique(networks[i][0])
        products = np.unique(networks[i][1])
        shared = np.intersect1d(substrates, products)
        #Calculate autofacilitation index for strain i
        auto_idx[i] = len(shared)/(comb(m, 2, exact = True)-(m-1))

    return np.mean(auto_idx)
    
def vector2tuple(substrate, product):
    '''
    Transform a matrix of substrates and products for a community into a 
    list of tuples of arrays

    Parameters:
        
        substrate (Pandas.Series): Column of substrates for all strains.
        product (Pandas.Series): Column of products for all strains.

    Returns:

        network (list of tuples): list of reaction networks of each strain in 
        tuples.
    '''
    nets = [tuple([substrate.loc[i], product.loc[i]]) 
                   for i in range(len(substrate))]

    return nets

def evol_richness(abundances_t):
    '''Get evolution of richness as time goes by'''
    it = len(abundances_t)
    richness = np.zeros(it)
    for i in range(it):
        richness[i] = sum(abundances_t.iloc[i,2:]>=1)
    return richness


def interaction_evolution(time_series, networks, m, t_points):
    '''Calculate matrix of interactions at every point in time'''
    s = len(networks)
    it = len(t_points)
    inter_evol = np.zeros(shape = (it, s))
    f_evol = np.zeros(shape = (it, s))
    c_evol = np.zeros(shape = (it, s))
    #Transform from string to numeric vector all elements in substrate and 
    #product colums
    networks['substrate'] = networks['substrate'].apply(lambda x: 
                            np.fromstring(x[1:-1], sep=' ', dtype = int))
    networks['product'] = networks['product'].apply(lambda x: 
                          np.fromstring(x[1:-1], sep=' ', dtype = int))
    #Calculate community interaction measures at extinction points. 
    j = 0
    for i in t_points:
        #Get present and extinct species, those that are greater than 1e-3
        present = np.where(time_series.iloc[i+1,2:]>=1)[0]
        extinct = np.where(time_series.iloc[i+1,2:]<1)[0]
        #Get networks from the present species
        networks_i = networks.iloc[present].reset_index(drop = True)
        #Get a tuple for reaction networks
        networks_tu = vector2tuple(networks_i['substrate'], 
                                   networks_i['product'])
        #Calculate interaction matrix for this set of reaction networks
        f_mat = facilitation_matrix(networks_tu, m)
        c_mat = competition_matrix(networks_tu, m)
        interactions = f_mat - c_mat
        f_evol[j,present] = np.mean(f_mat, axis = 1)
        c_evol[j,present] = np.mean(c_mat, axis = 1)
        inter_evol[j,present] = np.mean(interactions, axis = 1)
        j += 1
    return(f_evol, c_evol, inter_evol)






