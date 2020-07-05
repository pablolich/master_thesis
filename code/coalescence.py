#!/usr/bin/env python3

__appname__ = '[coalescence.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
from functions import *
from scipy.special import comb
from progressbar import ProgressBar
import pandas as pd

## CONSTANTS ##

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; DeltaGATP = 75e3 # J/mol
global T; T = 298 # K 
global E; E = 3e6 # J (Total energy contained in metabolites)
global seed; seed = 600034

## FUNCTIONS ##

def names(s, m):
    '''Create a vector of column names for results dataframe'''
    #Strains
    s_seq = [str(i) for i in range(s)]
    names_s = np.repeat('s', s)
    label_s = np.char.add(names_s, s_seq)
    m_seq = [str(i) for i in range(m)]
    names_m = np.repeat('m', m)
    label_m = np.char.add(names_m, m_seq)
    return np.concatenate([label_s, label_m])

def letter_number(letter, number):
    ret = [letter+str(i) for i in range(number)]
    return(ret)

def main(argv):
    '''Main function'''
    
    #Set fix parameters#
    ####################

    #Number of metabolites
    m = 10
    #Number of strains
    s = 75
    #Generate random uniformly distributed chemical potentials
    np.random.seed(seed)
    #mu = random_pack(E, m)
    #mu = decreasing_pack(E, m)
    mu = uniform_pack(E, m)
    #Set growth proportionality constant
    g = np.ones(s) 
    #Set supply and dilution rate of metabolites
    kappa = 2*np.ones(m) 
    gamma = 1*np.ones(m)
    #Average cost per reaction
    xi0 = 0.005
    epsilon = 1e-3
    #Choose Eta matrix for each strain
    Eta = np.zeros(shape = (s,m,m))
    triu_ind = np.triu_indices(m,1)
    tril_ind = np.tril_indices(m,-1)
    eta_values_u = (triu_ind[1]-triu_ind[0])/(m-1)
    eta_values_l = (tril_ind[0]-tril_ind[1])/(m-1)
    for i in range(s):
        Eta[i][triu_ind] = eta_values_u
        Eta[i][tril_ind] = eta_values_l
    #Choose q_max, ks and krn
    q_max = [np.ones(1) for i in range(s)]
    ks = [0.1*q_max[i] for i in range(s)]
    kr = [10*q_max[i] for i in range(s)]
    #Total number of posible reactions (only forward net reactions are taking 
    #place)
    N_reac = comb(m, 2, exact = True)
    #Perform n_simul simulations
    n_simul = 1
    #Generate reaction network for each strain
    tot = s*n_simul
    n_reac_s = np.zeros(tot, dtype = 'int')
    maintenance = np.zeros(tot)
    #Prealocate vector for storing initial resource surpluses
    surplus = np.zeros(tot)
    #Prealocate vector for storing cohesion indices
    cohesion = np.zeros(tot)
    #Prealocate vector for storing facilitation indices
    facilitation = np.zeros(tot)
    #Prealocate vector for storing facilitation indices
    providing = np.zeros(tot)
    #The corrected version
    cohesion_corrected = np.zeros(tot)
    #Prealocate vector for storing competition indices
    competition = np.zeros(tot)
    #Prealocate vector for storing auto_cohesion indices
    auto_cohesion = np.zeros(tot)
    #Prealocate vector for storing new cohesion and competition indices
    coh = np.zeros(tot)
    comp = np.zeros(tot)
    #Prealocate container of all networks 
    tot_reac_network = tot*[0]
    pbar = ProgressBar()
    print('Creating reaction networks')
    for i in pbar(range(tot)):
        #Randomly set number of reactions in the network
        num_reac = 0
        np.random.seed(i+seed)
        #while num_reac == 0: 
        #    num_reac = round(np.random.beta(1,5)*N_reac)
        num_reac = np.random.randint(1, N_reac+1)#m+1)
        #Get reaction network
        reac_network = network(m, num_reac)
        #Store reaction network 
        tot_reac_network[i] = reac_network
        #Store actual number of reactions
        n_reac_s[i] = len(tot_reac_network[i][0])
        #Set cost based on the number of reactions
        rand_cost = np.random.normal(0,1)
        while rand_cost < -1/epsilon:
            rand_cost = np.random.normal(0,1)
        #Free energy gap of reactions
        mult = sum(tot_reac_network[i][1]-tot_reac_network[i][0])
        maintenance[i] = xi0*mult#*(1+epsilon*rand_cost)

    #Integrate model#
    #################

    #Time span
    tspan = tuple([1, 1e6])
    #Initial conditions 
    N0 = s*[1]
    C0 = m*[2] 
    z0 = N0 + C0


    #Prealocate storing objects
    empty_z = [[] for i in range(m+s)]
    empty_deg = [[] for i in range(m)]
    zn = np.array(empty_z)
    zn_stable = np.zeros(shape = (n_simul, m+s))
    Fn_stable = np.zeros(n_simul)
    degree_n = np.array(empty_deg)
    betweenness_n = np.array(empty_deg)
    F_n = np.array([])
    tn = np.array([])
    n_label = np.array([])

    pbar = ProgressBar()
    for  n in pbar(range(n_simul)):

        #Network
        network_n = tot_reac_network[s*n:s*(n+1)]
        #Maintenance vector
        maintenance_n = maintenance[s*n:s*(n+1)]
        #Integrate the whole community
        #Calculate possible measures of performance of the network prior to 
        #its realization
        #Calculate initial surplus of species in the network
        #Get s matrices of zeros: reaction presence/absence for each strain
        rates = np.zeros(shape = (s,m,m))
        #Calculate initial surplus for each strain given environment C0
        for j in range(s):
            network_t = network_n[j];
            eta = Eta[j]
            q_m = q_max[j]
            k_s = ks[j]
            k_r = kr[j]
            rates[j] = rate_matrix(network_t, eta, q_m, k_s, k_r, np.array(C0),
                                   mu, m)
            surplus[s*n + j] = jgrow(rates[j], eta) - maintenance_n[j]
        #Calculate level of cohesion of each species
        cohesion_matrix = facilitation_matrix(network_n, m)
        #Providing index
        providing_index = np.sum(cohesion_matrix, axis = 0)
        #Facilitation index
        facilitation_index = np.sum(cohesion_matrix, axis = 1)
        #Auto cohesion indices
        auto_cohesion_indices = np.diag(cohesion_matrix)
        cohesion_indices = cooperation_index(cohesion_matrix)
        coh_ind_corr = cohesion_indices - n_reac_s[s*n:s*(n+1)] + auto_cohesion_indices
        #Calculate the level of competition of each species
        comp_mat = competition_matrix(network_n, m)
        comp_indices = sum(comp_mat)
        #Calculate cohesion-competition
        coh[s*n:s*(n+1)] = np.sum(cohesion_matrix, axis = 1)
        comp[s*n:s*(n+1)] = np.sum(comp_mat, axis = 1)
        #Store cohesion, autocohesion and competition indices
        cohesion[s*n:s*(n+1)] = cohesion_indices
        facilitation[s*n:s*(n+1)] = facilitation_index
        providing[s*n:s*(n+1)] = providing_index
        cohesion_corrected[s*n:s*(n+1)] = coh_ind_corr
        competition[s*n:s*(n+1)] = comp_indices
        auto_cohesion[s*n:s*(n+1)] = auto_cohesion_indices
        #Calculate net level of interactions species
        interactions = cohesion_matrix - comp_mat
        #Add lower elements to upper elements
        triu = np.triu(interactions)
        tril = np.tril(interactions)
        total_interactions = triu + tril.transpose()
        sum_tot_inter = [sum(total_interactions[i,:]) + 
                         sum(total_interactions[:,i])
                         for i in range(s)]
         
        sol = solve_ivp(lambda t,z: model(t,z,s,m,kappa,gamma,network_n,mu,Eta,
                        q_max,ks,kr,g,maintenance_n), tspan, z0, 
                        method = 'BDF', atol = 1e-3)
        #Extract solution and time evaluations
        z = sol.y
        t = sol.t
        t_points = len(t)
        #if n == 0:
        #    image = cohesion_matrix - comp_mat
        #    plt.imshow(image)
        #    plt.colorbar()
        #    plt.show()
        #    import ipdb; ipdb.set_trace(context = 20)
        
        #Create community network and evaluate centrality measures for each
        #timepoint
        #Prealocate storing elements
        degree_time = np.zeros(shape=(m, t_points))
        betweenness_time = np.zeros(shape=(m, t_points))
        F = np.zeros(t_points)
        #Don't worry about divergent cases
        if any(z[:, -1]>1e4):
            pass
        else:
            F = F_calculator(z, t, m, s, network_n)
            for i in range(t_points):
                #Abundance of each strain at timepoint i
                abundances = z[0:s,i]
                #Community networkx for those abundances
                community_network = networks2community(network_n, abundances, 
                                                       s, m) 
                G = community_network 
                layout = nx.circular_layout(G)
                #Calculate size of edges (importance of metabolites)
                degree = list(dict(G.degree(weight = 'weight')).values())
                betweenness = nx.betweenness_centrality(G, weight = 'weight')
                betweenness = list(dict(betweenness).values())
                degree_time[:,i] = degree
                betweenness_time[:,i] = betweenness
                #nx.draw_networkx(G, layout, node_size = degree)
                #plt.pause(0.01)
                #plt.close()

        #plt.plot(t, np.transpose(degree_time))
        #plt.xscale('log')
        #plt.legend(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9'))
        #plt.show()

        #Store
        zn = np.concatenate([zn, z], axis = 1)
        zn_stable[n,:] = zn[:,-1]
        F_n = np.concatenate([F_n, F])
        Fn_stable[n] = F_n[-1]
        degree_n = np.concatenate([degree_n, degree_time], axis = 1)
        betweenness_n = np.concatenate([betweenness_n, betweenness_time],
                                       axis = 1)
        tn = np.concatenate([tn, t])
        n_label = np.concatenate([n_label, n*np.ones(len(t), dtype = int)])
        

    #Create dataframe of networks
    df_network = pd.DataFrame(tot_reac_network, 
                              columns = ['substrate', 'product'])
    
    df_network = pd.concat([df_network, 
                            pd.DataFrame({'maintenance':maintenance, 
                                          'surplus':surplus,
                                          'n_reac':n_reac_s, 
                                          'cohesion':cohesion,
                                          'providing':providing,
                                          'facilitation':facilitation,
                                          'cohesion_corr':cohesion_corrected,
                                          'competition':competition,
                                          'autocohesion':auto_cohesion, 
                                          'coh':coh,
                                          'comp':comp})], 
                            axis = 1)
    df_network.insert(0, 'strain', np.tile(np.arange(s), n_simul))
    df_network.insert(0, 'n_simulation', np.repeat(np.arange(n_simul), s))
    #Create dataframe of time series solutions
    cols = names(s, m)
    names_met = [i for i in range(len(cols)) if cols[i].startswith('m')]
    df_sol = pd.DataFrame(np.transpose(zn), columns = cols)
    df_sol.iloc[:,names_met] = np.transpose(degree_n)
    #Names betweenness
    nam = letter_number('b', m)
    df_bet = pd.DataFrame(np.transpose(betweenness_n), columns = nam)
    df_sol = pd.concat([df_sol, df_bet], axis = 1)
    df_sol = pd.concat([df_sol, pd.DataFrame(F_n, columns = ['F'])], axis = 1)
    df_sol.insert(0, 't', tn)
    df_sol.insert(0, 'n_simulation', n_label)
    #Create dataframe of stable community compositions
    df_stable = pd.DataFrame(zn_stable, columns = cols)
    df_stable = pd.concat([df_stable, pd.DataFrame(Fn_stable, 
                                                   columns = ['F'])], 
                          axis = 1)
    df_stable.insert(0, 'n_simulation', np.arange(n_simul))
    #Save
    df_sol.to_csv('../data/coal_time_series.csv')
    df_network.to_csv('../data/coal_network.csv')
    df_stable.to_csv('../data/coal_composition.csv')

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status) 
