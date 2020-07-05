#!/usr/bin/env python3

__appname__ = '[coalescence_event.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp, odeint
from coalescence import uniform_pack
from functions import * 
import matplotlib.pylab as plt
from itertools import combinations
from progressbar import ProgressBar
from scipy.special import comb

## CONSTANTS ##

global R; R = 8.314462618 # J/(K mol)
global DeltaGATP; deltaGATP = 75e3 # J/mol
global T; T = 298 # k 
global E; E = 3e6 # j (Total energy contained in metabolites)

## FUNCTIONS ##


def main(argv):
    '''Main function'''

    #Load communities data 
    composition = pd.read_csv('../data/coal_composition.csv', index_col = 0)
    networks = pd.read_csv('../data/coal_network.csv', index_col = 0)
    time_series = pd.read_csv('../data/coal_time_series.csv', index_col = 0)
    all_data = pd.read_csv('../data/community_data.csv')

    #Transform from string to numeric vector all elements in substrate and 
    #product colums
    all_data['substrate'] = all_data['substrate'].apply(lambda x: 
                            np.fromstring(x[1:-1], sep=' ', dtype = int))
    all_data['product'] = all_data['product'].apply(lambda x: 
                          np.fromstring(x[1:-1], sep=' ', dtype = int))
    keys = list(composition.keys())
    #Number of metabolites
    m = len([i for i in keys if i.startswith('m')])
    s = len([i for i in keys if i.startswith('s')])
    #Vector of richnesses
    richness = np.unique(all_data.richness)[4:5]
    #Preallocate storing objects
    column_names = ["similarity", "delF", "richness", "delP", "delP2"]
    similarity_fitness = pd.DataFrame(columns = column_names)
    column_names = ['n_simulation', 'comm_number', 'strain', 'richness', 't',
                     'N']
    time_series = pd.DataFrame(columns = column_names)
    #Perform coalescence events between all posible pairs of communities with 
    #the same richness
    tn = np.array([])
    Nn = np.array([])
    richnessn = np.array([], dtype = int)
    n_simulationn = np.array([], dtype = int)
    strainsn = np.array([], dtype = int)
    c_numbern = np.array([], dtype = int)

    for r in richness:
        comp_long = all_data[all_data.richness == r]
        #How many communities do we have?
        com_simul = np.unique(comp_long.n_simulation)
        num_com = len(com_simul)
        all_pairs = list(combinations(com_simul, 2))
        it = len(all_pairs)
        #Initialize storing objects
        similarity = np.zeros(it)
        DF = np.zeros(it)
        DP = np.zeros(it)
        DP2 = np.zeros(it)
        rich = r*np.ones(it)
        pbar = ProgressBar()
        print('Coalescence of communities of richness:', r)
        for i in pbar(range(it)):
            #Pick communities c1 and c2
            c1 = all_pairs[i][0]
            c2 = all_pairs[i][1]
            comp_c1 = comp_long[comp_long.n_simulation == c1].reset_index(drop = True)
            comp_c2 = comp_long[comp_long.n_simulation == c2].reset_index(drop = True)
            #Obtain number of metabolites
            #Get number of strains in each community
            num_c1 = len(comp_c1)
            num_c2 = len(comp_c2)
            s = num_c1 + num_c2
            #Calculate fitness of initial community
            F_C1 = np.unique(comp_c1.F) 
            #Obtain reaction networks  of c1 as a list of tuples
            net_C1 = vector2tuple(comp_c1['substrate'],
                                  comp_c1['product'])
            #Facilitation matrix of initial community
            f_mat = facilitation_matrix(net_C1, m)
            #Competition matrix 
            c_mat = competition_matrix(net_C1, m)
            #Provided and facilitated indices
            providing_index = np.sum(f_mat, axis = 0)
            facilitation_index = np.sum(f_mat, axis = 1)
            #Competition index for each species
            competition_index = sum(c_mat)
            #Calculate fitness of each individual
            fitness = individual_fitness(_in = np.sum(f_mat, axis = 0),
                                         out_ = np.sum(f_mat, axis = 1), 
                                         competition = np.sum(c_mat, axis = 0))

            #Calclulate level of cooperation/cohesion
            #Predictor2
            P_c1 = np.mean((facilitation_index-providing_index))
            #Predictor3
            P2_c1 = sum(sum(f_mat-c_mat))
            #Perform coalescence event between the two communities
            t, z, nets = coalescence_event(C1 = comp_c1, 
                                           C2 = comp_c2, 
                                           m = m, 
                                           s = s)

            N = z[0:s]
            C = z[s:s+m]
            t_points = len(t)
            r_coal = len(N)
            tn = np.concatenate([tn, np.tile(t, r_coal)])
            Nn = np.concatenate([Nn, N.reshape(1, r_coal*t_points)[0]])
            richness = np.repeat(r, 2*r*t_points)
            richnessn = np.concatenate([richnessn, richness])
            c_number = np.concatenate([np.repeat(c1, r*t_points),
                                       np.repeat(c2, r*t_points)])
            c_numbern= np.concatenate([c_numbern, c_number])
            n_simulation = np.repeat(i, 2*r*t_points)
            n_simulationn = np.concatenate([n_simulationn, n_simulation])
            strains = np.repeat(np.arange(2*r), t_points)
            strainsn = np.concatenate([strainsn, strains])

            #Predictors of community similarity
            #Calculate fitness of new community
            F = F_calculator(z, t, m, s, nets)
            ##Check with plot
            #colors = 6*['red'] + 6*['blue']

            #for j in range(s):
            #    plt.plot(t, N[j], label= 'Strain'+str(F[j]), color = colors[j])
            #    
            ##plt.legend()
            #plt.xscale('log')
            #plt.show()
            #Fitness difference of both communities
            DF[i] = F[-1] - F_C1 
            #Get abundance vector of remaining species after coalescence.
            abundance_f = N[:, -1] 
            #Get indices of surviving species after coalescence
            survival = abundance_f > 1
            #Eliminate extinctions from abundance vector
            stable = abundance_f[survival]
            #Create dataframe of coalescence outcome
            outcome = pd.concat([comp_c1, comp_c2])
            outcome['stable.state'] = abundance_f
            outcome = outcome[outcome['stable.state'] > 1].reset_index(drop = True)
            #Calculate predictor 2 after coalescence
            #Obtain reaction networks  of c1 as a list of tuples
            net_outcome = vector2tuple(outcome['substrate'],
                                       outcome['product'])
            #Facilitation index of outcome community
            f_mat = facilitation_matrix(net_outcome, m)
            c_mat = competition_matrix(net_outcome, m)
            #Provided and facilitated indices
            providing_index = np.sum(f_mat, axis = 0)
            facilitation_index = np.sum(f_mat, axis = 1)
            #Competition index for each species
            competition_index = sum(c_mat)
            #Calclulate level of cooperation/cohesion
            #Predictor2
            P_outcome = np.mean((facilitation_index-providing_index))
            DP[i] = P_outcome - P_c1
            #Predictor3
            P2_outcome = sum(sum(f_mat-c_mat))
            DP2[i] = P2_outcome - P2_c1
            #Number of species present in community c1 originally
            abundance_0 = np.array(comp_c1['stable.state'])
            #Add as many 0 as species in community c2 to calculate similarity
            abundance_0 = np.concatenate([abundance_0, np.zeros(num_c2)])
            #Calculate similarity 
            similarity[i] = np.dot(abundance_0, abundance_f)/\
                            (np.sqrt(sum(abundance_0**2))*\
                            np.sqrt(sum(abundance_f**2)))

            print(comp_c1)
            print(comp_c2)
            print(outcome)
            print(similarity[i])
            print(P2_c1)
            print(P2_outcome)
            import ipdb; ipdb.set_trace(context = 20)
            
            #if (similarity[i] > 0.5) and (similarity[i]< 0.75):
            #    import ipdb; ipdb.set_trace(context = 20)
        
        time_ser = pd.DataFrame({'n_simulation':n_simulationn,
                                 'comm_number':c_numbern,
                                 'strain':strainsn,
                                 'richness':richnessn,
                                 't':tn,
                                 'N':Nn})
        time_series = pd.concat([time_series, time_ser])

        sim_fit = pd.DataFrame({'similarity':similarity,
                                'delF':DF, 
                                'richness':rich,
                                'delP':DP,
                                'delP2':DP2})
        similarity_fitness = pd.concat([similarity_fitness, sim_fit])



    similarity_fitness.to_csv('../data/similarity_fitness.csv')


    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
     
