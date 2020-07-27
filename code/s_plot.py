#!/usr/bin/env python3

__appname__ = '[s_plot.py]'
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
from itertools import combinations, product
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
    #Metabolite list
    met_list = [i for i in keys if i.startswith('m')]
    #Number of metabolites
    m = len(met_list)
    s = len([i for i in keys if i.startswith('s')])
    #Vector of richnesses
    richness = np.unique(all_data.richness)
    #Preallocate storing objects
    column_names = ["similarity", "cohesion"]
    similarity_fitness = pd.DataFrame(columns = column_names)
    tn = np.array([])
    Nn = np.array([])
    richnessn = np.array([], dtype = int)
    n_simulationn = np.array([], dtype = int)
    strainsn = np.array([], dtype = int)
    c_numbern = np.array([], dtype = int)

    #Get communities with richness 5
    r = 5
    comp_long = all_data[all_data.richness == r]
    #Perform coalescence events between all posible pairs of communities 
    #How many communities do we have?
    com_simul = np.unique(comp_long.n_simulation)
    num_com = len(com_simul)
    #Get all pairs of competing communities 
    all_pairs = list(combinations(com_simul, 2))
    #How many pairs?
    it = len(all_pairs)
    #Initialize storing objects
    similarity = np.zeros(it)
    DF = np.zeros(it)
    DP2 = np.zeros(it)
    pbar = ProgressBar()
    print('Coalescence of communities of richness:', r)
    for i in pbar(range(it)):
        #Pick communities numbers c1 and c2
        c1 = all_pairs[i][0]
        c2 = all_pairs[i][1]
        #Extract community information
        comp_c1 = comp_long[comp_long.n_simulation == c1].reset_index(drop = True)
        comp_c2 = comp_long[comp_long.n_simulation == c2].reset_index(drop = True)
        #Get number of strains in each community
        num_c1 = len(comp_c1)
        num_c2 = len(comp_c2)
        s = num_c1 + num_c2
        #Obtain reaction networks  of c1 as a list of tuples
        net_C1 = vector2tuple(comp_c1['substrate'],
                              comp_c1['product'])
        net_C2 = vector2tuple(comp_c2['substrate'],
                              comp_c2['product'])
        #Facilitation matrix of initial communities
        f_mat = facilitation_matrix(net_C1, m)
        f_mat2 = facilitation_matrix(net_C2, m)
        #Competition matrix 
        c_mat = competition_matrix(net_C1, m)
        c_mat2 = competition_matrix(net_C2, m)
        #Cohesion of each ccommunity
        P2_c1 = np.mean(np.sum(f_mat-c_mat, axis = 1))
        P2_c2 = np.mean(np.sum(f_mat2-c_mat2, axis = 1))
        #Calculate difference in cohesion
        DP2[i] = P2_c1 - P2_c2
        #Perform coalescence event between c1 and c2
        t, z, nets = coalescence_event(C1 = comp_c1, 
                                       C2 = comp_c2, 
                                       m = m, 
                                       s = s)

        #Get aboundance time series
        N = z[0:s]
        #Get abundance vector of species after coalescence at stable state
        abundance_f = N[:, -1] 
        #Create dataframe of coalescence outcome
        outcome = pd.concat([comp_c1, comp_c2])
        outcome['stable.state'] = abundance_f
        #Eliminate extinctions
        outcome = outcome[outcome['stable.state'] > 1].reset_index(drop = True)
        #Obtain reaction networks  of outcome as a list of tuples
        #net_outcome = vector2tuple(outcome['substrate'],
        #                           outcome['product'])
        #Number of species present in community c1 originally
        abundance_0 = np.array(comp_c1['stable.state'])
        #Add as many 0 as species in community c2 to calculate similarity
        abundance_0 = np.concatenate([abundance_0, np.zeros(num_c2)])
        #Calculate similarity 
        similarity[i] = np.dot(abundance_0, abundance_f)/\
                        (np.sqrt(sum(abundance_0**2))*\
                        np.sqrt(sum(abundance_f**2)))

    similarity_fitness = pd.DataFrame({'similarity':similarity,
                            'delP2':DP2})


    similarity_fitness.to_csv('../data/similarity_fitness.csv')
    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
 
 
