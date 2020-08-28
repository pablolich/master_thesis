#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import pandas as pd
import numpy as np
import networkx as nx
import re
import matplotlib.pylab as plt

## CONSTANTS ##

global m; m = 10

## FUNCTIONS ##

def main(argv):
    '''Main function'''

    #Load pool of communities
    communities = pd.read_csv('../data/community_data.csv')
    #Transform from string to numeric vector all elements in substrate and 
    #product colums
    communities['substrate'] = communities['substrate'].apply(lambda x: 
                               np.fromstring(x[1:-1], sep=' ', dtype = int))
    communities['product'] = communities['product'].apply(lambda x: 
                             np.fromstring(x[1:-1], sep=' ', dtype = int))

    #Create a community network for each simulation
    #Get number of simulations
    i_vec = np.unique(communities.n_simulation)
    n_simul = len(i_vec)
    #Initialize list were all networks will be stored
    networks = [0]*n_simul
    all_in_degree = np.zeros(shape=(n_simul, m))
    all_out_degree = np.zeros(shape=(n_simul, m))
    all_betweenness = np.zeros(shape=(n_simul, m))
    richness_com = np.zeros(n_simul)
    q = 0
    for i in i_vec :
        #Subset simulation i
        net = communities[communities['n_simulation'] == i]
        richness_com[q] = np.unique(net.loc[:,['richness']])
        #Reset index to start from 0
        net = net.reset_index(drop = True)
        #Initialize adjacency matrices for each strain network
        net_strains = np.zeros(shape = (len(net), m, m))
        #Transform dataframe to list of tuples in order to index matrices
        tuples = [tuple([net.loc[k, 'substrate'], 
                         net.loc[k, 'product']]) for k in 
                         range(len(net))]
        #Assign weighted adjacency matrix to each strain 
        abundances = net['stable.state']
        for j in range(len(net)):
            net_strains[j][tuples[j]] = abundances[j]

        #Add up all networks to form the community adjacency matrix
        community_adj = np.sum(net_strains, axis = 0)
        networks[q] = nx.from_numpy_matrix(community_adj, 
                                           create_using = nx.DiGraph)
        #Calculate network centrality measures
        in_degree = networks[q].in_degree(weight = 'weight')
        out_degree = networks[q].out_degree(weight = 'weight')
        betweenness = nx.betweenness_centrality(networks[q],
                                                weight = 'weight')
        #Store in matrices
        all_in_degree[q,:] = list(dict(in_degree).values())
        all_out_degree[q,:] = list(dict(out_degree).values())
        all_betweenness[q,:] = list(dict(betweenness).values())
        q += 1

        ##Check
        #G = networks[0]
        #layout = nx.spring_layout(G)
        #nx.draw_networkx(G, layout)

    ####################################################################
    #Take std of rows
    sigma_in = np.std(all_in_degree, axis = 1)
    sigma_out = np.std(all_out_degree, axis = 1)
    sigma_bet = np.std(all_betweenness, axis = 1)
    #Removing elements that have consecutive duplicates from richness vector
    df = pd.DataFrame({'sigma_in':sigma_in, 'sigma_out':sigma_out,
                       'sigma_bet':sigma_bet, 'richness':richness_com})
    df.to_csv('../data/scatterplot.csv')
    import ipdb; ipdb.set_trace(context = 20)

    plt.scatter(sigma_in, sigma_out)
    plt.show()



    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
     
