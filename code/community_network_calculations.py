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
    n_simul = max(communities.n_simulation)
    #Initialize list were all networks will be stored
    networks = [0]*n_simul
    for i in range(n_simul):
        #Subset simulation i
        net = communities[communities['n_simulation'] == i].loc[:,['substrate', 
                                                                   'product']]
        #Reset index to start from 0
        net = net.reset_index(drop = True)
        #Initialize adjacency matrices for each strain network
        net_strains = np.zeros(shape = (len(net), m, m))
        #Transform dataframe to list of tuples in order to index matrices
        tuples = [tuple([net.loc[k, 'substrate'], 
                         net.loc[k, 'product']]) for k in 
                         range(len(net))]
        #Assign weighted adjacency matrix to each strain 
        for j in range(len(net)):
            net_strains[j][tuples[j]] = 1

        #Add up all networks to form the community adjacency matrix
        community_adj = np.sum(net_strains, axis = 0)
        networks[i] = nx.from_numpy_matrix(community_adj, 
                                           create_using = nx.DiGraph)

    ####################################################################


    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
     
