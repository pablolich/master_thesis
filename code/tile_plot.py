
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
    #Number of metabolites
    m = len([i for i in keys if i.startswith('m')])
    s = len([i for i in keys if i.startswith('s')])
    #Vector of richnesses
    richness = np.unique(all_data.richness)
    #Preallocate storing objects
    column_names = ["x", "y", "su", "richness"]
    tile_plot = pd.DataFrame(columns = column_names)
    r = 4
    #Get communities of richnes r
    comp_long = all_data[all_data.richness == r]
    #Get communities in team 1
    communities_team_1 = comp_long[comp_long.team == 1]
    #Communities in team -1
    communities_team__1 = comp_long[comp_long.team == -1]
    #Get all pairs of competing communities between the two groups
    all_pairs = list(product(np.unique(communities_team_1.n_simulation), 
                             np.unique(communities_team__1.n_simulation)))

    #Number of iterations
    it = len(all_pairs)
    #Initialize storing objects for the plot
    x = np.repeat(np.arange(it), 2*r)
    y = np.tile(np.arange(2*r), it)
    su = np.zeros(len(x))
    inter_rank = np.zeros(len(x))
    pbar = ProgressBar()
    print('Coalescence of communities of richness:', r)
    for i in pbar(range(it)):
        #Pick communities indices c1 and c2
        c1 = all_pairs[i][0]
        c2 = all_pairs[i][1]
        #Get community composition
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

        #Perform coalescence event between the two communities
        t, z, nets = coalescence_event(C1 = comp_c1, 
                                       C2 = comp_c2, 
                                           m = m, 
                                           s = s)
        #Get abundances
        N = z[0:s]
        #Get abundance vector of remaining species after coalescence at equil.
        abundance_f = N[:, -1] 
        ##Create dataframe of coalescence outcome
        outcome = pd.concat([comp_c1, comp_c2])
        outcome['stable.state'] = abundance_f
        #Get the provenance of extinct species into vector su
        temp = np.zeros(2*r)
        inds_surv = np.where(outcome['stable.state']<1)[0]
        temp[inds_surv] = np.array(outcome['team'])[inds_surv]
        #Sort according to cohesion
        sort = np.argsort(np.array(outcome['coh'] - outcome['comp']))
        su[2*r*i:2*r*(i+1)] = temp[sort]

    tile_p = pd.DataFrame({'x':x, 'y':y, 'su':su, 
                           'richness':np.repeat(r, len(x))})

    #Sort according to number of extinctionss from one group
    xvec = np.unique(tile_p.x)
    it = len(xvec)
    proportion_1 = np.zeros(it)
    for i in range(it):
        column = tile_p[tile_p.x == i]
        proportion_1[i] = sum(column.su == -1) - sum(column.su == 1)

    sort_x =np.argsort(proportion_1)[::-1]
    #sort_x is a vector that would sort the vector of proportions
    #What is the index corresponing to the values of sort_x?
    index = np.zeros(len(tile_p.x), dtype = int)
    for i in range(it):
        index[2*r*i:2*r*(i+1)] = np.arange(2*r*sort_x[i], 
                                           2*r*sort_x[i]+2*r)

    reorder = np.array(tile_p.su)[index]
    tile_p.su = reorder
    tile_plot = pd.concat([tile_plot, tile_p])
    tile_plot.to_csv('../data/tile_plot.csv')

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     
     



