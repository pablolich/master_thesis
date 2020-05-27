#!/usr/bin/env python3

__appname__ = '[analysis_integration]'
__author__ = 'Pablo Lechon (pl1619@ic.ac.uk)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import matplotlib.pylab as plt
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages

## CONSTANTS ##


## FUNCTIONS ##


def main(argv):
    '''Main function'''
    #Types of data saved
    _type = ['classic', 'new']
    ind = 0

    #Load the data
    Rn = np.load('../data/reactions_' + _type[ind] + '.npy')
    Nn = np.load('../data/population_' + _type[ind] + '.npy')
    Cn = np.load('../data/concentration_' + _type[ind] + '.npy')
    
    #Number of simulations
    k_max = len(Rn)
    #Number of metabolites
    m = len(Cn[0])
    #Number of strains
    s = len(Nn[0])

    #Get last column of Nn for all strains and realizations 
    #Initialize quantities
    final_pop = np.zeros(shape = (k_max, s))
    survive_cum = np.zeros(s)
    richness = np.zeros(k_max)
    for i in range(k_max):
        stable = Nn[i,:,-1]
        survive = (stable > 1e-6) & (stable < 1)
        richness[i] = sum(survive)
        survive_cum = survive_cum + survive
        final_pop[i,:] = stable

    #Compare reaction networks of both sides
    Rn_c = Rn
    Rn_n = np.load('../data/reactions_' + _type[1] + '.npy')
    with PdfPages('../results/comparing_methods.pdf') as pdf:
        for i in range(20):
            R_c = Rn_c[i]
            R_n = Rn_n[i]
            R_tot = np.concatenate([R_c, R_n])
            fig = plt.figure()
            for j in np.arange(10)+1:
                plt.subplot(2,5,j)
                plt.imshow(R_tot[j-1])
            pdf.savefig(fig)
            plt.close()
        
        
    

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

