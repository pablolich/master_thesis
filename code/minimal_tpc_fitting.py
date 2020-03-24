#!/usr/bin/env python3

__appname__ = '[minimal_tpc_fitting.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from lmfit import Minimizer, Parameters

## CONSTANTS ##

global R
global Tref
R = 8.62e-5 #Boltzmann constant
Tref = 273.15 #Reference temperature - 0C


## FUNCTIONS ##

class schoolfield:
    
    def __init__(self, parameters):
        '''Initialize Parameter values'''
    
        self.parameters = parameters
        #                         Name, Start, Vary, Lower, Upper
        self.parameters.add_many(('B0',   0.022, True, -np.inf, np.inf ),
                                 ('E',    0.65 , True, 10E-3,   np.inf ), 
                                 ('E_D',  4.4,   True, 10E-10,  np.inf ),
                                 ('T_pk', 309,   True, 273.15,  np.inf))
        #Parameter values should come from an algorithm that calculates initial 
        #estimates. 
        return

    def math_expr(self, T):
        "Model implementation"

        #Unpack values of parameters.
        parameter_vals = self.parameters.valuesdict()
        #Basic metabolic rate
        B0 = parameter_vals['B0'] 
        #Activation energy of enzymes
        E = parameter_vals['E'] 
        #Inactivation energy of enzymes
        E_D = parameter_vals['E_D']         
        #Temperature at which peak response is observed
        T_pk = parameter_vals['T_pk']
        #Calculate theoretical predictions
        result  = B0 + \
                np.log(np.exp((-E / R) * ((1 / T) - (1 / Tref))) /\
                (1 + (E/(E_D - E)) * np.exp(E_D / R * (1 / T_pk - (1 / T))))) 
        return result


def residuals(params, expr, T, mu_max) :
    '''Calculate residuals between model and data'''

    theo_mu_max = expr(T)
    residuals = np.exp(theo_mu_max) - mu_max

    return residuals


def main(argv):
    '''Main function'''


    #LOAD DATA#
    ###########

    data = pd.read_csv('../data/minimal_example_data.csv')
    T = data.temp + 273
    mu_max = data.std_mu_max


    #MINIMIZE SCHOOLFIELD MODEL#
    ############################

    #Create the object model from the class schoolfield.
    model = schoolfield(Parameters())

    #Fit model
    
    params = model.parameters
    minner = Minimizer(userfcn = residuals, params = params, 
                       fcn_args = (model.math_expr, T, mu_max))
    fit_ = minner.minimize(method = 'leastsq')
    import ipdb; ipdb.set_trace(context = 20)

    plt.scatter(T, mu_max)
    plt.plot(T, np.exp(model.math_expr(T)))
    plt.show()

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

