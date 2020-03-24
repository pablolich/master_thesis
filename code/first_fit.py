#!/usr/bin/env python3

__appname__ = '[App_name_here]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from lmfit import Minimizer, Parameters, fit_report

## CONSTANTS ##

global R
global Tref
R = 8.62e-5 #Boltzmann constant
Tref = 273.15 #Reference temperature - 0C


## FUNCTIONS ##

class schoolfield():
    
    def set_parameters(self, params):
        
        #                         Name, Start, Vary, Lower, Upper
        params.add_many(('B0',   0.022, True, -np.inf, np.inf ),
                        ('E',    0.65 , True, 10E-3,   np.inf ), 
                        ('E_D',  4.4,   True, 10E-10,  np.inf ),
                        ('T_pk', 309,   True, 273.15,  np.inf))
        #Eventually the initial values will come from a function that will 
        #belong to the schoolfield class
        return params


    def model(self, params, T, mu_max):
        "Model implementation"

        #Unpack values of parameters.
        parameter_vals = params.valuesdict()
        #Basic metabolic rate
        B0 = parameter_vals['B0'] 
        #Activation energy of enzymes
        E = parameter_vals['E'] 
        #Inactivation energy of enzymes
        E_D = parameter_vals['E_D']         
        #Temperature at which peak response is observed
        T_pk = parameter_vals['T_pk']
        #Calculate theoretical predictions
        result = (B0 + \
                np.log(np.exp((-E / R) * ((1 / T) - (1 / Tref))) /\
                (1 + (E/(E_D - E)) * np.exp(E_D / R * (1 / T_pk - (1 / T))))) )

        return np.exp(result) - mu_max

def main(argv):
    '''Main function'''

    #Read data
    data = pd.read_csv('../data/minimal_example_data.csv')
    T = data.temp + 273
    mu_max = data.std_mu_max

    #Create model object from schoolfield class
    model = schoolfield()
    #Set parameters of the model
    params = Parameters()
    params = model.set_parameters(params)

    #Fit model
    minner = Minimizer(userfcn = model.model, params = params, 
                       fcn_args = (T, mu_max))
    fit_ = minner.minimize(method = 'leastsq')

    print(fit_report(fit_))
    #Plot results 
    plt.scatter(T, mu_max)
    plt.plot(T, model.model(fit_.params, T, mu_max) + mu_max, color = 'red')
    plt.savefig('../results/first_fit.pdf')

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

