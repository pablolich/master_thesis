from functions import interaction_evolution, evol_richness
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from progressbar import ProgressBar

time_series = pd.read_csv('../data/coal_time_series.csv',
                          index_col = 0)
networks = pd.read_csv('../data/coal_network.csv', 
                       index_col = 0)

#Get surplus of each species
surplus = networks.surplus
#Filter time_series
keys = list(time_series.keys())
#List of metabolites 
junk = [i for i in keys if not (i.startswith('s') or 
                                i.startswith('n') or 
                                i.startswith('t'))]
#Get rid of metabolites
time_series = time_series.drop(columns = junk)
#Get number of strains
s = len(time_series.keys())-2
#Get number of simulations
n_simul = len(np.unique(time_series.n_simulation))-1990
#Filter networks
keys = list(networks.keys())
junk = [i for i in keys if not (i.startswith('sub') or i.startswith('prod'))]
#Get rid of junk
networks = networks.drop(columns = junk)
#Initialize storing object
names_strains = ['s' + str(i) for i in range(s)]
column_names = ["n_simulation", "t"] + names_strains
inter_evol_df = pd.DataFrame(columns = column_names)
c_evol_df = pd.DataFrame(columns = column_names)
f_evol_df = pd.DataFrame(columns = column_names)
surv_vector = np.array([])
pbar = ProgressBar()
print('Calculating evolution of interactions during community assembly')
for i in (range(n_simul)):
    time_series_i = time_series[time_series.n_simulation == i].reset_index(drop = True)
    networks_i = networks.iloc[i*s:(i+1)*s, :].reset_index(drop = True)
    #Get evolution of richness during community assembly
    richness = evol_richness(abundances_t = time_series_i)
    #Get positions where extinctions happen
    indx = np.where(richness[:-1] != richness[1:])[0]
    #Get values of time where extinctions happen
    t_ext = time_series_i.t[indx]
    t_points = len(t_ext)
    t_ext_ind = np.array(t_ext.index)
    f_evol, c_evol, inter_evol = interaction_evolution(time_series_i, 
                                                       networks_i, m = 15, 
                                                       t_points = t_ext_ind)
    #Get vector of surviving and extinct species
    stable_state = time_series_i.iloc[-1, 2:]
    survivors = [int(i) for i in stable_state>1]
    surv_vector = np.concatenate((surv_vector, np.repeat(survivors, t_points)))
    #Create vector of simulation number to insert in dataframe
    n_sim_i = np.repeat(i, len(t_ext))
    #Create a pandas dataframe for interactions, competitions and facilitations
    inter_evol_i = pd.DataFrame(np.concatenate((n_sim_i[:, np.newaxis], 
                                       t_ext[:, np.newaxis], inter_evol), 
                                      axis = 1),
                       columns = column_names)
    c_evol_i = pd.DataFrame(np.concatenate((n_sim_i[:, np.newaxis], 
                                       t_ext[:, np.newaxis], c_evol), 
                                      axis = 1),
                       columns = column_names)
    f_evol_i = pd.DataFrame(np.concatenate((n_sim_i[:, np.newaxis], 
                                       t_ext[:, np.newaxis], f_evol), 
                                      axis = 1),
                       columns = column_names)
    #Add to our running total pandas datframes
    inter_evol_df = pd.concat([inter_evol_df, inter_evol_i])
    c_evol_df = pd.concat([c_evol_df, c_evol_i])
    f_evol_df = pd.concat([f_evol_df, f_evol_i])
    
#Get index of surviving species at end of simulation
inter_evol_df = pd.melt(inter_evol_df, 
                        id_vars = ['t', 'n_simulation'], 
                        var_name = ['strain'], 
                        value_name = 'interaction')
f_evol_df = pd.melt(f_evol_df, id_vars = ['t', 'n_simulation'], 
                    var_name = 'strain', 
                    value_name = 'facilitation')
c_evol_df = pd.melt(c_evol_df, id_vars = ['t', 'n_simulation'],
                    var_name = 'strain', 
                    value_name = 'competition')
f_evol_df = f_evol_df.drop(columns = ['strain', 'n_simulation', 't'])
c_evol_df = c_evol_df.drop(columns = ['strain', 'n_simulation', 't'])

inter_evol_df = pd.concat([inter_evol_df, f_evol_df, c_evol_df], axis = 1)
#Add column of survivors
inter_evol_df['survivor'] = surv_vector
import ipdb; ipdb.set_trace(context = 20)

#Get boolean where the value id (strain id) matches survivors
idx = inter_evol_df.strain.isin(survivors)
#Once I get survivors, assign a vector for colors
colors = np.zeros(len(inter_evol_df))
abundances = np.zeros(len(inter_evol_df))
abundance_survivors = np.repeat(time_series[survivors].iloc[-1], t_points)
colors[idx] = 1
abundances[idx] = abundance_survivors
#Add colors to dataframe
inter_evol_df['survivor'] = colors
inter_evol_df['abundance'] = abundances
inter_evol_df.insert(0, 't', np.tile(t_ext, s))
#Save dataframe
inter_evol_df.to_csv('../data/interaction_evolution.csv')
