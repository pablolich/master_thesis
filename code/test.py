from functions import interaction_evolution, evol_richness
import matplotlib.pylab as plt
import numpy as np
import pandas as pd

time_series = pd.read_csv('../data/coal_time_series_interaction_evolution.csv',
                          index_col = 0)
networks = pd.read_csv('../data/coal_network_interaction_evolution.csv', 
                       index_col = 0)

#Get number of strains
s = len(networks)
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
#Filter networks
keys = list(networks.keys())
junk = [i for i in keys if not (i.startswith('sub') or i.startswith('prod'))]
#Get rid of junk
networks = networks.drop(columns = junk)
richness = evol_richness(abundances_t = time_series)
#Get positions where extinctions happen
indx = np.where(richness[:-1] != richness[1:])[0]
#Get values of time where extinctions happen
t_ext = time_series.t[indx]
t_points = len(t_ext)
t_ext_ind = np.array(t_ext.index)
f_evol, c_evol, inter_evol = interaction_evolution(time_series, 
                                                   networks, m = 10, 
                                                   t_points = t_ext_ind)

#Get index of surviving species at end of simulation
col_names = ['s'+str(i) for i in range(s)]
df_inter = pd.DataFrame(inter_evol)
df_f = pd.DataFrame(f_evol)
df_c = pd.DataFrame(c_evol)
df_inter.columns = col_names
survivors = list(np.array(df_inter.iloc[-1][df_inter.iloc[-1]!=0].index))
inter_evol_df = pd.melt(df_inter, var_name = 'strain', 
                              value_name = 'interaction')
f_evol_df = pd.melt(df_f, var_name = 'strain', value_name = 'facilitation')
c_evol_df = pd.melt(df_c, var_name = 'strain', value_name = 'competition')
f_evol_df = f_evol_df.drop(columns = ['strain'])
c_evol_df = c_evol_df.drop(columns = ['strain'])

inter_evol_df = pd.concat([inter_evol_df, f_evol_df, c_evol_df], axis = 1)

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
