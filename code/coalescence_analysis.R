#setwd("~/Desktop/master_thesis/code")
#Load libraries and functions#
library(stringr)
require(ggplot2)
library(reshape2)
library(viridis)
library(grid)
library(gridExtra)
source('coal_analysis_functions.R')

#Load data#
network = read.table('../data/coal_network.csv', sep = ',', header = T, 
                     colClasses = c('character'), row.names = 1)
composition = read.table('../data/coal_composition.csv', sep = ',', header = T,
                         row.names = 1)
time_series = read.table('../data/coal_time_series.csv', sep = ',', header = T)

#Code#
#Get m and s
m = sum(startsWith(names(time_series), 'm'))
s = sum(startsWith(names(time_series), 's'))
n_sim = length(composition$n_simulation == 0)

#Get abundance and richness of each surviving strain
ind_strains = which(startsWith(names(composition), 's') |
                    startsWith(names(composition), 'n'))
strain_comp = composition[,ind_strains]
#Melt to combine
melt_names = names(strain_comp)[2:length(names(strain_comp))]
melted_str_comp = melt(strain_comp, measure.vars = melt_names,
                       variable.name = 'strain',
                       value.name = 'stable.state')
#Get richness at equilibrium for each simulation
#richness = rep(rowSums(strain_comp[2:length(strain_comp)] > 1),
#              each = s)
#melted_str_comp['richness'] = richness

#Eliminate extinctions and outliers
surv_strain = subset(melted_str_comp, melted_str_comp$stable.state>1 & melted_str_comp$stable.state<1e6)
#Sort conveniently
surv_strain = surv_strain[order(surv_strain$n_simulation, 
                                surv_strain$strain),]

richness_tab = table(surv_strain$n_simulation) 
richness = rep(richness_tab, richness_tab)
surv_strain['richness'] = richness
#Get networks of surviving species

#Plot time series along with community network
#Change strain vector in network so that it maches that of 
#surv_strain
network$strain = paste('s', network$strain, sep = '')

#Merge
all_data = merge(surv_strain, network, by = c('n_simulation', 'strain'))
#Save it to data
write.csv(all_data, '../data/community_data.csv', row.names = F)

# #Quick check
# v = all_data$richness
# #To remove contiguous duplicated elements only
# v[c(TRUE, !v[-length(v)] == v[-1])]
# hist(v)


#Plot time series to pdf
#First put time_series of strains into long format
ind_strains_t = which(startsWith(names(time_series), 's')|
                      startsWith(names(time_series), 'n')|
                      startsWith(names(time_series), 't'))

strains_time_s = time_series[,ind_strains_t]
#Melt
str_time_melt = melt(strains_time_s, measure.vars = melt_names,
                     variable.name = 'strain',
                     value.name = 'population')

#Plot and save time series
pdf('../results/coal_time_series.pdf', width = 9.5, height = 4)
pb = txtProgressBar(1, n_sim-1, style=3)
for (i in unique(str_time_melt$n_simulation)){
  
  #Get dataframe to plot (all strains of that simulation)
  data = str_time_melt[which(str_time_melt$n_simulation == i),]
  time = plot_time_series(data)
  
  #Plot and save community reaction networks
  #Get indices of networks corresponding to simulation i
  inds = which(all_data$n_simulation == i)
  #Get networks for that simulation
  sub_networks = all_data[inds,c('substrate', 'product')]
  abundances = all_data$stable.state[inds]
  #Get community network
  com_network = community_network(sub_networks, abundances, m)
  network = plot_network(com_network)
  
  #Join them in one plot
  grid.arrange(time, network, nrow = 1)
  #Print progress...
  setTxtProgressBar(pb, i)
}
dev.off()

