#setwd("~/Desktop/master_thesis/code")
rm(list=ls())
#Load libraries and functions#
library(stringr)
require(ggplot2)
library(reshape2)
library(viridis)
library(grid)
library(gridExtra)
library(colorspace)
library(plotwidgets)
source('coal_analysis_functions.R')

#Load data#
network = read.table('../data/coal_network.csv', sep = ',', header = T,
                     colClasses = c('n_simulation'='numeric','strain'='numeric','substrate'='character',
                                    'substrate'='character','maintenance'='numeric', 'surplus'='numeric',
                                    'n_reac'='numeric', 'cohesion'='numeric','competition'='numeric', 
                                    'autocohesion'='numeric'), row.names = 1)
#Change format of one column quickly
network$strain = paste('s', network$strain, sep = '')
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
                    startsWith(names(composition), 'n') |
                    startsWith(names(composition), 'F'))
strain_comp = composition[,ind_strains]
#Melt to combine
melt_names_str = names(strain_comp)[2:length(names(strain_comp))]
melted_str_comp = melt(strain_comp, id.vars = c('n_simulation', 'F'),
                       variable.name = 'strain',
                       value.name = 'stable.state'
                       )
#Get metabolite names


#Eliminate extinctions and outliers
out = melted_str_comp$n_simulation[which(melted_str_comp$stable.state > 1e4)]
surv_strain = melted_str_comp[which(!(melted_str_comp$n_simulation %in% out) & 
                                    melted_str_comp$stable.state > 1) , ]
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


#Merge
all_data = merge(surv_strain, network, by = c('n_simulation', 'strain'))
write.csv(all_data, '../data/community_data.csv', row.names = F)

#Quick check
v = all_data$richness
#To remove contiguous duplicated elements only
v = v[c(TRUE, !v[-length(v)] == v[-1])]
hist(v)

#Plot F for each simulation
ggplot(data = time_series, aes(x = t, y = F, colour = as.factor(n_simulation)))+ 
  geom_line()+
  theme(legend.position = 'none', legend.text = element_blank()) + 
  scale_x_continuous(trans = 'log10')


#Plot time series to pdf
#First put time_series of strains into long format
ind_strains_t = which(startsWith(names(time_series), 's')|
                      startsWith(names(time_series), 'n')|
                      startsWith(names(time_series), 't')|
                      startsWith(names(time_series), 'F'))
#Metabolites time series
ind_metabolite_t = which(startsWith(names(time_series), 'm')|
                         startsWith(names(time_series), 'n')|
                         startsWith(names(time_series), 't'))

#Metabolites, betweenness time series
ind_bet_t = which(startsWith(names(time_series), 'b')|
                           startsWith(names(time_series), 'n')|
                           startsWith(names(time_series), 't'))

strains_time_s = time_series[,ind_strains_t]
metabolite_time_s = time_series[,ind_metabolite_t]
bet_time_s = time_series[, ind_bet_t]

#Melt
str_time_melt = melt(strains_time_s,
                     variable.name = 'strain',
                     id.vars = c('n_simulation', 'F','t'),
                     value.name = 'population')
#get rid of extinctions? uncomment!
#str_time_melt=time_series_eliminator(str_time_melt, s)
met_time_melt = melt(metabolite_time_s, 
                     measure.vars = names(metabolite_time_s)[3:length((metabolite_time_s))],
                     variable.name = 'metabolite',
                     value.name = 'population')
bet_time_melt = melt(bet_time_s, 
                     measure.vars = names(bet_time_s)[3:length((bet_time_s))],
                     variable.name = 'betweenness',
                     value.name = 'population')

#Plot and save time series and community networks
pdf('../results/coal_time_series.pdf', width = 9.5, height = 4)
j = 1
pb = txtProgressBar(1, n_sim, style=3)
for (i in unique(sort(all_data$n_simulation))){

  #Get dataframe to plot (all strains of that simulation)
  data_str = str_time_melt[which(str_time_melt$n_simulation == i-1),]
  # data_met = met_time_melt[which(met_time_melt$n_simulation == i),]
  # data_bet = bet_time_melt[which(bet_time_melt$n_simulation == i),]
  #Initialite plot holder
  plot_holder = list()
  plot_c = 1
  time_str = plot_time_series(data_str, 
                              as.character(i),
                              data_str$strain)
  plot_holder[[plot_c]] = time_str
  plot_c = plot_c +1
  #Save it to data
  # time_met = plot_time_series(data_met, as.character(i),
  #                             data_met$metabolite)
  # time_bet = plot_time_series(data_bet, as.character(i),
  #                             data_bet$betweenness)
  
  #Plot and save community reaction networks
  #Get indices of networks corresponding to simulation i
  inds = which(all_data$n_simulation == i)
  #Get all data for that simulation
  all_data_sub = all_data[inds,]
  predictor = ggplot(data = all_data_sub, 
                     aes(x = (facilitation - providing)/competition, 
                         y = stable.state, 
                         colour = as.factor(n_simulation)))+
    geom_point(aes(shape = as.factor(richness), 
                   size = n_reac)) + 
    theme(legend.position = 'none')
  plot_holder[[plot_c]] = predictor
  plot_c = plot_c +1
  
  #Get networks for that simulation
  sub_networks = all_data[inds,c('substrate', 'product')]
  abundances = all_data$stable.state[inds]
  #Get community network
  com_network = community_network(sub_networks, abundances, m)
  network_p = plot_network(com_network)
  plot_holder[[plot_c]] = network_p
  plot_c = plot_c +1
  #Gather all networks from surviving species
  all_nets = all_networks(sub_networks, abundances, m)
  #Initialize plot holder
  #Plot them all
  for (u in seq(length(all_nets))+3){
    plot_holder[[u]] = plot_network(all_nets[[u-3]])
  }

  n = length(plot_holder)
  nCol = floor(sqrt(n))
  do.call('grid.arrange', c(plot_holder, nrow=nCol))
  #Join them in one plot
  #Print progress...
  j = j+1
  setTxtProgressBar(pb, j)
}
dev.off()


