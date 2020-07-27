#setwd("~/Desktop/master_thesis/code")
rm(list=ls())
#Load libraries and functions#
library(stringr)
require(ggplot2)
library(ggpubr)
library(igraph)
library(reshape2)
library(viridis)
library(grid)
library(gridExtra)
library(colorspace)
library(plotwidgets)
library(RColorBrewer)
library(wesanderson)


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
out = melted_str_comp$n_simulation[which(melted_str_comp$stable.state > 1e3)]
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

#Quick check
v = all_data$richness
#To remove contiguous duplicated elements only
v = v[c(TRUE, !v[-length(v)] == v[-1])]
ay = hist(v)
histogram_richness = ggplot(data.frame(v), aes(x=v)) + 
  geom_histogram(binwidth=1,
                 color = '#6BB1C1', 
                 fill = '#3B9AB2')+
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size = 1),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20),
        aspect.ratio = 1)+
  scale_x_continuous(expand = c(0, 0), 
                   breaks = c(3, 5, 7, 9), 
                   labels = c('3', '5', '7', '9'))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = 'Richness', y = 'Frequency')+
  ggsave(filename = '../results/histogram_richness.pdf', width = 4.3, height = 4)

#Another quick check
simulations = sample(seq(0,499), size = 10)
#get indicess that match sampled simulations
sampled_data = all_data[all_data$n_simulation %in% simulations,]
ggplot(sampled_data, aes(x = competition, y = stable.state,
                         size = n_reac, 
                         colour = as.factor(n_simulation)))+
  geom_point()

#Get rid of outliers
all_data = all_data[!(all_data$n_simulation == 270),]
levels_reac = length(unique(all_data$n_reac))
#Colors
paleta_reac = wes_palette("Zissou1", levels_reac, type = 'continuous')
#Plot median stable state for each reaction number
num_reac_fitness = aggregate(stable.state~n_reac, all_data, mean)
std_ =  sqrt(aggregate(stable.state~n_reac, all_data, var))/sqrt(table(all_data$n_reac))
num_reac_fitness['err'] = std_$stable.state
n_reac_stablestate = ggplot(num_reac_fitness, aes(x = n_reac, y = stable.state,
                             colour = n_reac))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin=stable.state-err, ymax=stable.state + err,
                    colour = n_reac), 
                width=.2,
                position=position_dodge(.9))+
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size = 1),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20),
        legend.position = c(0.85, 0.7),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.background = element_blank(),
        aspect.ratio = 1)+
  scale_x_continuous(breaks = c(1, 5, 10, 15), 
                     labels = c('1', '5', '10', '15'))+
  labs(colour = expression('n'[r]), x = expression('n'[r]), 
       y = expression(symbol("\341")*N[infinity]*symbol("\361")))+
  scale_color_gradientn(colors = paleta_reac,
                        breaks = c(1,5,10,15))+
  ggsave(filename = '../results/n_reac_stablestate.pdf', width = 4.3, height = 4)

#Plot F for each simulation
ggplot(data = time_series, aes(x = t, y = F, colour = as.factor(n_simulation)))+ 
  geom_line()+
  theme(legend.position = 'none', legend.text = element_blank()) + 
  scale_x_continuous(trans = 'log10')

#Plot interaction evolution
inter_evol = read.table('../data/interaction_evolution.csv', sep = ',', header = T,
                         row.names = 1)

#Get average curve for dying people, and average curve for surviving people, and compare.
extinct_before = inter_evol[inter_evol$survivor == 0,]
#extinct_before = extinct[extinct$interaction != 0,]
mean_extinct = aggregate(extinct_before[, 3], list(extinct_before$t), mean)
extant = inter_evol[inter_evol$survivor != 0,]
mean_extant = aggregate(extant[, 3], list(extant$t), mean)
#Get rid of 0 so tha the mean doesn't go down artificially
inter_evol_no0 = inter_evol[!(inter_evol$survivor == 0 & inter_evol$interaction == 0),]
mean_tot = aggregate(inter_evol_no0[,3], list(inter_evol_no0$t), mean)
plot(mean_tot$Group.1, mean_tot$x, log = 'x')
ggplot()+
  geom_point(data = inter_evol, aes(colour  =  log10(abundance+1), 
                                   alpha = as.factor(survivor), 
                                   x = t, 
                                   y = interaction, 
                                   group = strain))+
  theme_classic()+
  theme(legend.position = 'none')+

  scale_alpha_manual(values = c(0.1,0.5,1,1))+
  scale_x_continuous(trans = "log10")+
  scale_color_gradient(low = "blue", high = "red", na.value = NA)+
  geom_line(data = means[means$variable=='x.x',], aes(x = Group.1, y = value),
            color = 'blue', linetype = 'dashed', size = 3)+
  geom_line(data = means[means$variable=='x.y',], aes(x = Group.1, y = value),
            color = 'red', linetype = 'dashed', size = 3)

#Plot community facilitation vs community competition
facilitation = aggregate(all_data[,17], list(all_data$n_simulation), mean)
competition = aggregate(all_data[,18], list(all_data$n_simulation), mean)
n_reac_community = aggregate(all_data[,10], list(all_data$n_simulation), mean)
plot(competition$x, facilitation$x)
#Make dataframe with data on community-level information
colors = distance(competition$x, facilitation$x, 1.4, 0.2)
data_community = data.frame(competition$x, facilitation$x, colors)

#Add color data to all_data and save it
unique_richness = non_contiguous(all_data)
compete_vector = rep(colors, unique_richness)
#Sort all_data according to simulation before merging
all_data = all_data[order(all_data$n_simulation),]
#Merge to all_data
all_data['team'] = compete_vector 

#Save data
write.csv(all_data, '../data/community_data.csv', row.names = F)

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
  data_str = str_time_melt[which(str_time_melt$n_simulation == i),]
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
                     aes(x = (coh+surplus - comp), 
                         y = stable.state, 
                         colour = as.factor(n_simulation)))+
    ggtitle(as.character(mean(all_data_sub$competition)))+
    geom_point(aes(shape = as.factor(richness), 
                   size = n_reac)) + 
    theme(legend.position = 'none')
  plot_holder[[plot_c]] = predictor
  plot_c = plot_c +1
  
  #Get networks for that simulation
  sub_networks = all_data[inds,c('substrate', 'product', 'stable.state')]
  #Indexes of top 3 performing species
  sub_networks = sub_networks[order(-sub_networks$stable.state),]
  top_networks = sub_networks[1:3,]
  abundances = all_data$stable.state[inds]
  #Get community network
  com_network = community_network(sub_networks, abundances, m)
  network_p = plot_network(com_network)
  plot_holder[[plot_c]] = network_p
  plot_c = plot_c +1
  #Gather networks from top 3 surviving species
  all_nets = all_networks(top_networks, abundances[0:3], m)
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

#Plot time series for a nice one.
data_nice = str_time_melt[which(str_time_melt$n_simulation == 128),]
#Get number of time points evaluations
t_points = length(unique(data_nice$t))
#Add number of reactions to each strain
data_nice['n_reac'] = rep(network$n_reac[network$n_sim == 128], each = t_points)
#Colors
paleta = wes_palette("Zissou1", 6, type = 'continuous')
#Plot
model_example = ggplot(data_nice, aes(x = t, y = population, 
                      group = strain, 
                      colour = as.factor(n_reac)))+
  geom_line(size = 1.2) +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15, color = 'black'),
    legend.position = c(0.1, 0.75),
    legend.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size = 1),
    axis.title = element_text(size = 20), 
    legend.key = element_blank(), 
    axis.text = element_text(size = 15),
    #plot.margin=unit(c(2,0,-2,0),"cm")
    )+
  scale_x_continuous(limits = c(1, 5e4), trans = 'log10', expand = c(0, 0), 
                     labels = c('1', '1e2', '1e3')) + 
  scale_y_continuous(limits = c(-1, 35), expand = c(0, 0))+
  labs(x = 't [a.u.]',
       y = expression(paste('N'[alpha],' [a.u.]')))+
  guides(color=guide_legend(expression('n'[r]), 
                            title.hjust = 0.5))+
  scale_color_manual('', values=paleta)+
  ggsave(filename = '../results/model_example.pdf', width = 4.79, height = 3.53 )

#Plot community reaction network for this example.
#Get networks for simulation 128
ind_nice = which(all_data$n_simulation == 128)
sub_networks = all_data[ind_nice, c('substrate', 'product', 'stable.state')]
#Indexes of top 3 performing species
sub_networks = sub_networks[order(-sub_networks$stable.state),]
abundances = all_data$stable.state[ind_nice]
#Get community network
com_network = community_network(sub_networks, abundances, m)
#Plot a network, given the matrix of reactions
paint = melt(com_network)
#Set colors
levels = length(unique(paint$value))
pal <- wes_palette("Zissou1", levels, type = "continuous")
#Plot
community_reaction_network = ggplot(paint, aes(Var2, rev(Var1))) +
  geom_tile(aes(fill = value), color = '#3B9AB2')+
  theme(panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.ticks = element_line(),
        legend.position = c(0.2, 0.3),
        legend.background = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))+
  scale_x_continuous(expand = c(0, 0), 
                     breaks = c(1, 5, 10, 15), 
                     labels = c('1', '5', '10', '15'))+
  scale_y_continuous(expand = c(0, 0),
                     breaks = rev(c(1, 5, 10, 15)), 
                     labels = c('1', '5', '10', '15'))+
  scale_fill_gradientn(colours = paleta,
                       breaks = c(0, 40, 80))+
  labs(fill = expression(Omega[xy]), x = 'y (product)', y = 'x (substrate)')+
  coord_fixed(ratio = 1)+
  ggsave(filename = '../results/community_reaction_network.pdf', width = 4.3, height = 4)

#Arange all plots in a mega plot
tot = ggarrange(model_example, community_reaction_network,
          histogram_richness, n_reac_stablestate,
          labels = c("A", "B", "C", "D"),
          ncol =2, nrow = 2, 
          font.label = list(size = 25))

ggsave('../results/community_assembly_plots.pdf', plot = tot, height = 10, width = 10)


#Plot all communities in a competition facilitation scatter plot
ggplot(data = data_community, aes(x = competition$x, y = facilitation$x))+
  geom_point(aes(colour = as.factor(colors)))+
  scale_color_manual(values=brewer.pal(n = 3, name = 'RdYlBu'))

