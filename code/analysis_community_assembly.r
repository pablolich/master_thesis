rm(list=ls())
dev.off()
require(ggplot2)
library(reshape2)
library(stringr)
library(wesanderson)
require(cowplot)
library(ggpubr)
library(grid)
library(gridExtra)
library(igraph)
source('coal_analysis_functions.R')

#Compile all useful community data into one object
##############################################################################
#Load time series data and stable state community composition
composition = read.table('../data/coal_composition.csv', sep = ',', header = T,
                         row.names = 1)
time_series = read.table('../data/coal_time_series.csv', sep = ',', header = T)
#Get m and s
m = sum(startsWith(names(time_series), 'm'))
s = sum(startsWith(names(time_series), 's'))
n_sim = length(composition$n_simulation == 0)
#Get column names that start with s (for strain), and n (for n_simulation)
ind_strains = which(startsWith(names(composition), 's') |
                      startsWith(names(composition), 'n'))
#Filter those rows to get abundance and richness of each surviving strain
strain_comp = composition[,ind_strains]
#Get abundance and richness of each surviving strain
strain_comp_long = melt(strain_comp, id.vars = c('n_simulation'),
                        variable.name = 'strain',
                        value.name = 'stable.state')
#Eliminate extinctions and outliers
surv_comp_long = strain_comp_long[!(strain_comp_long$stable.state>1e3|
                                    strain_comp_long$stable.state<1),]
#Sort conveniently
surv_comp_long = surv_comp_long[order(surv_comp_long$n_simulation,
                                      surv_comp_long$strain),]
#Get a table of frequencies of richnesses
richness_tab = table(surv_comp_long$n_simulation)
#Create a vector to add to the dataframe
richness = rep(richness_tab, richness_tab)
#Add vector to dataframe
surv_comp_long['richness'] = richness

#Load data about reaction network of species#
network = read.table('../data/coal_network.csv', sep = ',', header = T,
                     colClasses = c('n_simulation'='numeric','strain'='numeric','substrate'='character',
                                    'substrate'='character','maintenance'='numeric', 'surplus'='numeric',
                                    'n_reac'='numeric', 'cohesion'='numeric','competition'='numeric', 
                                    'autocohesion'='numeric'), row.names = 1)
#Change format of one column quickly
network$strain = paste('s', network$strain, sep = '')

#Get proxy of individual fitness for all species in network
it = length(network$substrate)
#Prealocate
substrate_vectors = matrix(0, nrow = it, ncol = m)
fitness_proxy = rep(0, it)
for (i in seq(it)){
  #Get string of substrates
  substrates_i = network$substrate[i]
  #Extract numbers into a list
  matches <- regmatches(substrates_i, gregexpr("[[:digit:]]+", substrates_i))
  #Unlist to a vector
  vec_subs_i = as.numeric(unlist(matches))
  #Add vector to matrix
  substrate_vectors[i,1:length(vec_subs_i)] = vec_subs_i
  #Calculate fitness proxy
  fitness_proxy[i] = network$surplus[i] + 
    2*length(unique(vec_subs_i)) - length(vec_subs_i)
}
#Add fitness_proxy to the network dataframe
network['fitness'] = fitness_proxy

#############################################################################

#Merge with surv_comp_long
all_data = merge(surv_comp_long, network, by = c('n_simulation', 'strain'))
#Get facilitation and competition indices each community
facilitation = aggregate(all_data[,16], list(all_data$n_simulation), mean)
competition = aggregate(all_data[,17], list(all_data$n_simulation), mean)
#Assign each community to a group based on their comp and fac levels.
colors = distance(competition$x, facilitation$x, 1.4, 0.2)
#Get a unique value for richness per simulation
unique_richness = non_contiguous(all_data)
#Create a vector of teams for members in the community
compete_vector = rep(colors, unique_richness)
#Sort all_data according to simulation before merging
all_data = all_data[order(all_data$n_simulation),]
#Add the compete_vector to all_data
all_data['team'] = compete_vector 
#Save data
write.csv(all_data, '../data/community_data.csv', row.names = F)
##########################################################################################

###Plot histogram of richnesses###
##########################################################################################
#Get richness
v = all_data$richness
#To remove contiguous duplicated elements only
v = v[c(TRUE, !v[-length(v)] == v[-1])]
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
                     breaks = c(1,3, 5, 7, 9), 
                     labels = c('1','3', '5', '7', '9'))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = 'Richness', y = 'Counts')+
  ggsave(filename = '../results/histogram_richness.pdf', width = 4.3, height = 4)

################################################################################

###Plot one example of time series###
###############################################################################
#Load time series data
time_series = read.table('../data/coal_time_series.csv', sep = ',', header = T)
#Select rows of interest
ind_strains_t = which(startsWith(names(time_series), 's')|
                        startsWith(names(time_series), 'n')|
                        startsWith(names(time_series), 't'))
#Filter to get those rows
strains_time_s = time_series[,ind_strains_t]
#Put time_series of strains into long format
time_series_long = melt(strains_time_s,
                     variable.name = 'strain',
                     id.vars = c('n_simulation','t'),
                     value.name = 'population')



#Plot and save one example of community assembly for report.
data_nice = time_series_long[which(time_series_long$n_simulation == 128),]
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
    legend.position = c(0.1, 0.7),
    legend.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size = 1),
    axis.title = element_text(size = 20), 
    legend.key = element_blank(), 
    axis.text = element_text(size = 15),
    aspect.ratio = 1
    #plot.margin=unit(c(2,0,-2,0),"cm")
  )+
  scale_x_continuous(limits = c(1, 5e4), trans = 'log10', expand = c(0, 0), 
                     labels = c('1', '1e2', '1e3')) + 
  scale_y_continuous(limits = c(-1, 35), expand = c(0, 0))+
  labs(x = 't [a.u.]',
       y = expression(paste('Abundance ',italic(N[alpha]),' [a.u.]')))+
  guides(color=guide_legend(expression('n'[r]), 
                            title.hjust = 0.5))+
  scale_color_manual('', values=paleta)+
  ggsave(filename = '../results/model_example.pdf', width = 4.79, height = 3.53 )

#########################################################################################
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
        legend.position = c(0.2, 0.43),
        legend.background = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 23),
        legend.key.size = unit(10, 'mm'))+
  scale_x_continuous(expand = c(0, 0), 
                     breaks = c(1, 5, 10, 15), 
                     labels = c('1', '5', '10', '15'))+
  scale_y_continuous(expand = c(0, 0),
                     breaks = rev(c(1, 5, 10, 15)), 
                     labels = c('1', '5', '10', '15'))+
  scale_fill_gradientn(colours = paleta,
                       breaks = c(0, 35, 70))+
  labs(fill = expression(Omega[xy]), 
       x = expression(paste(italic(y), ' (product)')), 
       y = expression(paste(italic(x), ' (substrate)')))+
  coord_fixed(ratio = 1)+
  ggsave(filename = '../results/community_reaction_network.pdf', width = 4.3, height = 4)
##########################################################################################

#Plot scatter plot of mean abundance at equilibrium and proportion of survivors
#as a function of the number of reactions.
#Prealocate matrix of results
nrows = max(all_data$n_reac)
ncols = length(unique(all_data$n_simulation))
reactions_community = rep(0, nrows)
for (i in seq(ncols)){
  #Get a vector of number of reactions of each surviving strain
  n_reactions_i = all_data$n_reac[all_data$n_simulation == i]
  #Get frequencies of each reaction number
  n_reac_freq = table(n_reactions_i)
  #Add the number of strains with n_reac reactions to each corresponding bin
  ind = as.numeric(names(n_reac_freq))
  reactions_community[ind] = reactions_community[ind] + as.numeric(n_reac_freq)
}

#Get vector of counts in number of reactions before community assembly
h = hist(network$n_reac, breaks = seq(16)-0.5)$counts
#Get proportion of survivors after community assembly
alpha = reactions_community/h
###Plot mean abundance and alpha as a function of reaction number###
##########################################################################################
#Get number of different reaction numbers
levels_reac = length(unique(all_data$n_reac))
#Get colors pallette
paleta_reac = wes_palette("Zissou1", levels_reac, type = 'continuous')
#Create dataframe for plotting
#Dealing with mean abundance
num_reac_fitness = aggregate(stable.state~n_reac, all_data, mean)
#Need to divide by square root of samples
#Get number of elements in each reaction number
len = aggregate(stable.state~n_reac, all_data, length)$stable.state
std_ =  aggregate(stable.state~n_reac, all_data, sd)
num_reac_fitness['err'] = std_$stable.state/(sqrt(len)*max(num_reac_fitness$stable.state))
#Normalize to 1 in order to plot in the same place
num_reac_fitness$stable.state = num_reac_fitness$stable.state/max(num_reac_fitness$stable.state)
#################################################
#Calculate standard error of binomial distribution for survival probabilities
se_binom = sqrt(alpha*(1-alpha)/len)
#This is hovigs section
num_reac_fitness[16:30,] =  c(seq(max(num_reac_fitness$n_reac)), 
                             alpha, 
                             se_binom)
#Add vector specifying group of points that we are in
num_reac_fitness['carrcap_survprop'] = c(rep(0, max(num_reac_fitness$n_reac)),
                                         rep(1, max(num_reac_fitness$n_reac)))
#################################################
#Add errorbars to red points by calculating the standard error of survival proportion 
#across simulations. I think this can be done using a binomial distribution.
#Add second axis to put the name of the red points as a variable
#Plot median stable state for each reaction number
n_reac_stablestate = ggplot(num_reac_fitness, aes(x = n_reac, 
                                                  y = stable.state,
                                                  colour = as.factor(carrcap_survprop)))+
  geom_errorbar(aes(ymin=stable.state-err, ymax=stable.state + err), 
                width=.2)+
  geom_point(size = 3)+
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size = 1),
        axis.text = element_text(size = 15), 
        axis.title.x = element_text(size = 17),
        axis.title.y.left = element_text(size = 15, color = "#3B9AB2"),
        axis.title.y.right = element_text(size = 17, color = "#F21A00"),
        legend.position = 'none',
        aspect.ratio = 1,
        axis.text.y.right = element_text(colour = "#F21A00"),
        axis.line.y.right = element_line(colour = "#F21A00"),
        axis.text.y.left = element_text(colour = "#3B9AB2"),
        axis.line.y.left = element_line(colour = "#3B9AB2"),
        axis.ticks.y.left = element_line(colour = "#3B9AB2", size = 1.25),
        axis.ticks.y.right = element_line(colour = "#F21A00", size = 1.25))+
  scale_x_continuous(breaks = c(1, 5, 10, 15), 
                     labels = c('1', '5', '10', '15'))+
  scale_y_continuous(breaks = c(0.3, 0.65, 1),
                     labels = c('0.3', '0.65', '1'),
                     sec.axis = sec_axis(trans = ~ . ,
                                         breaks = c(0.3, 0.65, 1),
                                         labels = c(),
                                         name = expression(paste("Survival probability ", 
                                         italic(p)))))+
  scale_color_manual(values = c("#3B9AB2", "#F21A00"))+
  labs(x = expression(paste('Number reactions ',italic(n[r]))), 
       y = expression(paste('Rescaled mean abundance ',
                            italic(symbol("\341")*n[infinity]*symbol("\361")))))+
  ggsave(filename = '../results/n_reac_stablestate.pdf', width = 4.3, height = 4)

########################################################################################
  bar_plot = ggplot(data=data_barplot, aes(x=n_reac, 
                                y=n_strains, 
                                fill=as.factor(random_assembled))) +
    geom_bar(stat="identity", color="black", position=position_dodge(),
             width = 0.7)+
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size = 1),
          axis.text = element_text(size = 15), 
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 19),
          legend.position = c(0.35, 0.15),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          legend.key.size = unit(5, 'mm'),
          aspect.ratio = 1,
          axis.text.y.right = element_text(colour = "#F21A00"),
          axis.line.y.right = element_line(colour = "#F21A00"),
          axis.text.y.left = element_text(colour = "#3B9AB2"),
          axis.line.y.left = element_line(colour = "#3B9AB2"))+
    scale_fill_manual(values = c("#3B9AB2", "#F21A00"),
                      labels = c('Random', 'Assembled'))+
    scale_x_continuous(limits = c(0.5, 15.5),
                       breaks = c(1,5,10,15))+
    scale_y_continuous(breaks = c(0, 0.5, 1),
                       labels = c('0', '750', '1500'),
                       expand = c(0,0),
                       sec.axis = sec_axis(trans = ~ . /2,
                                           breaks = c(0, 0.25, 0.5),
                                           labels = c('0','375', '750')))+
    labs(x = expression(italic(n[r])),
         y = 'Frequency')
    
    
##########################################################################################
#Put all plots together
  tot_assembly = ggdraw() +
    draw_plot(histogram_richness, x = 0, y = 0, width = 0.5, height = 0.5) +
    draw_plot(n_reac_stablestate, x = 0.5, y = 0, width = 0.5, height = 0.5) +
    draw_plot(model_example, x = 0, y = 0.5, width = 0.5, height = 0.5)+
    draw_plot(community_reaction_network, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
    draw_plot_label(label = c("A", "B", "C", "D"), size = 25,
                    x = c(0,0.49,0,0.49), y = c(1,1,0.5,0.5))+
    ggsave(filename = '../results/community_assembly_plots.pdf', width = 10, height = 8)


##########################################################################################

#Evolution of interactions during community coalescence

#Load data
inter_evol = read.table('../data/interaction_evolution.csv', sep = ',', header = T,
                        row.names = 1)


#Get rid of 0 and high times for nice plot
#Get rid of 0 so tha the mean doesn't go down artificially
inter_evol_no0 = inter_evol[!(inter_evol$survivor == 0 & 
                              inter_evol$interaction == 0 |
                              inter_evol$t > 500),]

data = aggregate(inter_evol_no0[,4:6], 
                 list(inter_evol_no0$t,
                      inter_evol_no0$n_simulation), 
                 sum)

#Bin time for nice plotting
h_t = hist(data$Group.1, 
           breaks=seq(floor_dec(min(data$Group.1, 1)),
                      ceiling_dec(max(data$Group.1), 1), 
                      length.out = 70))
#Find the interval containing each element of delP2 in breaks.
ind_av = findInterval(data$Group.1, h_t$breaks)
#Bin similar data together
time = rep(0, length(unique(ind_av)))
interaction = rep(0, length(unique(ind_av)))
facilitation = rep(0, length(unique(ind_av)))
competition = rep(0, length(unique(ind_av)))
interaction_error = rep(0, length(unique(ind_av)))
j = 1
for (i in (unique(ind_av))){
  #Average equivalent measurements by bin
  av_data = data[which(ind_av == i),]
  time[j] = median(av_data$Group.1)
  interaction[j] = median(av_data$interaction)
  facilitation[j] = median(av_data$facilitation)
  competition[j] = median(av_data$competition)
  interaction_error[j] = sd(av_data$interaction)
  j = j+1
}
averaged_data = data.frame('time'= time, 
                           'interaction' = interaction,
                           'competition' = competition,
                           'facilitation' = facilitation)
plot(time, interaction)
plot(time, competition, col = 'red')
points(time, facilitation, col = 'green')

interaction_plot = ggplot(data = averaged_data, aes(x = time, y = interaction))+
  geom_point(colour = "#F21A00", size = 4)+
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size = 2),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 25),
        aspect.ratio = 1)+
  scale_x_continuous(breaks = c(0, 250, 500))+
  scale_y_continuous(breaks = c(-1,-0.55,-0.1),
                     labels = c('-1','-0.5', '0'))+
  labs(x = 't [a.u.]', 
       y = expression(paste(symbol("\341"),Theta,symbol("\361"))))
  
#Plot inset
subplot_data = melt(averaged_data, id.vars = c('time'))
subplot_data = subplot_data[subplot_data$variable != 'interaction',]
inset_plot = ggplot(data = subplot_data, aes(x = time, y = value,
                                colour = as.factor(variable)))+
  geom_point(size = 2)+
  theme_minimal()+
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size = 1),
        axis.title = element_blank(),
        axis.text = element_text(size = 15), 
        axis.ticks = element_line(),
        legend.position = c(0.75, 0.7),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.spacing.x = unit(0, 'mm'),
        aspect.ratio = 1)+
   guides(color = guide_legend(override.aes = list(size=5)))+
   scale_x_continuous(breaks = c(0, 250, 500),
                      limits = c(0,500),
                      expand = c(0.02,0.02))+
   scale_y_continuous(breaks = c(1.5,2.5,3.5))+
  scale_colour_manual(values = c("#EBCC2A", "#3B9AB2"),
                      name = "", labels = c(expression(paste(symbol("\341"),
                                                             italic(bar(C)),
                                                             symbol("\361"))), 
                                             expression(paste(symbol("\341"),
                                                              italic(bar(F)),
                                                              symbol("\361")))))+
  labs(fill = '')
  

###############################################################################################

#Plot network of bacteria before and after the assembly

#Load data
inter_random_df = read.table('../data/interaction_random.csv', sep = ',',
                          header = T)
inter_assembled_df = read.table('../data/interaction_assembled.csv', sep = ',',
                             header = T)
#Initialize vector of survival score
surv_score = rep(0, ncol(inter_random_df)-1)
surv_score_weighted = rep(0, ncol(inter_random_df)-1)
surv_score_fitness = rep(0, ncol(inter_random_df)-1)
#Initialize matrix of number of reactions
reac_number_ranked = matrix(0, ncol(inter_random_df)-1, 
                            length(unique(inter_random_df$n_simulation)))
#Initialize vector of extinctions
number_extinctions = rep(0,  length(unique(inter_random_df$n_simulation)))

#Plot random and assembled networks for all simulations
pb = txtProgressBar(1, length(unique(inter_random_df$n_simulation)), style=3)
pdf("../results/interactions.pdf") 
for (i in unique(inter_random_df$n_simulation)){
  #Transform to matrix
  inter_random_i = data.matrix(inter_random_df[inter_random_df$n_simulation==i,2:ncol(inter_random_df)])
  #Get number of survivors 
  n_surv = length(all_data$strain[all_data$n_simulation == i])
  #Get rid of 0 in the assembled matrix
  inter_assembled_i = data.matrix(inter_assembled_df[inter_assembled_df$n_simulation == i, 2:ncol(inter_assembled_df)])[0:n_surv,0:n_surv]
  if (length(inter_assembled_i) == 0 | length(inter_assembled_i) == 1){
    next
  }
  #Get rid of column names
  colnames(inter_random_i) = NULL
  colnames(inter_assembled_i) = NULL
  #Create graph
  g_random = graph_from_adjacency_matrix(inter_random_i, mode = 'undirected',
                                         weighted = T)
  #Get a pallette with as many colours as different weights
  n_weights = length(E(g_random)$weight)
  colours = wes_palette('Zissou1', type = 'continuous', n = n_weights)
  #Get indices that would order the array of weitghts
  ordered_index = order(order(E(g_random)$weight, decreasing = F))
  #Create a color dataframe
  color_df = data.frame(weight = E(g_random)$weight, 
                        color = colours[ordered_index])
  #Get colors for edges
  E(g_random)$color = colours[ordered_index]
  #Add width to edges based in the weights
  E(g_random)$width = 5*range01(E(g_random)$weight)
  #Color for nodes
  #Get a pallette with as many colours as different node_strengths
  degree_random = strength(g_random)
  n_nodes = length(degree_random)
  colours_nodes = wes_palette('Zissou1', type = 'continuous', n = n_nodes)
  #Get indices that would order the array of node strengths
  #Get size of nodes.
  reac_number = network$n_reac[network$n_simulation == i]
  V(g_random)$size = 3*(5+reac_number)
  ordered_index = order(order(degree_random, decreasing = F))
  V(g_random)$color = colours_nodes[ordered_index]
  #Add labels to nodes
  label_random = colnames(inter_random_df)[2:ncol(inter_random_df)]
  #Get individual fitness of species in the community
  fitness = network$fitness[network$n_simulation == i]
  cohesion_rank = data.frame('label' = label_random, 
                             'cohesion_rank' = order(order(degree_random, 
                                                           decreasing = T)),
                             'fitness_rank' = order(order(fitness,
                                                          decreasing = T)),
                             'reac_number' = reac_number)
  V(g_random)$label = label_random
  #Plot
  par(mfrow=c(2,1),
      mar = c(2, 0, 0, 0),
      oma = c(0,0,0,0),
      xpd = NA)
  plot(g_random, layout = layout_in_circle(g_random, order = order(degree_random,
                                                                   decreasing = T)),
       vertex.label.cex = 1.8)
  
  #Create assembled grapph
  g_assembled = graph_from_adjacency_matrix(inter_assembled_i, mode = 'undirected',
                                            weighted = T)
  #Add width to edges based in the weights
  E(g_assembled)$width = 5*range01(E(g_assembled)$weight)
  #Get colors for edges
  indices = match(E(g_assembled)$weight, color_df$weight)
  E(g_assembled)$color = as.character(color_df$color[indices])
  #Color for nodes
  #Get a pallette with as many colours as different node_strengths
  degree_assembled = strength(g_assembled)
  n_nodes = length(degree_assembled)
  colours_nodes = wes_palette('Zissou1', type = 'continuous', n = n_nodes)
  #Get indices that would order the array of node strengths
  ordered_index = order(order(degree_assembled, decreasing = F))
  #Assign colours
  V(g_assembled)$color = colours_nodes[ordered_index]
  #Get size of nodes.
  reac_number = all_data$n_reac[all_data$n_simulation == i]
  V(g_assembled)$size = 3*(6+reac_number)
  #Get labels of surviving sepecies
  surv_labels = as.character(all_data$strain[all_data$n_simulation == i])
  number_extinctions[i] = length(label_random) - length(surv_labels)
  #Add labels to nodes
  V(g_assembled)$label = surv_labels
  #Get original cohesion rank from surviving species
  rank_surv = cohesion_rank$cohesion_rank[cohesion_rank$label %in% surv_labels]
  rank_fitness = cohesion_rank$fitness_rank[cohesion_rank$label %in% surv_labels]
  #Add to vector of survival proportion
  surv_score[rank_surv] = surv_score[rank_surv] + 1
  surv_score_fitness[rank_fitness] = surv_score_fitness[rank_fitness] + 1 
  #Get rank in assembled community
  rank_assembled = order(order(degree_assembled))
  #Add to vector of weighted survival population
  cohesion_rank_assembled = data.frame('label' = surv_labels, 
                                       'cohesion_rank' = rank_assembled)
  #Get weights
  weights = cohesion_rank_assembled$cohesion_rank/max(cohesion_rank_assembled$cohesion_rank)
  surv_score_weighted[rank_surv] = surv_score_weighted[rank_surv] + weights
  
  #Get number of reactions of surviving species
  reac_number_surv = cohesion_rank$reac_number[cohesion_rank$label %in%
                                                 surv_labels]
  #Add to the matrix of ranked number of reactions
  reac_number_ranked[rank_surv, i] = reac_number_surv  
  #Plot
  plot(g_assembled, 
       layout = layout_in_circle(g_assembled, 
                                 order = order(degree_assembled, decreasing = T)),
       vertex.label.cex = 1.8)
  setTxtProgressBar(pb, i)
}
dev.off()

###############################################################################

#Plot presence proportion in the assembled community as a function of the cohesion rank 
#in the random community. Blue curve is the weighted version of the plot by the rank in the
#assembled community.

#Prepare data
surv_score_unit = surv_score/length(unique(inter_random_df$n_simulation))
surv_score_weighted_unit = surv_score_weighted/length(unique(inter_random_df$n_simulation))
surv_score_fitness_unit = surv_score_fitness/length(unique(inter_random_df$n_simulation))
#Prealocate vector of reaction number median
number_reac_median_ranked = rep(0, 10)
n_samples = rep(0,10)
std_samples = rep(0,10)
#Get median number of reactions at each rank
for (i in seq(nrow(reac_number_ranked))){
  #Get vector of reaction number at rank position i
  reac_number_i = reac_number_ranked[i,]
  #Eliminate all 0
  reac_number_i = reac_number_i[reac_number_i != 0]
  #Get number of samples
  n_samples[i] = length(reac_number_i)
  #Get standar deviations of each rank
  std_samples[i] = sd(reac_number_i)
  #Calculate median
  number_reac_median_ranked[i] = mean(reac_number_i)
}
#Calculate vector of errors
k = length(unique(inter_random_df$n_simulation))
err_cohesion = sqrt((surv_score/k*(1 - surv_score/k))/k)
err_fitness = sqrt((surv_score_fitness/k*(1 - surv_score_fitness/k))/k)
#Create data for plot
presence_proportion = data.frame('x' = rep(seq(10), 2),
                                 'y' = c(surv_score_fitness_unit, surv_score_weighted_unit),
                                 'group' = c(rep(0, 10), rep(1, 10)),
                                 'number_reac' = rep(number_reac_median_ranked,2),
                                 'err' = c(err_cohesion, err_fitness))
#Get number of different reaction numbers
levels_reac = length(unique(presence_proportion$number_reac))
#Get colors pallette
paleta_reac = wes_palette("Zissou1", levels_reac, type = 'continuous')
survival_proportion = ggplot(data = presence_proportion, 
       aes(x = x, y = y))+
  geom_hline(yintercept=median(number_extinctions)/s, linetype="dashed", size=1)+
  geom_point(aes(shape = as.factor(group),
                 fill = number_reac), size = 4)+
  geom_errorbar(aes(ymin = y-err, ymax= y + err), 
                width=.2)+
  annotate(geom="text", x=8, y=0.56, label=expression(paste(symbol("\341"),
                                                          italic(p),
                                                          symbol("\361"))),
            size = 8, parse = TRUE)+
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size = 2),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20),
        aspect.ratio = 1,
        legend.position = c(0.27, 0.25),
        legend.background = element_blank(),
        legend.margin = margin(t = -15, r = 0, b = -15, l = 0, unit = "pt"),
        legend.key = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17))+
  scale_fill_gradientn(colors = paleta_reac,
                       breaks = c(6,7,8,9),
                       guide = guide_colourbar(direction = "horizontal", 
                                                title.position = "top"))+
  scale_shape_manual(values = c(24,21),
                     labels = c('performance', expression(paste('cohesion ',
                                                                italic(s[alpha])))),
                     guide = guide_legend())+
  # scale_colour_manual(values = c("#F21A00", "#3B9AB2"),
  #                     labels = c('weighted', 'unweighted'))+
  scale_x_continuous(breaks = c(1, 10))+
  scale_y_continuous(limits = c(0,0.7),
                     breaks = c(0,0.35,0.7),
                     labels = c('0', '0.35', '0.7'))+
  guides(fill = guide_colourbar(order = 1, 
                                direction = 'horizontal',
                                title.position = "top"), 
         shape = guide_legend(order = 2))+
  labs(fill = expression(paste(symbol("\341"),italic(n[r]),symbol("\361"))), 
       shape = '',
       x = 'Rank',
       y = expression(paste('Survival probability ', italic(p))))+
  ggsave('../results/survival_proportion.pdf')


ggdraw() +
  draw_plot(interaction_plot, y = 0, height = 0.5) +
  draw_plot(inset_plot, x = .45, y = .09, width = .4, height = .25) +
  draw_plot(survival_proportion, y = 0.5, height = 0.5)+
  draw_plot_label(label = c("C", "D"), size = 25,
                  x = c(0, 0), y = c(1,0.5))+
  ggsave(filename = '../results/interaction_evolution_plots.pdf', 
         height = 8,
         width = 5)

##############################################################################
#Draw a plot of correlation between abundance at equilibrium and my proxy for

#Load data
fitness_abundance = read.table('../data/abundance_fitness.csv', sep = ',', header = T,
                         row.names = 1)
#Get rid of outlier
fitness_abundance = fitness_abundance[fitness_abundance$K < 1000 ,]

#Get table of frequencies of values
tab_K = table(fitness_abundance$K)
#Get very repeated indices
rep_inds = which(as.vector(table(fitness_abundance$K)) >= 6)
#Find values for those indices
rep_values = as.numeric(names(tab_K)[rep_inds])
#Get rid of them
fitness_abundance_norep = fitness_abundance[!(fitness_abundance$K >= min(rep_values)-0.0001&
                                          fitness_abundance$K <= max(rep_values)+0.0001),]

fitness_abundance = fitness_abundance_norep
#Bin and average plot
#Average plot
h = hist(fitness_abundance$fitness, 
         breaks=seq(floor_dec(min(fitness_abundance$fitness, 1)),
                    ceiling_dec(max(fitness_abundance$fitness), 1), 
                    length.out = 10))
#Find the interval containing each element of delP2 in breaks.
ind_av = findInterval(fitness_abundance$fitness, h$breaks)
#Bin similar data together
fit_av = rep(0, length(unique(ind_av)))
K_av = rep(0, length(unique(ind_av)))
K_std = rep(0, length(unique(ind_av)))
j = 1
for (i in (unique(ind_av))){
  #Average equivalent measurements by bin
  av_data = fitness_abundance[which(ind_av == i),]
  fit_av[j] = median(av_data$fitness)
  K_av[j] = median(av_data$K)
  K_std[j] = sd(av_data$K)/sqrt(length(av_data$K))
  j = j+1
}
#Create dataframe with averaged measures
average_fitness_abundance = data.frame('fitness' = fit_av, 
                                       'K' = K_av, 
                                       'K_err' = K_std)
#Plot correlation

fitness_abundance = ggplot(average_fitness_abundance, aes(x = fitness, y = K_av))+
  geom_errorbar(aes(ymin=K_av-K_err, ymax=K_av + K_err), 
                width=.2,
                position=position_dodge(.9),
                linetype = 'dashed')+
  geom_point(size = 3)+
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size = 1),
        axis.text = element_text(size = 15), 
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 19),
        legend.position = c(0.85, 0.7),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.background = element_blank(),
        legend.key.size = unit(10, 'mm'),
        aspect.ratio = 1)+
  # scale_x_continuous(breaks = c(1, 5, 10, 15), 
  #                    labels = c('1', '5', '10', '15'))+
  labs(x = 'Individual Fitness', 
       y = expression(italic(symbol("\341")*N[infinity]^0*symbol("\361"))))+
  ggsave(filename = '../results/fitness_abundance.pdf', width = 4.3, height = 4)

