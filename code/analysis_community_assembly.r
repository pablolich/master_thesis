require(ggplot2)
library(reshape2)
library(wesanderson)
require(cowplot)
source('coal_analysis_functions.R')
rm(list=ls())
dev.off()
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

###Plot mean abundance as a function of reaction number###
##########################################################################################
#Get number of different reaction numbers
levels_reac = length(unique(all_data$n_reac))
#Get colors pallette
paleta_reac = wes_palette("Zissou1", levels_reac, type = 'continuous')
#Plot median stable state for each reaction number
num_reac_fitness = aggregate(stable.state~n_reac, all_data, mean)
std_ =  aggregate(stable.state~n_reac, all_data, sd)
num_reac_fitness['err'] = std_$stable.state
n_reac_stablestate = ggplot(num_reac_fitness, aes(x = n_reac, y = stable.state,
                                                  colour = n_reac))+
  geom_point(size = 5)+
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
        legend.key.size = unit(10, 'mm'),
        aspect.ratio = 1)+
  scale_x_continuous(breaks = c(1, 5, 10, 15), 
                     labels = c('1', '5', '10', '15'))+
  labs(colour = expression('n'[r]), x = expression('n'[r]), 
       y = expression(symbol("\341")*N[infinity]*symbol("\361")))+
  scale_color_gradientn(colors = paleta_reac,
                        breaks = c(1,5,10,15))+
  ggsave(filename = '../results/n_reac_stablestate.pdf', width = 4.3, height = 4)

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
  labs(x = 'Richness', y = 'Frequency')+
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
        legend.position = c(0.2, 0.4),
        legend.background = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 25),
        legend.key.size = unit(10, 'mm'))+
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
##########################################################################################

#Put all plots together
#Arange all plots in a mega plot
tot = ggarrange(model_example, community_reaction_network,
                histogram_richness, n_reac_stablestate,
                labels = c("A", "B", "C", "D"),
                ncol =2, nrow = 2, 
                font.label = list(size = 25))

ggsave('../results/community_assembly_plots.pdf', plot = tot, height = 10, width = 10)

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
  geom_point(colour = "#F21A00", size = 6)+
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size = 2),
        axis.text = element_text(size = 25), 
        axis.title = element_text(size = 30),
        aspect.ratio = 0.7)+
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
  geom_point(size = 3)+
  theme_minimal()+
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size = 1),
        axis.title = element_blank(),
        axis.text = element_text(size = 20), 
        axis.ticks = element_line(),
        legend.position = c(0.8, 0.7),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.spacing.x = unit(0, 'mm'),
        aspect.ratio = 0.7)+
   guides(color = guide_legend(override.aes = list(size=6)))+
   scale_x_continuous(breaks = c(0, 250, 500),
                      limits = c(0,500),
                      expand = c(0.02,0.02))+
   scale_y_continuous(breaks = c(1.5,2.5,3.5))+
  scale_colour_manual(values = c("#EBCC2A", "#3B9AB2"),
                      name = "", labels = c(expression(paste(symbol("\341"),
                                                             italic(C),
                                                             symbol("\361"))), 
                                             expression(paste(symbol("\341"),
                                                              italic(F),
                                                              symbol("\361")))))+
  labs(fill = '')

ggdraw() +
  draw_plot(interaction_plot) +
  draw_plot(inset_plot, x = .34, y = .22, width = .75, height = .31) +
  ggsave(filename = '../results/interaction_evolution_plots.pdf', height = 8, width = 10 )
  

###############################################################################################

#Plot network of bacteria before and after the assembly


#Things I want to do:
#1.   Order vertices as a function of their strength
#2.   Size of the circle represents the number of reactions that it has.


#Load data
inter_random_df = read.table('../data/interaction_random.csv', sep = ',',
                          header = T)
inter_assembled_df = read.table('../data/interaction_assembled.csv', sep = ',',
                             header = T)

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
  if (length(inter_assembled_i) == 0){
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
  degree_ = strength(g_random)
  n_nodes = length(degree_)
  colours_nodes = wes_palette('Zissou1', type = 'continuous', n = n_nodes)
  #Get indices that would order the array of node strengths
  #Get size of nodes.
  reac_number = network$n_reac[network$n_simulation == i]
  V(g_random)$size = 3*(5+reac_number)
  ordered_index = order(order(degree_, decreasing = F))
  V(g_random)$color = colours_nodes[ordered_index]
  #Add labels to nodes
  V(g_random)$label = colnames(inter_random_df)[2:ncol(inter_random_df)]
  #Plot
  par(mfrow=c(2,1),
      mar = c(2, 0, 0, 0),
      oma = c(0,0,0,0),
      xpd = NA)
  plot(g_random, layout = layout_in_circle(g_random, order = order(degree_,
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
  degree_ = strength(g_assembled)
  n_nodes = length(degree_)
  colours_nodes = wes_palette('Zissou1', type = 'continuous', n = n_nodes)
  #Get indices that would order the array of node strengths
  ordered_index = order(order(degree_, decreasing = F))
  #Assign colours
  V(g_assembled)$color = colours_nodes[ordered_index]
  #Get size of nodes.
  reac_number = all_data$n_reac[all_data$n_simulation == i]
  V(g_assembled)$size = 3*(5+reac_number)
  #Get labels of surviving sepecies
  surv_labels = as.character(all_data$strain[all_data$n_simulation == i])
  #Add labels to nodes
  V(g_assembled)$label = surv_labels
  #Plot
  plot(g_assembled, 
       layout = layout_in_circle(g_assembled, 
                                 order = order(degree_, decreasing = T)),
       vertex.label.cex = 1.8)
  setTxtProgressBar(pb, i)
}
dev.off()
