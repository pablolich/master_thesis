require(ggplot2)
library(reshape2)
library(wesanderson)
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