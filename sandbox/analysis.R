#setwd("~/Desktop/master_thesis/code")
#Load libraries#
library(stringr)
require(ggplot2)
library(reshape2)
library(viridis)
library(grid)
library(gridExtra)

source('analysis_functions.R')

#Load data#
network = read.table('../data/network.csv', sep = ',', header = T, 
                           colClasses = c('character'), row.names = 1)
composition = read.table('../data/composition.csv', sep = ',', header = T,
                         row.names = 1)
time_series = read.table('../data/time_series.csv', sep = ',', header = T)
performance = read.table('../data/performance.csv', sep = ',', header = T,
                         row.names = 1)


#Functions#

#Code#
#Get number of strains and metabolites
network$X = label_it(network, c('simulation', 'strain'))
lab = network$X
s = sum(startsWith(lab, '0'))
m = length(composition) - s
n_sim = dim(composition)[1]

#Abundance of each surviving strains
strains_composition = composition[,1:s]
survivor_label = lab_survivor(strains_composition)
#What are the characteristics of survivors?
#Winners networks
ind_win = which(network$X %in% survivor_label)
win_network = network[ind_win,]
#Loosers networks
ind_loose = which(!network$X %in% survivor_label)
loose_network = network[ind_loose,]
#Compare mean harvested energies
# harvest_win_mean = mean(harvested_energy(win_network))
# harvest_loose_mean = mean(harvested_energy(loose_network))
#Found that harvest_win_mean is a bit higher, but not significantly.

#Lets keep going!
#Calculate individual performance. A nice measure of the 
#individual performance is the maximum growth rate.
#Drop the last m columns of the stable composition of dataframe
stable_strains = composition[1:(length(composition)-m)]
names(stable_strains) = seq(0,s-1)
#Put them together
one = melt(stable_strains, measure.vars = names(stable_strains),
           variable.name = 'strain',
           value.name = 'stable.state')
names(performance) = seq(0,s-1)
two = melt(performance, measure.vars = names(performance), 
           variable.name = 'strain',
           value.name = 'performance')
#Get richness at equilibrium for each simulation
richness = rowSums(strains_composition > 1)
#Repeat s times
richness = rep(richness, s)
#Set simulation number
sim_number = rep(seq(0,n_sim-1), s)

prediction = data.frame('Strain' = one$strain,
                        'Simulation' = sim_number,
                        'Performance' = two$performance, 
                        'Stable.State' = one$stable.state,
                        'Richness' = richness, 
                        row.names = seq(1,n_sim*200))

#Eliminate extinctions and outliers
prediction = prediction[!(prediction$Stable.State<1),]
prediction = prediction[which(predict_mean$Stable.State<1e5),]
prediction = prediction[(prediction$Performance<1000),]
prediction = prediction[order(prediction$Simulation, prediction$Strain),]

p = ggplot(data = prediction, aes(x = Performance, y = Stable.State))+
    theme_bw()+
    geom_point(aes(colour = Richness))+
    scale_color_gradientn(colours = rainbow(6))

#Plot the median instead of all of them
predict_mean = aggregate(prediction[,3:5], list(prediction$Simulation), mean)
p2 = ggplot(data = predict_mean, aes(x = Performance, 
                                       y = Stable.State, 
                                       color = as.factor(Richness)))+
     theme_minimal()+
     theme(legend.title = element_text(),
           legend.key = element_blank(),
           panel.border = element_rect(fill = NA))+
     geom_point()+
     geom_smooth(method = lm, se = F)+
     #scale_color_brewer(palette="Set1")+
     labs(color = 'Richness')

print(p2)

#Analyzing reaction networks
#Get reaction network of surviving species for each simulation
#Get label of surviving species
labels = paste(prediction$Simulation, prediction$Strain, sep = '_')
surv_networks = network[which(network$X %in% labels),]

#Plot time series along with community network
plot_time_series(time_series, surv_networks, prediction, m)
