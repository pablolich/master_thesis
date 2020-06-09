#setwd("~/Desktop/master_thesis/code")
#Load libraries#
library(stringr)
require(ggplot2)
library(reshape2)
library(viridis)
library(grid)
library(gridExtra)

#Load data#
network = read.table('../data/network.csv', sep = ',', header = T, 
                           colClasses = c('character'), row.names = 1)
composition = read.table('../data/composition.csv', sep = ',', header = T,
                         row.names = 1)
time_series = read.table('../data/time_series.csv', sep = ',', header = T)
performance = read.table('../data/performance.csv', sep = ',', header = T,
                         row.names = 1)


#Functions#
extract_vector = function(string){
  #Pattern to extract: numbers of one, two, or three digits
  pattern = '\\d{1,3}'
  #Extract vector of integers from a string
  vector = as.numeric(str_extract_all(string, pattern)[[1]])
  return(vector)
}

lab_survivor = function(composition){
  n_sim = dim(composition)[1]
  lab = c()
  for (i in seq(n_sim)){
    surv = as.character(which(composition[i,] > 1))
    lab = c(lab, paste(as.character(i-1), surv, sep = '_'))
  }
  return(lab)
}

label_it = function(data, cols){
  #Mergers contents in columns to create a single id column
  new_col = apply(data[, cols], 1, paste, collapse = "_")
  return(new_col)
}
pop_time_series = function(time_series, s){
  #Plot population time series for a single simulation
  population = time_series[,0:s+2]
  #Make long format
  dat = melt(population, id.vars = c('t'))
  p = ggplot(dat, aes(x = t, y = value, color = variable))+
      geom_line() +
      theme(legend.title = element_blank(),
            legend.position = "none") + 
      scale_x_continuous(trans = "log10") 
  return(p)
}

plot_network = function(mat){
  paint = melt(mat)
  p = ggplot(paint, aes(Var2, Var1)) +
      geom_tile(aes(fill = value), color = NA)+
      scale_y_reverse()+
      theme(panel.background = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            legend.margin = margin(0,0,0,0),
            legend.box.margin = margin(0, 0, -20,0),
            legend.title = element_blank())
  return(p)
}

harvested_energy = function(networks){
  #Get network length
  l = dim(networks)[1]
  #Prealocate vector of harvested energy per reaction network
  harvest = rep(0,l)
  for (i in seq(l)){
    #Extract network into a couple of vectors
    substrate = extract_vector(networks[i,2])
    product = extract_vector(networks[i,3])
    #Compute the energy difference per element
    harvest[i] = sum(product - substrate)
  }
  return(harvest)
}

extract_network = function(network){
  #Put one reaction network into a dataframe
  substrate = extract_vector(network[2])+1
  product = extract_vector(network[3])+1
  #Note new index notation because we are in R
  return(cbind(substrate, product))
}

plot_time_series = function(time_series, surv_networks, prediction, m){
  #Plot the all simulations and save to pdf#
  ##########################################
  sections = c(0,which(time_series$t == max(time_series$t)))
  pdf('../results/populations.pdf', width = 9.5, height = 4)
  #Set progress bar
  pb = txtProgressBar(1, n_sim, style=3)
  for (i in seq(n_sim)){
    #Trim time_series per simulation
    series_i = time_series[(sections[i]+1):sections[i+1],]
    p = pop_time_series(time_series = series_i, s)
    
    #Get reaction networks of surviving species in simulation i
    if (i == 123){
      browser()
    }
    inds = which(surv_networks$simulation == i-1)
    sub_network = surv_networks[inds,]
    abundances = prediction$Stable.State[inds]
    #Get community network
    com_network = community_network(sub_network, abundances, m)
    q = plot_network(com_network)
    #Plot next to eachother
    grid.arrange(p, q, nrow = 1)
    #Print progress...
    setTxtProgressBar(pb, i)
  }
  dev.off()
}

get_network = function(substrates, products){
  #Given a string of substrates and a string of products, yield
  #the corresponding reaction matrix
  sub_vec = extract_vector(substrates)
  prod_vec = extract_vector(products)
  ind = cbind(sub_vec+1, prod_vec+1)
  mat = matrix(0, m, m)
  mat[ind] = 1
  return(mat)
}

community_network = function(networks, abundances, m){
  #Add up all reaction networks weighted by abundance
  #Preallocate community network matrix
  community_network = matrix(0, m, m)
  for (i in seq(length(abundances))){
    mat = get_network(networks[i,3], networks[i,4])
    w_mat = abundances[i]*mat
    community_network = community_network + w_mat
  }
  return(community_network)
}
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
                        row.names = seq(1,150*200))

#Eliminate extinctions and outliers
prediction = prediction[!(prediction$Stable.State<1),]
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
labels = paste(prediction$Simulation, prediction$Strain, sep = '')
surv_networks = network[which(network$X %in% labels),]

#Plot time series along with community network
plot_time_series(time_series, surv_networks, prediction, m)
