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
    # if (i == 123){
    #   browser()
    # }
    #
    #Get indices of networks corresponding to simulation i-1
    inds = which(surv_networks$simulation == i-1)
    #Get networks for that simulation
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
  #Add up all reaction networks in a community weighted by abundance
  #Preallocate community network matrix
  community_network = matrix(0, m, m)
  for (i in seq(length(abundances))){
    mat = get_network(networks[i,3], networks[i,4])
    w_mat = abundances[i]*mat
    community_network = community_network + w_mat
  }
  return(community_network)
}

coal_plot_time_series = function(time_series, surv_networks, prediction, m){
  #Plot the all simulations and save to pdf#
  ##########################################
  sections = c(0,which(time_series$t == max(time_series$t)))
  pdf('../results/coal_populations.pdf', width = 9.5, height = 4)
  #Set progress bar
  pb = txtProgressBar(1, n_sim, style=3)
  for (i in seq(n_sim)){
    #Trim time_series per simulation
    series_i = time_series[(sections[i]+1):sections[i+1],]
    p = pop_time_series(time_series = series_i, s)
    
    #Get reaction networks of surviving species in simulation i
    # if (i == 123){
    #   browser()
    # }
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