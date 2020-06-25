
plot_time_series = function(data, sim, col){
  #Plot time series data
  p = ggplot(data, aes(x = t, y = population, color = col))+
    ggtitle(sim)+
    geom_line() +
     theme(
          #legend.title = element_text(),
          legend.position = "none") + 
    scale_x_continuous(trans = "log10") 
  return(p)
}

plot_network_colors = function(mat, color){
  #Plot a network, given the matrix of reactions
  paint = cbind(melt(mat), tot_colors)
  #Concatenate column of colors
  p = ggplot(paint, aes(Var2, Var1)) +
    geom_tile(aes(fill = tot_colors, colour = tot_colors, alpha = try(log(value))), color = NA)+
    scale_y_reverse()+
    theme(panel.background = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(0, 0, -20,0),
          legend.title = element_blank())+
    scale_fill_manual(values=unique(tot_colors))
  return(p)
}

plot_network = function(mat){
  #Plot a network, given the matrix of reactions
  paint = melt(mat)
  #Concatenate column of colors
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

extract_vector = function(string){
  #Pattern to extract: numbers of one, two, or three digits
  pattern = '\\d{1,3}'
  #Extract vector of integers from a string
  vector = as.numeric(str_extract_all(string, pattern)[[1]])
  return(vector)
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
    mat = get_network(networks[i,1], networks[i,2])
    w_mat = abundances[i]*mat
    community_network = community_network + w_mat
  }
  return(community_network)
}

rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)

all_networks = function(networks, abundances, m){
  #Preallocate community network matrix
  store_here = list()
  community_network = matrix(0, m, m)
  for (i in seq(length(abundances))){
    mat = get_network(networks[i,1], networks[i,2])
    w_mat = abundances[i]*mat
    store_here[[i]] = w_mat
  }
  return(store_here)
}

color_mixer = function(networks, abundances, m){
  store_here = all_networks(networks, abundances, m)
  tot = melt(store_here)
  #Get as many colors as species from the rainbow
  colors = rainbow(length(abundances))
  #Transform colors to RGB system
  colors_rgb = col2rgb(colors)
  #Transform colors to RGB system
  #colors_hsl = col2hsl(colors)
  #Get normalized abundances
  norm_abundanceas = abundances / sqrt(sum(abundances^2))
  #Prealocate matrix of colors 
  colors_mat = matrix('0', m, m)
  #Traverse the matrix
  for (i in seq(m)){
    for (j in seq(m)){
      #Get element i j from all strains
      elements_ij = tot[which(tot$Var1 == i & tot$Var2 == j),]$value
      #Normalize elements_ij
      abundances_norm = elements_ij / sum(elements_ij)
      #Mix colors of the rainbow weighted by the abundance at element_ij
      if (is.finite(abundances_norm[1])){
        weighted_colors = sweep(colors_rgb, MARGIN = 2, abundances_norm, '*')
      } else {
        #Elements are all 0
        weighted_colors = sweep(colors_rgb, MARGIN = 2, elements_ij, '*')
      }
      #Mix colors
      mixed = rowSums(weighted_colors)
      #Transform to hex and assign to that matrix element
      #colors_mat[i,j] = hsl2col(t(t(mixed)))
      colors_mat[i,j] = rgb2hex(mixed['red'], mixed['green'], mixed['blue'])
      
    }
  }
  tot_colors = melt(colors_mat)
  return(tot_colors$value)
}

time_series_eliminator = function(time_series_melt, s){
  #Eliminates time series that end up extinct.
  #Get length of time 
  t_points = length(unique(time_series_melt$t))
  #Get number of simulations
  n_simul = length(unique(time_series_melt$n_simulation))
  for (i in seq(n_simul)){
    #Loop over strains
    for (j in seq(s)){
      #Get data for strain j and simulation i
      #Get data for that simulation
      str = paste0('s',j-1)
      inds = (time_series_melt$n_simul == i-1) & 
             (time_series_melt$strain == str)
      sub_data = time_series_melt[inds,]
      if (sub_data$population[t_points] < 1){
        #If species went extinct remove the whole time series
        time_series_melt = time_series_melt[!inds,]
      }
    }
  }
  return(time_series_melt)
}
