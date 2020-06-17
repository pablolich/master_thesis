
plot_time_series = function(data, sim, col){
  #Plot time series data
  p = ggplot(data, aes(x = t, y = population, color = col))+
    ggtitle(sim)+
    geom_line() +
     theme(
          legend.title = element_text(),
          legend.position = "bottom") + 
    scale_x_continuous(trans = "log10") 
  return(p)
}

plot_network = function(mat){
  #Plot a network, given the matrix of reactions
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