#Quick plot
require(ggplot2)
scatterr = read.table('../data/scatterplot.csv', sep = ',', header = T, row.names = 1)
similarity_fitness = read.table('../data/similarity_fitness.csv', sep = ',',
                                header = T)
extinct_plot = read.table('../data/tile_plot.csv', sep = ',',
                                header = T)
subdata = extinct_plot[extinct_plot$richness == 4,]
ggplot(subdata, aes(x+1, y+1)) +
  geom_tile(aes(fill = as.factor(su)))+
  scale_fill_manual(values = c('cyan', 'black', 'magenta'))

 #Get rid of outliers
# similarity_fitness = similarity_fitness[abs(similarity_fitness$delF)<10,]
#Initialize similarity and p vectors
similarity = c()
fitness_dif = c()
richness_p = c()
for (j in unique(similarity_fitness$richness)[2:length(unique(similarity_fitness$richness))]){
  sub_data = similarity_fitness[similarity_fitness$richness == j,]
  #Get rid of 0
  data = sub_data[sub_data$similarity != 0,]
  #Average plot
  h = hist(data$delP2, breaks = seq(floor_dec(min(data$delP2, 1)),
                                    ceiling_dec(max(data$delP2), 1),  0.1)
  )
  ind_av = findInterval(data$delP2, h$breaks)
  #Bin similar data together
  P_av = rep(0, length(unique(ind_av)))
  sim_av = rep(0, length(unique(ind_av)))
  for (i in seq(length(unique(ind_av)))){
    #Average equivalent measurements by temperature
    av_data = data[which(ind_av == i),]
    P_av[i] = median(av_data$delP2)
    sim_av[i] = median(av_data$similarity)
  }
  similarity = c(similarity, sim_av)
  fitness_dif = c(fitness_dif, P_av)
  richness_p = c(richness_p, rep(j, length(sim_av)))
}

plot_data = data.frame(similarity, fitness_dif, richness_p)
ggplot(data = plot_data, aes(x = fitness_dif, y = similarity))+
  geom_point(aes(colour = as.factor(richness_p),
                 shape = as.factor(richness_p)))+
  theme(legend.title = element_text())+
  guides(color=guide_legend("my title"), shape = guide_legend('hi'))+
  scale_colour_discrete(labels = c('1', '2', '3', '4', '5', '6', '7'))
ggplot(data = sub_data, aes(x = delP2, y = similarity)) +
  geom_point(aes(colour = richness))+
  scale_color_gradient(low="black", high="red")
rich = unique(scatterr$richness)
plots = rep(NA, length(rich))
j = 1
 ggplot(scatterr, aes(x = sigma_in, y = sigma_out))+geom_point(aes(colour=sigma_bet,
                                                                   shape = as.factor(richness)))+
  scale_color_gradient(low="black", high="red")+
   scale_shape_manual(values=c(15,16,17,18,19,7,8,9,10))
pdf('../results/heatmap.pdf',  width = 5, height = 4)
for (i in rich){
  #Subset by richness
  data = subset(scatterr, richness == i)
  p = ggplot(data, aes(x = sigma_in, y = sigma_out))+geom_point(aes(colour=sigma_bet))+
    scale_color_gradient(low="black", high="red")+ 
    ggtitle(as.character(i))
  print(p)
  j = j+1
}
dev.off()



                       