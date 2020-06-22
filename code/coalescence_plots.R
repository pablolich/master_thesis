#Quick plot
require(ggplot2)
scatterr = read.table('../data/scatterplot.csv', sep = ',', header = T, row.names = 1)
similarity_fitness = read.table('../data/similarity_fitness.csv', sep = ',',
                                header = T)
ggplot(data = similarity_fitness, aes(x = delF, y = similarity)) +
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



                       