#Quick plot
require(ggplot2)


scatterr = read.table('../data/scatterplot.csv', sep = ',', header = T, row.names = 1)
data = read.table('../data/similarity_fitness.csv', sep = ',',
                                header = T)
extinct_plot = read.table('../data/tile_plot.csv', sep = ',',
                                header = T)
subdata = extinct_plot[extinct_plot$richness == 4,]
ggplot(subdata, aes(x+1, y+1)) +
  geom_tile(aes(fill = as.factor(su)))+
  scale_fill_manual(values = c('cyan', 'black', 'magenta'))

 #Get rid of outliers
data = data[similarity_fitness$similarity>0,]
#Initialize similarity and p vectors
# similarity = c()
# fitness_dif = c()
# richness_p = c()
# data = similarity_fitness[(similarity_fitness$richness > 3)&
#                           (similarity_fitness$richness < 7),]
plot_data = data.frame(similarity = numeric(0), 
                       cohesion = numeric(0),
                       error = numeric(0),
                       richness = numeric(0))

#Average plot
h = hist(data$delP2, 
         breaks=seq(floor_dec(min(data$delP2, 1)),
                    ceiling_dec(max(data$delP2), 1), 
                    length.out = 100))
#Find the interval containing each element of delP2 in breaks.
ind_av = findInterval(data$delP2, h$breaks)
#Bin similar data together
P_av = rep(0, length(unique(ind_av)))
sim_av = rep(0, length(unique(ind_av)))
sim_std = rep(0, length(unique(ind_av)))
j = 1
for (i in (unique(ind_av))){
  #Average equivalent measurements by bin
  av_data = data[which(ind_av == i),]
  P_av[j] = median(av_data$delP2)
  sim_av[j] = median(av_data$similarity)
  sim_std[j] = sd(av_data$similarity)
  j = j+1
}
richness = rep(r, length(sim_av))
plot_data = rbind(plot_data, 
                  data.frame(sim_av, P_av, 
                             sim_std,richness))



# #Get counts with less than 10 points:
# trash = which(h$counts<20)
# #Find the interval containing each element of delP2 in breaks.
# pre_ind_av = findInterval(data$delP2, h$breaks)
# #Find indices of shit elements
# ind_shit = which(pre_ind_av %in% trash)
# ind_av = pre_ind_av[!pre_ind_av %in% trash]
# #Eliminate them from data
# data = data[-ind_shit,]
#Bin similar data together

#richness_p = c(richness_p, rep(j, length(sim_av)))

plot_data = plot_data[plot_data$richness == 5,]
plot_data = plot_data[complete.cases(plot_data),]
# plot_data$sim_std[plot_data$sim_av - plot_data$sim_std<0] = 
#   plot_data$sim_av[plot_data$sim_av - plot_data$sim_std<0]
# plot_data[plot_data$sim_av + plot_data$sim_std>1] = 1
ggplot(data = plot_data, aes(x = P_av, y = sim_av))+
  geom_line(aes(colour = as.factor(richness)), size = 1)+
  geom_ribbon(aes(ymin = sim_av-sim_std, ymax = sim_av+sim_std,
                  x = P_av, fill = 'band', alpha = 0.3))+
  theme(panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.ticks = element_line(),
        panel.border = element_rect(colour = "black", fill=NA, size = 1),
        legend.position = 'none',
        aspect.ratio = 1)+
  scale_color_manual("Richness",values="#F21A00")+
  scale_x_continuous(expand = c(-0.06, -0.05),
                     limits = c(-1.2,1.2),
                     breaks = c(-1, 0, 1),
                     labels = c('-1', '0', '1')) +
  scale_y_continuous(expand = c(0, -0.05),
                     limits = c(-0.15, 1.11),
                     breaks = c(0, 0.5, 1),
                     labels = c('0','0.5', '1'))+
  labs(x = expression(Delta~Theta), 
       y = expression(paste(italic(S),
                            '(',italic(C[P]),','
                            ,italic(C[R]),')')))+
  ggsave(filename = '../results/cohesiondiff_similarity.pdf', width = 4.3, height = 4)
# ggplot(data = sub_data, aes(x = delP2, y = similarity)) +
#   geom_point(aes(colour = richness))+
#   scale_color_gradient(low="black", high="red")
# rich = unique(scatterr$richness)
# plots = rep(NA, length(rich))
# j = 1
#  ggplot(scatterr, aes(x = sigma_in, y = sigma_out))+geom_point(aes(colour=sigma_bet,
#                                                                    shape = as.factor(richness)))+
#   scale_color_gradient(low="black", high="red")+
#    scale_shape_manual(values=c(15,16,17,18,19,7,8,9,10))
# pdf('../results/heatmap.pdf',  width = 5, height = 4)
# for (i in rich){
#   #Subset by richness
#   data = subset(scatterr, richness == i)
#   p = ggplot(data, aes(x = sigma_in, y = sigma_out))+geom_point(aes(colour=sigma_bet))+
#     scale_color_gradient(low="black", high="red")+ 
#     ggtitle(as.character(i))
#   print(p)
#   j = j+1
# }
# dev.off()
# 
#                        