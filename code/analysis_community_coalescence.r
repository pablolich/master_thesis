#setwd("~/Desktop/master_thesis/code")
require(ggplot2)
require(cowplot)
rm(list=ls())
source('coal_analysis_functions.R')


data = read.table('../data/similarity_fitness.csv', sep = ',',
                  header = T)
extinct_plot = read.table('../data/tile_plot.csv', sep = ',',
                          header = T)

all_data = read.table('../data/community_data.csv', sep = ',',
                      header = T)

data_hist = read.table('../data/similarity_fitness copy.csv', sep = ',',
                       header = T)

#Get facilitation and competition indices each community
facilitation = aggregate(all_data[,16], list(all_data$n_simulation), mean)
competition = aggregate(all_data[,17], list(all_data$n_simulation), mean)
#Assign each community to a group based on their comp and fac levels.
colors = distance(competition$x, facilitation$x, 1.4, 0.2)
data_community = data.frame(competition$x, facilitation$x, colors)

#Plot all communities in a competition facilitation scatter plot
##############################################################################
comp_fac_scatter = ggplot(data = data_community, 
                          aes(x = competition$x,
                              y = facilitation$x))+
  geom_point(aes(colour = as.factor(colors), 
                 size = as.factor(colors)))+
  scale_color_manual(values=c("#F21A00", "#EBCC2A", "#3B9AB2"))+
  scale_size_manual(values = c(2,0.5,2))+
  theme(panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.ticks = element_line(),
        legend.position = 'none',
        panel.border = element_rect(colour = "black", fill=NA, size = 1),
        legend.key.size = unit(10, 'mm'),
        aspect.ratio = 1)+
  scale_x_continuous(expand = c(0,0),
                     limits = c(2,6)) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(2, 5))+
  geom_abline(intercept = 0, slope = 1,
              size = 2, linetype = 'dashed')+
  labs(x = expression(paste('Competition ', symbol("\341"), italic(C), symbol("\361"))),
       y = expression(paste('Facilitation ', symbol("\341"), italic(F), symbol("\361"))))+
  annotate('text', x = 2.7, y = 4.5, 
           label = " Theta>0 ",parse = TRUE,size=10) +
  annotate('text', x = 5.2, y = 2.4, 
           label = " Theta<0 ",parse = TRUE,size=10) +
  ggsave(filename = '../results/comp_fac_scatter.pdf', height = 4.3, width = 4)
##############################################################################

#Plot tile plot with extinctions and provenance
extinction_provenance = ggplot(extinct_plot, aes(x+1, y+1)) +
  geom_tile(aes(fill = as.factor(su)), colour = NA)+
  theme(panel.background = element_blank(),
        legend.position = 'none',
        panel.grid.major = element_blank(), #Get rid of grids
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(angle = 0, vjust = 0.5,
                                    size = 20,
                                    margin = margin(r = -20)),
        axis.title.x = element_text(size = 20),
        axis.ticks.y = element_blank())+
  scale_x_continuous(expand = c(0,0),
                     breaks = c(1,100,200,300, 400),
                     labels = c('1','','200','', '400')) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = c(1,8),
                     labels = c('low','high'))+
  scale_fill_manual(values = c("#F21A00", "#EBCC2A", "#3B9AB2"))+
  labs(x = 'Coalescence instances (ordered by dominant color)',
       y = expression(Theta))+
  ggsave(filename = '../results/extinctions_provenance.pdf', width = 8, height = 2)

#####################################################################################

#Get rid of outliers
r = 5
data = data[data$similarity>0,]
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

plot_data = plot_data[complete.cases(plot_data),]

cohesion_similarity = ggplot(data = plot_data, aes(x = P_av, y = sim_av))+
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
  scale_color_manual("Richness",values="red")+
  scale_x_continuous(expand = c(-0.06, -0.05),
                     limits = c(-1.2,1.2),
                     breaks = c(-1, 0, 1),
                     labels = c('-1', '0', '1')) +
  scale_y_continuous(expand = c(0, -0.05),
                     limits = c(-0.15, 1.11),
                     breaks = c(0, 0.5, 1),
                     labels = c('0','0.5', '1'))+
  labs(x = expression(paste('Fitness difference ', Delta~Theta)), 
       y = expression(paste('Similarity ', italic(S),
                            '(',italic(C[P]),','
                            ,italic(C[R]),')')))+
  ggsave(filename = '../results/cohesiondiff_similarity.pdf', width = 4.3, height = 4)
#########################################################################################
data_hist = data_hist[(data_hist$richness==4 & data_hist$similarity>0),]
num_elem = length(data_hist$similarity)

 similarity_hist = ggplot(data_hist, aes(x = similarity)) + 
  geom_histogram(breaks = seq(0, 1, 0.1),
                 color = '#6BB1C1', 
                 fill = '#3B9AB2')+
  theme_minimal()+
  theme(panel.background= element_rect(),
        plot.title = element_text(hjust = 0.5,
                                  margin = margin(b = 0)),
        panel.grid = element_blank(),
        panel.border = element_rect(color = NA, 
                                    fill = NA),
        axis.text.x = element_text(size = 12,
                                 margin = margin(b = 0)),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 0.7,
        rect = element_rect(fill = "transparent"))+
  scale_x_continuous(breaks = c(0.05, 0.5, 0.95),
                     labels = c('0', '0.5', '1'))+
  scale_y_continuous(labels = c())+
  labs(title = expression(paste(italic(S),
                            '(',italic(C[P]),','
                            ,italic(C[R]),')')),
       y = '')+
  ggsave(filename = '../results/histogram_similarity.pdf', width = 4.3, height = 4)


tot_coalescence = ggdraw(ylim = c(0,0.87), xlim = c(0,1)) +
    draw_plot(comp_fac_scatter, x = 0, y = .31, width = .5, height = .5) +
    draw_plot(cohesion_similarity, x = .5, y = .31, width = .5, height = .5) +
    draw_plot(similarity_hist, x = 0.81, y = 0.358, width = 0.17, height = 0.2)+
    draw_plot(extinction_provenance, x = 0.02, y = 0, width = 0.97, height = 0.3) +
    draw_plot_label(label = c("A", "B", "C"), size = 25,
                    x = c(0.05, 0.5, 0.05), y = c(0.87, .87, 0.35)) + 
    ggsave(filename = '../results/community_coalescence_plots.pdf', height = 8, width = 10 )

  