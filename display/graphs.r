rm(list=ls())
library('ggplot2')
library('ggtern')
library('plotly')
library('colorspace')
library('plotrix')
library('plyr')
# Sys.setenv("plotly_username"="blac")
# Sys.setenv("plotly_api_key"="25rk3TlEakY37FkC5F2p")

#--------- convenience function ------------------------------------------------

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#--------- convenience function ------------------------------------------------

D = read.table('history.csv', sep='\t', header=T)
D = D[ D$kept == "True" & D$repetition == 0, ] # Eliminate not kept sequences

fitness = D$fitness
time = D$time
event = D$event
kept = D$kept
gene_ratio = D$gene_ratio
plasmid_size = D$plasmid_size
up_down_ratio = D$up_down_ratio
mean_space = D$mean_space
repetition = D$repetition

# Fitness plot
seg = data.frame(x0=time[1:(length(time)-1)],
                 y0=fitness[1:(length(fitness)-1)],
                 x1=time[2:length(time)],
                 y1=fitness[2:length(fitness)],
                 events=event[2:length(event)])

plotA = ggplot() + 
  geom_segment(aes(x=x0,
                   y=y0,
                   xend=x1,
                   yend=y1,
                   colour=events),
               data=seg,
               size=1.5) +
  geom_point(aes(x=time, y=fitness), size=1) + 
  xlab('Time') + 
  ylab('Fitness') + 
  ggtitle('Fitness evolution and events') + 
  theme(plot.title = element_text(lineheight=1, face="bold"))

pdf('fitness.pdf', width=10, height=6)
plotA
dev.off()

# Multiplot macro-statistics

seg2 = data.frame(x0=time[1:(length(time)-1)],
                 y0=plasmid_size[1:(length(plasmid_size)-1)],
                 x1=time[2:length(time)],
                 y1=plasmid_size[2:length(plasmid_size)],
                 events=event[2:length(event)])

pA =  ggplot() + 
      geom_segment(aes(x=x0,
                       y=y0,
                       xend=x1,
                       yend=y1,
                       colour=events),
                   data=seg2,
                   size=1) +
      geom_point(aes(x=time, y=plasmid_size), size=0.5) +
      xlab('Time') + 
      ylab('Plasmid size') + 
      ggtitle('Plasmid size evolution') + 
      theme(plot.title = element_text(lineheight=1, face="bold"), 
            legend.position='none')

seg3 = data.frame(x0=time[1:(length(time)-1)],
                  y0=up_down_ratio[1:(length(up_down_ratio)-1)],
                  x1=time[2:length(time)],
                  y1=up_down_ratio[2:length(up_down_ratio)],
                  events=event[2:length(event)])

pB =  ggplot() + 
  geom_segment(aes(x=x0,
                   y=y0,
                   xend=x1,
                   yend=y1,
                   colour=events),
               data=seg3,
               size=1) +
  geom_point(aes(x=time, y=up_down_ratio), size=0.5) +
  xlab('Time') + 
  ylab('Orientation ratio (+/-)') + 
  ggtitle('Orientation ratio evolution') + 
  theme(plot.title = element_text(lineheight=1, face="bold"), 
        legend.position='none')

seg4 = data.frame(x0=time[1:(length(time)-1)],
                  y0=mean_space[1:(length(mean_space)-1)],
                  x1=time[2:length(time)],
                  y1=mean_space[2:length(mean_space)],
                  events=event[2:length(event)])

pC =  ggplot() + 
  geom_segment(aes(x=x0,
                   y=y0,
                   xend=x1,
                   yend=y1,
                   colour=events),
               data=seg4,
               size=1) +
  geom_point(aes(x=time, y=mean_space), size=0.5) +
  xlab('Time') + 
  ylab('Mean space between genes') + 
  ggtitle('Mean space between genes evolution') + 
  theme(plot.title = element_text(lineheight=1, face="bold"), 
        legend.position='none')  

seg5 = data.frame(x0=time[1:(length(time)-1)],
                  y0=gene_ratio[1:(length(gene_ratio)-1)],
                  x1=time[2:length(time)],
                  y1=gene_ratio[2:length(gene_ratio)],
                  events=event[2:length(event)])

pD =  ggplot() + 
  geom_segment(aes(x=x0,
                   y=y0,
                   xend=x1,
                   yend=y1,
                   colour=events),
               data=seg5,
               size=1) +
  geom_point(aes(x=time, y=gene_ratio), size=0.5) + 
  xlab('Time') + 
  ylab('Gene proportion in plasmid') + 
  ggtitle('Gene proportion evolution') + 
  theme(plot.title = element_text(lineheight=1, face="bold"), 
        legend.position='none')  


pdf('stats.pdf', width=10, height=6)
multiplot(pA,pB,pC,pD,cols=2)
dev.off()

# Macro-stats to events

fitness = D$fitness
time = D$time
event = D$event
kept = D$kept
gene_ratio = D$gene_ratio
plasmid_size = D$plasmid_size
up_down_ratio = D$up_down_ratio
mean_space = D$mean_space
repetition = D$repetition

Devents = event[2:length(event)]
Dfitness = fitness[2:length(fitness)] - fitness[1:length(fitness)-1]
Dgratio = gene_ratio[2:length(gene_ratio)] - gene_ratio[1:length(gene_ratio)-1]
Dsize = plasmid_size[2:length(plasmid_size)] - plasmid_size[1:length(plasmid_size)-1]
Dudratio = up_down_ratio[2:length(up_down_ratio)] - up_down_ratio[1:length(up_down_ratio)-1]
Dmspace = mean_space[2:length(mean_space)] - mean_space[1:length(mean_space)-1]

dfd = data.frame(Devents, Dfitness, Dsize, Dgratio, Dudratio, Dmspace)
colnames(dfd) = c('event', 'fitness', 'size', 'gratio', 'udratio', 'mspace')

Pfit = ggplot(dfd) +
  geom_histogram(aes(fitness, fill=event), bins=25, color=1) +
  xlab('Fitness') + 
  ylab('Count') +
  theme(legend.position='none')  

Pgratio = ggplot(dfd) +
  geom_histogram(aes(gratio, fill=event), bins=25, color=1) +
  xlab('Gene ratio') + 
  ylab('Count') +
  theme(legend.position='none')  

Psize = ggplot(dfd) +
  geom_histogram(aes(size, fill=event), bins=25, color=1) +
  xlab('Plasmid size') + 
  ylab('Count') +
  theme(legend.position='none') 

Pudratio = ggplot(dfd) +
  geom_histogram(aes(udratio, fill=event), bins=25, color=1) + 
  xlab('+/- ratio') + 
  ylab('Count') +
  theme(legend.position='none') 

Pmspace = ggplot(dfd) +
  geom_histogram(aes(mspace, fill=event), bins=25, color=1) + 
  xlab('Mean space between genes') + 
  ylab('Count') +
  theme(legend.position='none') 

pdf('Dstats.pdf', width=6, height=6)
multiplot(Psize, Pudratio, Pmspace, Pgratio,cols=2)
dev.off()

#--------------------------PLASMID PLOT ---------------------------

P = read.table('plasmid.csv', header=T, sep='\t')
P = P[P$repetition == 0,]
last_sim_time = max(P$time)
sub_P = P[P$repetition== 0 & P$time==last_sim_time,]
sub_P = sub_P[order(sub_P$location),]

s1 = sub_P$location - 1
s2 = sub_P$location + sub_P$length - 1
ty = sub_P$type
or = sub_P$strand
co = paste(sub_P$type, sub_P$strand)

dfp = data.frame(s1, s2, ty, or, co)
colnames(dfp) = row.names=c('start', 'stop', 'type', 'orient', 'color')
p_size = max(sapply(dfp[,c('start', 'stop')], max))
dfp$start = 2*pi*(1-dfp$start/p_size)
dfp$stop  = 2*pi*(1-dfp$stop/p_size)

offset = 90/360*pi*2

pdf('plasmid.pdf', width=6, height=6)
plot(0, 0, 
     col= NA, 
     xlim=c(-1,1), 
     ylim=c(-1, 1), 
     main=NA,
     xlab=NA, 
     ylab=NA,
     asp=1,
     axes=F)
for(i in 1:length(dfp$start)){
  
  xs = cos(offset + seq(dfp$stop[i], dfp$start[i], 0.01))
  ys = sin(offset + seq(dfp$stop[i], dfp$start[i], 0.01))
  
  if(dfp$color[i] == 'G 1') acol = rgb(0, 0.4, 1.0)
  else if (dfp$color[i] == 'G -1') acol = rgb(0, 1.0, 0.4)
  else if (dfp$color[i] == 'V 0') acol = rgb(0, 0, 0)
  else acol = rgb(0.5, 0.5, 0.5)
  
  if( dfp$type[i] != 'P'){
    lines(xs, ys , col = acol , lwd=10, lend=1, ljoin=0)
  }
  else{
    points(cos(offset + dfp$stop[i]), 
           sin(offset + dfp$start[i]), 
           pch=20, 
           col=2, 
           cex=3)
  }
    
}

legend(-0.25, 0.6, 
       c('G+', 'G-', 'barrier'), 
       fill=c(rgb(0,0.4,1.0),
              rgb(0,1.0,0.4),
              NA),
       col = c(NA, NA, 2), pch = c(NA ,NA, 20), border=c(1,1,NA),
       cex=0.9)

text(0, 0, paste('time:', last_sim_time, '\n',
                 'size:', p_size, '\n'
                 ), offset = 0.5,  cex = 1, col = 1)

dev.off()


