rm(list=ls())

library('ggplot2')
library('colorspace')
library('plyr')
library('rmarkdown')

args <- commandArgs(trailingOnly = TRUE)

path_to_files = args[1]
history_file = paste(path_to_files, 'history.csv', sep='')
plasmid_file = paste(path_to_files, 'plasmid.csv', sep='')

#--- Convenience function ------------------------------------------------------

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

# Summary

DATA = read.table(history_file, sep='\t', header=T)
#DATA = read.table('history.csv', sep='\t', header=T)

D = DATA[ DATA$kept == "True", ] # Eliminate not kept sequences

fitness = D$fitness
time = D$time
event = D$event
kept = D$kept
gene_ratio = D$gene_ratio
plasmid_size = D$plasmid_size
up_down_ratio = D$up_down_ratio
mean_space = D$mean_space
repetition = D$repetition

means = aggregate(fitness ~ repetition, D, mean);
maxs = aggregate(fitness ~ repetition, D, max);
vars = aggregate(fitness ~ repetition, D, var);

GMEAN = mean(means$fitness);
GMAX = max(maxs$fitness);

D = DATA[DATA$event != 'Beg',]

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank())

rep = ggplot(data=D,
             aes(x='',
                 fill=event)) +
  stat_count(width=1) +
  coord_polar('y') +
  blank_theme

pdf(paste(path_to_files, sep='', 'plt_events.pdf'), width=5, height=5)
rep
dev.off()

# Fitness evolution

D = DATA[ DATA$kept == "True", ] # Eliminate not kept sequences

fitness = D$fitness
time = D$time
event = D$event
kept = D$kept
gene_ratio = D$gene_ratio
plasmid_size = D$plasmid_size
up_down_ratio = D$up_down_ratio
mean_space = D$mean_space
repetition = D$repetition

plot_evol =  ggplot(D) + 
  geom_line(aes(x=time, y=fitness), 
            size=1, 
            color = rainbow_hcl(5)[repetition+1] ) +
  facet_grid(repetition ~ . ) +
  ggtitle('Fitness evolution') +
  xlab('Time') +
  ylab(time)

pdf(paste(path_to_files, sep='', 'plt_evol.pdf'), width=8, height=6)
plot_evol
dev.off()


plot_means = ggplot(D, aes(repetition, fitness, group=repetition)) + 
  geom_violin(aes(fill=factor(repetition)), draw_quantiles = 0.5) +
  ggtitle('Fitness repartition') +
  xlab('Repetition') +
  ylab('Fitness') +
  theme(legend.position='none')

pdf(paste(path_to_files, sep='', 'plt_means.pdf'), width=7, height=6)
plot_means
dev.off()

# Fitness plot

D = DATA[ DATA$kept == "True" & DATA$repetition == 0, ]

fitness = D$fitness
time = D$time
event = D$event
kept = D$kept
gene_ratio = D$gene_ratio
plasmid_size = D$plasmid_size
up_down_ratio = D$up_down_ratio
mean_space = D$mean_space
repetition = D$repetition

seg = data.frame(x0=time[1:(length(time)-1)],
                 y0=fitness[1:(length(fitness)-1)],
                 x1=time[2:length(time)],
                 y1=fitness[2:length(fitness)],
                 events=event[2:length(event)])

plot_fit = ggplot() + 
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

pdf(paste(path_to_files, sep='', 'plt_fitness.pdf'), width=10, height=6)
plot_fit
dev.off()

# Statistics

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

pdf(paste(path_to_files, sep='', 'plt_stats.pdf'), width=10, height=6)
multiplot(pA,pB,pC,pD,cols=2)
dev.off()

# Stats and events

Devents = event[2:length(event)]
Dfitness = fitness[2:length(fitness)] - fitness[1:length(fitness)-1]
Dgratio = gene_ratio[2:length(gene_ratio)] - gene_ratio[1:length(gene_ratio)-1]
Dsize = plasmid_size[2:length(plasmid_size)] - plasmid_size[1:length(plasmid_size)-1]
Dudratio = up_down_ratio[2:length(up_down_ratio)] - up_down_ratio[1:length(up_down_ratio)-1]
Dmspace = mean_space[2:length(mean_space)] - mean_space[1:length(mean_space)-1]

dfd = data.frame(Devents, Dfitness, Dsize, Dgratio, Dudratio, Dmspace)
colnames(dfd) = c('event', 'fitness', 'size', 'gratio', 'udratio', 'mspace')

L = length(dfd$event)
L2 = as.integer(0.33*L)
L3 = as.integer(0.66*L)

R = cbind(c(rep('0%-33%',L2),
            rep('34%-66%',L3-L2),
            rep('67%-100%',L-L3)))
dfd['type'] = as.factor(R)

perevent_start = ggplot(dfd, aes(x=event, y=fitness, fill=event)) +
  geom_violin(draw_quantiles=c(0.5)) +
  facet_wrap(~type) +
  xlab('Event') +
  ylab('Effet en fitness') +
  theme(legend.position='none') +
  ggtitle('Events effet on fitness at different times')

pdf(paste(path_to_files, sep='', 'plt_stats_var.pdf'), width=9, height=5)
perevent_start
dev.off()

#--- OTHERS --------------------------------------------------------------------

D = DATA[ DATA$kept == "True", ] # Eliminate not kept sequences

DATA = DATA[DATA$event != 'Beg',]
Dkept = DATA[ DATA$kept == "True", ]
Dnot = DATA[ DATA$kept == "False", ]
Dall = DATA

cts_all = count(Dall, 'event')
cts_all$freq = cts_all$freq/sum(cts_all$freq)

cts_kept = count(Dkept, 'event')
cts_kept$freq = cts_kept$freq/sum(cts_kept$freq)

pall = ggplot(cts_all, aes(x='', weight=freq, fill=event)) +
  geom_bar(width=0.5) +
  xlab('All mutations') +
  ylab('%') +
  theme(legend.position='none')

pkept = ggplot(cts_kept, aes(x='', weight=freq, fill=event)) +
  geom_bar(width=0.5) +
  xlab('Kept mutations') +
  ylab('') +
  theme(legend.position='none')

pdf(paste(path_to_files, sep='','plt_diffevt.pdf'), width=5, height=8)
multiplot(pall, pkept, cols=2)
dev.off()

# Additional stats

D = DATA[DATA$event != 'Beg',]

dat2 = ddply(D, 
             c('repetition'), 
             function(df) c(mean(df$fitness),
                            sd(df$fitness),
                            max(df$fitness)))
colnames(dat2) = c('rep', 'mean', 'sd', 'max')

MEAN = mean(dat2$mean)
SD = sd(dat2$mean)
MAX = mean(dat2$max)

# Plasmid plot

PLAS = read.table(plasmid_file, header=T, sep='\t')

P = PLAS[PLAS$repetition == 0,]
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

pdf(paste(path_to_files, sep='', 'plt_plasmid.pdf'), width=6, height=6)
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

text(0, 0, paste('time:', last_sim_time, '\n',
                 'size:', p_size, '\n'
), offset = 0.5,  cex = 1, col = 1)

dev.off()

#--- RMD FILE ------------------------------------------------------------------

render('display/exporter.rmd', 
       output_format='pdf_document',
       output_file='summary.pdf',
       output_dir=path_to_files,
       quiet=T)
