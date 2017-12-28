rm(list=ls())
library('ggplot2')
#setwd("~/Documents/INSA/BC/Projet/Evolution-Simulation")

D = read.table('history.csv', sep='\t', header=T)

fitness = D$fitness
time = D$time
event = D$event
kept = D$kept


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
  geom_point(aes(x=time, y=fitness, shape=factor(kept)), size=2) + 
  xlab('Time') + 
  ylab('Fitness') + 
  ggtitle('Fitness evolution and events') + 
  theme(plot.title = element_text(lineheight=1, face="bold"))

pdf('graph.pdf', width=10, height=6)
plotA
dev.off()
