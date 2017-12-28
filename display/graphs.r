rm(list=ls())
library('ggplot2')
library('ggtern')
library('plotly')

D = read.table('history.csv', sep='\t', header=T)

fitness = D$fitness
time = D$time
event = D$event
kept = D$kept
gene_ratio = D$gene_ratio
plasmid_size = D$plasmid_size
up_down_ratio = D$up_down_ratio
mean_space = D$mean_space

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

#-----------

pA = ggplot() + geom_line(aes(x=time, y=plasmid_size)) +
  geom_point(aes(x=time, y=plasmid_size, colour=event, shape=kept), size=3)
pB = ggplot() + geom_line(aes(x=time, y=up_down_ratio)) +
  geom_point(aes(x=time, y=up_down_ratio, colour=event, shape=kept), size=3)
pC = ggplot() + geom_line(aes(x=time, y=mean_space)) +
  geom_point(aes(x=time, y=mean_space, colour=event, shape=kept), size=3)
pD = ggplot() + geom_line(aes(x=time, y=gene_ratio)) +
  geom_point(aes(x=time, y=gene_ratio, colour=event, shape=kept), size=3)

multiplot(pA,pB,pC,pD,cols=2)


ggplot() + geom_line(aes(x=time, y=plasmid_size)) +
  geom_point(aes(x=time, y=plasmid_size, colour=event, shape=kept), size=3)








#--------------------


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