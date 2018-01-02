rm(list=ls())
library('ggplot2')
library('ggtern')
library('plotly')
library('colorspace')
library('plotrix')
library('plyr')
Sys.setenv("plotly_username"="blac")
Sys.setenv("plotly_api_key"="25rk3TlEakY37FkC5F2p")

D = read.table('history.csv', sep='\t', header=T)
D = D[ D$kept == "True", ]

fitness = D$fitness
time = D$time
event = D$event
kept = D$kept
gene_ratio = D$gene_ratio
plasmid_size = D$plasmid_size
up_down_ratio = D$up_down_ratio
mean_space = D$mean_space
repetition = D$repetition


p = plot_ly(D, x = ~time) %>%
  add_trace(y = ~fitness, name = 'fitness', mode = 'lines', type='scatter') %>%
  add_trace(y = ~fitness, color= ~event, mode='markers', type='scatter')
p

api_create(p, filename = "midwest-boxplots")



ggplot() + geom_line(aes(x=time, y=fitness, colour=factor(repetition))) + geom_point(aes(x=time, y=fitness, shape=event))

seg = data.frame(x0=time[1:(length(time)-1)],
                 y0=fitness[1:(length(fitness)-1)],
                 x1=time[2:length(time)],
                 y1=fitness[2:length(fitness)],
                 events=event[2:length(event)])


p = plot_ly(seg) %>%
  add_segments(
    data = seg,
    x = ~x0, xend = ~x1,
    y = ~y0, yend = ~y1, #split=~id,
    alpha = 1, size = I(3),
    color= ~events
  ) %>%
  add_markers(
    data=D,
    x = ~time, y = ~fitness,
    name="Fitness"
    )

p

plotA = ggplot() + 
  geom_segment(aes(x=x0,
                   y=y0,
                   xend=x1,
                   yend=y1,
                   colour=events),
               data=seg,
               size=1.5) +
  geom_point(aes(x=time, y=fitness), size=2) + 
  xlab('Time') + 
  ylab('Fitness') + 
  ggtitle('Fitness evolution and events') + 
  theme(plot.title = element_text(lineheight=1, face="bold"))

#pdf('graph.pdf', width=10, height=6)
plotA
#dev.off()

#-----------

pA = ggplot() + geom_line(aes(x=time, y=plasmid_size)) +
  geom_point(aes(x=time, y=plasmid_size, colour=event), size=2)
pB = ggplot() + geom_line(aes(x=time, y=up_down_ratio)) +
  geom_point(aes(x=time, y=up_down_ratio, colour=event), size=2)
pC = ggplot() + geom_line(aes(x=time, y=mean_space)) +
  geom_point(aes(x=time, y=mean_space, colour=event), size=2)
pD = ggplot() + geom_line(aes(x=time, y=gene_ratio)) +
  geom_point(aes(x=time, y=gene_ratio, colour=event), size=2)

multiplot(pA,pB,pC,pD,cols=2)




ggplot() + geom_line(aes(x=time, y=plasmid_size)) +
  geom_point(aes(x=time, y=plasmid_size, colour=event, shape=kept), size=3)


#--------------------------PLASMID---------------------------

P = read.table('plasmid.csv', header=T, sep='\t')
sub_P = P[P$repetition== 0 & P$time== 238,]
sub_P = sub_P[order(sub_P$location),]

palette_t = c(G='blue', V='red', P='grey')
palette_s = c('blue', 'red', 'grey', 'yellow')
cols = palette_s[sub_P$strand+1]
cols= palette_t[sub_P$type]

p = plot_ly(sub_P) %>%
  add_trace(hole=0.9,
            sort=F,
            labels=~location,
            type='pie',
            values=~length,
            marker=list(
              colors=cols),
            hoverinfo='text',
            hovertext=~length,
            textinfo='text',
            text=~strand,
            insidetextfont=list(color='white'))

p

#--------------------

radius = rep(1, length(sub_P))

x = sub_P$location/max(sub_P$location+sub_P$length)*360
x1 = (sub_P$location + sub_P$length) /max(sub_P$length)*360


r = rep(10, length(x))

polar.plot(lengths = r, polar.pos = x, rp.type='s', 
           clockwise=T, point.symbols = 20, show.grid.labels = F,
           cex=1, point.col = sub_P$type, radial.lim=c(0,11), start=0)

polar.plot(lengths = r, polar.pos = x + x1, rp.type='s', start=0,
           clockwise=T, point.symbols = 21, show.grid.labels = F,
           cex=1, point.col = sub_P$type, radial.lim=c(0,11), add=T)



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
