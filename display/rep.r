rm(list=ls())
library('ggplot2')
library('ggtern')
library('plotly')
library('colorspace')
library('plotrix')
library('plyr')

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
D = D[ D$kept == "True", ] # Eliminate not kept sequences

fitness = D$fitness
time = D$time
event = D$event
kept = D$kept
gene_ratio = D$gene_ratio
plasmid_size = D$plasmid_size
up_down_ratio = D$up_down_ratio
mean_space = D$mean_space
repetition = D$repetition

pdf('evol.pdf', width=8, height=6)
pA =  ggplot(D) + 
  geom_line(aes(x=time, y=fitness), 
            size=1, 
            color = rainbow_hcl(5)[repetition+1] ) +
  facet_grid(repetition ~ . )
pA
dev.off()
