rm(list=ls())

library('plyr')
library('ggplot2')


DA = read.table('alphac.csv', sep='\t', header=T)

xy1 = ggplot(DA, aes(x=alphac)) +
  geom_line(aes(y=max, 
                colour='Max'),
            size=1.5) +
  geom_line(aes(y=mean,
                colour='Mean'),
            size=1.5) +
  theme(legend.position='right') +
  xlab('Alpha_C') +
  ylab('Fitness') +
  ggtitle('Mean and max fitness vs. Alpha_C') +
  scale_color_discrete(name='Legend')

pdf('plt_meanmax.pdf', width=8, height=4)
xy1
dev.off()
