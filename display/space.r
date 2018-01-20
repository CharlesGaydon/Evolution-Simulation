rm(list=ls())

library('ggplot2')
library('ggtern')
library('colorspace')

DATA = read.table('space.csv', sep='\t', header=T)

# to avoid errors in interpolation and display issues
DATA$pinv[DATA$pinv == 0] = 0.02
DATA$pdel[DATA$pdel == 0] = 0.02
DATA$pins[DATA$pins == 0] = 0.02

# normalizing values
DATA$max = (DATA$max - min(DATA$max))
DATA$max = (DATA$max)/(max(DATA$max))

DATA$mean = DATA$mean - min(DATA$mean)
DATA$mean = DATA$mean/max(DATA$mean)

DATA$sd = (DATA$sd-min(DATA$sd))/(max(DATA$sd)-min(DATA$sd))

plt1 = ggtern(DATA,
             aes(pinv,
                 pdel,
                 pins,
                 value=max)) + 
        stat_interpolate_tern(geom="polygon",
                              formula=value ~ x+y,
                              method='loess',
                              aes(fill=..level..),
                              show.legend=F,
                              bins=40) +
        geom_point(aes(color=max),
                   size=3) +
        theme_latex() +
        theme_arrowlarge() +
        Tlab('$P_{del}$') +
        Llab('$P_{inv}$') +
        Rlab('$P_{ins}$') +
        scale_fill_gradientn(colours=diverge_hcl(40, h=c(255, 330)), 
                             name='Max value (%)')
        

pdf('triangle_l1.pdf', width=6, height=6)
plt1
dev.off()


plt2 = ggtern(DATA,
              aes(pinv,
                  pdel,
                  pins,
                  value=mean)) + 
  stat_interpolate_tern(geom="polygon",
                        formula=value~x+y,
                        method='glm',
                        aes(fill=..level..),
                        show.legend=F,
                        bins=30) +
  geom_point(aes(color=mean),
             size=3) +
  theme_latex() +
  theme_arrowlarge() +
  Tlab('$P_{del}$') +
  Llab('$P_{inv}$') +
  Rlab('$P_{ins}$') +
  scale_color_continuous(name='Max value (%)')

pdf('triangle_l2.pdf', width=6, height=6)
plt2
dev.off()


plt3 = ggtern(DATA,
              aes(pinv,
                  pdel,
                  pins,
                  value=sd)) + 
  stat_interpolate_tern(geom="polygon",
                        formula=value ~ x+y,
                        aes(fill=..level..),
                        show.legend=F) +
  geom_point(aes(color=sd),
             size=3) +
  theme_latex() +
  theme_arrowlarge() +
  Tlab('$P_{del}$') +
  Llab('$P_{inv}$') +
  Rlab('$P_{ins}$') +
  scale_color_continuous(name='Max value (%)')

pdf('triangle_l3.pdf', width=6, height=6)
plt3
dev.off()



