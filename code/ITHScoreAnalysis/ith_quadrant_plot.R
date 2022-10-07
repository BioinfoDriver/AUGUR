
# load ith score
load(file='/data/intra.inter.ith.score.RData')

QuadrantPlot <- function(ith.score, file.name){
 
 ##  input: rna-ith scores
 colnames(ith.score) <- c('intra.score', 'inter.score')
 
 within.ith <- ifelse(ith.score$intra.score > mean(ith.score$intra.score), 'top', 'bottom')
 between.ith <- ifelse(ith.score$inter.score > mean(ith.score$inter.score), 'right', 'left')
 ith.score$quadrant <- paste(within.ith, between.ith, sep='_')
 
  ## Make new quadrant type figure - 2 density plots, one with four to scale circles
  a <- table(factor(ith.score[,'quadrant'],levels=c("top_left", "bottom_left", "top_right", "bottom_right")))
  b <- a/min(a)
  
  pdf(paste0(file.name, '.pdf'))
  # scale circles
  plot(c(100,100,400,400),c(400,100,400,100), pch=16, cex=b*2, ylim=c(0,500),xlim=c(0,500), 
   col=c("firebrick1", "darkorchid2", "gold1", "turquoise2"),axes=F,ylab='',xlab='')
   
  abline(h=250,v=250,col=1,lty=2,lwd=4)
  text(c(100,100,400,400),c(450,200,450,200), paste(a,'genes'))
  box()
  
  # density plots
  hist(ith.score[,'inter.score'], breaks=100, prob=TRUE,col='#6baed699',axes=F,ylab='', xlab='',main='Between tumour, LIHC')
  lines(density(ith.score[,'inter.score']), lwd = 3,col = '#08519c')
  abline(v=mean(ith.score[,'inter.score']), col=1,lwd=4, lty=2)
  
  
  hist(ith.score[,'intra.score'], breaks=100, prob=TRUE,col='#fd8d3c99',axes=F,ylab='', xlab='',main='Within tumour, LIHC')
  lines(density(ith.score[,'intra.score']), lwd = 3,col = '#a63603')
  abline(v=mean(ith.score[,'intra.score']), col=1,lwd=4, lty=2)
  
  dev.off()
}


setwd('/result/Section2/')
QuadrantPlot(losic.ith.score[, c('losic.intra.sd', 'losic.inter.sd')], 'losic_quadrant_plot_sd')



