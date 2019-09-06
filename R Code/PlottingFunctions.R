library(lattice)

plot.IterationMatrix <- function(IM, L, U, n=5, cex=1, ...){
  
  # absolutely no idea how to get the axis text to change size, cex.axis 
  # doesn't work:
  xlist=list(at = seq(1, ncol(IM), l = n), labels = seq(L, U, l = n), cex=cex)
  ylist=list(at = seq(1, ncol(IM), l = n), labels = seq(U, L, l = n), cex=cex)
  IM %>% `class<-`('matrix') %>% rotate %>%
    levelplot(scale=list(x=xlist, y=ylist),
              xlab = 'Size at previous census / size of parent',
              ylab='Size at current census', main='IPM Kernel', ...)
}

rotate <- function(M) M %>% apply(2, rev) %>% t
normalise <- function(x) x/sum(x)

plotSizeCounts <- function(dataMatrix, palette="Blues", breaks = NULL){
  # purpose : Produces a stacked area plot for size count distributions over
  #           time.
  # inputs  : dataMatrix - The matrix which contains the counts. Each row is a 
  #                        size class, and each column is a census.
  #           palette    - The character name of the brewer colour palette which
  #                        should be used to create the plot.
  
  require(ggplot2)
  
  # To avoid unecessary typing:
  nC <- ncol(dataMatrix) ; nR <- nrow(dataMatrix)
  
  # get the names from the breaks :
  if (is.null(breaks)) names <- as.character(1:nC)
  else names <- cut(breaks[-1] - diff(breaks)/2, breaks) %>% as.character
  
  # reorganise the data for ggplot:
  DF <- data.frame(sizeClass = factor(rep(1:nR, nC)),
                            census = rep(1:nC, rep(nR, nC)),
                            count = as.vector(dataMatrix[1:nR,]))
  
  # produce the plot:
  ggplot(DF, aes(x = census, y = count, fill = sizeClass, col = sizeClass)) +
    geom_area(aes(colour = sizeClass, fill = sizeClass), position = 'stack') +
    scale_fill_brewer(palette = palette, labels = names) + ylab("Count") +
    xlab("Census") + scale_color_brewer(palette = palette, labels = names) +
    labs(fill = "Size class", col = "Size class") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14))
}
