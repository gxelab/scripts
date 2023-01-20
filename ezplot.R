# author: zh(mt1022)

# color pallete
library(RColorBrewer)
brewer.pal.set1 <- brewer.pal(9, name = 'Set1')

# graphic parameter for publication
op <- par(las = 1, family = 'Liberation Sans')

# plot multiple ecdf in a plot
PlotEcdf <-function(ll, line.col = 1:length(ll), axes = T, ...){
    plot(ecdf(ll[[1]]), col = line.col[1], axes = axes, ...)
    for(i in 2:length(ll)){
        plot(ecdf(ll[[i]]), col = line.col[i], add = T, ...)
    }
}

# general plot
Plot <- function(x, y = NULL, ...){
	plot(x, y, main = NA, xlab = NA, ylab = NA, axes = F, lwd = 4, cex = 2, ...)
	axis(1, lwd = 4, cex.axis = 2)
	axis(1, lwd = 4, cex.axis = 2)
}

# plot reversed ecdf
# example:
#	x <- 1:10
#	y = rcdf(x)
#	plot(y[, 1], y[, 2], type = 's')
Rcdf <- function(x){
    n <- length(x)
    x.min <- min(x)
    x.max <- max(x)
    y <- sapply(0:n/n * (x.max - x.min) + x.min, function(i) sum(x >= i)/n)
    res <- data.frame(x1 = 0:n/n * (x.max - x.min) + x.min, x2 = y)
    return(res)
}

# given vector of p.value, return significance code
SigCode <- function(x){
    res <- sapply(x, function(i){
        if(i < 0.001){
            return('***')
        }else if(i < 0.01){
            return('**')
        }else if(i < 0.05){
            return('*')
        }else{
            return('')
        }
    })
    return(res)
}

SigCode.2 <- function(x){
    res <- sapply(x, function(i){
        if(i < 0.001){
            return('***')
        }else if(i < 0.01){
            return('**')
        }else if(i < 0.05){
            return('')
        }else{
            return('')
        }
    })
    return(res)
}

# given relative coordiante, return actual coordiante;
RelPos <- function(x, y, r = par('usr')){
    res <- c(
        x = r[1] + (r[2] - r[1]) * x,
        y = r[3] + (r[4] - r[3]) * y
    )
}

# return plot expression in scientific format for p.value
# @param x input p.value
# @return res expression
SciExp <- function(x, prefix = 'P = '){
    sci <- strsplit(sprintf('%.2e', x), 'e')[[1]]
    sci.text <- paste("'", prefix, sci[1], "x10'^", as.integer(sci[2]), sep = '')
    res <- parse(text = sci.text)
    return(res)
}
