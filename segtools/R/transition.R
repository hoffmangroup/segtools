## requires common.R
## requires track_statistics.R

library(plyr)
library(reshape)
library(cluster)
library(lattice)
library(RColorBrewer)
library(latticeExtra)

## Relabel and reorder matrix row/colnames with mnemonics
relabel.matrix <- function(mat, mnemonics = NULL, 
                           relabel.cols = TRUE, relabel.rows = TRUE)
{
  if (length(mnemonics) > 0) {
    mnemonics.frame <- data.frame(old = as.integer(as.character(mnemonics[,1])),
                                  new = as.character(mnemonics[,2]))
    row.order <- 1:nrow(mat)
    col.order <- 1:ncol(mat)
    if (relabel.rows) {
      ## Substitute label names (assume default labels are named)
      rownames(mat) <- mnemonics.frame$new[match(0:(nrow(mat)-1), 
                                                 mnemonics.frame$old)]
      row.order <- mnemonics.frame$old + 1
    }
    if (relabel.cols) {
      colnames(mat) <- mnemonics.frame$new[match(0:(ncol(mat)-1), 
                                                 mnemonics.frame$old)]
      col.order <- mnemonics.frame$old + 1
    }
    ## Reorder
    mat <- mat[row.order, col.order]
  } else {
    ## Default names
    rownames(mat) <- as.character(0:(nrow(mat)-1))
    colnames(mat) <- as.character(0:(ncol(mat)-1))
  }

  mat
}

read.transition <- function(filename, mnemonics = NULL, ..., header = FALSE)
{
  res <- read.delim(filename, ..., header = FALSE)
  res.labeled <-  relabel.matrix(res, mnemonics = mnemonics)

  res.labeled
}

read.gmtk.transition <- function(filename, mnemonics = NULL, ...)
{
  lines <- readLines(filename)
  
  start <- grep("^seg_seg", lines) + 3
  con <- textConnection(lines[start - 1])
  dims <- as.numeric(scan(con, what = "numeric", n = 2, quiet = TRUE))
  close(con)
  end <- start + dims[1] - 1
  
  con <- textConnection(lines[start:end])
  res <- read.table(con)
  close(con)
  res.labeled <- relabel.matrix(res, mnemonics = mnemonics)

  res.labeled
}


matrix.find_quantile <- function(x, q) 
{
  v = as.vector(x)
  quantile(x, q, names=FALSE)
}

matrix.asymmetry <- function(x)
{
  x[x < 1e-5] = 0  # Threshold low values
  x = x/t(x)
  x[!is.finite(x)] = 0  # Kill weird areas
  x = log2(x)
  x[!is.finite(x)] = 0  # Kill weird areas
  x[abs(x) < 0.5] = 0  # Threshold low values

  x
}

## Generates levelplot of data in given file
## Returns data used to generate plot
## mnemonics: an array where [,1] is old labels and [,2] is new labels
levelplot.transition <-
  function(data,
           ddgram = FALSE, 
           asymmetry = FALSE,
           aspect = "iso",
           scales = list(x = list(rot = 90)),
           legend = if (ddgram) ddgram.legend(dd.row, dd.col, row.ord, col.ord)
                    else list(),
           palette = colorRampPalette(rev(brewer.pal(11, "PiYG")),
             interpolate = "spline", 
             space = "Lab")(100),
           ...)
{
  ## Looking at reciprocal probabilities for this run
  if (asymmetry) {
    data <- matrix.asymmetry(data)
  }
  
  if (ddgram) {
    dd.row <- as.dendrogram(hclust(dist(data)))
    row.ord <- order.dendrogram(dd.row)
    dd.col <- as.dendrogram(hclust(dist(t(data))))
    col.ord <- order.dendrogram(dd.col)
  } else {
    row.ord <- nrow(data):1
    col.ord <- 1:ncol(data)
  }
  
  colorkey <-
    if (ddgram) {
      list(space = "left")
    } else {
      list(space = "right")
    }
  
  par(oma=c(1, 1, 1, 1))  # Add a margin
  levelplot(t(data[row.ord, col.ord]),
            aspect = aspect,
            scales = scales,
            xlab = "End label", 
            ylab = "Start label",
            cuts = 99,
            col.regions = palette,
            colorkey = colorkey,
            legend = legend,
            ...)
}

plot.transition <- function(filename, mnemonics = NULL, gmtk = FALSE, ...) 
{
  if (gmtk) {
    data <- read.gmtk.transition(filename, mnemonics = mnemonics)
  } else {
    data <- read.transition(filename, mnemonics = mnemonics)
  }
  levelplot.transition(data, ...)
}
