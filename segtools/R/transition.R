## requires common.R

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

ddgram.legend <- function(dd.row, dd.col, row.ord, col.ord) 
{
  legend <-
    list(right = list(fun = dendrogramGrob,
           args = list(x = dd.row, 
             ord = row.ord,
             side = "right",
             size = 10)),
         top = list(fun = dendrogramGrob,
           args = list(x = dd.col, 
             ord = col.ord,
             side = "top")))
  legend
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

plot.transition <- function(filename, mnemonics = NULL, ...) 
{
  data <- read.transition(filename, mnemonics = mnemonics)
  levelplot.transition(data, ...)
}


################## GMTK params ####################

## File should have fields: label, trackname, mean, sd, ...
read.params <- function(filename, mnemonics = NULL, ...)
{
  params <- read.delim(filename, ...)
  params$label <- factor(params$label)
  params$trackname <- factor(params$trackname)
  
  if (length(mnemonics) == 0) {
    ## Generate our own mnemonics
    mnemonics <- generate.param.mnemonics(params)
  }
  data$label <- relabel.factor(data$label, mnemonics)
}

num.gmtk.labels <- function(filename)
{
  lines <- readLines(filename)
  
  start <- grep("^seg_seg", lines)
  con <- textConnection(lines[start + 2])
  num_labels <- as.numeric(scan(con, what = "numeric", n = 1, quiet = TRUE))
  close(con)
  
  num_labels
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

read.gmtk.params <- function(filename, normalize = TRUE, mnemonics = NULL, 
                             cov = FALSE, ...) {
  data <- parse.track.params(filename)  
  data$label <- factor(data$label)
  
  if (length(mnemonics) > 0) {
    data$label <- relabel.factor(data$label, mnemonics)
  }
  
  data <- rename.tracks(data)
  data.melt <- melt.params(data)
  
  res <- covar2sd(data.melt)
  if (normalize) res <- normalize.params(res, cov = cov) 
  
  res
}

## Read gmtk params and either use given mnemonics or generate some.
read.labeled.params <- function(filename, mnemonics = NULL, ...)
{
  params <- read.gmtk.params(filename, mnemonics = mnemonics, ...)
  if (is.null(mnemonics)) {
    ## Generate our own mnemonics
    mnemonics <- generate.param.mnemonics(params)
    ## Copied from relabel.matrix (didn't work because this is 3D array)
    mnemonics.frame <- data.frame(old = as.integer(as.character(mnemonics[,1])),
                                  new = as.character(mnemonics[,2]))
    params.labeled <- params
    colnames(params.labeled) <- mnemonics.frame$new[match(0:(ncol(params)-1), 
                                                          mnemonics.frame$old)]
    col.order <- mnemonics.frame$old + 1
    params.labeled[,col.order,]
  } else {
    params
  }
}

read.observed.params <- function(filename, normalize = TRUE, mnemonics = NULL, 
                                 cov = FALSE, ...) {
  data <- parse.track.params(filename)  
  data$label <- factor(data$label)
  
  if (length(mnemonics) > 0) {
    data$label <- relabel.factor(data$label, mnemonics)
  }
  
  data <- rename.tracks(data)
  data.melt <- melt.params(data)
  
  res <- covar2sd(data.melt)
  if (normalize) res <- normalize.params(res, cov = cov) 
  
  res
}

hclust.params <- function(params) {
  params.mean <- t(params[,,"mean"])
  
  hclust(dist(params.mean))
}

seq.0based <- function(x) {
  seq(0, x - 1)
}

## If filename is specified, outputs mnemonics to a file and returns filename
## Else, returns mnemonic data.frame
generate.param.mnemonics <- function(params, filename = NULL)
{
  hclust.col <- hclust.params(params)
  
  stems <- (cutree(hclust.col, h = 1.0) - 1)[hclust.col$order]
  stems.reorder <- as.integer(factor(stems, levels = unique(stems))) - 1
  
  stem.starts <- c(0, which(diff(stems.reorder) == 1), length(stems.reorder))
  leaves <- do.call(c, lapply(diff(stem.starts), seq.0based))
  stems.leaves <- paste(stems.reorder, leaves, sep = ".")

  ## Before: hclust.col$labels[hclust.col$order], After: stems.leaves
  mnemonics <- data.frame(index = with(hclust.col, labels[order]),
                          mnemonic = stems.leaves, stringsAsFactors = FALSE)

  if (!is.null(filename))  # Write mnemonics to output file
  {
    write.table(mnemonics, file = filename, quote = FALSE, 
                col.names = TRUE, row.names = FALSE, sep = "\t")
    filename
  } else {
    mnemonics
  }
}


abbr.tracknames <- function(tracknames) {
  ## try eliminating everything
  assaynames <- gsub("[^_.]+[_.]", "", tracknames)
  duplicated.assaynames <- assaynames[duplicated(assaynames)]
  which.duplicated <- which(assaynames %in% duplicated.assaynames)
  assaynames[which.duplicated] <- tracknames[which.duplicated]

  ## for failures, try eliminating cell type (middle) only
  ## XXX: duplicative
  assaynames <- gsub("[_.][^_.]+[_.]", ".", assaynames)
  duplicated.assaynames <- assaynames[duplicated(assaynames)]
  which.duplicated <- which(assaynames %in% duplicated.assaynames)
  assaynames[which.duplicated] <- tracknames[which.duplicated]

  assaynames
}

rename.tracks <- function(params) {
  levels(params$trackname) <- gsub("_", ".", levels(params$trackname))
  levels(params$trackname) <- abbr.tracknames(levels(params$trackname))

  params
}

covar2sd <- function(params) {
  ## Convert covar to sd
  params[, , "covar"] <- sqrt(params[, , "covar"])
  dimnames(params)$param <- sub("covar", "sd", dimnames(params)$param)

  params
}

normalize.params <- function(params, cov = FALSE) {
  ## Normalize mean
  mean <- params[, , "mean"]
  mean.range <- t(apply(mean, 1, range))
  mean.min <- mean.range[, 1]
  mean.max <- mean.range[, 2]
  params[, , "mean"] <- (mean - mean.min) / (mean.max - mean.min) 

  ## Make sd into coefficient of variation (capped at 1)
  sds <- params[, , "sd"]
  if (cov) { 
    sds <- sds / rowMeans(mean)
    sds[sds > 1] <- 1
  } else {  # Normalize same as mean
    sds <- sds / (mean.max - mean.min)
    # If any are over 1, scale all down
    if (any(sds > 1)) {
      sds <- sds / max(sds)
    }
  }
  params[, , "sd"] <- sds

  params
}

melt.params <- function(data) {
  params.melted <- melt(data, id.vars = c("label", "trackname", "param"))
  params.cast <- cast(params.melted, trackname ~ label ~ param)
}

parse.track.params <- function(filename) {
  lines <- readLines(filename)

  start <- grep("% means", lines) + 1
  end <- grep("% Components", lines) - 1

  lines.norm <- lines[start:end]

  anonfile <- file()

  lines.interesting.mean <-
    lines.norm[grep("^[^_ ]+_seg[^_ ]+_[^ ]+ 1 .* $", lines.norm)]

  reformulated <- gsub("^([^_ ]+)_seg([^_ ]+)_([^ ]+) 1 (.*)",
                       "\\1\t\\2\t\\3\t\\4", lines.interesting.mean,
                       perl = TRUE)

  writeLines(reformulated, anonfile)

  lines.interesting.covar <-
    lines.norm[grep("^covar_[^ ]+ 1 .* $", lines.norm)]

  reformulated <- gsub("^(covar)_([^ ]+) 1 (.*)",
                      "\\1\t0\t\\2\t\\3", lines.interesting.covar, perl = TRUE)
  writeLines(reformulated, anonfile)

  params <- read.delim(anonfile, header = FALSE,
                       col.names = c("param", "label", "trackname", "val"))
  close(anonfile)

  ## replicate covar for other seg labels
  params.covar <- subset(params, param == "covar")
  copy.covar <- function(label) {
    res <- params.covar
    res$label <- label
    res
  }

  res <- do.call(rbind, c(list(params), 
                          lapply(1:max(params$label), copy.covar)))

  res
}

panel.params <-
  function(x, y, z, subscripts, at = pretty(z), params.sd, 
           ncolors = 2, threshold = FALSE, sd.shape = "box", 
           panel.fill = "mean", box.fill = "gradient", 
           sd.box.size = 0.4, sd.line.size = 0.1,
           panel.outline = FALSE, horizontal.sd = TRUE, ...)
{
  require("grid", quietly = TRUE)
  x <- as.numeric(x)[subscripts]
  y <- as.numeric(y)[subscripts]
  z <- as.numeric(z)[subscripts]
  sds <- as.numeric(params.sd)[subscripts]

  z.low <- z - sds
  z.high <- z + sds

  if (threshold) {
    z.low[z.low < 0] <- 0
    z.high[z.high > 1] <- 1
  }

  ##zcol.low <- level.colors(z.low, at = at, ...)
  ##zcol.high <- level.colors(z.high, at = at, ...)
  zcol <- level.colors(z, at = at, ...)
  for (i in seq(along = z))
  {
    sd.size <- sds[i]
    col.mean <- zcol[i]
    z.range <- seq(from = z.low[i], to = z.high[i], length = ncolors)
    col.gradient <- level.colors(z.range, at = at, ...)
    
    panel.offsets <- seq(from = - 0.5, by = 1 / ncolors, length = ncolors)
    panel.grad.size <- 1 / ncolors
    box.offsets <- seq(from = - sd.size / 2, by = sd.size / ncolors, 
                       length = ncolors)
    box.grad.size <- sd.size / ncolors

    if (horizontal.sd) {
      xs <- x[i] + panel.offsets
      ys <- y[i]
      box.xs <- x[i] + box.offsets
      box.ys <- y[i]
      box.width <- sd.size
      box.height <- sd.box.size
      box.grad.width <- box.grad.size
      box.grad.height <- box.height
      grad.just <- "left"
      panel.grad.width <- panel.grad.size
      panel.grad.height <- 1
      line.width <- sd.size
      line.height <- sd.line.size
    } else {
      xs <- x[i]
      ys <- y[i] + panel.offsets
      box.xs <- x[i]
      box.ys <- y[i] + box.offsets
      box.width <- sd.box.size
      box.height <- sd.size
      box.grad.width <- box.width
      box.grad.height <- box.grad.size
      grad.just <- "bottom"
      panel.grad.width <- 1
      panel.grad.height <- panel.grad.size
      line.width <- sd.size
      line.height <- sd.line.size
    }

    if (panel.fill == "mean") {
      grid.rect(x = x[i], y = y[i], height = 1, width = 1,
                default.units = "native",
                gp = gpar(col = NA, fill = col.mean))
    } else if (panel.fill == "gradient") {
      grid.rect(x = xs, y = ys, height = panel.grad.height, 
                width = panel.grad.width,
                just = grad.just, default.units = "native",
                gp = gpar(col = NA, fill = col.gradient))
    }

    if (!is.null(box.fill)) {
      if (box.fill == "mean") {
        grid.rect(x = x[i], y = y[i], height = 1, width = 1,
                  default.units = "native",
                  gp = gpar(col = NA, fill = col.mean))
      } else if (box.fill == "gradient") {
        grid.rect(x = box.xs, y = box.ys, height = box.grad.height, 
                  width = box.grad.width, just = grad.just, 
                  default.units = "native", 
                  gp = gpar(col = NA, fill = col.gradient))
      }
    }

    if (!is.null(sd.shape)) {
      if (sd.shape == "box") {
        grid.rect(x = x[i], y = y[i], height = box.height, width = box.width,
                  default.units = "native",
                  gp = gpar(col = "black", fill = NA))
      } else if (sd.shape == "line") {
        grid.rect(x = x[i], y = y[i], height = line.height, width = line.width,
                  default.units = "native",
                  gp = gpar(col = NA, fill = "black"))
      }
    }

    if (panel.outline) {
      grid.rect(x = x[i], y = y[i], height = 1, width = 1,
                default.units = "native",
                gp = gpar(col = "black", fill = NA))
    }

  }
}


## params should be a 3D array with [,,"mean"] and [,,"sd"]
levelplot.params <-
  function(params, 
           axis.cex = 1.0,
           scale.cex = 1.0,
           xlab = list("Segment label", cex = axis.cex),
           ylab = list("Input track", cex = axis.cex),
           aspect = "iso",
           scales = list(x = list(rot = 90), cex = scale.cex),
           panel = panel.params,
           threshold = FALSE,
           legend = ddgram.legend(dd.row, dd.col, row.ord, col.ord),
           colorkey = list(space = "left", at = colorkey.at),
           palette = colorRampPalette(rev(brewer.pal(11, "RdYlBu")),
                                      interpolate = "spline", 
                                      space = "Lab")(100),
           ...)
{
  params.means <- params[, , "mean"]
  params.sd <- params[, , "sd"]

  if (threshold) {
    z.range <- range(params.means)
  } else {
    z.range <- c(min(params.means - params.sd), max(params.means + params.sd))
  }
  colorkey.at <- seq(from = z.range[1], to = z.range[2], length = 101)

  dd.row <- as.dendrogram(hclust(dist(params.means)))
  row.ord <- order.dendrogram(dd.row)

  dd.col <- as.dendrogram(hclust(dist(t(params.means))))
  col.ord <- order.dendrogram(dd.col)

  par(oma = c(1, 1, 1, 1))  # Add a margin
  levelplot(t(params.means[row.ord, col.ord]),
            params.sd = t(params.sd[row.ord, col.ord]),
            aspect = aspect,
            scales = scales,
            panel = panel,
            threshold = threshold,
            xlab = xlab, 
            ylab = ylab,
            at = colorkey.at,
            col.regions = palette,
            colorkey = colorkey,
            legend = legend,
            ...)
}

plot.gmtk.params <-
  function(filename, mnemonics = NULL, ...)
{
  params <- read.labeled.params(filename, mnemonics = mnemonics, ...)

  levelplot.params(params, ...)
}


plot.gmtk.transition <-
  function(filename, mnemonics = NULL, ...)
{
  probs <- read.gmtk.transition(filename, mnemonics = mnemonics, ...)

  levelplot.transition(probs, ...)
}

make.gmtk.mnemonic.file <- function(gmtk_file, filename, ...)
{
  params <- read.gmtk.params(gmtk_file, ...)
  mnemonicfile <- generate.param.mnemonics(params, filename = filename, ...)

  mnemonicfile
}

make.mnemonic.file <- function(tabfilename, filename, ...)
{
  params <- read.observed.params(gmtk_file, ...)
  mnemonicfile <- generate.param.mnemonics(params, filename = filename, ...)

  mnemonicfile
}
