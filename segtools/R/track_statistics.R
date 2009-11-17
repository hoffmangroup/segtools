library(plyr)
library(reshape)
library(cluster)

COLNAMES <- c("label", "trackname", "mean", "sd")

## File should have fields: label, trackname, mean, sd, ...
read.track.stats <- function(filename, mnemonics = NULL, ...,
                        colClasses = list(label = "character"))
{
  stats <- read.delim(filename, colClasses = colClasses)
  if (!all(names(stats) == COLNAMES)) {
    stop("Unrecognized track statistic file format")
  }
  
  stats$label <- factor(stats$label)

  if (length(mnemonics) > 0) {
    stats$label <- relabel.factor(stats$label, mnemonics)
  }
  stats.renamed <- rename.tracks(stats, ...)
  stats.melted <- melt.track.stats(stats.renamed)
  
  stats.melted
}

read.gmtk.track.stats <- function(filename, normalize = TRUE, mnemonics = NULL, 
                                  cov = FALSE, ...) {
  data <- parse.gmtk.track.stats(filename)
  data$label <- factor(data$label)
  
  if (length(mnemonics) > 0) {
    data$label <- relabel.factor(data$label, mnemonics)
  }
  data <- rename.tracks(data, ...)
  data.shaped <- cast(data, trackname ~ label ~ variable)
  
  res <- covar2sd(data.shaped)
  if (normalize) res <- normalize.track.stats(res, cov = cov) 
  
  res
}

covar2sd <- function(stats) {
  ## Convert covar to sd
  stats[, , "covar"] <- sqrt(stats[, , "covar"])
  dimnames(stats)[[3]] <- sub("covar", "sd", dimnames(stats)[[3]])

  stats
}

parse.gmtk.track.stats <- function(filename) {
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
  
  stats <- read.delim(anonfile, header = FALSE,
                      col.names = c("variable", "label", "trackname", "value"))
  close(anonfile)

  ## replicate covar for other seg labels
  stats.covar <- subset(stats, variable == "covar")
  copy.covar <- function(label) {
    res <- stats.covar
    res$label <- label
    res
  }
  res <- do.call(rbind, c(list(stats),
                          lapply(1:max(stats$label), copy.covar)))

  res
}

hclust.track.stats <- function(stats) {
  stats.mean <- t(stats[,,"mean"])
  hclust(dist(stats.mean))
}

seq.0based <- function(x) {
  seq(0, x - 1)
}

## Creates mnemonics from 3d stats array
## You should probably normalize with normalize.track.stats() first
generate.mnemonics <- function(stats)
{
  hclust.col <- hclust.track.stats(stats)
  cut.height <- median(hclust.col$height)
  stems <- (cutree(hclust.col, h = cut.height) - 1)[hclust.col$order]
  stems.reorder <- as.integer(factor(stems, levels = unique(stems))) - 1
  
  stem.starts <- c(0, which(diff(stems.reorder) == 1), length(stems.reorder))
  leaves <- do.call(c, lapply(diff(stem.starts), seq.0based))
  stems.leaves <- paste(stems.reorder, leaves, sep = ".")

  ## Before: hclust.col$labels[hclust.col$order], After: stems.leaves
  mnemonics <- data.frame(index = with(hclust.col, labels[order]),
                          mnemonic = stems.leaves, stringsAsFactors = FALSE)
  mnemonics
}

## Can rename tracks according to a regex pattern/replacement
## or a translation table (a Nx2 array, where a trackname found at [i,1]
## are replaced by the string at [i,2])
## If both a translation table and regex replacements are specified, the
## transition table is applied first
rename.tracks <- function(stats, patterns = "_", replacements = ".", translation = NULL, ...) {
  tracknames <- levels(stats$trackname)
  # Apply translation table substitutions
  if (!is.null(translation)) {
    indices <- match(tracknames, translation[, 1])
    if (any(is.finite(indices))) {
      tracknames[is.finite(indices)] <- translation[indices[is.finite(indices)], 2]
    }
  }
  # Apply regex substitutions
  if (!is.null(patterns) || !is.null(replacements)) {
    if (length(patterns) != length(replacements)) {
      stop("Must have an equal number of patterns and replacements")
    }
    for (i in 1:length(patterns)) {
      tracknames <- gsub(patterns[[i]], replacements[[i]], tracknames)
    }
  }
  levels(stats$trackname) <- tracknames
  
  stats
}

normalize.track.stats <- function(stats, cov = FALSE) {
  ## Normalize mean
  mean <- stats[, , "mean"]
  mean.range <- t(apply(mean, 1, range))
  mean.min <- mean.range[, 1]
  mean.max <- mean.range[, 2]
  stats[, , "mean"] <- (mean - mean.min) / (mean.max - mean.min) 

  if ("sd" %in% dimnames(stats)[[3]]) {
    sds <- stats[, , "sd"]
    if (cov) {  # Make sd into coefficient of variation (capped at 1)
      sds <- sds / rowMeans(mean)
      sds[sds > 1] <- 1
    } else {  # Normalize same as mean
      sds <- sds / (mean.max - mean.min)
      ## If any are over 1, scale all down
      if (any(sds > 1)) {
        sds <- sds / max(sds, na.rm = TRUE)
      }
    }
    stats[, , "sd"] <- sds
  }

  stats
}

melt.track.stats <- function(data) {
  stats.melted <- melt(data, id.vars = c("label", "trackname"))
  stats.cast <- cast(stats.melted, trackname ~ label ~ variable)
  
  stats.cast
}

panel.track.stats <-
  function(x, y, z, subscripts, at = pretty(z), sds = NULL, 
           ncolors = 2, threshold = FALSE, sd.shape = "box", 
           panel.fill = "mean", box.fill = "gradient", 
           sd.box.size = 0.4, sd.scale = 1, sd.line.size = 0.1,
           panel.outline = FALSE, horizontal.sd = TRUE, ...)
{
  require("grid", quietly = TRUE)
  x <- as.numeric(x)[subscripts]
  y <- as.numeric(y)[subscripts]
  z <- as.numeric(z)[subscripts]
  if (is.null(sds)) {
    sds <- 0
  } else {
    sds <- as.numeric(sds)[subscripts]
  }
  
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
    if (is.null(sds) || !is.finite(sds)) {
      sd.size <- 0
    } else {
      sd.size <- sds[i] * sd.scale
    }
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

    if (!is.null(sd.shape) && sd.size > 0) {
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

levelplot.track.stats <-
  function(track.stats,
           axis.cex = 1.0,
           scale.cex = 1.0,
           xlab = list("Segment label", cex = axis.cex),
           ylab = list("Signal track", cex = axis.cex),
           aspect = "iso",
           scales = list(x = list(rot = 90), cex = scale.cex),
           panel = panel.track.stats,
           threshold = FALSE,
           legend = ddgram.legend(dd.row, dd.col, row.ord, col.ord),
           colorkey = list(space = "left", at = colorkey.at),
           palette = colorRampPalette(rev(brewer.pal(11, "RdYlBu")),
                                      interpolate = "spline", 
                                      space = "Lab")(100),
           ...)
{
  stats.norm <- normalize.track.stats(track.stats)
  means <- stats.norm[,,"mean"]
  sds <- stats.norm[,,"sd"]

  if (!any(is.finite(means))) {
    stop("No finite mean values found. Nothing to plot!")
  } else if (!any(is.finite(sds))) {
    ## Pretent no sds were specified
    sds <- NULL
  }

  mask.rows <- apply(!is.finite(means), 1, any)  # Any row with an NaN
  if (any(mask.rows)) {
    warning("Found infinite or NaN mean values. Ignoring those signal tracks")
    means <- means[!mask.rows,]
    sds <- sds[!mask.rows,]
  }
  
  if (threshold || is.null(sds)) {
    z.range <- range(means, na.rm = TRUE)
  } else {
    z.range <- c(min(means, means - sds), max(means + sds))
  }
  colorkey.at <- seq(from = z.range[1], to = z.range[2], length = 101)

  dd.row <- as.dendrogram(hclust(dist(means)))
  row.ord <- order.dendrogram(dd.row)

  dd.col <- as.dendrogram(hclust(dist(t(means))))
  col.ord <- order.dendrogram(dd.col)

  par(oma = c(1, 1, 1, 1))  # Add a margin

  sds.ordered = NULL
  if (! is.null(sds)) {
    sds.ordered = t(sds[row.ord, col.ord])
  }
  levelplot(t(means[row.ord, col.ord]),
            sds = sds.ordered,
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

plot.track.stats <- function(filename, mnemonics = NULL, gmtk = FALSE, ...)
{
  if (gmtk) {
    stats <- read.gmtk.track.stats(filename, mnemonics = mnemonics, ...)
  } else {
    stats <- read.track.stats(filename, mnemonics = mnemonics, ...)
  }
  
  levelplot.track.stats(stats, ...)
}

## infilename: stats data.frame tab file
## outfilename: name of mnemonic file to create
make.mnemonic.file <- function(infilename, outfilename, gmtk = FALSE, ...)
{
  if (gmtk) {
    stats <- read.gmtk.track.stats(infilename, ...)
  } else {
    stats <- read.track.stats(infilename, ...)
  }
  stats.norm <- normalize.track.stats(stats)
  mnemonics <- generate.mnemonics(stats.norm, ...)
  write.table(mnemonics, file = outfilename, quote = FALSE, 
              col.names = TRUE, row.names = FALSE, sep = "\t")
  
  outfilename
}
