library(plyr)
library(reshape)
library(cluster)

COLNAMES <- c("label", "trackname", "mean", "sd")

## File should have fields: label, trackname, mean, sd, ...
read.track.stats <- function(filename, mnemonics = NULL,
                             ...,
                             colClasses = list(label = "factor",
                               trackname = "factor",
                               mean = "numeric",
                               sd = "numeric"))
{
  stats <- read.delim(filename, colClasses = colClasses)
  if (!all(names(stats) == COLNAMES)) {
    stop("Unrecognized track statistic file format")
  }
  stats$label <- relevel.mnemonics(stats$label, mnemonics)
  stats <- rename.tracks(stats, ...)
  
  stats
}

read.gmtk.track.stats <- function(filename, mnemonics = NULL,
                                  cov = FALSE, ...) {
  means <- parse.gmtk.means(filename)
  covars <- parse.gmtk.covars(filename)
  data <- merge.gmtk.track.stats(means, covars)
  
  data$label <- relevel.mnemonics(data$label, mnemonics)
  data <- rename.tracks(data, ...)
  data <- covar2sd(data)
  
  data
}


merge.gmtk.track.stats <- function(means, covars) {
  if (!all(means$trackname[1:nrow(covars)] == covars$trackname)) {
    stop("Track ordering different between means and covars!")
  } else {
    ## Implicit replication across all segs and subsegs
    stats <- subset(means, select=-value)
    stats$mean <- means$value
    stats$covar <- covars$value
    
    stats
  }
}

is.hierarchical.gmtk <- function(filename) {
  lines <- readLines(filename)
  hierarchical <- length(grep("^seg_subseg $", lines)) > 0

  hierarchical
}

parse.gmtk.means <- function(filename, hierarchical =
                             is.hierarchical.gmtk(filename)) {
  ## hierarchical:
  ##   TRUE: parse the params file expecting a hierarchical format
  ##   FALSE: parse the params file as a normal format
  
  lines <- readLines(filename)
  
  if (hierarchical) {
    pattern <- "^mean_seg([^_ ]+)_subseg([^_ ]+)_([^ ]+) 1 (.*)"
    replacement <- "\\1_\\2\t\\3\t\\4"
  } else {
    pattern <- "^mean_seg([^_ ]+)_([^ ]+) 1 (.*)"
    replacement <- "\\1\t\\2\t\\3"
  }
  
  start <- grep("^% means$", lines) + 1
  lines.norm <- lines[start:length(lines)]

  anonfile <- file()
  lines.interesting <-
    lines.norm[grep(pattern, lines.norm)]

  reformulated <- gsub(pattern, replacement, lines.interesting, perl = TRUE)
  writeLines(reformulated, anonfile)
  
  means <- read.delim(anonfile, header = FALSE,
                      col.names = c("label", "trackname", "value"),
                      colClasses = c("factor", "factor", "numeric"))
  close(anonfile)

  means
}


parse.gmtk.covars <- function(filename) {
  lines <- readLines(filename)
  start <- grep("^% diagonal covariance matrices$", lines) + 1
  lines.norm <- lines[start:length(lines)]

  anonfile <- file()
  lines.interesting <-
    lines.norm[grep("^covar_[^ ]+ 1 .* $", lines.norm)]
  reformulated <- gsub("^covar_([^ ]+) 1 (.*)",
                      "\\1\t\\2", lines.interesting, perl = TRUE)
  writeLines(reformulated, anonfile)
  
  covars <- read.delim(anonfile, header = FALSE,
                       col.names = c("trackname", "value"),
                       colClasses = c("factor", "numeric"))
  close(anonfile)

  covars
}

## You should probably normalize with normalize.track.stats() first
generate.mnemonics <- function(stats, hierarchical = FALSE, clust.func = hclust, ...) {
  if (!is.array(stats)) {
    stats <- as.stats.array(stats)
  }
  
  stats.mean <- t(stats[, , "mean"])
  dist.func <- if (hierarchical) hierarchical.dist else dist
  hclust.col <- clust.func(dist.func(stats.mean))
  
  cut.height <- mean(range(hclust.col$height))
  stems <- (cutree(hclust.col, h = cut.height) - 1)[hclust.col$order]
  stems.reorder <- as.integer(factor(stems, levels = unique(stems))) - 1

  seq.0based <- function(x) {
    seq(0, x - 1)
  }
  stem.starts <- c(0, which(diff(stems.reorder) == 1), length(stems.reorder))
  leaves <- do.call(c, lapply(diff(stem.starts), seq.0based))
  stems.leaves <- paste(stems.reorder, leaves, sep = ".")

  ## Before: hclust.col$labels[hclust.col$order], After: stems.leaves
  mnemonics <- data.frame(old = with(hclust.col, labels[order]),
                          new = stems.leaves, stringsAsFactors = FALSE)
  mnemonics
}

## Can rename tracks according to a regex pattern/replacement
## or a translation table (a Nx2 array, where a trackname found at [i,1]
## are replaced by the string at [i,2])
## If both a translation table and regex replacements are specified, the
## translation table is applied first
rename.tracks <- function(stats, patterns = "_", replacements = ".",
                          translation = NULL, ...) {
  if (!is.data.frame(stats)) {
    stop("rename.tracks expected stats as data.frame")
  }
  
  tracknames <- levels(stats$trackname)
  # Apply translation table substitutions
  if (!is.null(translation)) {
    indices <- match(tracknames, translation[, 1])
    if (any(is.finite(indices))) {
      tracknames[is.finite(indices)] <-
        translation[indices[is.finite(indices)], 2]
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

as.stats.array <- function(data) {
  if (!is.data.frame(data)) {
    stop("as.stats.array expected data as data.frame")
  }
  data.melted <- melt(data, id.vars = c("label", "trackname"))
  cast(data.melted, trackname ~ label ~ variable)
}

as.stats.data.frame <- function(stats) {
  if (!is.array(data)) {
    stop("as.stats.data.frame expected data as array")
  }
  stats.melted <- melt(stats)
  cast(stats.melted, label + trackname ~ variable)
}

covar2sd <- function(stats) {
  ## Convert covar to sd
  stats$sd <- sqrt(stats$covar)
  stats$covar <- NULL

  stats
}

normalize.track.stats <- function(stats, cov = FALSE, sd.scale = TRUE) {
  if (!is.array(stats)) {
    stop("normalize.track.stats expected stats in array format!")
  }
  
  ## Normalize mean
  means <- stats[, , "mean"]
  means.range <- t(apply(means, 1, range))
  means.min <- means.range[, 1]
  means.max <- means.range[, 2]
  stats[, , "mean"] <- (means - means.min) / (means.max - means.min) 

  if ("sd" %in% dimnames(stats)[[3]]) {
    sds <- stats[, , "sd"]
    if (cov) {  # Make sd into coefficient of variation (capped at 1)
      sds <- sds / rowMeans(means)
      sds[sds > 1] <- 1
    } else {  # Normalize same as mean
      sds <- sds / (means.max - means.min)
      ## If any are over 1, scale all down
      if (sd.scale && any(is.finite(sds)) && any(sds > 1)) {
        sds <- sds / max(sds, na.rm = TRUE)
      }
    }
    stats[, , "sd"] <- sds
  }

  stats
}

normalize.binary.track.stats <- function(stats, cov = FALSE) {
  if (!is.array(stats)) {
    stop("normalize.binary.track.stats expected stats in array format!")
  }
  if (! "sd" %in% dimnames(stats)[[3]]) {
    stop("normalize.binary.track.stats expected standard deviations!")
  }
  ## Compute normalization statistic
  means <- stats[, , "mean"]
  sds <- stats[, , "sd"]
  if (dim(means)[2] != 2) {
    stop("normalize.binary.track.stats expected 2-label segmentation results!")
  }

  means.diff <- apply(means, 1, diff)
  stats[, , "mean"] <- cbind(-means.diff, means.diff) / sds
  stats[, , "sd"] <- NA

  stats
}

hierarchical.dist <- function(mat) {
  ## Returns a distance matrix, accomodating a hierarchical model
  
  ## Hierarchy is denoted in the rownames, with "a_b" specifying a's child, b
  ## All nodes of a given level must have the same number of children
  levels.tokens <- strsplit(rownames(mat), "_", fixed = TRUE)
  depth <- length(levels.tokens[[1]])
  for (level.row in levels.tokens) {
    if (length(level.row) != depth) {
      stop("Hierarchy tree does not have a uniform depth")
    }
  }
  if (depth == 1) {
    dist(mat)
  } else {
    ## Convert the list of arrays to a single array,
    ##   since all depths are the same
    levels.mat <- matrix(unlist(levels.tokens), ncol = depth, byrow = TRUE)
    
    mat.dist <- as.matrix(dist(mat))
    dist.thresh <- max(mat.dist) + 1  # Separate hierarchy levels this far
    
    for (row.i in 1:nrow(mat)) {
      diff.rows <- rep(FALSE, nrow(mat))
      for (depth.i in 1:(depth - 1)) {
        cur.level <- levels.mat[row.i, depth.i]
        ## Which rows are different at this level
        diff.level <- levels.mat[, depth.i] != cur.level
        ## Add to distance of rows that are different at this level or a
        ##   higher one
        diff.rows <- diff.rows | diff.level
        mat.dist[diff.rows, row.i] <- mat.dist[diff.rows, row.i] + dist.thresh
      }
    }
    as.dist(mat.dist)
  }
}

panel.track.stats <-
  function(x, y, z, subscripts, at = pretty(z), sds = NULL, 
           ncolors = 2, threshold = FALSE, sd.shape = "box", 
           panel.fill = "mean", box.fill = "gradient", 
           sd.box.size = 0.4, sd.scale = 1, sd.line.size = 0.18,
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
    if (is.finite(sds[i])) {
      sd.size <- sds[i] * sd.scale
    } else {
      sd.size <- 0
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
      line.width <- sd.line.size
      line.height <- sd.size
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
           hclust.label = TRUE,
           hclust.track = TRUE,
           clust.func = hclust,
           hierarchical = FALSE,
           aspect = "iso",
           sd.shape = "line",
           box.fill = NULL,
           symmetric = FALSE,
           scales = list(x = list(rot = 90), cex = scale.cex),
           panel = panel.track.stats,
           threshold = FALSE,
           legend = ddgram.legend(dd.row,  row.ord, dd.col,col.ord),
           colorkey = list(space = "left", at = colorkey.at),
           palette = colorRampPalette(rev(brewer.pal(11, "RdYlBu")),
                                      interpolate = "spline", 
                                      space = "Lab")(100),
           ...)
{
  if (!is.array(track.stats)) {
    track.stats <- as.stats.array(track.stats)
  }

  if (dim(track.stats[, , "mean"])[2] == 2) {
    norm.func <- normalize.binary.track.stats
    symmetric <- TRUE
  } else {
    norm.func <- normalize.track.stats
  }
  stats.norm <- norm.func(track.stats)
  
  means <- stats.norm[, , "mean"]
  sds <- stats.norm[, , "sd"]

  if (!any(is.finite(means))) {
    warning("No finite mean values found. Nothing to plot!")
    return(NULL)
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
  
  if (threshold || is.null(sds) || sd.shape == "line") {
    z.range <- range(means, na.rm = TRUE)
  } else {
    z.range <- c(min(means, means - sds), max(means + sds))
  }
  
  if (symmetric) {
    z.max <- max(abs(z.range))
    z.range <- c(-z.max, z.max)
  }
  colorkey.at <- seq(from = z.range[1], to = z.range[2], length = 101)

  # Set up dendrogram
  dd.row <- NULL
  dd.col <- NULL
  row.ord <- 1:nrow(means)
  col.ord <- 1:ncol(means)
  if (hclust.track) {
    dd.row <- as.dendrogram(clust.func(dist(means)))
    row.ord <- order.dendrogram(dd.row)
  }
  if (hclust.label) {
    dist.func <- if (hierarchical) hierarchical.dist else dist
    dd.col <- as.dendrogram(clust.func(dist.func(t(means))))
    col.ord <- order.dendrogram(dd.col)
  }
  par(oma = c(1, 1, 1, 1))  # Add a margin

  sds.ordered = NULL
  if (! is.null(sds)) {
    sds.ordered = t(sds[row.ord, col.ord])
  }
  levelplot(t(means[row.ord, col.ord, drop = FALSE]),
            sds = sds.ordered,
            aspect = aspect,
            scales = scales,
            panel = panel,
            sd.shape = sd.shape,
            box.fill = box.fill,
            threshold = threshold,
            xlab = xlab, 
            ylab = ylab,
            at = colorkey.at,
            col.regions = palette,
            colorkey = colorkey,
            legend = legend,
            ...)
}

plot.track.stats <- function(filename, mnemonics = NULL,
                             translation_file = NULL,
                             gmtk = FALSE, symmetric = FALSE, ...) {
  translation <-
    if (length(translation_file) > 0 && nchar(translation_file) > 0) {
      read.mnemonics(translation_file)
    } else {
      NULL
    }
  
  if (gmtk) {
    stats <- read.gmtk.track.stats(filename, mnemonics = mnemonics,
                                   translation = translation, ...)
    hierarchical <- is.hierarchical.gmtk(filename)
  } else {
    stats <- read.track.stats(filename, mnemonics = mnemonics,
                              translation = translation, ...)
    hierarchical <- FALSE
  }
  
  levelplot.track.stats(stats, hierarchical = hierarchical,
                        symmetric = symmetric, ...)
}

## infilename: stats data.frame tab file
## outfilename: name of mnemonic file to create
make.mnemonic.file <- function(infilename, outfilename, gmtk = FALSE, ...)
{
  if (gmtk) {
    stats <- read.gmtk.track.stats(infilename, ...)
    hierarchical <- is.hierarchical.gmtk(infilename)
  } else {
    stats <- read.track.stats(infilename, ...)
    hierarchical <- FALSE
  }
  stats <- as.stats.array(stats)
  stats.norm <- normalize.track.stats(stats)
  mnemonics <- generate.mnemonics(stats.norm, hierarchical = hierarchical, ...)
  write.table(mnemonics, file = outfilename, quote = FALSE, 
              col.names = TRUE, row.names = FALSE, sep = "\t")
  
  outfilename
}
