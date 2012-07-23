library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(plyr)
library(reshape2)

############### OVERLAP ANALYSIS ##############

read.overlap <- function(filename, mnemonics = NULL, col_mnemonics = NULL,
                         ..., check.names = FALSE, colClasses = "character",
                         comment.char = "#")
{
  res <- read.delim(filename, ..., check.names = check.names,
                    colClasses = colClasses, comment.char = comment.char)

  # Substitute colnames and reorder cols (except first and last two)
  set.cols <- c(1, ncol(res) - 1, ncol(res))
  colname.map <- map.mnemonics(colnames(res)[-set.cols], col_mnemonics)
  colnames(res) <- c("label", colname.map$labels,
                     colnames(res)[set.cols[-1]])
  res <- res[, c(set.cols[1], colname.map$order + 1, set.cols[-1])]

  # Substitute rownames and reorder rows
  label.map <- map.mnemonics(res$label, mnemonics)
  res$label <- label.map$labels
  res <- res[label.map$order,]
  
  if (ncol(res) == 2) {
    res[, 2] <- as.numeric(res[, 2])
  } else{
    res[, 2:ncol(res)] <- apply(res[, 2:ncol(res)], 2, as.numeric)
  }
  
  metadata <- read.metadata(filename, comment.char = comment.char)

  list(data = res, metadata = metadata)
}

# Given a list of labels, returns a list of colors to give those labels
label.colors <-
  function(labels)
{
  # Determine label groups
  labels <- labels.classify(labels)
  label.groups <- unique(labels$group)

  # Assign a different color to each group
  color.groups <- rainbow(length(label.groups))
  label.colors <- vector(mode = "character", length = nrow(labels))
  
  for (i in 1:length(label.groups)) {
    # Adjust color darker for each increasing number in same group
    label.group <- subset(labels, group == label.groups[i])
    group.ordered <- label.group[order(label.group$index),]

    # Transpose and scale to 0-1 for rgb()
    color.group <- t(col2rgb(color.groups[i]) / 255)
    group.size <- nrow(label.group)
    for (j in 1:group.size) {
      # Subtrack up to 1/3
      color.rgb <- color.group * (1 - 0.33 * (j - 1) / group.size)
      color.rgb[color.rgb < 0] <- 0
      color <- rgb(red = color.rgb[1], 
                   green = color.rgb[2], 
                   blue = color.rgb[3])

      # Insert color back in at appropriate point in color vector
      label <- group.ordered[j,]
      label.index <- (labels$group == label$group &
                      labels$index == label$index)
      label.colors[label.index] <- color
    }
  }

  label.colors
}

panel.overlap <-
  function(x, y, groups, subscripts, labels, colors,
           plot.text = TRUE, plot.points = TRUE, cumulative = TRUE, ...)
{
  if (plot.points) {
    panel.xyplot(x, y, pch = 20, col = "black",
                 subscripts = subscripts, ...)
  }

  if (plot.text) {
    if (plot.points) {
      adj <- c(-0.25, -0.25)
    } else {
      adj <- NULL
    }
    labels <- as.character(labels)
    if (cumulative) {
      ## Prepend a "+" to all but the first label
      x.ord <- order(x)
      indices <- x.ord[2:length(x.ord)]
      labels[indices] <- paste("+", labels[indices], sep = "")
    }
    panel.text(x, y, 
               labels = labels,
               col = colors,
               adj = adj,
               #cex = c(0.2, 0.9),
               ...)
  }
}

# Returns a data table of tp/fp/fn for the given counts
overlap.stats <- function(counts, cumulative = TRUE) {
  labels <- counts$label
  total.counts <- counts$total
  none.counts <- counts$none
  feature.counts <- subset(counts, select = -c(label, total, none))
  feature.sums <- colSums(feature.counts)
#  total.sum <- sum(as.numeric(total.counts))

  tp <- feature.counts
  fp <- total.counts - feature.counts
  fn <- t(feature.sums - t(feature.counts))
#  tn <- total.sum - total.counts - fn

  if (cumulative) {
    precision <- tp / (tp + fp)
    ## Replace columns with cumulative data (sorted by precision)
    for (col.i in 1:ncol(precision)) {
      label.ord <- order(precision[, col.i], decreasing = TRUE)
      ## Calculate cumulative statistics
      total.cum <- cumsum(total.counts[label.ord])
      tp.cum <- cumsum(feature.counts[label.ord, col.i])
      fp.cum <- total.cum - tp.cum
      fn.cum <- feature.sums[col.i] - tp.cum
      tp[label.ord, col.i] <- tp.cum
      fp[label.ord, col.i] <- fp.cum
      fn[label.ord, col.i] <- fn.cum
    }
  }
  
  res <- list(label = labels,
              tp = tp,
              fp = fp,
              fn = fn)
                     
  res
}

# Returns a data frame of p/r for the given counts
pr.stats <- function(data, ...) {
  stats <- overlap.stats(data$data, ...)

  precision <- suppressMessages(melt(with(stats, tp / (tp + fp))))
  recall <- suppressMessages(melt(with(stats, tp / (tp + fn))))

  res <- data.frame(label = stats$label, group = precision$variable,
                    precision = precision$value, recall = recall$value)
  
  list(data = res, metadata = data$metadata)
}

# Write an overlap data frame to a file
write.stats <- function(stats, namebase, dirpath = ".", clobber = FALSE) {
  data <- stats$data
  tabfilename <- .extpaste(namebase, "tab")
  tabfilepath <- file.path(dirpath, tabfilename)

  if (file.exists(tabfilepath) && !clobber) {
    stop("Found stats table: ", tabfilepath,
         ". Use --clobber to overwrite!", sep = "")
  }
  row.ord = order(data$group, data$precision, decreasing = TRUE)
  data.formatted <- format(data[row.ord,], digits = 3)
  write.table(data.formatted, tabfilepath, quote = FALSE,
              sep = "\t", row.names = FALSE)
}

xyplot.overlap <-
  function(overlap.data = NULL,
           data = overlap.data$data,
           metadata = overlap.data$metadata,
           mode = metadata$mode,
           x = precision ~ recall | group, 
           small.cex = 1.0,
           large.cex = 1.0,
           as.table = TRUE,
           auto.key = FALSE, #list(space = "right"),
           xlab = list("Recall (TP / (TP + FN))", cex = large.cex),
           ylab = list("Precision (TP / (TP + FP))", cex = large.cex),
           cumulative = TRUE,
           sub = if (cumulative) "Cumulative P-R (left-to-right, by precision)"
                 else NULL,
           x.lim = c(0, 1),
           y.lim = NULL,
           scales = list(cex = small.cex,
             x = list(limits = x.lim),
             y = list(limits = y.lim)),
           panel = panel.overlap,
           labels = data$label,
           colors = label.colors(labels),
           par.strip.text = list(cex = small.cex),
           ...)
{
  if (mode != "bases") {
    warning("P-R performance plot only mathematically valid for base-overlap")
  }
  if (is.null(y.lim)) {
    y.max <- max(data$precision)
    y.lim <-
      if (y.max >= 0.5) {
        c(0, 1)
      } else {
        c(0, y.max * 1.05)
      }
  }
  
  xyplot(x, data, groups = label,
         as.table = as.table, 
         auto.key = auto.key,
         xlab = xlab,
         ylab = ylab,
         sub = sub,
         scales = scales,
         cumulative = cumulative,
         panel = panel,
         labels = labels,
         colors = colors,
         par.strip.text = par.strip.text,
         ...)
}

plot.overlap.performance <-  function(tabfile, mnemonics = NULL,
                                      col_mnemonics = NULL,
                                      cumulative = TRUE,
                                      comment.char = "#", ...) {
  ## Plot the predictive ability of each segment label for each feature class
  ##   in PR space
  ##
  ## tabfile: a tab file containing overlap data with segment labels on the rows
  ##   and feature classes on the columns and the overlap "count" at the
  ##   intersection. Row and columns should have labels.
  
  data <- read.overlap(tabfile, mnemonics = mnemonics,
                       col_mnemonics = col_mnemonics,
                       comment.char = comment.char)
  stats <- pr.stats(data, cumulative = cumulative)
  
  xyplot.overlap(stats, cumulative = cumulative, ...)
}

save.overlap.performance <- function(dirpath, namebase, tabfilename,
                                     mnemonic_file = NULL,
                                     col_mnemonic_file = NULL,
                                     clobber = FALSE,
                                     cumulative = TRUE,
                                     panel.size = 300,  # px
                                     comment.char = "#",
                                     ...) {
  mnemonics <- read.mnemonics(mnemonic_file)
  col.mnemonics <- read.mnemonics(col_mnemonic_file)
  data <- read.overlap(tabfilename, mnemonics = mnemonics,
                       col_mnemonics = col.mnemonics,
                       comment.char = comment.char)
  stats <- pr.stats(data, cumulative = cumulative)
  write.stats(stats, namebase, dirpath, clobber = clobber)
  
  panels.sqrt <- max(c(sqrt(nlevels(stats$data$group)), 1))
  width <- 100 + panel.size * ceiling(panels.sqrt)
  height <- 200 + panel.size * floor(panels.sqrt)
  save.plots(dirpath, namebase,
              xyplot.overlap(stats, cumulative = cumulative, ...),
              width = width,
              height = height,
              clobber = clobber)
}

############### OVERLAP HEATMAP ############

levelplot.overlap <- function(overlap.data = NULL,
                              data = overlap.data$data,
                              metadata = overlap.data$metadata,
                              mode = metadata[["mode"]],
                              row.normalize = TRUE,
                              y.mode = if (row.normalize) "Fraction"
                                       else "Count",
                              sub = paste(y.mode, "of", mode, "in subject",
                                "label that overlap at least one in query",
                                "label"),
                              xlab = "label in query file",
                              ylab = "label in subject file",
                              none.col = TRUE,  # Include 'none' column
                              cluster = FALSE,  # Cluster both dimensions
                              cluster.rows = cluster,
                              cluster.cols = cluster,
                              num.colors = 100,
                              aspect = "fill",
                              max.contrast = FALSE,  # Saturate color range
                              col.range = if (max.contrast) NULL else c(0, 1),
                              scales = list(x = list(rot = 90)),
                              palette = colorRampPalette(
                                rev(brewer.pal(11, "RdYlBu")),
                                interpolate = "spline",
                                space = "Lab")(num.colors),
                              cuts = num.colors - 1,
                              ...)
{
  ## Create a levelplot showing overlap proportions
  ##
  ## data: data frame with subject labels on rows, query labels on cols, and
  ##   proportion of coverage at intersection
  ## mode: "segments", "bases" or whatever the units of overlap are
  ## col.range: NULL sets the colorscale to the range of the data, else
  ##   it should be a vector or list of two integers which specify the
  ##   lower and upper bounds of the color scale. Overrides max.contrast

  if (is.null(mode)) stop("Overlap file missing 'mode' metadata.")
  
  ## Convert to matrix
  mat <- subset(data, select = -c(label, total))
  mat.nonone <- subset(mat, select = -c(none))
  if (!none.col) {
    mat <- mat.nonone
  }
  mat <- as.matrix(mat)
  if (row.normalize) {
    mat <- mat / data$total
  }
  rownames(mat) <- data$label
  
  ## Order rows and columns
  row.ord <- nrow(mat.nonone):1
  col.ord <- 1:ncol(mat.nonone)
  
  ## Cluster, holding out none col
  ## This "clustering", or more appropiately, reording of rows and columns
  ## to increase density along the diagonal is done using multidimensional
  ## scaling techniques to find a linear order. This MDS approach to
  ## ordering a matrix is discussed in the "seriation" package, under
  ## the seriate.dist documentation.
  if (cluster.rows) {
    row.ord <- order(cmdscale(dist(mat.nonone), k = 1))
  }
  if (cluster.cols) {
    col.ord <- rev(order(cmdscale(dist(t(mat.nonone)), k = 1)))
  }
  if (none.col) {  # Re-add none col to end
    col.ord <- c(col.ord, length(col.ord) + 1)
  }
  
  if (is.null(col.range)) {
    col.range <- range(mat)
  } else if (length(col.range) != 2) {
    stop("Invalid value of col.range")
  }

  colorkey.at <- seq(col.range[[1]], col.range[[2]], length = num.colors - 1)
  levelplot(t(mat[row.ord, col.ord]),
            sub = sub,
            xlab = xlab,
            ylab = ylab,
            cuts = cuts,
            scales = scales,
            aspect = aspect,
            at = colorkey.at,
            col.regions = palette,
            ...)
}

plot.overlap.heatmap <- function(filename, mnemonics = NULL,
                                 col_mnemonics = NULL,
                                 comment.char = "#",
                                 ...) {
  ## Plot a heatmap from overlap data
  ##
  ## filename: overlap table file
  ## mnemonics, col_mnemonics: mnemonic list (as per read.mnemonics)
  
  data <- read.overlap(filename, mnemonics = mnemonics,
                       col_mnemonics = col_mnemonics,
                       comment.char = comment.char)

  levelplot.overlap(data, ...)
}

save.overlap.heatmap <- function(dirpath, namebase, tabfilename,
                                 mnemonic_file = NULL,
                                 col_mnemonic_file = NULL,
                                 clobber = FALSE,
                                 max_contrast = FALSE,
                                 panel.size = 30,  # px
                                 comment.char = "#",
                                 ...) {
  mnemonics <- read.mnemonics(mnemonic_file)
  col.mnemonics <- read.mnemonics(col_mnemonic_file)
  data <- read.overlap(tabfilename, mnemonics = mnemonics,
                       col_mnemonics = col.mnemonics,
                       comment.char = comment.char)

  height <- 400 + panel.size * nrow(data$data)
  width <- 400 + panel.size * ncol(data$data)
  save.plots(dirpath, namebase,
              levelplot.overlap(data, ..., max.contrast = max_contrast),
              height = height,
              width = width,
              clobber = clobber)
}
