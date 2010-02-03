library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(plyr)
library(reshape)


############### OVERLAP ANALYSIS ##############

read.overlap <- function(filename, mnemonics = NULL, col_mnemonics = NULL,
                         ..., check.names = FALSE, colClasses = "character")
{
  res <- read.delim(filename, ..., check.names = check.names,
                    colClasses = colClasses)

  # Substitute colnames and reorder cols
  colname.map <- map.mnemonics(colnames(res)[-1], col_mnemonics)
  colnames(res) <- c("label", colname.map$labels)
  res <- res[, c(1, colname.map$order + 1)]

  # Substitute rownames and reorder rows
  label.map <- map.mnemonics(res$label, mnemonics)
  res$label <- label.map$labels
  res <- res[label.map$order,]
  
  if (ncol(res) == 2) {
    res[, 2] <- as.numeric(res[, 2])
  } else{
    res[, 2:ncol(res)] <- apply(res[, 2:ncol(res)], 2, as.numeric)
  }
  
  res
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
  function(x, y, groups, subscripts, labels, colors, reference, 
           plot.text = TRUE, plot.points = FALSE, significances = NULL, ...)
{
  # Plot x == y reference line
  panel.abline(c(0, 1), reference = TRUE, ...)

  if (plot.points) {
    panel.xyplot(x, y, groups = groups, col = colors, 
                 subscripts = subscripts, ...)
  }

  if (plot.text) {
    if (plot.points) {
      pos <- 1
      offset <- 1
    } else {
      pos <- NULL
      offset <- NULL
    }
    panel.text(x, y, 
               labels = labels,
               col = colors,
               pos = pos,
               offset = offset,
               #cex = c(0.2, 0.9),
               ...)
  }

  # Plot significances next to points
  if (! is.null(significances)) {
    significant <- significances < 0.05
    panel.text(x[significant], y[significant],
               labels = format(significances[significant], digits = 1),
               col = "black",
               pos = 4,
               #cex = 0.8,
               offset = 1,
               ...)
  }
}

# Returns a data table of tp/tn/fp/fn for the given counts
overlap.stats <-
  function(counts)
{
  total.counts <- counts$total
  feature.counts <- subset(counts, select = -c(label, total, none))
  total.sum <- sum(as.numeric(total.counts))
  feature.sums <- colSums(feature.counts)
  tp <- feature.counts
  fp <- total.counts - feature.counts
  fn <- t(feature.sums - t(feature.counts))
  tn <- total.sum - total.counts - fn

  res <- list()
  res$label <- counts$label
  res$tp <- tp
  res$tn <- tn
  res$fp <- fp
  res$fn <- fn
  res$col.names <- c("label",
                     paste("tp", colnames(tp), sep = "."),
                     paste("tn", colnames(tn), sep = "."),
                     paste("fp", colnames(fp), sep = "."),
                     paste("fn", colnames(fn), sep = "."))
                     
  res
}

# Convert the stats to a single data frame (suitable for write.stats)
stat.data.frame <-
  function(stats)
{
  stat.df <- data.frame(stats$label, stats$tp, stats$tn, stats$fp, stats$fn)
  colnames(stat.df) <- stats$col.names
  
  stat.df
}

# Write a stat data frame to a file
write.stats <-
  function(tabbasename, dirpath = ".", stat.df)
{
  tabfilename <- extpaste(tabbasename, "tab")
  tabfilepath <- file.path(dirpath, tabfilename)
  
  write.table(stat.df, tabfilepath, quote = FALSE, sep = "\t", row.names = FALSE)
}

xyplot.overlap <-
  function(x = tpr ~ fpr | factor, 
           data,
           small.cex = 1.0,
           large.cex = 1.0,
           as.table = TRUE,
           aspect = "iso",
           auto.key = FALSE, #list(space = "right"),
           xlab = list("False positive rate (FP / (FP + TN))", cex = large.cex),
           ylab = list("True positive rate (TP / (TP + FN))", cex = large.cex),
           scales = list(cex = small.cex),
           panel = panel.overlap,
           labels = data$label,
           colors = label.colors(labels),
           par.strip.text = list(cex = small.cex),
           ...)
{
  tpr <- suppressMessages(melt(data$tp / (data$tp + data$fn)))
  fpr <- suppressMessages(melt(1 - data$tn / (data$tn + data$fp)))
  data.merged <- data.frame(label = labels, factor = tpr$variable,
                            tpr = tpr$value, fpr = fpr$value)

  xyplot(x, data.merged, groups = label,
         as.table = as.table, 
         aspect = aspect,
         auto.key = auto.key,
         xlab = xlab, ylab = ylab,
         scales = scales,
         panel = panel,
         labels = labels,
         colors = colors,
         par.strip.text = par.strip.text,
         ...)
}

plot.overlap <-  function(tabfile, mnemonics = NULL, col_mnemonics = NULL, ...) 
{
  ## Plot the predictive ability of each segment label for each feature class
  ##   in ROC space
  ##
  ## tabfile: a tab file containing overlap data with segment labels on the rows
  ##   and feature classes on the columns and the overlap "count" at the
  ##   intersection. Row and columns should have labels.
  
  counts <- read.overlap(tabfile, mnemonics = mnemonics,
                         col_mnemonics = col_mnemonics)
  stats <- overlap.stats(counts)

  #stat.df <- stat.data.frame(stats)
  #if (!is.null(basename)) {
  #  write.stats(tabbasename = basename, dirpath = dirpath, stat.df)
  #}

  xyplot.overlap(data = stats, ...)
}


############### P-VALUE ANALYSIS ############

read.pvalues <-  function(...) {
  read.overlap(...)
}

panel.pvalues <- function(x, y, subscripts = NULL, groups = NULL,
                          reference = 0.01, col = NULL, ...)
{
  ref.x <- log10(reference)
  panel.refline(v = ref.x)
  out.of.bounds = is.infinite(x)
  x[out.of.bounds] <- min(x[!out.of.bounds]) * 2

  panel.barchart(x, y, subscripts = subscripts, groups = groups,
                 stacked = FALSE, col = col[y], ...)
}

barchart.pvalues <- function(data,
                             x = reorder(label, -value) ~ value,
                             groups = if (ngroups > 1) variable else NULL,
                             as.table = TRUE,
                             main = "Approximate significance of overlap",
                             panel = panel.pvalues,
                             xlab = "p-value",
                             ylab = "Segment label",
                             reference = 0.01,
                             origin = 0,
                             auto.key = if (ngroups > 1) list(points = FALSE)
                                        else FALSE,
                             colors = label.colors(data.melted$label),
                             scales = list(x = list(log = TRUE)),
                             ...)
{
  ## Create a barchart from overlap pvalue data
  ##
  ## data: a data frame containing pvalue data
  
  data.melted <- melt(data, id.vars = "label")
  ngroups = nlevels(data.melted$variable)
  
  colors.ordered <- colors[order(-data.melted$value)]
  xyplot(x = x,
         data = data.melted,
         groups = groups,
         panel = panel,
         as.table = as.table,
         reference = reference,
         main = main,
         col = colors.ordered,
         scales = scales,
         origin = origin,
         auto.key = auto.key,
         xlab = xlab,
         ylab = ylab,
         #par.settings = list(clip = list(panel = "off")),
         ...)
}

plot.overlap.pvalues <- function(tabfile, mnemonics = NULL,
                                 col_mnemonics = NULL, ...) {
  ## Plot the overlap pvalue data
  pvalues <- read.pvalues(tabfile, mnemonics = mnemonics,
                          col_mnemonics = col_mnemonics)
  barchart.pvalues(pvalues, ...)
}


############### OVERLAP HEATMAP ############

levelplot.overlap <- function(mat,
                              mode = "segments",
                              y.mode = "Fraction",
                              sub = paste(y.mode, "of", mode, "in subject",
                                "label that overlap at least one in query",
                                "label"),
                              xlab = "label in query file",
                              ylab = "label in subject file",
                              num.colors = 100,
                              col.range = NULL,
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
  ## mat: matrix with subject labels on rows, query labels on cols, and
  ##   proportion of coverage at intersection
  ## mode: "segments", "bases" or whatever the units of overlap are
  ## col.range: NULL sets the colorscale to the range of the data, else
  ##   it should be a vector or list of two integers which specify the
  ##   lower and upper bounds of the color scale
  row.ord <- nrow(mat):1
  col.ord <- 1:ncol(mat)

  
  
  plot.levelplot <- function(...) {
    levelplot(t(mat[row.ord, col.ord, drop = FALSE]),
              sub = sub,
              xlab = xlab,
              ylab = ylab,
              cuts = cuts,
              scales = scales,
              col.regions = palette,
              ...)
  }
  
  if (!is.null(col.range)) {
    if (length(col.range) != 2) stop("Invalid value of col.range")
    at <- seq(col.range[[1]], col.range[[2]], length = num.colors - 1)

    plot.levelplot(at = at)
  } else {
    plot.levelplot()
  }
}

plot.overlap.heatmap <- function(filename, mnemonics = NULL,
                                 col_mnemonics = NULL,
                                 none_col = FALSE,
                                 row_normalize = TRUE,
                                 max_contrast = TRUE,
                                 ...) {
  ## Plot a heatmap from overlap data
  ##
  ## filename: overlap table file
  ## mnemonics, col_mnemonics: mnemonic list (as per read.mnemonics)
  ## none_col: include the none column (TRUE) or not (FALSE)
  ## row_normalize: always (TRUE), never (FALSE), or only when there are
  ##   at least two columns (default).
  ## max_contrast: set color range to range of data instead of [0, 1]
  
  data <- read.overlap(filename, mnemonics = mnemonics,
                       col_mnemonics = col_mnemonics)

  # Convert to matrix
  data.mat <- subset(data, select = -c(label, total))
  if (!none_col) {
    data.mat <- subset(data.mat, select = -c(none))
  }
  data.mat <- as.matrix(data.mat)
  
  if (row_normalize) {
    data.mat <- data.mat / data$total
    y.mode <- "Fraction"
  } else {
    y.mode <- "Number"
  }
  rownames(data.mat) <- data$label
  
  col.range <- if (max_contrast) NULL else c(0, 1)
  levelplot.overlap(data.mat, col.range = col.range,
                    y.mode = y.mode, ...)
}
