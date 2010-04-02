library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(plyr)
library(reshape)

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
  ## Plot x == y reference line
  ## panel.abline(c(0, 1), reference = TRUE, ...)

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
  none.counts <- counts$none
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
  function(data,
           metadata = NULL,
#           x = tpr ~ fpr | factor, 
           x = precision ~ recall | factor, 
           small.cex = 1.0,
           large.cex = 1.0,
           as.table = TRUE,
           aspect = "fill",
           auto.key = FALSE, #list(space = "right"),
#           xlab = list("False positive rate (FP / (FP + TN))", cex = large.cex),
#           ylab = list("True positive rate (TP / (TP + FN))", cex = large.cex),
           xlab = list("Recall (TP / (TP + FN))", cex = large.cex),
           ylab = list("Precision (TP / (TP + FP))", cex = large.cex),
           x.lim = c(0, 1),
           y.lim = c(0, 1),
           scales = list(cex = small.cex,
             x = list(limits = x.lim),
             y = list(limits = y.lim)),
           panel = panel.overlap,
           labels = data$label,
           colors = label.colors(labels),
           par.strip.text = list(cex = small.cex),
           ...)
{
  stats <- overlap.stats(data)
  precision <- suppressMessages(melt(stats$tp / (stats$tp + stats$fp)))
  recall <- suppressMessages(melt(stats$tp / (stats$tp + stats$fn)))
  stats.merged <- data.frame(label = labels, factor = precision$variable,
                            precision = precision$value, recall = recall$value)

  xyplot(x, stats.merged, groups = label,
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

plot.overlap <-  function(tabfile, mnemonics = NULL, col_mnemonics = NULL,
                          comment.char = "#", ...) {
  ## Plot the predictive ability of each segment label for each feature class
  ##   in ROC space
  ##
  ## tabfile: a tab file containing overlap data with segment labels on the rows
  ##   and feature classes on the columns and the overlap "count" at the
  ##   intersection. Row and columns should have labels.
  
  counts <- read.overlap(tabfile, mnemonics = mnemonics,
                         col_mnemonics = col_mnemonics,
                         comment.char = comment.char)
  
  #stat.df <- stat.data.frame(stats)
  #if (!is.null(basename)) {
  #  write.stats(tabbasename = basename, dirpath = dirpath, stat.df)
  #}

  xyplot.overlap(stats, ...)
}

save.overlap <- function(dirpath, namebase, tabfilename,
                         mnemonic_file = NULL,
                         col_mnemonic_file = NULL,
                         clobber = FALSE,
                         panel.size = 300,  # px
                         comment.char = "#",
                         ...) {
  mnemonics <- read.mnemonics(mnemonic_file)
  col.mnemonics <- read.mnemonics(col_mnemonic_file)
  data <- read.overlap(tabfilename, mnemonics = mnemonics,
                       col_mnemonics = col.mnemonics,
                       comment.char = comment.char)
  metadata <- read.metadata(tabfilename, comment.char = comment.char)

  panels.sqrt <- max(c(sqrt(nlevels(data$label)), 1))
  width <- 100 + panel.size * ceiling(panels.sqrt)
  height <- 200 + panel.size * floor(panels.sqrt)
  save.images(dirpath, namebase,
              xyplot.overlap(data, metadata = metadata, ...),
              width = width,
              height = height,
              clobber = clobber)
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
                                 col_mnemonics = NULL,
                                 comment.char = "#", ...) {
  ## Plot the overlap pvalue data
  pvalues <- read.pvalues(tabfile, mnemonics = mnemonics,
                          col_mnemonics = col_mnemonics,
                          comment.char = comment.char)
  barchart.pvalues(pvalues, ...)
}


############### OVERLAP HEATMAP ############

levelplot.overlap <- function(data,
                              metadata = NULL,
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
  metadata <- read.metadata(filename, comment.char = comment.char)

  levelplot.overlap(data, metadata = metadata, ...)
}

save.overlap.heatmap <- function(dirpath, namebase, tabfilename,
                                 mnemonic_file = NULL,
                                 col_mnemonic_file = NULL,
                                 clobber = FALSE,
                                 panel.size = 30,  # px
                                 comment.char = "#",
                                 ...) {
  mnemonics <- read.mnemonics(mnemonic_file)
  col.mnemonics <- read.mnemonics(col_mnemonic_file)
  data <- read.overlap(tabfilename, mnemonics = mnemonics,
                       col_mnemonics = col.mnemonics,
                       comment.char = comment.char)
  metadata <- read.metadata(tabfilename, comment.char = comment.char)

  height <- 400 + panel.size * nrow(data)
  width <- 400 + panel.size * ncol(data)
  save.images(dirpath, namebase,
              levelplot.overlap(data,  metadata = metadata, ...),
              height = height,
              width = width,
              clobber = clobber)
}
