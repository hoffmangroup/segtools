library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(plyr)
library(reshape)


############### OVERLAP ANALYSIS ##############

read.overlap <- function(filename, mnemonics = NULL, ...,
                         check.names = FALSE)
{
  res <- read.delim(filename, ..., check.names = check.names)
  colnames(res)[1] = "label"
  res$label <- factor(res$label, levels = unique(res$label))

  if (length(mnemonics) > 0) {
    res$label <- relabel.factor(res$label, mnemonics)
  } else {
    res$label <- smart.factor.reorder(res$label)
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
               #cex = 0.9,
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
  total.sum <- sum(total.counts)
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
           small.cex = 2.0,
           large.cex = 2.4,
           as.table = TRUE,
           aspect = "fill",
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

plot.overlap <- 
  function(tabfile, dirpath = ".", basename = NULL, mnemonics = NULL, ...) 
{
  counts <- read.overlap(tabfile, mnemonics = mnemonics)
  stats <- overlap.stats(counts)

  stat.df <- stat.data.frame(stats)
  if (!is.null(basename)) {
    write.stats(tabbasename = basename, dirpath = dirpath, stat.df)
  }

  xyplot.overlap(data = stats, ...)
}


############### P-VALUE ANALYSIS ############

read.pvalues <- 
  function(filename, mnemonics = NULL,
    ..., header = TRUE, check.names = FALSE)
{
  res <- read.delim(filename, ..., header = header, check.names = check.names)
  colnames(res)[1] = "label"

  # Order label factor by order in label col
  res$label <- factor(res$label, levels = unique(res$label))

  if (length(mnemonics) > 0) {
    res$label <- relabel.factor(res$label, mnemonics)
  } else {
    res$label <- smart.factor.reorder(res$label)
  }

  res
}

panel.pvalues <-
  function(x, y, groups, subscripts, reference, ...)
{
  # Plot p-value == 0.01
  panel.abline(c(0.01, 0), reference = TRUE, ...)

  panel.levelplot(x, y, groups = groups, subscripts = subscripts, ...)
}

levelplot.pvalues <-
  function(x = log10(value) ~ variable + label,
           data, 
           min.val = 1e-6,
           scales = list(x = list(rot = 90), 
                         z = list(log = TRUE)),
           aspect = "iso",
           as.table = TRUE,
           main = "Overlap significance [log10(p-value)]",
           panel = panel.pvalues,
           xlab = "Feature class",
           ylab = "Segment label",
           ..., 
           palette = colorRampPalette(rev(brewer.pal(11, "RdBu")),
             interpolate = "spline", space = "Lab")(100))
{
  data.melted <- melt(data, id.vars="label")

  # Reverse label ordering for plotting
  data.melted$label <- factor(data.melted$label, 
                              levels=rev(levels(data.melted$label)))

  # Truncate small values
  data.melted$value[data.melted$value < min.val] <- min.val

  levelplot(x = x,
            data = data.melted,
            scales = scales,
            aspect = aspect,
            as.table = as.table,
            cuts = 99,
            col.regions = palette,
            main = main,
            xlab = xlab,
            ylab = ylab)
}

plot.overlap.pvalues <- 
  function(tabfile, mnemonics = NULL, ...)
{
  pvalues <- read.pvalues(tabfile, mnemonics = mnemonics)
 
  levelplot.pvalues(data = pvalues, ...)
}