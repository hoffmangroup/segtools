library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(plyr)
library(reshape)


############### OVERLAP ANALYSIS ##############

read.overlap <- function(filename, mnemonics = NULL, ...,
                         colClasses = list(X = "character"))  # label column
{
  res <- read.delim(filename, ..., colClasses = colClasses)
  colnames(res)[1] <- "label"
  # Order by file order
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
           small.cex = 1.8,
           large.cex = 2.0,
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
           ..., header = TRUE,
           colClasses = list(X = "character")) # Label column
{
  res <- read.delim(filename, ..., header = header, colClasses = colClasses)
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
  function(x, y, reference = 0.01, colors = NULL, ...)
{

  ref.x <- log10(reference)
  panel.refline(v = ref.x)
  #panel.axis(side = "bottom", at = c(ref.x), outside = TRUE,
  #           labels = c(paste("p = ", as.character(reference))))
  special = is.infinite(x)
  if (!is.null(colors)) {
    col = colors[special]
  } else {
    col = NA
  }
  x[special] <- min(x[!special]) * 2
  #panel.abline(x = c(-10, 1), h = y[special], col = col, lwd = 2, lty = 3, ...)
  
  panel.barchart(x, y, stacked = FALSE, col = colors, ...)
}

barchart.pvalues <-
  function(data,
           x = reorder(label, -value) ~ value,
           as.table = TRUE,
           main = "Approximate significance of overlap",
           panel = panel.pvalues,
           xlab = "p-value",
           ylab = "Segment label",
           reference = 0.01,
           origin = 0,
           auto.key = if (nlevels(data$variable) > 1) { TRUE } else { FALSE },
           colors = label.colors(data$label),
           scales = list(x = list(log = TRUE)),
           ...)
{
  data.melted <- melt(data, id.vars="label")
  xyplot(x = x,
         data = data.melted,
         group = variable,
         panel = panel,
         as.table = as.table,
         reference = reference,
         main = main,
         colors = colors,
         scales = scales,
         origin = origin,
         auto.key = auto.key,
         xlab = xlab,
         ylab = ylab,
         #par.settings = list(clip = list(panel = "off")),
         ...)
}

plot.overlap.pvalues <- 
  function(tabfile, mnemonics = NULL, ...)
{
  pvalues <- read.pvalues(tabfile, mnemonics = mnemonics)
  barchart.pvalues(pvalues, ...)
}
