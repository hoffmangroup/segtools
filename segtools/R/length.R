read.length <- function(filename, mnemonics = NULL, ...) 
{
  res <- read.delim(filename, ...)

  res$label <- factor(res$label)
  
  if (length(mnemonics) > 0) {
    res$label <- relabel.factor(res$label, mnemonics)
  } else {
    res$label <- smart.factor.reorder(res$label)
  }

  # Reverse rows for plot
  res$label <- factor(res$label, levels = rev(levels(res$label)))

  res
}

expected.len <- function(prob.loop) {
  prob.loop / (1 - prob.loop)
}

panel.violinplot <-
  function(x, y, ..., box.ratio, mark.prior = NULL,
           marks.posterior = NULL, rug = FALSE)
{
  if (!is.null(mark.prior)) {
    log10.mark.prior <- log10(mark.prior)
    panel.refline(v = log10.mark.prior, h = 0)
  }

  # In order to do a violin-plot, each label must have at least
  # two data points. Ignore any labels with fewer data points.
  keep_levels = c()
  for (lev in levels(y)) {
    if (sum(y == lev) > 1) {
      # Keep entries for level if enough data points
      keep_levels = c(keep_levels, lev)
    }
  }
  keep = is.element(y, keep_levels)

  # Violin-plot of non-unique label entries
  panel.violin(x[keep], y[keep], ..., box.ratio = box.ratio)

  # Add remaining plots
  panel.bwplot(x, y, ..., box.ratio = 0)
  if (rug) {
    panel.rug(x, NULL, ..., end = 0.01)
  }
  if (!is.null(marks.posterior)) {
    log10.marks.posterior <- log10(marks.posterior)
    panel.points(log10.marks.posterior, seq_along(marks.posterior),
                 pch = "|", col = "black", cex = 3)
  }
}

# data: a data.frame with a row for each segment
#   label: (factor) the segment label
#   length: (integer) the length of the segment
# mnemonics: a matrix where [,1] is old labels and [,2] is new labels
violinplot.length <-
  function(data, 
           x = label ~ length,
           nint = 100,
           xlab = "Length", ylab = "Label",
           scales = list(x = list(log = TRUE,
                                  at = at.log(),
                                  labels = labels.log() )),
           panel = panel.violinplot, 
           expected.lens = NULL)
{
  marks.posterior <-
    if (is.null(expected.lens)) {
      NULL
    } else {
      expected.len(expected.lens)
    }

  bwplot(x, data, nint = nint, 
         xlab = xlab, ylab = ylab,
         scales = scales, 
         panel = panel,
         marks.posterior = marks.posterior)
}

plot.length <- function(filename, mnemonics = NULL) {
  data <- read.length(filename, mnemonics = mnemonics)

  violinplot.length(data = data)
}


############### SEGMENT SIZES ###############

read.segment.sizes <-
  function(filename, ..., check.names = FALSE)
{
  res.all <- read.delim(filename, ..., check.names = check.names)

  res <- res.all[res.all$label != "all",]

  res$label <- factor(res$label)
  res$frac.segs <- res$num.segs / sum(res$num.segs)

  res
}

panel.segment.sizes <- 
  function(...)
{
  panel.barchart(..., origin = 0)
}

barchart.segment.sizes <- 
  function(data, 
           x = reorder(label, -(frac.segs + frac.bp)) ~ frac.bp + frac.segs,
           main = NULL, #"Segment sizes",
           xlab = "Fraction of segmentation",
           ylab = "Segment label",
           panel = panel.segment.sizes,
           as.table = TRUE,
           xlim = c(0, max.frac),
           col = c("red", "blue"),
           auto.key = list(space = "top", 
                           columns = 2, 
                           text = c("Bases", "Segments"),
                           rectangles = TRUE,
                           points = FALSE),
           par.settings = list(superpose.polygon = list(col = col)),
           ...)
{
  
  max.frac <- min(c(1, 1.05 * max(c(data$frac.bp, data$frac.segs))))

  barchart(x = x, 
           data = data,
           main = main,
           xlab = xlab, 
           ylab = ylab,
           panel = panel,
           as.table = as.table,
           xlim = xlim,
           auto.key = auto.key,
           par.settings = par.settings,
           ...)
}

plot.segment.sizes <-
  function(filename, ...)
{
  data <- read.segment.sizes(filename)
  barchart.segment.sizes(data = data, ...)
}
