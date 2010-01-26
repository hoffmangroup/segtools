# requires common.R

library(lattice)
library(latticeExtra)
library(grid)

read.signal <- function(filename, mnemonics = NULL, comment.char = "#", ...,
                        colClasses = list(label = "character")) {
  comment.lines <- read.comment.lines(filename, comment.char = comment.char)
  res <- read.delim(filename, ..., colClasses = colClasses)
  res$label <- relevel.mnemonics(factor(res$label), mnemonics)

  res
}

panel.histogram.ecdf <- function(x, y, ref = FALSE, ...) {
  panel.ecdfplot(y, ref = FALSE, ...)
}

panel.histogram.precomp <- function(x, y, ...) {
  box.width <- diff(x[1:2])
  x <- x + (box.width/2)
  panel.barchart(x, y, box.width = box.width, ...)
}

panel.histogram.normalized <- function(x, y, type = NULL, ...) {
  keep <- (is.finite(y) & y > 0)
  panel.xyplot(x[keep], y[keep], type = c("l", "h"), ...)
}

strip.highlight.segtracks <- function(segtracks, horizontal.outer = FALSE, ...) {
  # Define custom strip function that highlights strips
  # when the trackname is in 'segtracks'
  strip.highlighter <- function(which.given,
                                which.panel,
                                factor.levels,
                                ...,
                                horizontal = NULL,  # Ignore
                                bg = NULL) {

    if (factor.levels[which.panel[which.given]] %in% segtracks) {
      strip.default(which.given = which.given, which.panel = which.panel,
	            factor.levels = factor.levels, ..., 
                    horizontal = horizontal.outer, bg = "yellow")
    } else {
      strip.default(which.given = which.given, which.panel = which.panel,
      		    horizontal = horizontal.outer, 
                    factor.levels = factor.levels, ...)
    }
  }

  strip.highlighter
}

## Generates pretty scales for the data.
## layout is a 2-element vector: c(num_rows, num_cols) for the xyplot
## tolerance - fraction of extreme value allowed to not be included in plot.
panel.scales <- function(data, layout, log.y = TRUE, mode = "histogram",
                         tolerance = 0.05, ...) {

  scale.y = list(log = log.y,
    tck = c(0, 1),
    alternating = c(2, 0))
  
  ## Get the x and y axis limits depending on plot mode
  if (mode == "ecdf") {
    limits.x <- list()
    at.x <- list()
    for (track in levels(data$trackname)) {
      track.subset <- subset(data, trackname == track)
      max.x <- max(track.subset$lower_edge, na.rm = TRUE)
      limits.x <- c(limits.x, list(c(0, max.x)))
      at.x <- c(at.x, at.pretty(from = 0, to = max.x, length = 5, largest = TRUE))
    }

    scale.y$labels <- c("0", "1")
    scale.y$limits <- c(0, 1)
    scale.y$at <- c(0, 1)
    log.y <- FALSE
  } else if (mode == "histogram") {
    limits.x <- list()
    at.x <- list()
    for (track in levels(data$trackname)) {
      track.subset <- subset(data, trackname == track)
      max.x <- max(track.subset$lower_edge, na.rm = TRUE)
      limits.x <- c(limits.x, list(c(0, max.x)))
      at.x <- c(at.x, at.pretty(from = 0, to = max.x, length = 5, largest = TRUE))
    }
    
    max.y <- max(data$count, na.rm = TRUE)
    high.exp <- floor(log10(max.y))
    min.y <- if (log.y) 1 else 0
    at.y <- c(min.y, 10^high.exp)
    
    scale.y$labels <- c(as.character(min.y), label.log(10^high.exp))
    scale.y$limits <- c(head(at.y, n = 1), tail(at.y, n = 1))
    scale.y$at <- at.y
  } else if (mode == "normalized") {
    limits.x <- list()
    at.x <- list()
    for (track in levels(data$trackname)) {
      track.subset <- subset(data, trackname == track)
      count.subset <- subset(track.subset, is.finite(count) & count > 0)
      if (nrow(count.subset) == 0) {
        min.x <- min(track.subset$lower_edge)
        max.x <- max(track.subset$lower_edge)
      } else {
        ##min.y <- quantile(count.subset$count, probs = tolerance, na.rm = TRUE)
        min.y <- max(count.subset$count, na.rm = TRUE) * tolerance
        within.tolerance <- subset(count.subset, count >= min.y)
        min.x <- min(within.tolerance$lower_edge)
        max.x <- max(within.tolerance$lower_edge)
      }
      limits.x <- c(limits.x, list(c(min.x, max.x)))
      at.x.pretty <- unlist(at.pretty(from = min.x, to = max.x, length = 10))
      at.x <- c(at.x, list(c(0, tail(at.x.pretty, n = 1))))
      ## print(paste(max(count.subset$lower_edge), " ---> ", max.x))
    }
    scale.y$log <- FALSE
  }

  scales = list(x = list(relation = "free",  # Allow x axes to vary
                  at = c(rep(list(NULL), layout[2] * (layout[1] - 1)), at.x),
                  limits = limits.x,
                  rot = 90),
    y = scale.y,
    cex = 1,
    ...)
  
  scales
}

sum.nona <- function(..., na.rm = TRUE) {
  sum(..., na.rm = na.rm)
}

# Generates a histogram from the precomputed histogram in 'data'
xyplot.signal <-
  function(data,
           x = count ~ lower_edge | trackname + label,  # formula, explained in the lattice book
                # ~ means "is modeled by", as in y ~ x
                # So count is the y axis,  lower_edge is the x axis
                # | means conditioned on
                # So it makes a plot for every label and trackname
           segtracks = NULL,  # List of tracks in segmentation
           mode = "histogram",  # Mode determines display mode
                                # Choice of histogram, ecdf, normalized
           text.cex = 1,
           par.settings = list(add.text = list(cex = text.cex),
                               layout.heights = list(axis.panel = axis.panel,
                                                     strip = strips.heights),
                               layout.widths = list(strip.left = strips.widths)
                               ),
           xlab = "Signal Value",
           ylab = if (mode == "ecdf") "Empirical cumulative density"
                  else if (mode == "normalized") "Fraction of bases"
                  else if (log.y) "Number of bases (log10)"
                  else "Number of bases",
           log.y = FALSE,
           panel = if (mode == "histogram") panel.histogram.precomp
                   else if (mode == "ecdf") panel.histogram.ecdf
                   else if (mode == "normalized") panel.histogram.normalized
                   else stop(paste("Unrecognized mode: ", mode)),
           strip = strip.custom(horizontal = FALSE),
           strip.left = strip.custom(horizontal = TRUE),
           strip.height = 14,
           strip.width = 3,
           border = "transparent",
           col = trellis.par.get("superpose.symbol")$col[1],
           ...)
{
  # defining a new lattice function: see D. Sarkar, _Lattice:
  # Multivariate Data Visualization with R,_ sec. 14.3.1;
  # https://stat.ethz.ch/pipermail/r-help/2008-December/182370.html

  if (length(segtracks) > 0) {
    # Put all seg tracks on one side
    tracks.levels <- levels(data$trackname)
    tracks.highlight <- tracks.levels %in% segtracks
    tracks.ordered <- c(tracks.levels[tracks.highlight],
                        tracks.levels[ ! tracks.highlight])
    data$trackname <- factor(data$trackname, levels = tracks.ordered)

    # Highlight seg track strips instead
    strip <- strip.highlight.segtracks(segtracks)
  }

  if (mode == "normalized") {
    group.sums <- with(data,
                       aggregate(count, list(trackname, label), sum.nona))$x
    group.size <- with(data,
                       aggregate(count, list(trackname, label), length))$x
    data.divisor <- rep(group.sums, times = group.size)
    zero.out <- cumsum(c(1, group.size[-length(group.size)]))
    data.divisor[zero.out] <- NA
    data$count <- data$count / data.divisor
  }

  # this computes the sum of the counts, for the "all" label
  if (mode == "histogram" || mode == "ecdf") {
    labelsums <- with(data, aggregate(count, list(trackname, lower_edge), sum))
    names(labelsums) <- c("trackname", "lower_edge", "count")
    labelsums$label <- factor("all")
    data.full <- rbind(data, labelsums)
    data.full$label <- relevel(data.full$label, ref = "all") 
  } else {
    data.full <- data
  }

  num_rows <- length(levels(data.full$label))
  num_cols <- length(levels(data.full$trackname))
  layout <- c(num_rows, num_cols)
  # Remove space between panels (where axes were)
  axis.panel <- rep(c(0, 1), c(num_rows - 1, 1))
  # Make strips wider and longer to fit full text
  strips.heights <- rep(c(strip.height, 0), c(1, num_rows - 1))
  strips.widths <- rep(c(strip.width, 0), c(1, num_cols - 1))
  
  scales <- panel.scales(data.full, layout, log.y = log.y, mode = mode)

  trellis <- xyplot(x = x,
                    data = data.full,
                    scales = scales,      
                    as.table = TRUE,
                    xlab = xlab,
                    ylab = ylab,
                    panel = panel,
                    border = border,
                    horizontal = FALSE,
                    col = col,
                    ...)

  trellis.outer <- useOuterStrips(trellis, strip = strip,
                                  strip.left = strip.left)

  update(trellis.outer, par.settings = par.settings, evaluate = FALSE)
}


plot.signal <- function(filename, segtracks = NULL, mnemonics = NULL, ...) {
  data <- read.signal(filename, mnemonics = mnemonics)

  xyplot.signal(data = data, segtrack = segtracks, ...)
}
