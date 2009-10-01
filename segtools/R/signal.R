# requires common.R

library(lattice)
library(latticeExtra)

read.signal <- function(filename, mnemonics=NULL, ...) {
  res <- read.delim(filename, ...)
  res$label <- factor(res$label)    # adds a label attribute

  # Replace labels with mnemonics, if supplied
  if (length(mnemonics) > 0) {
    res$label <- relabel.factor(res$label, mnemonics)
  } else {
    res$label <- smart.factor.reorder(res$label)
  }

  res
}

panel.histogram.precomp <- function(x, y, ecdf = FALSE, ref = FALSE, ...) {
  if (ecdf) {
    panel.ecdfplot(y, ref = FALSE, ...)
  } else {
    box.width <- diff(x[1:2])
    x <- x + (box.width/2)
    panel.barchart(x, y, box.width = box.width, ...)
  }
}

strip.highlight.segtracks <- function(segtracks, horizontal.outer=FALSE, ...) {
  # Define custom strip function that highlights strips
  # when the trackname is in 'segtracks'
  strip.highlighter <- function(which.given,
                                which.panel,
                                factor.levels,
                                ...,
                                horizontal = NULL,  # Ignore
                                bg = NULL) {

    if (factor.levels[which.panel[which.given]] %in% segtracks) {
      strip.default(which.given=which.given, which.panel=which.panel,
	            factor.levels=factor.levels, ..., 
                    horizontal=horizontal.outer, bg="yellow")
    } else {
      strip.default(which.given=which.given, which.panel=which.panel,
      		    horizontal=horizontal.outer, 
                    factor.levels=factor.levels, ...)
    }
  }

  strip.highlighter
}


## Generates pretty scales for the data.
## layout is a 2-element vector: c(num_rows, num_cols) for the xyplot
panel.scales <- function(data, layout, log.y = TRUE, ecdf = FALSE, ...) {
  # Get the max x for each track
  limits.x <- list()
  at.x <- list()
  for (track in levels(data$trackname)) {
    track_subset <- subset(data, data$trackname == track)
    max.x <- max(track_subset$lower_edge, na.rm=TRUE)
    limits.x <- c(limits.x, list(c(0, max.x)))
    at.x <- c(at.x, at.pretty(from=0, to=max.x, length=5, largest=TRUE))
  }

  if (ecdf) {
    at.y <- c(0, 1)
    labels.y <- c("0", "1")
    log.y <- FALSE
    limits.y <- c(0, 1)
  } else {
    max.y <- max(data$count, na.rm=TRUE)
    high.exp <- floor(log10(max.y))
    min.y <- if (log.y) 1 else 0
    at.y <- c(min.y, 10^high.exp)
    labels.y <- c(as.character(min.y), label.log(10^high.exp))
    limits.y <- c(head(at.y, n = 1), tail(at.y, n = 1))
  }

  scales = list(x = list(relation = "free",  # Allow x axes to vary
                         limits = limits.x,
                         at = c(rep(list(NULL), layout[2] * (layout[1] - 1)), 
                                at.x),
                         rot = 90
                         ),
                y = list(log = log.y,
                         tck = c(0, 1),
                         labels = labels.y,
                         limits = limits.y,
                         at = at.y,
                         alternating = c(2, 0)
                         ),
                cex = 1,
                ...)
  
  scales
}


# Generates a histogram from the precomputed histogram in 'data'
xyplot.signal <-
  function(x = count ~ lower_edge | trackname + label,  # formula, explained in the lattice book
                # ~ means "is modeled by", as in y ~ x
                # So count is the y axis,  lower_edge is the x axis
                # | means conditioned on
                # So it makes a plot for every label and trackname
           data,
           segtracks = NULL,  # List of tracks in segmentation
           ecdf = FALSE,  # Plot an ecdf for each packet instead of a histogram
           text.cex = 1,
           par.settings = list(add.text = list(cex = text.cex),
                               layout.heights = list(axis.panel = axis.panel,
                                                     strip = strips.heights),
                               layout.widths = list(strip.left = strips.widths)
                               ),
           xlab = "Signal Value",
           ylab = if (!ecdf) { if (log.y) "Counts (log10)" else "Counts" }
                  else "Empirical cumulative density",
           log.y = TRUE,
           panel = panel.histogram.precomp, 
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

  # this computes the sum of the counts, for the "all" label
  labelsums <- with(data, aggregate(count, list(trackname = trackname,
                                                lower_edge = lower_edge),
                                    sum))
  names(labelsums) <- sub("^x$", "count", names(labelsums))
  labelsums$label <- "all"

  #data.full <- data
  data.full <- rbind(data, labelsums)
  data.full$label <- relevel(data.full$label, ref = "all") 

  num_rows <- length(levels(data.full$label))
  num_cols <- length(levels(data.full$trackname))
  layout <- c(num_rows, num_cols)

  # Can use just labelsums, since it should cover the same range
  scales <- panel.scales(labelsums, layout, log.y = log.y, ecdf = ecdf)
  #scales <- panel.scales(data.full, layout, log.y = log.y)

  # Remove space between panels (where axes were)
  axis.panel <- rep(c(0, 1), c(num_rows - 1, 1))
  # Make strips wider and longer to fit full text
  strips.heights <- rep(c(strip.height, 0), c(1, num_rows - 1))
  strips.widths <- rep(c(strip.width, 0), c(1, num_cols - 1))

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
                    ecdf = ecdf, 
                    ...)

  trellis.outer <- useOuterStrips(trellis, strip = strip,
                                  strip.left = strip.left)

  update(trellis.outer, par.settings = par.settings, evaluate = FALSE)
}


plot.signal <- function(filename, segtracks=NULL, mnemonics=NULL, ...) {
  data <- read.signal(filename, mnemonics=mnemonics)

  xyplot.signal(data=data, segtrack=segtracks, ...)
}