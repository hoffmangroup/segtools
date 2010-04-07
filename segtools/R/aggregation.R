library(lattice)
library(RColorBrewer)
library(latticeExtra)
library(plyr)
library(reshape)

COLNAMES <- c("group", "component", "offset")
lattice.options(panel.error="stop")

read.aggregation <- function(filename, mnemonics = NULL, ..., 
                             comment.char = "#",
                             check.names = FALSE) {
  if (!file.exists(filename)) {
    stop(paste("Error: could not find aggregation file:", filename))
  }
  data.raw <- read.delim(filename, ..., 
                         comment.char=comment.char, 
                         check.names=check.names)
  colnames(data.raw)[1] <- "group"
  if (!all(colnames(data.raw)[1:3] == COLNAMES)) {
    stop(paste("Error: unsupported aggregation file format.",
               "Expected first three columns to be:",
               paste(COLNAMES, collapse=", ")))
  }

  data <- melt.aggregation(data.raw)
  
  data$label <- relevel.mnemonics(data$label, mnemonics)

  ## Order components by the order observed in the file
  data$component <- factor(data$component, levels=unique(data$component))
  
  data
}


## Generates pretty scales for the data.
## layout is a 2-element vector: c(num_rows, num_cols) for the xyplot
## num_panels is the number of panels/packets in the trellis
panel.scales <- function(data, layout, num_panels, x.axis = FALSE) {
  components <- levels(data$component)
  num_components <- length(components)

  ## Avoid overlapping scales if there is not an even row at the bottom
  remove.extra.scales <- (layout[1] * layout[2] != num_panels) * num_components

  ## Figure out x axis labels (should be same within component)
  at.x <- list()
  limits.x <- list()
  for (cur_component in components) {
    component_subset <- subset(data, component == cur_component)
    min.x <- min(component_subset$offset, na.rm=TRUE)
    max.x <- max(component_subset$offset, na.rm=TRUE)
    at.x.pretty <- at.pretty(from=min.x, to=max.x, length=5, largest=TRUE)
    at.x <- c(at.x, at.x.pretty)
    ## Cludgy shrink of axes to compensate for automatic trellis expansion
    ## We want the x axis to go to the edge exactly, and axs="i" doesn't
    ## do it.
    limits.x[[length(limits.x) + 1]] <-
      extendrange(c(min.x, max.x), f = -0.05)
  }

  at.x.full <-
    if (x.axis) {
      ## Remove internal axes and space where axes were
      at.x.nonnull.full <- rep(at.x,
                               as.integer((layout[1] - remove.extra.scales) /
                                          num_components))
      c(rep(list(NULL), num_panels - layout[1] + remove.extra.scales),
        at.x.nonnull.full)
    } else {
      NULL
    }

  range.y <- range(data$overlap, na.rm=TRUE)
  min.y <- min(range.y[1], 0)
  max.y <- max(range.y[2], 0)
  limits.y <- extendrange(c(min.y, max.y))
  at.y <- unique(round(c(min.y, 0, max.y), digits = 2))
  scales <- list(x = list(relation = "free",
                          tck = c(1, 0),
                          at = at.x.full,
                          rot = 90,
                          limits = limits.x),
                 y = list(alternating = c(2, 0),
                          tck = c(0, 1),
                          limits = limits.y,
                          at = at.y))

  scales
}

transpose.levels <- function(data.group, dim.length) {
  lev <- levels(data.group)
  nlev <- nlevels(data.group)
  res <- vector(mode = "character", length = nlev)
  num.added <- 0
  for (i in seq(from = 1, length = dim.length)) {
    dest.indices <- seq(from = i, to = nlev, by = dim.length)
    source.indices <- seq(from = num.added + 1, along = dest.indices)
    res[dest.indices] <- lev[source.indices]
    num.added <- num.added + length(dest.indices)
  }
  
  factor(data.group, levels = res)
}

melt.aggregation <- function(data.raw) {
  id.vars <- colnames(data.raw)[1:3]
  data <- melt(data.raw, id.vars = id.vars)
  colnames(data)[(colnames(data) == "variable")] <- "label"
  colnames(data)[(colnames(data) == "value")] <- "count"

  data$group <- factor(data$group)
  data$component <- factor(data$component)
  data$label <- factor(data$label)

  data
}

cast.aggregation <- function(data.df) {
  cast(data.df, group + component + offset ~ label, value = "count")
}

calc.signif <- function(count, total, random.prob) {
  ## Calculate "significance" of count and total given
  ## current random.prob.
  if (!is.finite(count) || !is.finite(total)) {
    stop("Found non-finite count (", count,
         ") or total (", total,
         ") for label: ", label,
         sep = "")
  }
  
  expected <- total * random.prob
  sep <- abs(expected - count)
  
  signif <- 
    if (total == 0) {
      ## No significance when there are no observed overlaps
      1
    } else if (count > 100 && expected < 10) {
      ## Calculate two-tailed probability 
      ppois(expected - sep, expected) + 
        ppois(expected + sep - 1, expected, lower.tail = FALSE)
    } else {
      ## Binomial is symmetric, so probability is two-tailed with * 2
      if (count < expected) {
        pbinom(count, total, random.prob) * 2
      } else {
        pbinom(count - 1, total, random.prob, lower.tail = FALSE) * 2
      }
    }
  signif <- min(signif, 1)
##   signif <- -log10(signif)
##   signif[!is.finite(signif)] <- Inf
##   signif[signif > 20] <- 20
##   signif <- if(count < expected) -signif else signif
  
  signif
}

## Given an aggregation data frame, process the counts by optionally
## normalizing the counts over all labels and by the sizes of the labels
## as well as calculating a significance level for each row.
process.counts <- function(data, label.sizes, pseudocount = 1,
                           normalize = TRUE, pval.thresh = 0.001) {
  if (!is.vector(label.sizes)) stop("Expected vector of label sizes")
  total.size <- sum(label.sizes)
  if (length(total.size) != 1) stop("Error summing label size vector")
  if (!is.finite(total.size)) stop("Unexpected sum of sizes (not finite)")

  ##data.mat <- cast.aggregation(data)
  
  ## Ideally, a multinomial test would be applied once for each bin
  ## (testing whether the distribution of overlaps by label is as expected),
  ## but the label-wise binomial seems a decent approximation.
  labels.sum <- with(data, aggregate(count, list(offset, component),
                               sum))$x
  for (label in levels(data$label)) {
    random.prob <- label.sizes[label] / total.size
    cur.rows <- data$label == label
    cur.counts <- data$count[cur.rows]

    calc.row.signif <- function(row) {
      count <- row[1]
      total <- row[2]
      calc.signif(count, total, random.prob)
    }

    if (normalize) {
      data$count[cur.rows] <- log2((cur.counts / labels.sum + 1) /
                                   (random.prob + 1))
    }
    pvals <- apply(cbind(cur.counts, labels.sum), 1, calc.row.signif)
    data$significant[cur.rows] <- (is.finite(pvals) & (pvals < pval.thresh))
    ##data$count[cur.rows] <- -log10(pvals)
  }
  
  data
}


panel.aggregation <- function(x, y, significant, ngroups, groups = NULL,
                              subscripts = NULL, font = NULL, col = NULL,
                              col.line = NULL, pch = NULL,
                              group.number = NULL, ...) {
  ## hide 'font' from panel.segments, since it doesn't like having
  ## font and fontface
  panel.refline(h = 0)

  significant <- as.logical(significant)[subscripts]
  significant[!is.finite(significant)] <- FALSE
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (any(significant)) {
    ## Only shade region for first.
    fill.col <- "black"
    #fill.col <- if (ngroups == 1) "black" else col.line
    if (ngroups == 1) {
      y.sig <- y
      y.sig[!significant] <- 0
      panel.polygon(cbind(c(min(x), x, max(x)), c(0, y.sig, 0)),
                    col = fill.col)
    } else {
      x.sig <- x[significant]
      y.sig <- y[significant]
      panel.points(x.sig, y.sig, col = fill.col, pch = "*")
    }
  }
##  panel.xyplot(x, y, groups = groups, subscripts = subscripts, ...)
  panel.xyplot(x, y, font = font, col = col, col.line = col.line,
               pch = pch, ...)
}

get.metadata.label.sizes <- function(metadata, data) {
  label.sizes <- NULL
  for (label in levels(data$label)) {
    label.size.raw <- metadata[[as.character(label)]]
    label.size <- as.numeric(label.size.raw)
    if (length(label.size) == 0 || !is.finite(label.size)) {
      stop(paste("Error: encountered invalid size for label:", label,
                 paste("(", label.size.raw,")", sep = "")))
    }
    label.sizes[as.character(label)] <- label.size
  }
  label.sizes
}

## Plots overlap vs position for each label
##   data: a data frame with fields: overlap, offset, label
##   spacers should be a vector of indices, where a spacer will be placed after
##     that component (e.g. c(1, 3) will place spacers after the first and third
##     components
xyplot.aggregation <- function(data,
    metadata = NULL,
    x = overlap ~ offset | component * label,
    spacers = metadata[["spacers"]],
    normalize = TRUE,
    x.axis = FALSE,  # Show x axis
    pval.thresh = 0.001,
    text.cex = 1,
    spacing.x = 0.4,
    spacing.y = 0.4,
    ngroups = nlevels(data$group),      
    panel = panel.aggregation,
    par.settings = list(add.text = list(cex = text.cex),
                        layout.heights = list(axis.panel = axis.panel,
                                              strip = strips.heights)),
    auto.key = if (ngroups < 2) FALSE
               else list(points = FALSE, lines = TRUE),
    strip = strip.custom(horizontal = FALSE),
    strip.height = 10,
    xlab = NULL,
    ylab = if (normalize) "Enrichment {log2[(fObs + 1)/(fRand + 1)]}"
           else "Count",
    sub = if (normalize) paste("Black regions are significant with p<",
                               pval.thresh, sep = "")
          else NULL,
    ...)
{
  metadata <- as.list(metadata)
  ## Normalize and/or calculate significance if metadata exists
  if (length(metadata) > 0) {
    label.sizes <- get.metadata.label.sizes(metadata, data)
    data <- process.counts(data, label.sizes, pval.thresh = pval.thresh,
                             normalize = normalize)
  } else {
    stop("Cannot normalize/calculate significance without metadata")
  }

  colnames(data) <- gsub("^count$", "overlap", colnames(data))
  data$overlap[!is.finite(data$overlap)] <- 0

  ## Determine panel layout
  num_levels <- nlevels(data$label)
  num_components <- nlevels(data$component)
  num_rows <- num_components
  num_cols <- num_levels
  num_panels <- num_rows * num_cols

  ## Rework layout to optimize organization
  while (num_cols > num_rows) {
    num_cols <- num_cols / 2
    num_rows <- num_rows * 2
  }
  num_cols <- ceiling(num_cols)
  layout <- c(num_rows, num_cols)

  ## Reorder labels so they are in order downward in panels
  data$label <- transpose.levels(data$label, num_rows / num_components)

  ## Separate distinct groups
  spaces.x <- rep(0, num_components - 1)
  if (is.numeric(spacers) && length(spacers) > 0) {
    if (any(spacers < 1 | spacers >= num_components)) {
      stop("Spacer vector should only contain values in the range [1, ",
           num_components - 1,
           "] since there are ", num_components, " components")
    }
    spaces.x[spacers] <- spacing.x
  }
  between <- list(x = c(spaces.x, spacing.x), 
                  y = spacing.y)

  scales <- panel.scales(data, layout, num_panels, x.axis = x.axis)
  axis.panel <- rep(c(0, 1), c(num_cols - 1, 1))

  ## Make top strips longer
  strips.heights <- rep(c(strip.height, 0), c(1, num_cols - 1))

  args <- list(x, data = data, type = "l", groups = quote(group),
               auto.key = auto.key, as.table = TRUE, strip = strip,
               xlab = xlab, ylab = ylab, sub = sub,
               significant = data$significant, ngroups = ngroups,
               panel = "panel.superpose", panel.groups = panel, ...)

  trellis.raw <- do.call(xyplot, args)

  trellis <- useOuterStrips(trellis.raw, strip = strip)

  update(trellis, layout = layout, between = between, scales = scales,
         par.settings = par.settings)
}

plot.aggregation <- function(filename, mnemonics = NULL, ...,
                             comment.char = "#") {
  data <- read.aggregation(filename, mnemonics = mnemonics)
  metadata <- read.metadata(filename, comment.char = comment.char)
  # Rename metadata keys with mnemonics
  names(metadata) <- map.mnemonics(names(metadata), mnemonics)$labels

  xyplot.aggregation(data = data, metadata = metadata, ...)
}

save.aggregation <- function(dirpath, namebase, tabfilename,
                             mnemonic_file = NULL,
                             clobber = FALSE,
                             panel.size = 200,  # px
                             comment.char = "#",
                             ...) {
  mnemonics <- read.mnemonics(mnemonic_file)
  data <- read.aggregation(tabfilename, mnemonics = mnemonics)
  metadata <- read.metadata(tabfilename, comment.char = comment.char)
  # Rename metadata keys with mnemonics
  names(metadata) <- map.mnemonics(names(metadata), mnemonics)$labels
  
  image.size <- panel.size * ceiling(sqrt(nlevels(data$label) / 2))
  save.images(dirpath, namebase,
              xyplot.aggregation(data = data, metadata = metadata, ...),
              width = image.size,
              height = image.size,
              clobber = clobber)
}
