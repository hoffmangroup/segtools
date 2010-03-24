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
  for (cur_component in components) {
    component_subset <- subset(data, component == cur_component)
    min.x <- min(component_subset$offset, na.rm=TRUE)
    max.x <- max(component_subset$offset, na.rm=TRUE)
    at.x.pretty <- at.pretty(from=min.x, to=max.x, length=5, largest=TRUE)
    at.x <- c(at.x, at.x.pretty)
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
                          axs = "i"
                          ),
                 y = list(alternating = c(2, 0),
                          tck = c(0, 1),
                          limits = limits.y,
                          at = at.y
                          )
                 #cex = 0.7
                 )

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

## Given an aggregation data frame, normalize the counts over all labels
## and by the sizes of the labels (to calculate enrichment
normalize.counts <- function(data, label.sizes, pseudocount = 1,
                             pval.thresh = 0.01) {
  if (!is.vector(label.sizes)) stop("Expected vector of label sizes")
  total.size <- sum(label.sizes)
  if (length(total.size) != 1) stop("Error summing label size vector")
  if (!is.finite(total.size)) stop("Unexpected sum of sizes (not finite)")

  data.mat <- cast.aggregation(data)
  labels.sum <- with(data, aggregate(count, list(offset, component),
                               sum))$x
  for (label in levels(data$label)) {
    random.prob <- label.sizes[label] / total.size
    cur.rows <- data$label == label
    cur.counts <- data$count[cur.rows]

    calc.pval <- function(row) {
      count <- row[1]
      total <- row[2]
      lambda <- total * random.prob
      if (count > 100 && lambda < 10) {
        ppois(count, lambda)
      } else {
        pbinom(count, total, random.prob)
      }
    }
    
    pvals <- apply(cbind(cur.counts, labels.sum), 1, calc.pval)
    data$count[cur.rows] <- log2((cur.counts / labels.sum + 1) /
                                 (random.prob + 1))
    data$significant[cur.rows] <- (pvals < pval.thresh | pvals > (1 - pval.thresh))
  }
  data
##   ## Convert data to matrix, row-normalize, and then convert back
##   data.mat <- cast.aggregation(data)
##   data.values <- data.mat[, 4:ncol(data.mat)]
##   group.counts <- rowSums(data.values)
##   ## Calculate enrichment above random probability
##   for (col in 4:ncol(data.mat)) {
##     label.size <- label.sizes[colnames(data.mat)[col]]
##     random.prob <- label.size / total.size
##     probs <- pbinom(data.mat[, col], group.counts, random.prob)
##     significant <- probs < pval.thresh | probs > (1 - pval.thresh)
##     data.mat[, col] <- log2((data.mat[, col] + pseudocount) /
##                             (random.prob + pseudocount))
##   }
##   ## Normalize over labels (cols)
##   data.mat[, 4:ncol(data.mat)] <- data.values / group.counts
##   data <- melt.aggregation(data.mat)
}


panel.aggregation <- function(x, y, significant, ngroups, groups = NULL, subscripts = NULL,
                              font = NULL, col = NULL, col.line = NULL, ...) {
  ## hide 'font' from panel.segments, since it doesn't like having
  ## font and fontface
  panel.refline(h = 0)

  significant <- as.logical(significant)[subscripts]
  ## Only shade region for first
  if (any(significant)) {
    fill.col <-
      if (ngroups == 1) "black"
      else col.line
    x.sig <- as.numeric(x)[significant]
    y.sig <- as.numeric(y)[significant]
    panel.segments(x.sig, y.sig, x.sig, 0, col = fill.col, ...)
  }
#  panel.xyplot(x, y, groups = groups, subscripts = subscripts, ...)
  panel.xyplot(x, y, font = font, col = col, col.line = col.line, ...)
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
    sub = if (normalize) paste("Shaded regions are significant with p<",
                               pval.thresh, sep = "")
          else NULL,
    ...)
{
  metadata <- as.list(metadata)
  ## Normalize by number of segments in each labels
  if (normalize) {
    if (length(metadata) > 0) {
      label.sizes <- get.metadata.label.sizes(metadata, data)
      data <- normalize.counts(data, label.sizes, pval.thresh = pval.thresh)
    } else {
      stop("Cannot normalize y-axis without metadata")
    }
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
