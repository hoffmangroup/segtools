library(lattice)
library(RColorBrewer)
library(latticeExtra)

read.aggregation <- function(filename, mnemonics=NULL, ..., 
                             comment.char = "#",
                             check.names=FALSE) {
  data <- read.delim(filename, ..., 
                     comment.char=comment.char, 
                     check.names=check.names)
  data$label <- factor(data$label)
  data$factor <- factor(data$factor)
  data$component <- factor(data$component)
  
  if (length(mnemonics) > 0) {
    data$label <- relabel.factor(data$label, mnemonics)
  } else {
    data$label <- smart.factor.reorder(data$label)
  }

  # Order components by the order observed in the file
  data$component <- factor(data$component, levels=unique(data$component))
  
  data
}


## Generates pretty scales for the data.
## layout is a 2-element vector: c(num_rows, num_cols) for the xyplot
## num_panels is the number of panels/packets in the trellis
panel.scales <- function(data, layout, num_panels) 
{
  components <- levels(data$component)
  num_components <- length(components)

  # Avoid overlapping scales if there is not an even row at the bottom
  remove.extra.scales <- (layout[1] * layout[2] != num_panels) * num_components

  # Figure out x axis labels (should be same within component)
  at.x <- list()
  for (cur_component in components) {
    component_subset <- subset(data, component == cur_component)
    min.x <- min(component_subset$offset, na.rm=TRUE)
    max.x <- max(component_subset$offset, na.rm=TRUE)
    at.x.pretty <- at.pretty(from=min.x, to=max.x, length=5, largest=TRUE)
    at.x <- c(at.x, at.x.pretty)
  }

  at.x.nonnull.full <- rep(at.x, 
    as.integer((layout[1] - remove.extra.scales) / num_components))

  # Remove internal axes and space where axes were
  at.x.full <- c(rep(list(NULL), num_panels - layout[1] + remove.extra.scales),
                 at.x.nonnull.full)

  max.y <- max(data$frequency, na.rm=TRUE)
  limits.y <- c(-0.001, max.y*1.1)
  at.y.full <- unlist(at.pretty(from=0, to=max.y, length=5))
  # Take only 0 and the last
  at.y <- c(0, tail(at.y.full, n=1))
  scales <- list(x = list(relation = "free",
                          tck = c(1, 0),
                          #at = at.x.full,
                          at = NULL,
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

transpose.levels <-
  function(data.factor, dim.length)
{
  lev <- levels(data.factor)
  nlev <- nlevels(data.factor)
  res <- vector(mode = "character", length = nlev)
  num.added <- 0
  for (i in seq(from = 1, length = dim.length)) {
    dest.indices <- seq(from = i, to = nlev, by = dim.length)
    source.indices <- seq(from = num.added + 1, along = dest.indices)
    res[dest.indices] <- lev[source.indices]
    num.added <- num.added + length(dest.indices)
  }
  
  factor(data.factor, levels = res)
}

# Plots frequency vs position for each label
#   data: a data frame with fields: frequency, offset, label
#   spacers should be a vector of indices, where a spacer will be placed after
#     that component (e.g. c(1, 3) will place spacers after the first and third
#     components
xyplot.aggregation <- 
  function(data,
    x = frequency ~ offset | component * label,
    spacers = NULL,
    genemode = FALSE,
    normalize_labels = FALSE,
    text.cex = 1,
    spacing.x = 0.4,
    spacing.y = 0.4,
    par.settings = list(add.text = list(cex = text.cex),
                        layout.heights = list(axis.panel = axis.panel,
                                              strip = strips.heights)),
    auto.key = if (nlevels(data$factor) < 2) FALSE
               else list(points = FALSE, lines = TRUE),
    strip = strip.custom(horizontal = FALSE),
    strip.height = 10,
    xlab = NULL,
    ylab = "Frequency",
    ...) 
{

  # Calculate frequencies from counts
  data$frequency <- 
    if (normalize_labels) {
      data$count / data$label_count
    } else {
      data$count / data$component_count
    }
  data$frequency[!is.finite(data$frequency)] <- 0


  # Determine panel layout
  num_levels <- nlevels(data$label)
  num_components <- nlevels(data$component)
  num_rows <- num_components
  num_cols <- num_levels
  num_panels <- num_rows * num_cols

  # Rework layout to optimize organization
  while (num_cols > num_rows) {
    num_cols <- num_cols / 2
    num_rows <- num_rows * 2
  }
  num_cols <- ceiling(num_cols)
  layout <- c(num_rows, num_cols)

  # Reorder labels so they are in order downward in panels
  data$label <- transpose.levels(data$label, num_rows / num_components)

  # Separate distinct groups
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


  scales <- panel.scales(data, layout, num_panels)
  axis.panel <- rep(c(0, 1), c(num_cols - 1, 1))

  # Make top strips longer
  strips.heights <- rep(c(strip.height, 0), c(1, num_cols - 1))

#  if (genemode) {
#    par.settings$strip.border <- list(col="transparent")
#    trellis <- resizePanels(trellis, w = panel.widths(
#  }

  trellis.raw <- xyplot(x, 
                        data = data, 
                        type = "l",
                        groups = factor, 
                        auto.key = auto.key,
                        as.table = TRUE,
                        strip = strip,
                        xlab = xlab,
                        ylab = ylab,
                        ...)

  trellis <- useOuterStrips(trellis.raw, strip = strip)

  update(trellis, layout = layout, between = between, scales = scales,
         par.settings = par.settings)
}

plot.aggregation <- function(filename, mnemonics=NULL, ..., verbose=FALSE) {
  data <- read.aggregation(filename, mnemonics=mnemonics)

  if (verbose) {
    xyplot.aggregation(data=data, ...)
  } else {
    suppressWarnings(xyplot.aggregation(data=data, ...))
  }
}
