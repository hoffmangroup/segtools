library(lattice)
library(RColorBrewer)
library(latticeExtra)
library(plyr)
library(reshape2)

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

  ## Read metadata and mnemonic-ize it as well
  metadata <- read.metadata(filename, comment.char = comment.char)
  names(metadata) <- map.mnemonics(names(metadata), mnemonics)$labels

  list(data = data, metadata = metadata)
}


## Generates pretty scales for the data.
## layout is a 2-element vector: c(num_rows, num_cols) for the xyplot
## num_panels is the number of panels/packets in the trellis
panel.scales.aggregation <- function(data, layout, num_panels, x.axis = FALSE,
                         significance = FALSE) {
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
    if (!is.finite(min.x))
      print(component_subset)
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

  ## Extend range for rug
  ngroups <- nlevels(data$group)
  rug.height <-
    if (significance && ngroups > 1)
      diff(range.y) / (1 - get.rug.start(ngroups + 2)) - diff(range.y)
    else
      0

  min.y <- min(range.y[1], 0) - rug.height
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
                          at = at.y)
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

calc.pvals <- function(count, total, random.prob) {
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

  signif
}

## Given an aggregation data frame, process the counts by optionally
## normalizing the counts over all labels and by the sizes of the labels
## as well as calculating a significance level for each row.
process.counts <- function(data, label.sizes, pseudocount = 1,
                           normalize = TRUE) {
  if (!is.vector(label.sizes)) stop("Expected vector of label sizes")
  total.size <- sum(label.sizes)
  if (length(total.size) != 1) stop("Error summing label size vector")
  if (!is.finite(total.size)) stop("Unexpected sum of sizes (not finite)")

  ## Ideally, a multinomial test would be applied once for each bin
  ## (testing whether the distribution of overlaps by label is as expected),
  ## but the label-wise binomial seems a decent approximation.

  ## Get sum over labels for each component and offset
  labels.sum <- with(data, aggregate(count, list(offset = offset,
                                                 component = component,
                                                 group = group), sum))
  names(labels.sum) <- gsub("^x$", "sum", names(labels.sum))

  for (label in levels(data$label)) {
    random.prob <- label.sizes[label] / total.size
##    if (label == "NONE") cat(paste("\n\nNONE:\nrandom.prob:", random.prob))
    cur.rows <- data$label == label
    cur.data <- data[cur.rows,]
    ## Add order row to keep track of order extracted
    ## (so we can reorder for re-insertion)
    cur.data$order <- 1:nrow(cur.data)
    ## Add sum over labels as data frame column
    cur.data <- merge(cur.data, labels.sum)
    ## Reorder merged data back to the order extracted
    cur.data <- cur.data[order(cur.data$order),]
    calc.row.pvals <- function(row) {
      count <- row[1]
      total <- row[2]
      calc.pvals(count, total, random.prob)
    }

##    if (label == "NONE") print(cur.data[1:1,])
    data$pval[cur.rows] <- apply(subset(cur.data, select = c(count, sum)), 1,
                                 calc.row.pvals)

    if (normalize) {
      cur.enrichments <- log2((cur.data$count / cur.data$sum + pseudocount) /
                              (random.prob + pseudocount))
      ## Any NaN regions (from 0/0) are masked
      cur.enrichments[is.nan(cur.enrichments)] <- 0

      data$count[cur.rows] = cur.enrichments
    }
  }
  data
}

panel.significance <- function(x, y, height, signif, col = plot.line$col,
                               lty = plot.line$lty, lwd = plot.line$lwd, ...) {
  ## signif should be a boolean vector as long as x
  n <- length(x)
  if (n != length(signif)) stop()
  require("grid", quietly = TRUE)
  x.units <- "native"
  y.units <- "npc"
  plot.line <- trellis.par.get("plot.line")

  ## Identify contiguous blocks of significance
  starts <- x[c(signif[1], !signif[-n] & signif[-1])]
  ends <- x[c(signif[-n] & !signif[-1], signif[n])]

  gp <- gpar(col = col, lty = 1, lwd = 0, fill = col)
  x0 <- unit(starts - 0.5, x.units)
  #x1 <- unit(ends + 0.5, x.units)
  width <- unit(ends + 0.5 - starts, x.units)
  y0 <- unit(y, y.units)
  height <- unit(height, y.units)
#  grid.segments(x0, y0, x1, y0, gp = gp)
  grid.rect(x0, y0, width, height, just = c("left", "bottom"), gp = gp)
}

get.rug.start <- function(group.number, rug.height = 0.03, rug.spacing = 0.015,
                          rug.offset = 0.025, ...) {
  rug.offset + (rug.spacing + rug.height) * (group.number - 1)
}

panel.aggregation <- function(x, y, ngroups, signif = NULL, groups = NULL,
                              subscripts = NULL, font = NULL, col = NULL,
                              col.line = NULL, pch = NULL,
                              col.high = rgb(.773,.098,.133),
                              col.low = rgb(.224,.365,.655),
                              group.number = NULL, rug.height = 0.03, ...) {
  ## hide 'font' from panel.segments, since it doesn't like having
  ## font and fontface
  panel.refline(h = 0)

  x <- as.numeric(x)
  y <- as.numeric(y)

  if (!is.null(signif)) {
    signif <- as.logical(signif)[subscripts]
    signif[!is.finite(signif)] <- FALSE

    if (any(signif)) {
      ## Only shade region for first.
      if (ngroups == 1) {
        fill.col <- "black"
        y.sig <- y
        y.sig[!signif] <- 0
        panel.polygon(cbind(c(min(x), x, max(x)), c(0, y.sig, 0)),
                      col = fill.col)
      } else {
        ##x.sig <- x[signif]
        ##y.sig <- y[signif]
        ##panel.points(x.sig, y.sig, col = fill.col, pch = "*")
        rug.start <- get.rug.start(group.number, rug.height = rug.height, ...)
        panel.significance(x, rug.start, rug.height, signif = signif,
                           col = col.line, ...)
      }
    }
  }

  y.fill <- y
  y.fill[y < 0] <- 0
  panel.polygon(cbind(c(min(x), x, max(x)), c(0, y.fill, 0)),
                col=col.high, border=NA)
  y.fill <- y
  y.fill[y > 0] <- 0
  panel.polygon(cbind(c(min(x), x, max(x)), c(0, y.fill, 0)),
                col=col.low, border=NA)

  ##  panel.xyplot(x, y, groups = groups, subscripts = subscripts, ...)
  #panel.xyplot(x, y, font = font, col = NULL, col.line = NULL,
  #             pch = pch, ...)
}

get.label.sizes <- function(data, metadata) {
  if (length(metadata) == 0) stop("Missing aggregation metadata")
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

calculate.signif <- function(data, fdr.level = 0.05) {
  require("qvalue", quietly = TRUE)
  qobj <- qvalue(p = data$pval, pi0.method = "bootstrap", fdr.level = fdr.level)
  qobj$significant
}

## Plots overlap vs position for each label
##   data: a data frame with fields: count, offset, label
##   spacers should be a vector of indices, where a spacer will be placed after
##     that component (e.g. c(1, 3) will place spacers after the first and third
##     components
xyplot.aggregation <-
  function(agg.data = NULL,
           data = agg.data$data,
           x = overlap ~ offset | component * label,
           metadata = agg.data$metadata,
           spacers = metadata$spacers,
           normalize = TRUE,
           significance = FALSE,
           pseudocount = 1,
           x.axis = FALSE,  # Show x axis
           fdr.level = 0.05,
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
           strip.height = 11,
           par.strip.text = list(cex = 0.85),
           xlab = NULL,
           ylab = if (normalize) {
             if (pseudocount != 0)
               substitute(expression(paste("Enrichment: ",
                   log[2], group("[", (f[plain(obs)] + PSEUDOCOUNT)/(f[plain(rand)] + PSEUDOCOUNT),
                                 "]"))),
                          list(PSEUDOCOUNT = pseudocount)
                          )
             else expression(paste("Enrichment: ", log[2],
                 group("[", f[plain(obs)]/f[plain(rand)], "]")))
           }
           else "Count",
           sub = if (significance) {
             substitute(expression(paste(FORMAT, " regions are significant with ",
                 italic(q) <= FDR)),
                        list(FORMAT = if (ngroups > 1) "Rug" else "Black",
                             FDR = fdr.level))
           } else {
             NULL
          },
    ...)
{
  label.sizes <- get.label.sizes(data, metadata)
  ## Normalize and/or calculate significance if metadata exists
  data <- process.counts(data, label.sizes, pseudocount = pseudocount, normalize = normalize)
  if (significance) {
    data$signif <- calculate.signif(data, fdr.level)
  }

  names(data) <- gsub("^count$", "overlap", names(data))
  data$overlap[!is.finite(data$overlap)] <- 0

  ## Determine panel layout
  num_levels <- nlevels(data$label)
  num_components <- nlevels(data$component)
  num_rows <- num_components
  num_cols <- num_levels
  num_panels <- num_rows * num_cols

  num_cols <- ceiling(num_cols)
  layout <- c(num_rows, num_cols)

  ## Reorder labels so they are in order downward in panels
  data$label <- transpose.levels(data$label, num_rows / num_components)

  ## Separate distinct groups
  spaces.x <- rep(0, num_components - 1)
  if (!is.null(spacers) && length(spacers) > 0) {
    spacers <- as.integer(spacers)
    if (any(spacers < 1 | spacers >= num_components)) {
      stop("Spacer vector should only contain values in the range [1, ",
           num_components - 1,
           "] since there are ", num_components, " components")
    }
    spaces.x[spacers] <- spacing.x
  }
  between <- list(x = c(spaces.x, spacing.x),
                  y = spacing.y)

  scales <- panel.scales.aggregation(data, layout, num_panels, x.axis = x.axis,
                                     significance = significance)
  axis.panel <- rep(c(0, 1), c(num_cols - 1, 1))

  ## Make top strips longer
  strips.heights <- rep(c(strip.height, 0), c(1, num_cols - 1))

  args <- list(x, data = data, type = "l", groups = quote(group),
               auto.key = auto.key, as.table = TRUE, strip = strip,
               xlab = xlab, ylab = ylab, sub = sub,
               signif = data$signif, ngroups = ngroups,
               panel = "panel.superpose", panel.groups = panel,
               par.strip.text = par.strip.text,
               ...)

  trellis.raw <- do.call(xyplot, args)

  trellis <- useOuterStrips(trellis.raw, strip = strip)

  update(trellis, layout = layout, between = between, scales = scales,
         par.settings = par.settings)
}

plot.aggregation <- function(filename, mnemonics = NULL, ...,
                             comment.char = "#") {
  data <- read.aggregation(filename, mnemonics = mnemonics,
                           comment.char = comment.char)

  xyplot.aggregation(data, ...)
}

save.aggregation <- function(dirpath, namebase, tabfilename,
                             mnemonic_file = NULL,
                             clobber = FALSE,
                             panel.size = 200,  # px
                             margin.size = 200,  #px
                             comment.char = "#",
                             ...) {
  mnemonics <- read.mnemonics(mnemonic_file)
  data <- read.aggregation(tabfilename, mnemonics = mnemonics,
                           comment.char = comment.char)

  image.size <- margin.size +
    panel.size * ceiling(sqrt(nlevels(data$data$label) / 2))
  save.plots(dirpath, namebase,
              xyplot.aggregation(data, ...),
              width = image.size,
              height = image.size,
              clobber = clobber)
}

save.gene.aggregations.data <-
  function(data, dirpath, namebase1 = NULL, namebase2 = NULL, clobber = FALSE,
           panel.size = 200,  # px
           margin.size = 200,  #px
           width = image.size, height = image.size,
           width.pdf = image.size / 72, height.pdf = image.size / 72,
           ...)
{
  image.size <- margin.size +
    panel.size * ceiling(sqrt(nlevels(data$data$label) / 2))
  ngroup1 <- as.integer(data$metadata[["spacers"]])
  if (ngroup1 <= 0 || ngroup1 >= nlevels(data$data$component))
    stop("Invalid metadata for gene aggregation plots")

  plot.one <- function(namebase, components.sub, ...) {
    if (is.null(namebase))
        return()

    data.sub <- subset(data$data, component %in% components.sub, drop = TRUE)
    data.sub$component <- factor(data.sub$component, levels=components.sub)
    save.plots(dirpath, namebase,
                xyplot.aggregation(data = data.sub, metadata = data$metadata,
                                   spacers = NULL, ...),
                width = width,
                height = height,
                width.pdf = width.pdf,
                height.pdf = height.pdf,
                clobber = clobber)
  }

  components <- levels(data$data$component)
  ## Transcriptional components
  plot.one(namebase1, head(components, ngroup1), ...)
  ## Translational components (plus flankings again)
  plot.one(namebase2, c(components[1],
                        tail(components, -ngroup1),
                        components[ngroup1]), ...)
}

save.gene.aggregations <-
  function(dirpath, namebase1 = NULL, namebase2 = NULL, tabfilename, mnemonic_file = NULL,
           comment.char = "#", ...)
{
  mnemonics <- read.mnemonics(mnemonic_file)
  data <- read.aggregation(tabfilename, mnemonics = mnemonics,
                           comment.char = comment.char)

  save.gene.aggregations.data(data, dirpath, namebase1, namebase2, ...)
}
