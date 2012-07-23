read.distances <-
  function(filename, mnemonics = NULL, ..., check.names = FALSE)
{
  res <- read.delim(filename, ..., check.names = check.names)
  res$label <- relevel.mnemonics(factor(res$label), mnemonics)
  res
}

panel.distances <-
  function(...)
{
  panel.barchart(..., origin = 0)
}

figure.distances <- 
  function(data, 
           x = count ~ distance | label,
           xlab = "Distance from feature (0 if overlapping)",
           ylab = "Number of segments",
           panel = panel.distances,
           col = "blue",
           auto.key = as.logical(nlevels(data$group) > 1),
           label.count = 11,
           ...)
{
  layout <- c(1, nlevels(data$label))
  distance.vals <- sort(unique(data$distance))
  
  finite.distances <- data$distance[is.finite(distance.vals)]
  distance.max <- max(finite.distances)
  distance.min <- min(finite.distances)
  
  distance.labels <- format(distance.vals, scientific = FALSE)

  distance.labels[distance.vals == Inf] <-
    paste(">", format(distance.max, scientific = FALSE))
  distance.labels[distance.vals == -Inf] <-
    paste("<", format(distance.min, scientific = FALSE))
  keep.label.indices <-
    if (distance.min == 0) {
      seq(1, length(distance.labels), length.out=label.count)
    } else {
      zero.index <- which(distance.vals == 0)
      c(seq(1, zero.index, length.out=as.integer(label.count / 2 + 1)),
        tail(seq(zero.index, length(distance.labels),
                 length.out=as.integer(label.count / 2 + 1)), -1))
    }


  keep.labels <- vector(mode = "logical", length = length(distance.labels))
  keep.labels[keep.label.indices] <- TRUE
  distance.labels[!keep.labels] <- ""
  
  # Replace infinities with bounded vals so everything sorts correctly
  data$distance[data$distance == Inf] <- distance.max + 1
  data$distance[data$distance == -Inf] <- distance.min - 1
  barchart(x = x, 
           data = data,
           group = group,
           xlab = xlab, 
           ylab = ylab,
           panel = panel,
           col = col,
           auto.key = auto.key,
           horizontal = FALSE,
           stack = TRUE,
           layout = layout,
           scales = list(x = list(rot = 90, labels=distance.labels)),
           ...)
}

plot.distances <- function(filename, ...) {
  data <- read.distances(filename)
  figure.distances(data = data, ...)
}

save.distances <- function(dirpath, namebase, tabfilename,
                           clobber = FALSE, mnemonic_file = NULL,
                           height = 200 * rows, width = 600,
                           ...) {
  mnemonics <- read.mnemonics(mnemonic_file)
  data <- read.distances(tabfilename, mnemonics = mnemonics)
  rows <- nlevels(data$label)
  height <- 200 * rows
  save.plots(dirpath, namebase,
              figure.distances(data = data, ...),
              height = height,
              width = width,
              clobber = clobber)
}
  
