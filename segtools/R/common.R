library(grDevices)
library(lattice)
library(RColorBrewer)

# Set default theme
lattice.options(default.theme = "theme.dark2",
                legend.bbox = "full")


# Generate log labels
at.log <- function(low.exp = 0, high.exp = 9) {
  10^seq(from = low.exp, to = high.exp)
}


## Generate axis label list for a reasonable axis
at.pretty <- function(from = 0, to = 10, length = 4,
                      largest = FALSE, trim = TRUE) {
  if (from == to) {
    list(from)
  } else {
    at.seq <- pretty(c(from, to), n = length, min.n = length)
    if (trim) {
      # Trim sequence from 0-end
      trimmer <- if (abs(from) < abs(to)) head else tail
      at.seq <- trimmer(at.seq, n = length)
    }
    trimmer <- if (abs(from) < abs(to)) tail else head
    if (largest) {
      # Get only largest value (by abs)
      trimmer(at.seq, n = 1)
    } else {
      at.seq
    }
  }
}

label.log <- function(x) {
  if (x <= 100) {
    x
  } else {
    label <- substitute(10^X, list(X=log10(x)))
    do.call(expression, list(label))
  }
}

labels.log <- function(...) {
  sapply(at.log(...), label.log)
}

read.mnemonics <- function(filename, stringsAsFactors = NULL, colClasses = NULL,
                           ...) {
  read.delim(filename, stringsAsFactors = FALSE, colClasses = "character",
             ...)[, 1:2]
}

# Given a list of labels (e.g. levels(data$label)), returns a dataframe
# with group and index fields, corresponding to the character and numeric 
# components of each label.
labels.classify <-
  function(labels)
{
  labels.str <- as.character(labels)
  labels.split <- matrix(nrow = length(labels.str), ncol = 2)

  # First, split on "."
  labels.parts <- strsplit(labels.str, ".", fixed = TRUE)

  for (i in 1:length(labels.parts)) {
    label <- labels.parts[[i]]
    if (length(label) > 1) {
      ## Splitting worked!
      labels.split[i, 1:2] <- label[1:2]
    } else {
      ## Splitting didn't work. Try treating like normal mnemonics (e.g. TF0)
      match <- regexpr("[0-9]+$", label)
      match.len <- attributes(match)$match.length
      if (match.len > 0) {
        labels.split[i, 1] <- substring(label, 1, match - 1)
        labels.split[i, 2] <- substring(label, match, nchar(label))
      } else {
        labels.split[i, 1] <- label
        labels.split[i, 2] <- "0"
      }
    }
  }

  data.frame(group=factor(labels.split[,1]), index=factor(labels.split[,2]),
             stringsAsFactors = FALSE)
}

# Substitutes names for mnemonics and reorders to match mnemonic ordering
# Factor should have integer labels
# Usage: x$factor <- relabel.factor(x$factor, mnemonics)
relabel.factor <- function(field, mnemonics) {
  # Order by mnemonics
  mnemonics <- data.frame(old=as.integer(mnemonics[,1]),
                          new=as.character(mnemonics[,2]),
                          stringsAsFactors=FALSE)
  # Substitute label names
  levels.raw <- levels(field)
  mapping.rows <- match(levels.raw, mnemonics$old)
  mapping.valid <- is.finite(mapping.rows)
  levels.mapped <- mnemonics$new[mapping.rows[mapping.valid]]
  levels(field)[mapping.valid] <- levels.mapped

  # Fix level ordering
  levels.ordered <- c(levels(reorder(levels.mapped, 
                                     mapping.rows[mapping.valid])),
                      levels.raw[!mapping.valid])
  factor(field, levels=levels.ordered)
}

# Reorders factor first numerically, and then lexicographically
# Usage: x$factor <- smart.factor.reorder(x$factor)
smart.factor.reorder <- function(field) {
  levels.raw <- levels(field)
  suppressWarnings(levels.int <- as.integer(levels.raw))
  levels.valid <- (is.finite(levels.int) &
                   as.character(levels.int) == as.character(levels.raw))
  levels.full <- c(sort(levels.int[levels.valid]), 
                   sort(levels.raw[!levels.valid]))

  # Reordered factor
  factor(field, levels=levels.full)
}

theme.dark2 <- function() {
  brew.dark2.8 <- brewer.pal(8, "Dark2")
  brew.paired.12 <- brewer.pal(12, "Paired")
  brew.paired.light.6 <- brew.paired.12[seq(1, 12, 2)]
  brew.paired.dark.6 <- brew.paired.12[seq(2, 12, 2)]

  res <- col.whitebg()

  res$reference.line$col <- gray(0.7) # 30% grey
  res$axis.line$col <- gray(0.7)
  res$dot.line$col <- gray(0.7)
  res$axis.line$lwd <- 2

  res$superpose.line$col <- brew.dark2.8
  res$superpose.symbol$col <- brew.dark2.8

  res$plot.symbol$col <- brew.dark2.8[1]
  res$plot.line$col <- brew.dark2.8[1]

  res$dot.symbol$col <- brew.dark2.8[1]

  res$box.rectangle$col <- brew.dark2.8[1]
  res$box.umbrella$col <- brew.dark2.8[1]

  res$strip.background$col <- brew.paired.light.6
  res$strip.border$col <- gray(0.7)
  res$strip.border$lty <- 1
  res$strip.border$lwd <- 2

  res$strip.shingle$col <- brew.paired.dark.6

  res
}

theme.slide <- function() {
 res <- theme.dark2()

 additions <-
   list(axis.text = list(cex = 1.3),
        par.main.text = list(cex = 2),
        par.sub.text = list(cex = 2),
        layout.heights = list(axis.top = 1.5, axis.bottom = 1.25),
        layout.widths = list(axis.left = 1.25),
        add.text = list(cex = 1.5),
        par.xlab.text = list(cex = 1.5),
        par.ylab.text = list(cex = 1.5),
        par.zlab.text = list(cex = 1.5))

 modifyList(res, additions)
}

lattice.optionlist.slide <-
 list(layout.heights = list(top.padding = list(x = 0.05, units = "snpc"),
        bottom.padding = list(x = 0.05, units = "snpc")), # strip = list(x = 2, units = "lines")),
      layout.widths = list(left.padding = list(x = 0.05,
                             units = "snpc"),
        right.padding = list(x = 0.05, units = "snpc")),
      default.theme = theme.slide)

extpaste <- function(...) {
  paste(..., sep=".")
}

as.slide <- function() {
 update(trellis.last.object(),
        lattice.options = lattice.optionlist.slide,
        par.settings = theme.slide())
}

calc.slide.res <- function(width, height, 
                           screen.width = 1600, 
                           screen.height = 1200,
                           screen.res = 72) {
  ratios <- c(screen.width / width, screen.height / height) / 1.5
  if (max(ratios) < screen.res) {
    max(ratios)
  } else {
    100
  }
}

calc.slide.scale <- function(width, height, pivot.dim = 1000) {
  ratios <- c(width / pivot.dim, height / pivot.dim)
  if (any(ratios > 1)) {
    scaler <- 1 + ratios[which(ratios > 1)[1]] / 5
  } else {
    scaler <- 1
  }
  125 * scaler
}

save.image <- 
  function(basename, ext, dirname, device, as.slide = FALSE, ...) 
{
  filename.ext <- extpaste(basename, ext)
  filename.fq <- file.path(dirname, filename.ext)

  if (as.slide) {
    device(filename.fq, ...)
    print(as.slide())
    dev.off()
  } else {
    dev.print(device=device, file=filename.fq, ...)
  }

  filename.ext
}

dev.print.images <- function(basename, dirname,
                             width = 800, height = 800,
                             width.slide = 1280, height.slide = 1024,
                             width.pdf = 11, height.pdf = 8.5,
                             device.png = png,
                             device.pdf = pdf, 
                             downsample = FALSE,
                             make.png = FALSE,  # PNG made in python by default
                             make.slide = TRUE,
                             make.pdf = TRUE,
                             make.thumb = TRUE,
                             ...)
{
  # No need for PNG plot since it is done python-side
  if (make.png) {
    filename.raster <-
      save.image(basename, "png", dirname, device.png, 
                 width = width, height = height, units = "px", ...)
  }

  if (downsample && make.slide) {
    # Useful for tweaking text size on large plots without replotting
    slide.scale <- calc.slide.scale(width, height)
    width.slide <- width/slide.scale + height/width
    height.slide <- height/slide.scale + width/height
    res.slide <- calc.slide.res(width = width.slide, height = height.slide)
    filename.slide <-
      save.image(basename, "slide.png", dirname, device.png, 
                 width = width.slide, height = height.slide, 
                 units = "in", res = res.slide, ...)
  } else if (make.slide) {
    filename.slide <-
      save.image(basename, "slide.png", dirname, device.png, 
                 width = width.slide, height = height.slide, 
                 units = "px", as.slide = TRUE, ...)
  }

  if (make.pdf) {
    filename.pdf <-
      save.image(basename, "pdf", dirname, device.pdf,
                 width = width.pdf, height = height.pdf,
                 useDingbats = FALSE, as.slide = TRUE, ...)
  }

  if (make.thumb) {
    # Suppress warnings about minimum font size
    filename.thumb <-
      suppressWarnings(
        save.image(basename, "thumb.png", dirname, device.png, 
                   width = 10, height = 10, 
                   units = "in", res = 10, ...))
  }
} 
