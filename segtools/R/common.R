library(grDevices)
library(lattice)
library(RColorBrewer)
library(latticeExtra)
library(plyr)

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

read.mnemonics <- function(filename, stringsAsFactors = NULL,
                           colClasses = NULL, na.strings = NULL, ...) {
  filename <- as.character(filename)
  if (length(filename) > 0 && nchar(filename) > 0) {
    if (!file.exists(filename)) {
      stop(paste("Error: could not find mnemonic file:", filename))
    }
    read.delim(filename, stringsAsFactors = FALSE, colClasses = "character",
               na.strings = NULL, ...)[, 1:2]
  } else {
    NULL
  }
}

read.metadata <- function(filename, ...) {
  ## Read metadata of the form: # <variable>=<value> <variable>=<value>
  if (!file.exists(filename)) {
    stop(paste("Error reading metadata: file not found:", filename))
  }
  lines <- read.comment.lines(filename, ...)
  metadata <- list()
  for (line in lines) {
    tokens.split <- strsplit(line, "=")
    for (tokens in tokens.split) {
      if (length(tokens) != 2) stop("Error occurred processing file metadata")
      variable <- tokens[1]
      value <- tokens[2]
      metadata[[variable]] <- value
    }
  }
  metadata
}

read.comment.lines <- function(filename, comment.char = "#", ...) {
  ## Comment lines must begin with:
  # ...
  conn <- file(filename, "r")
  lines <- list()
  line <- readLines(conn, n = 1)
  while (length(line) > 0) {
    if (substr(line, 1, 1) == comment.char) {
      lineConn <- textConnection(substr(line, 2, nchar(line)))
      tokens <- scan(lineConn, what = character(0), quiet = TRUE)
      close(lineConn)
      lines[[length(lines) + 1]] <- tokens
      line <- readLines(conn, n = 1)
    } else {
      break
    }
  }
  close(conn)
  lines
}

read.matrix <- function(filename, nrow = NULL, ncol = NULL,
                        comment.char = "#") {
  comment.lines <- read.comment.lines(filename, comment.char = comment.char)
  skip <- length(comment.lines)
  
  if (is.null(ncol)) {
    ## Read first line to determine number of columns
    first.line <- scan(filename, quiet = TRUE, skip = skip, nlines = 1)
    ncol <- length(first.line)
  }
  if (!is.null(nrow)) {
    n <- as.integer(nrow * ncol)
  } else {
    n <- -1
  }
  ## Read data
  data <- scan(filename, n = n, quiet = TRUE, skip = skip)
  if (is.null(nrow)) {
    n <- length(data)
    nrow <- as.integer(n / ncol)
  }
  matrix(data, nrow, ncol, byrow = TRUE)
}

# Given a list of labels (e.g. levels(data$label)), returns a dataframe
# with group and index fields, corresponding to the character and numeric 
# components of each label.
labels.classify <- function(labels) {
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
      match <- regexpr("^[a-zA-Z]+", label)
      match.len <- attributes(match)$match.length
      if (match.len > 0) {
        labels.split[i, 1] <- substring(label, 1, match.len)
        if (nchar(label) > match.len) {
          labels.split[i, 2] <- substring(label, match.len + 1, nchar(label))
        } else {
          labels.split[i, 2] <- "0"
        }
      } else {
        labels.split[i, 1] <- label
        labels.split[i, 2] <- "0"
      }
    }
  }

  data.frame(group=factor(labels.split[,1]), index=factor(labels.split[,2]),
             stringsAsFactors = FALSE)
}

relevel.mnemonics <- function(labels, mnemonics = NULL) {
  if (!is.factor(labels)) stop("Expected labels as factor")
  label.map <- map.mnemonics(levels(labels), mnemonics = mnemonics)
  if (nunique(label.map$labels) != length(label.map$labels)) {
    stop(paste("Error: mnemonic replacement labels are not unique:",
               paste(label.map$labels, collapse = ", ")))
  }
  levels(labels) <- label.map$labels
  factor(labels, levels(labels)[label.map$order])
}

map.mnemonics <- function(labels, mnemonics = NULL) {
  ## Given a vector of string labels and a mnemonics data.frame,
  ## Returns a list with two elements:
  ##   labels: a vector of new string labels
  ##   order: the order these new string labels should then be arranged in
  labels.raw <- as.character(labels)
  labels.new <- labels.raw
  labels.order <- 1:length(labels)
  if (length(mnemonics) > 0) {
    if (ncol(mnemonics) < 2) stop("Mnemonics should have at least 2 columns")
    
    ## Order by mnemonics
    mnemonics <- data.frame(old=as.character(mnemonics[,1]),
                            new=as.character(mnemonics[,2]),
                            stringsAsFactors=FALSE)
    
    ## Map the labels that we can
    mapping.rows <- match(labels.raw, mnemonics$old)
    mapping.valid <- is.finite(mapping.rows)
    labels.mapped <- mnemonics$new[mapping.rows[mapping.valid]]
    labels.new[mapping.valid] <- labels.mapped
    ## Leave the labels we can't
    
    ## Order the mapped labels
    labels.order[mapping.valid] <- order(mapping.rows[mapping.valid])
  } else {
    mapping.valid <- vector(length = length(labels))
  }
  
  ## Order the unmapped labels, first numerically
  nmapped <- sum(mapping.valid)
  labels.unmapped <- labels.new[!mapping.valid]
  labels.num <- suppressWarnings(as.numeric(labels.unmapped))
  labels.num.valid <- is.finite(labels.num)
  labels.order[!mapping.valid][labels.num.valid] <-
    order(labels.num[labels.num.valid]) + nmapped
  ## ... then lexicographically
  nmapped <- nmapped + sum(labels.num.valid)
  labels.order[!mapping.valid][!labels.num.valid] <-
    order(labels.unmapped[!labels.num.valid]) + nmapped
  
  list(labels = labels.new, order = labels.order)
}

## Create a dendrogram legend
ddgram.legend <- function(dd.row = NULL, row.ord = NULL,
                          dd.col = NULL, col.ord = NULL) {
  ## if dd.row or dd.col is NULL, that dendrogram is not drawn
  legend <- list()

  if (!is.null(dd.row)) {
    legend$right <- list(fun = dendrogramGrob,
                         args = list(x = dd.row, 
                           ord = row.ord,
                           side = "right",
                           size = 10))
  }
  if (!is.null(dd.col)) {
    legend$top <- list(fun = dendrogramGrob,
                       args = list(x = dd.col, 
                         ord = col.ord,
                         side = "top"))
  }
  
  legend
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
  function(basename, ext, dirname, device, as.slide = FALSE, ...,
           clobber = FALSE) 
{
  filename.ext <- extpaste(basename, ext)
  filename.fq <- file.path(dirname, filename.ext)

  if (!clobber && file.exists(filename.fq)) {
    stop(paste("Error:", filename.fq, "already exists!",
               "Specify --clobber to overwrite."))
  }
  
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
                             make.png = TRUE,
                             make.slide = TRUE,
                             make.pdf = TRUE,
                             make.thumb = TRUE,
                             pdf.as.slide = TRUE,
                             clobber = FALSE,
                             ...) {
  # No need for PNG plot since it is done python-side
  if (make.png) {
    filename.raster <-
      save.image(basename, "png", dirname, device.png, 
                 width = width, height = height, units = "px",
                 ..., clobber = clobber)
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
                 units = "in", res = res.slide, ..., clobber = clobber)
  } else if (make.slide) {
    filename.slide <-
      save.image(basename, "slide.png", dirname, device.png, 
                 width = width.slide, height = height.slide, 
                 units = "px", as.slide = pdf.as.slide,
                 ..., clobber = clobber)
  }

  if (make.pdf) {
    filename.pdf <-
      save.image(basename, "pdf", dirname, device.pdf,
                 width = width.pdf, height = height.pdf,
                 useDingbats = FALSE, as.slide = TRUE, ...,
                 clobber = clobber)
  }

  if (make.thumb) {
    # Suppress warnings about minimum font size
    filename.thumb <-
      suppressWarnings(
        save.image(basename, "thumb.png", dirname, device.png, 
                   width = 10, height = 10, 
                   units = "in", res = 10, ..., clobber = clobber))
  }
}

save.images <- function(dirpath, basename, image, ..., clobber = FALSE) {
  cat("Saving images...", file=stderr())
  plot(image)
  dev.print.images(basename, dirpath, ..., clobber = clobber)
  cat("... done saving images.", file=stderr())
}
