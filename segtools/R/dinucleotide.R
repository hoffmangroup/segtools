library(plyr)
library(reshape)
library(RColorBrewer)
library(lattice)
library(latticeExtra)

read.nuc <- function(filename = "nucleotide.tab", mnemonics = NULL, ...,
                     check.names = FALSE)
{
  res <- read.delim(filename, ..., check.names = check.names)
  res$label <- relevel.mnemonics(factor(res$label), mnemonics)

  res
}

normalize.dinuc <- function(data) {
  data.df <- data.frame(data, check.names = FALSE)
  dinuc <- subset(data.df, select = -c(`?`, `??`, `A|T`, `C|G`))

  # Compensate for symmetric dinucleotide
  dinuc$`AT|AT` <- 2*dinuc$AT
  dinuc$`CG|CG` <- 2*dinuc$CG
  dinuc$`GC|GC` <- 2*dinuc$GC
  dinuc$`TA|TA` <- 2*dinuc$TA
  dinuc <- subset(dinuc, select = -c(AT, CG, GC, TA))

  # Normalize each row to one
  label.col <- match("label", colnames(dinuc))
  dinuc.nolabel <- dinuc[, -label.col]
  dinuc[, -label.col] <- dinuc.nolabel / rowSums(dinuc.nolabel)

  # order columns
  dinuc[,c(label.col, order(colnames(dinuc)[-label.col]) + 1)]
}

levelplot.dinuc <-
  function(x = value ~ variable + label,
           data, 
           ..., 
           deselect = c("C|G", "A|T", "?", "??"),
           aspect = "fill",
           palette = colorRampPalette(rev(brewer.pal(11, "RdYlBu")),
             interpolate = "spline", space = "Lab")(100))
{
  data.melted <- melt(data, id.vars="label")
  data.subset <- subset(data.melted, !variable %in% deselect)

  # Reverse label ordering for plotting
  data.subset$label <- factor(data.subset$label, 
                              levels=rev(levels(data.subset$label)))

  levelplot(x = x,
            data = data.subset,
            aspect = aspect,
            cuts = 99,
            col.regions = palette,
            xlab = "Dinucleotide",
            ylab = "Segment label")
}

plot.dinuc <- function(tabfile, mnemonics = NULL) {
  nucs <- read.nuc(tabfile, mnemonics = mnemonics)
  dinuc <- normalize.dinuc(nucs)

  levelplot.dinuc(data = dinuc)
}


save.dinuc <- function(dirpath, namebase, tabfilename,
                       mnemonic_file = NULL,
                       clobber = FALSE,
                       image.size = 600,  # px
                       ...) {
  mnemonics <- read.mnemonics(mnemonic_file)
  nucs <- read.nuc(tabfilename, mnemonics = mnemonics)
  dinuc <- normalize.dinuc(nucs)
  
  save.images(dirpath, namebase,
              levelplot.dinuc(data = dinuc, ...),
              height = image.size,
              width = image.size,
              clobber = clobber)
}
  
