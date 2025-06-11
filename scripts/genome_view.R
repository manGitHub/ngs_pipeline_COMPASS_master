#!/usr/bin/env Rscript

# script adapted from sequenza plot
# https://rdrr.io/cran/sequenza/src/R/graphics.R

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
in_file <- args[1]
out_file <- args[2]
in_file = read.delim(in_file,sep='\t')
# chrM is not needed
in_file = in_file[in_file$chromosome!='chrM',]

genome.view <- function(seg.cn, info.type = "AB", ...) {
    # adapt cnvkit input 
    seg.cn = seg.cn[,c('chromosome','start','end','cn1','cn2','cn')]
    colnames(seg.cn) = c("chromosome", "start.pos", "end.pos", "A", "B", "CNt")
    #
    chr.order <- unique(seg.cn$chromosome)
    seg.list <- split(x = seg.cn[,
        c("chromosome", "start.pos", "end.pos", "A", "B", "CNt")],
        f = seg.cn$chromosome)
    
    seg.list <- seg.list[order(order(chr.order))]
    seg.max <- lapply(X = seg.list, FUN = function(x) x[nrow(x), "end.pos" ])
    seg.pos <- lapply(seg.list, "[", TRUE, c("start.pos", "end.pos"))
    seg.max <- cumsum(as.numeric(do.call(rbind, seg.max)))
    chr.offset <- 0
    for (i in 1:length(seg.pos)){
        seg.pos[[i]] <- seg.pos[[i]] + chr.offset
        colnames(seg.pos[[i]]) <- c("abs.start", "abs.end")
        chr.offset <- seg.max[i]
    }
    seg.max <- sapply(X = seg.pos, FUN = function(x) x[nrow(x), "abs.end" ])
    abs.list <- mapply(cbind, seg.list, seg.pos, SIMPLIFY = FALSE)
    abs.segments <- do.call(rbind, abs.list)
    if (info.type == "AB") {
        na_As <- is.na(abs.segments$A)
        max_A <- max(abs.segments$A, na.rm = TRUE)
        abs.segments$A[na_As] <- abs.segments$CNt[na_As]
        plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
            y = c(-0.1, (max_A + 0.1)), type = "n",
            ylab = "Copy number", xlab = "Position (Mb)",
            xaxt = "n",  yaxt = "n", xaxs = "i", ...)
        axis(labels = 0:max_A, at = 0:max_A, side = 2, line = 0, las = 1)
	# adjusted spacing between lines
        segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
            y0 = (abs.segments$B), y1 = (abs.segments$B),
            col = "blue", lwd = 5, lend = 1)
        segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
            y0 = (abs.segments$A + 0.05), y1 = (abs.segments$A + 0.05),
            col = "red", lwd = 5, lend = 1)
    } else {
        min_CNt <- min(abs.segments$CNt, na.rm = TRUE)
        max_CNt <- max(abs.segments$CNt, na.rm = TRUE)
        plot(x = c(min(abs.segments$abs.start), max(abs.segments$abs.end)),
            y = c(min_CNt, max_CNt), type = "n",
            ylab = "Copy number", xlab = "Position (Mb)",
            xaxt = "n", yaxt = "n", xaxs = "i", ...)
        axis(labels = min_CNt:max_CNt,
            at = min_CNt:max_CNt,
            side = 2, line = 0, las = 1)
        segments(x0 = abs.segments$abs.start, x1 = abs.segments$abs.end,
            y0 = abs.segments$CNt, y1 = abs.segments$CNt, col = "red",
            lwd = 5, lend = 1)
    }
    abline(v = c(0, seg.max), lty = 3)
    for (i in 1:length(abs.list)){
        max.pos <- nrow(abs.list[[i]])
        mtext(chr.order[i], side = 3, line = 0,
            at = sum(abs.list[[i]]$abs.start[1],
                abs.list[[i]]$abs.end[max.pos]) / 2)
    }
    axis(labels = as.character(round(seq(abs.list[[1]]$start.pos[1] / 1e6,
        abs.list[[1]]$end.pos[nrow(abs.list[[1]])] / 1e6, by = 50), 0)),
        at = seq(abs.list[[1]]$abs.start[1],
            abs.list[[1]]$abs.end[nrow(abs.list[[1]])], by = 5e7),
        outer = FALSE, cex = par("cex.axis") * par("cex"), side = 1,
        line = 1)
}
# output results
pdf(file = out_file,height = 5, width = 15)
genome.view(in_file,info.type = "CNt")
genome.view(in_file)
dev.off()
