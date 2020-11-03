\name{plot.EarlyOrLate}
\alias{plot.EarlyOrLate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot Early or Late 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.EarlyOrLate(seg.mat.patient, TCGA.earlyLate, TCGA.purity, TCGA.barcode, max.cpn = 5, min.probes = 10, sub.clonal = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg.mat.patient}{
%%     ~~Describe \code{seg.mat.patient} here~~
}
  \item{TCGA.earlyLate}{
%%     ~~Describe \code{TCGA.earlyLate} here~~
}
  \item{TCGA.purity}{
%%     ~~Describe \code{TCGA.purity} here~~
}
  \item{TCGA.barcode}{
%%     ~~Describe \code{TCGA.barcode} here~~
}
  \item{max.cpn}{
%%     ~~Describe \code{max.cpn} here~~
}
  \item{min.probes}{
%%     ~~Describe \code{min.probes} here~~
}
  \item{sub.clonal}{
%%     ~~Describe \code{sub.clonal} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (seg.mat.patient, TCGA.earlyLate, TCGA.purity, TCGA.barcode, 
    max.cpn = 5, min.probes = 10, sub.clonal = 1) 
{
    seg.mat.patient$Chromosome <- as.numeric(as.character(seg.mat.patient$Chromosome))
    seg.mat.patient$StartPosition <- as.numeric(as.character(seg.mat.patient$StartPosition))
    seg.mat.patient$EndPosition <- as.numeric(as.character(seg.mat.patient$EndPosition))
    seg.mat.patient$nr.probes <- as.numeric(as.character(seg.mat.patient$nr.probes))
    chrom.length.copy <- fun.chrom.length(seg.mat.patient[, 2:4])
    chrom.segs <- fun.add.chrom(seg.mat.patient[, 2:4], chrom.length.copy)
    seg.mat.plot <- seg.mat.patient
    seg.mat.plot[, 2:4] <- chrom.segs
    major <- seg.mat.plot$Copy.Number - seg.mat.plot$min.allele
    seg.mat.plot <- cbind(seg.mat.plot, major)
    seg.mat.plot <- seg.mat.plot[seg.mat.plot$nr.probes >= min.probes, 
        ]
    TCGA.plot <- TCGA.earlyLate
    Chromosome <- as.numeric(do.call(rbind, strsplit(unlist(TCGA.earlyLate$mutation_id), 
        split = ":"))[, 2])
    Start_pos <- as.numeric(do.call(rbind, strsplit(unlist(TCGA.earlyLate$mutation_id), 
        split = ":"))[, 3])
    TCGA.plot <- cbind(TCGA.plot, Chromosome, Start_pos)
    TCGA.plot <- data.frame(apply(TCGA.plot, 2, unlist), stringsAsFactors = FALSE)
    min.x <- 0
    min.y <- -0.25
    max.y <- max.cpn + 1
    max.x <- as.numeric(max(seg.mat.plot$EndPosition))
    seg.mat.plot$major <- ifelse(seg.mat.plot$major > max.cpn, 
        max.cpn, seg.mat.plot$major)
    seg.mat.plot$min.allele <- ifelse(seg.mat.plot$min.allele > 
        max.cpn, max.cpn, seg.mat.plot$min.allele)
    TCGA.plot$mut.multi <- ifelse(TCGA.plot$mut.multi > max.cpn, 
        max.cpn, TCGA.plot$mut.multi)
    plot(1, xlim = c(min.x, max.x), ylim = c(min.y, max.y), type = "n", 
        xaxt = "n", yaxs = "i", yaxt = "n", xaxs = "i", ylab = "", 
        xlab = "", lwd = 2, bty = "n")
    box(lwd = 0.5)
    axis(side = 2, at = seq(0, max.cpn, by = 1), labels = c(seq(0, 
        max.cpn - 1, by = 1), paste(">", max.cpn, sep = "")), 
        las = 2, cex.axis = 0.7, lwd = 0.5)
    mtext(text = TCGA.barcode, side = 3, line = 2)
    fun.plot.segment <- function(x, seg.mat.plot) {
        x0 <- as.numeric(seg.mat.plot[x, "StartPosition"])
        x1 <- as.numeric(seg.mat.plot[x, "EndPosition"])
        y0 <- seg.mat.plot[x, "major"] + 0.05
        y1 <- y0
        segments(x0, y0, x1, y1, col = "black", lwd = 2.4)
        x0 <- as.numeric(seg.mat.plot[x, "StartPosition"])
        x1 <- as.numeric(seg.mat.plot[x, "EndPosition"])
        y0 <- as.numeric(seg.mat.plot[x, "min.allele"]) - 0.05
        y1 <- y0
        segments(x0, y0, x1, y1, col = "#009E73", lwd = 2.4)
    }
    sapply(1:nrow(seg.mat.plot), fun.plot.segment, seg.mat.plot)
    for (i in sort(unique(as.numeric(seg.mat.plot$Chromosome)))) {
        abline(v = max(as.numeric(seg.mat.plot[seg.mat.plot$Chromosome == 
            i, "EndPosition"])), lwd = 0.5)
    }
    TCGA.plot <- TCGA.plot[order(as.numeric(TCGA.plot$Chromosome), 
        as.numeric(TCGA.plot$Start_pos)), ]
    TCGA.plot$Start_pos <- as.numeric(fun.add.chrom(cbind(TCGA.plot$Chromosome, 
        TCGA.plot$Start_pos, TCGA.plot$Start_pos), chrom.length.copy)[, 
        2])
    early.muts <- TCGA.plot[TCGA.plot$timing \%in\% c("early"), 
        ]
    early.clonal <- early.muts[early.muts$absoltue.ccf.0.95 >= 
        sub.clonal, ]
    early.subclonal <- early.muts[early.muts$absoltue.ccf.0.95 < 
        sub.clonal, ]
    late.muts <- TCGA.plot[TCGA.plot$timing \%in\% c("late"), ]
    late.clonal <- late.muts[late.muts$absoltue.ccf.0.95 >= sub.clonal, 
        ]
    late.subclonal <- late.muts[late.muts$absoltue.ccf.0.95 < 
        sub.clonal, ]
    nontimed.muts <- TCGA.plot[!TCGA.plot$mutation_id \%in\% c(early.muts$mutation_id, 
        late.muts$mutation_id), ]
    nontimed.clonal <- nontimed.muts[nontimed.muts$absoltue.ccf.0.95 >= 
        sub.clonal, ]
    nontimed.subclonal <- nontimed.muts[nontimed.muts$absoltue.ccf.0.95 < 
        sub.clonal, ]
    plot.earlyORlateCol <- function(x, timed.muts, nontimed = FALSE) {
        mut <- timed.muts[x, , drop = FALSE]
        mut.multiplicity <- mut$mut.multi
        if (nontimed) {
            if (as.numeric(mut$mut.multi.bstr.0.95) < sub.clonal) {
                points(mut$Start_pos, mut.multiplicity, cex = 0.7, 
                  pch = 17, col = "#99999965")
            }
            if (as.numeric(mut$mut.multi.bstr.0.95) >= sub.clonal) {
                points(mut$Start_pos, mut.multiplicity, cex = 0.8, 
                  pch = 16, col = "#99999995")
            }
        }
        if (!nontimed) {
            if (mut$timing == "early") {
                if (as.numeric(mut$mut.multi.bstr.0.95) < sub.clonal) {
                  points(mut$Start_pos, mut.multiplicity, cex = 0.7, 
                    pch = 17, col = "#0072B245")
                }
                if (as.numeric(mut$mut.multi.bstr.0.95) >= sub.clonal) {
                  points(mut$Start_pos, mut.multiplicity, cex = 0.9, 
                    pch = 16, col = "#0072B295")
                }
            }
            if (mut$timing == "late") {
                if (as.numeric(mut$mut.multi.bstr.0.95) < sub.clonal) {
                  points(mut$Start_pos, mut.multiplicity, cex = 0.7, 
                    pch = 17, col = "#D55E0045")
                }
                if (as.numeric(mut$mut.multi.bstr.0.95) >= sub.clonal) {
                  points(mut$Start_pos, mut.multiplicity, cex = 0.9, 
                    pch = 16, col = "#D55E0095")
                }
            }
        }
    }
    if (nrow(early.muts) >= 1) {
        sapply(1:nrow(early.muts), plot.earlyORlateCol, early.muts)
    }
    if (nrow(late.muts) >= 1) {
        sapply(1:nrow(late.muts), plot.earlyORlateCol, late.muts)
    }
    if (nrow(nontimed.muts) >= 1) {
        sapply(1:nrow(nontimed.muts), plot.earlyORlateCol, nontimed.muts, 
            nontimed = TRUE)
    }
    mtext(side = 3, at = fun.chrom.mean(seg.mat.plot[, 2:4]), 
        text = sort(unique(seg.mat.plot[, 2])), cex = seq(0.6, 
            0.4, length.out = length(unique(seg.mat.plot[, 2]))), 
        line = -1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
