\name{plot.TCGA.sample.ccfs}
\alias{plot.TCGA.sample.ccfs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot TCGA sample ccf
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.TCGA.sample.ccfs(TCGA.earlyLate, clonal.cut.off = 1, spec.genes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TCGA.earlyLate}{
%%     ~~Describe \code{TCGA.earlyLate} here~~
}
  \item{clonal.cut.off}{
%%     ~~Describe \code{clonal.cut.off} here~~
}
  \item{spec.genes}{
%%     ~~Describe \code{spec.genes} here~~
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
function (TCGA.earlyLate, clonal.cut.off = 1, spec.genes = NULL) 
{
    cancer.patient <- data.frame(TCGA.earlyLate)
    layout(rbind(1, 2))
    par(mar = c(0.5, 5, 5, 5))
    plot(density(as.numeric(cancer.patient$ccf), bw = "SJ"), 
        xlim = c(0, 1.4), xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "Density (a.u.)", main = cancer.patient$patient[1], 
        lwd = 3)
    par(mar = c(5, 5, 0.5, 5))
    plot(as.numeric(cancer.patient$ccf), log(as.numeric(cancer.patient$var_counts) + 
        as.numeric(cancer.patient$ref_counts)), ylim = c(log(5), 
        log(1000)), xlim = c(0, 1.6), xaxt = "n", yaxt = "n", 
        xlab = "Adjusted VAF", ylab = "Tumour Coverage", main = "", 
        lwd = 3, cex = 1.5, pch = 21, col = ifelse(as.numeric(cancer.patient$absolute.ccf.0.95) >= 
            1, "#de2d2699", "#0072B299"), bg = ifelse(cancer.patient$absolute.ccf.0.95 >= 
            1, "#de2d2699", "#0072B299"))
    if (!is.null(spec.genes)) {
        for (gene in spec.genes) {
            if (!gene \%in\% cancer.patient$Hugo_Symbol) {
                next
            }
            gene.mut <- cancer.patient[cancer.patient$Hugo_Symbol \%in\% 
                gene, , drop = FALSE]
            gene.mut <- gene.mut[gene.mut$Variant_Classification \%in\% 
                c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site", 
                  "Frame_Shift_Del", "Translation_Start_Site", 
                  "Nonstop_Mutation"), , drop = FALSE]
            if (nrow(gene.mut) == 0) {
                next
            }
            points(gene.mut$ccf, log(gene.mut$var_counts + gene.mut$ref_counts), 
                cex = 1.8, lwd = 3)
            text(gene.mut$ccf, log(gene.mut$var_counts + gene.mut$ref_counts), 
                labels = gene, pos = 4, cex = 1, lwd = 3, font = 2)
        }
    }
    axis(side = 1, at = seq(0, 1.6, by = 0.2))
    cov.seq <- c(5, 10, 20, 50, 100, 200, 500, 1000)
    axis(side = 2, at = log(cov.seq), labels = cov.seq, las = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
