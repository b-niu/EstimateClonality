\name{earlyORlate.strict}
\alias{earlyORlate.strict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
strinct early or late function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
earlyORlate.strict(patient, complete.mutation.table, purity)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{patient}{
%%     ~~Describe \code{patient} here~~
}
  \item{complete.mutation.table}{
%%     ~~Describe \code{complete.mutation.table} here~~
}
  \item{purity}{
%%     ~~Describe \code{purity} here~~
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
##--  or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (patient, complete.mutation.table, purity) 
{
    suppressPackageStartupMessages(library(sequenza))
    suppressPackageStartupMessages(library(bootstrap))
    suppressPackageStartupMessages(library(boot))
    mut.table <- complete.mutation.table[complete.mutation.table$patient == 
        patient, ]
    cellularity <- as.numeric(purity)
    major.cn <- unlist(mut.table$major_cn)
    abs.cn <- unlist(mut.table$minor_cn) + unlist(mut.table$major_cn)
    minor.cn <- unlist(mut.table$minor_cn)
    depth.t <- unlist(mut.table$ref_counts) + unlist(mut.table$var_counts)
    max.cn <- max(abs.cn)
    VAF <- unlist(mut.table$var_counts)/(unlist(mut.table$var_counts) + 
        unlist(mut.table$ref_counts))
    types <- types.matrix(CNt.min = 1, CNt.max = max.cn, CNn = 2)
    types.xy <- types.matrix(CNt.min = 1, CNt.max = max.cn, CNn = 1)
    types <- rbind(types, types.xy)
    types <- types[types$Mt >= 1, ]
    types$F <- 0
    for (i in 1:nrow(types)) {
        types$F[i] <- theoretical.mufreq(cellularity = cellularity, 
            CNn = types$CNn[i], CNt = types$CNt[i], Mt = types$Mt[i])
    }
    get.Mt <- function(F, depth.t, types, CNt, CNn, Mt) {
        types <- types[types$CNn == CNn, ]
        l <- sequenza:::mufreq.dpois(mufreq = F, types$F[types$CNt == 
            CNt & types$Mt <= Mt], depth.t = depth.t)
        l <- l/sum(l)
        L <- data.frame(l = l, Mt = types$Mt[types$CNt == CNt & 
            types$Mt <= Mt])
    }
    get.conf <- function(F, depth.t) {
        conf.int <- cbind(prop.test(round(F * depth.t, 0), depth.t)$conf[1], 
            prop.test(round(F * depth.t, 0), depth.t)$conf[2])
        return(conf.int)
    }
    bootstrap.cf <- function(Vaf, cellularity, CNn, CNt, depth.t) {
        if (Vaf == 1) {
            conf.int <- cbind(prop.test(round(Vaf * depth.t, 
                0), depth.t)$conf[1], prop.test(round(Vaf * depth.t, 
                0), depth.t)$conf[2])
            lower <- get.mut.mult(Vaf = conf.int[1], cellularity = cellularity, 
                CNt = CNt, CNn = CNn)
            higher <- get.mut.mult(Vaf = conf.int[2], cellularity = cellularity, 
                CNt = CNt, CNn = CNn)
            conf.int <- cbind(lower, higher)
            return(conf.int)
        }
        x <- c(rep(1, round(Vaf * depth.t, 0)), rep(0, (depth.t - 
            round(Vaf * depth.t, 0))))
        theta <- function(x, i) {
            data <- x[i]
            est <- sum(data)/length(data)
            mut.multi <- (est * 1/cellularity) * ((cellularity * 
                CNt) + CNn * (1 - cellularity))
            return(mut.multi)
        }
        bt.res <- boot(x, theta, R = 1000)
        bt.ci <- boot.ci(bt.res, type = "norm")
        out <- c(bt.ci$normal[2], bt.ci$normal[3])
        return(out)
    }
    get.mut.mult <- function(CNt, Vaf, cellularity, CNn) {
        return((Vaf * 1/cellularity) * ((cellularity * CNt) + 
            CNn * (1 - cellularity)))
    }
    get.cancer.cell.fraction <- function(Max.Likelihood, mut.mult) {
        predicted.Mtn <- Max.Likelihood[, "Mt"]
        ccf <- mut.mult/predicted.Mtn
        return(ccf)
    }
    get.all.mut.info <- function(i) {
        obs.VAF <- VAF[i]
        mut.conf.0.05 <- get.conf(F = VAF[i], depth.t = depth.t[i])[1]
        mut.conf.0.95 <- get.conf(F = VAF[i], depth.t = depth.t[i])[2]
        if (abs.cn[i] == 0) {
            output <- cbind(obs.VAF, mut.conf.0.05, mut.conf.0.95, 
                mut.multi = NA, mut.multi.0.05 = NA, mut.multi.bstr.0.05 = NA, 
                mut.multi.0.95 = NA, mut.multi.bstr.0.95 = NA, 
                Exp.Cpn = NA, Exp.Cpn.Likelihood = NA, ccf = NA, 
                ccf.0.05 = NA, ccf.btstr.0.05 = NA, ccf.0.95 = NA, 
                ccf.btstr.0.95 = NA, timing = "Not.Poss")
            return(output)
        }
        L <- get.Mt(F = VAF[i], depth.t = depth.t[i], CNt = abs.cn[i], 
            types = types, CNn = unlist(mut.table$normal_cn[i]), 
            Mt = major.cn[i])
        mut.multi <- get.mut.mult(CNt = abs.cn[i], Vaf = VAF[i], 
            cellularity = cellularity, CNn = unlist(mut.table$normal_cn[i]))
        mut.multi.0.05 <- get.mut.mult(CNt = abs.cn[i], Vaf = mut.conf.0.05, 
            cellularity = cellularity, CNn = unlist(mut.table$normal_cn[i]))
        mut.multi.0.95 <- get.mut.mult(CNt = abs.cn[i], Vaf = mut.conf.0.95, 
            cellularity = cellularity, CNn = unlist(mut.table$normal_cn[i]))
        mut.multi.bstr <- bootstrap.cf(Vaf = VAF[i], cellularity = cellularity, 
            CNn = unlist(mut.table$normal_cn[i]), CNt = abs.cn[i], 
            depth.t = depth.t[i])
        mut.multi.bstr.0.05 <- mut.multi.bstr[1]
        mut.multi.bstr.0.95 <- mut.multi.bstr[2]
        if (is.na(L$l)[1]) {
            output <- cbind(obs.VAF, mut.conf.0.05, mut.conf.0.95, 
                mut.multi, mut.multi.0.05, mut.multi.bstr.0.05, 
                mut.multi.0.95, mut.multi.bstr.0.95, Exp.Cpn = NA, 
                Exp.Cpn.Likelihood = NA, ccf = NA, ccf.0.05 = NA, 
                ccf.btstr.0.05 = NA, ccf.0.95 = NA, ccf.btstr.0.95 = NA, 
                timing = "Not.Poss")
            return(output)
        }
        Max.Likelihood <- L[which.max(L$l), ]
        ccf <- get.cancer.cell.fraction(Max.Likelihood, mut.multi)
        ccf.0.05 <- get.cancer.cell.fraction(Max.Likelihood, 
            mut.multi.0.05)
        ccf.btstr.0.05 <- get.cancer.cell.fraction(Max.Likelihood, 
            mut.multi.bstr.0.05)
        ccf.0.95 <- get.cancer.cell.fraction(Max.Likelihood, 
            mut.multi.0.95)
        ccf.btstr.0.95 <- get.cancer.cell.fraction(Max.Likelihood, 
            mut.multi.bstr.0.95)
        expected.copy <- Max.Likelihood[2]
        timing <- "Not.Poss"
        if (Max.Likelihood$Mt > 1) {
            timing <- "early"
        }
        if (Max.Likelihood$Mt <= 1 & minor.cn[i] != 1) {
            timing <- "late"
        }
        if (major.cn[i] <= 1) {
            timing <- "Not.Poss"
        }
        output <- cbind(obs.VAF, mut.conf.0.05, mut.conf.0.95, 
            mut.multi, mut.multi.0.05, mut.multi.bstr.0.05, mut.multi.0.95, 
            mut.multi.bstr.0.95, Exp.Cpn = Max.Likelihood$Mt, 
            Exp.Cpn.Likelihood = Max.Likelihood$l, ccf, ccf.0.05, 
            ccf.btstr.0.05, ccf.0.95, ccf.btstr.0.95, timing)
        return(output)
    }
    output <- t(sapply(1:nrow(mut.table), get.all.mut.info))
    colnames(output) <- c("obs.VAF", "mut.conf.0.05", "mut.conf.0.95", 
        "mut.multi", "mut.multi.0.05", "mut.multi.bstr.0.05", 
        "mut.multi.0.95", "mut.multi.bstr.0.95", "Exp.Cpn", "Exp.Cpn.Likelihood", 
        "ccf", "ccf.0.05", "ccf.btstr.0.05", "ccf.0.95", "ccf.btstr.0.95", 
        "timing")
    out <- cbind(mut.table, output)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
