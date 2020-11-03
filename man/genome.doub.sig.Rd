\name{genome.doub.sig}
\alias{genome.doub.sig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genome Doubling function
}
\description{
The following function estimates the probability that a genome doubling has occured at some point during the evolutionary history of a tumour. 
}
\usage{
genome.doub.sig(sample, seg.mat.minor, seg.mat.copy, number.of.sim = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample}{
TCGA sample identifier
}
  \item{seg.mat.minor}{
segmented minor allele copy number matrix 
}
  \item{seg.mat.copy}{
segmented total copy number matrix
}
  \item{number.of.sim}{
The number of simulations to assess genome doubling likelihood, defaults to 10,000
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
%%  Nichlas McGranahan
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
function (sample, seg.mat.minor, seg.mat.copy, number.of.sim = 10000) 
{
    print(sample)
    require(limma)
    sub.minor <- subset(seg.mat.minor, seg.mat.minor[, 1] == 
        sample)
    sub.major <- subset(seg.mat.copy, seg.mat.copy[, 1] == sample)
    chr.names <- c("1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", 
        "5", "5.5", "6", "6.5", "7", "7.5", "8", "8.5", "9", 
        "9.5", "10", "10.5", "11", "11.5", "12", "12.5", "13.5", 
        "14.5", "15.5", "16", "16.5", "17", "17.5", "18", "18.5", 
        "19", "19.5", "20", "20.5", "21.5", "22.5")
    if (!identical(unique(sub.minor[, 2]), chr.names)) {
        print("expected chr.names:")
        print(chr.names)
        stop(c(("seg.mat.minor chr.names!= expected chr.names")))
    }
    if (!identical(unique(sub.major[, 2]), chr.names)) {
        print("expected chr.names:")
        print(chr.names)
        stop(c(("seg.mat.copy chr.names!= expected chr.names")))
    }
    chr.arm.ploidy.minor <- c()
    for (chr.arm in unique(sub.minor[, 2])) {
        sub.chr.minor <- rbind(subset(sub.minor, sub.minor[, 
            2] == chr.arm)[, 2:6])
        sub.chr.minor <- apply(sub.chr.minor, 2, as.numeric)
        sub.chr.minor <- rbind(sub.chr.minor)
        if (length(unique(sub.chr.minor[, 5])) == 1) {
            arm.ploidy <- unique(sub.chr.minor[, 5])
        }
        else if (length(unique(sub.chr.minor[, 5])) > 1) {
            arm.ploidy <- weighted.median(sub.chr.minor[, 5], 
                w = sub.chr.minor[, 3] - sub.chr.minor[, 2], 
                na.rm = T)
        }
        chr.arm.ploidy.minor <- c(chr.arm.ploidy.minor, arm.ploidy)
    }
    names(chr.arm.ploidy.minor) <- unique(sub.minor[, 2])
    chr.arm.ploidy.major <- c()
    for (chr.arm in unique(sub.major[, 2])) {
        sub.chr.major <- rbind(subset(sub.major, sub.major[, 
            2] == chr.arm)[, 2:6])
        sub.chr.major <- apply(sub.chr.major, 2, as.numeric)
        sub.chr.major <- rbind(sub.chr.major)
        if (length(unique(sub.chr.major[, 5])) == 1) {
            arm.ploidy <- unique(sub.chr.major[, 5])
        }
        else if (length(unique(sub.chr.major[, 5])) > 1) {
            arm.ploidy <- weighted.median(sub.chr.major[, 5], 
                w = sub.chr.major[, 3] - sub.chr.major[, 2], 
                na.rm = T)
        }
        chr.arm.ploidy.major <- c(chr.arm.ploidy.major, arm.ploidy)
    }
    names(chr.arm.ploidy.major) <- unique(sub.major[, 2])
    chr.arm.ploidy.major <- chr.arm.ploidy.major - chr.arm.ploidy.minor
    total.aber <- sum(abs(chr.arm.ploidy.minor - 1)) + sum(abs(chr.arm.ploidy.major - 
        1))
    chr.probs <- c()
    if (sum(total.aber) == 0) {
        p.val.genome.doubl <- c(1)
    }
    if (sum(total.aber) != 0) {
        for (chr.arm in chr.names) {
            chr.prob.A <- (chr.arm.ploidy.major[chr.arm] - 1)
            chr.prob.A <- c(sum(chr.prob.A[chr.prob.A > 0])/total.aber, 
                abs(sum(chr.prob.A[chr.prob.A < 0])/total.aber))
            names(chr.prob.A) <- c(paste(chr.arm, "_Again", sep = ""), 
                paste(chr.arm, "_Aloss", sep = ""))
            chr.prob.B <- (chr.arm.ploidy.minor[chr.arm] - 1)
            chr.prob.B <- c(sum(chr.prob.B[chr.prob.B > 0])/total.aber, 
                abs(sum(chr.prob.B[chr.prob.B < 0])/total.aber))
            names(chr.prob.B) <- c(paste(chr.arm, "_Bgain", sep = ""), 
                paste(chr.arm, "_Bloss", sep = ""))
            chr.prob <- c(chr.prob.A, chr.prob.B)
            chr.probs <- c(chr.probs, chr.prob)
        }
        prop.major.even.obs <- length(which(chr.arm.ploidy.major >= 
            2))/length(chr.arm.ploidy.major)
        prop.major.even.sim <- c()
        k <- 1
        while (k <= number.of.sim) {
            chr.sim <- table(sample(names(chr.probs), total.aber, 
                prob = chr.probs, replace = T))
            chr.sim.table <- chr.probs
            chr.sim.table <- rep(0, length(chr.probs))
            names(chr.sim.table) <- names(chr.probs)
            chr.sim.table[names(chr.sim)] <- chr.sim
            chr.sim.table <- cbind(chr.sim.table[seq(1, length(chr.sim.table), 
                by = 4)], chr.sim.table[seq(2, length(chr.sim.table), 
                by = 4)], chr.sim.table[seq(3, length(chr.sim.table), 
                by = 4)], chr.sim.table[seq(4, length(chr.sim.table), 
                by = 4)])
            rownames(chr.sim.table) <- chr.names
            colnames(chr.sim.table) <- c("gain.A", "loss.A", 
                "gain.B", "loss.B")
            chr.sim.major <- apply(cbind(c(1 + chr.sim.table[, 
                1] - chr.sim.table[, 2]), c(1 + chr.sim.table[, 
                3] - chr.sim.table[, 4])), 1, max)
            prop.major.even <- length(which(chr.sim.major >= 
                2))/length(chr.sim.major)
            prop.major.even.sim <- c(prop.major.even.sim, prop.major.even)
            k <- k + 1
        }
        p.val.genome.doubl <- length(which(prop.major.even.sim >= 
            prop.major.even.obs))/number.of.sim
        if (prop.major.even.obs == 1) {
            p.val.genome.doubl = 0
        }
    }
    return(p.val.genome.doubl)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
