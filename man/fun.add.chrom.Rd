\name{fun.add.chrom}
\alias{fun.add.chrom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
add chromosomes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fun.add.chrom(seg, chrom.length)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
%%     ~~Describe \code{seg} here~~
}
  \item{chrom.length}{
%%     ~~Describe \code{chrom.length} here~~
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
function (seg, chrom.length) 
{
    if (1 \%in\% seg[, 1]) 
        seg.tmp <- subset(seg, seg[, 1] == 1)
    else seg.tmp <- matrix(NA, nr = 0, nc = ncol(seg))
    for (i in 2:max(max(as.numeric(seg[, 1])), 3)) {
        if (!(i \%in\% seg[, 1])) 
            (next)()
        sub.chrom <- subset(seg, seg[, 1] == i)
        sub.chrom[, 2] <- as.numeric(sub.chrom[, 2]) + sum(as.numeric(chrom.length[1:(i - 
            1)]))
        sub.chrom[, 3] <- as.numeric(sub.chrom[, 3]) + sum(as.numeric(chrom.length[1:(i - 
            1)]))
        seg.tmp <- rbind(seg.tmp, sub.chrom)
    }
    return(seg.tmp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
