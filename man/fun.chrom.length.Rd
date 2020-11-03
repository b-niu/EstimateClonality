\name{fun.chrom.length}
\alias{fun.chrom.length}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
check chromosome length
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fun.chrom.length(seg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seg}{
%%     ~~Describe \code{seg} here~~
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
function (seg) 
{
    chrom.length <- c()
    for (i in sort(as.numeric(unique(seg[, 1])))) {
        sub.chrom <- subset(seg, seg[, 1] == i)
        chrom.length <- c(chrom.length, sub.chrom[nrow(sub.chrom), 
            3])
    }
    return(chrom.length)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
