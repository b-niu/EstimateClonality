\name{fun.GD.status}
\alias{fun.GD.status}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genome Doubling function (part 2)
}
\description{
The following function interprets the p-value from genome.doub.sig, to give a binary variable, GD (genome-doubled), or nGD (non-genome-doubled)
}
\usage{
fun.GD.status(GD.pval, ploidy.val)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GD.pval}{
The p-value obtained from genome.doub.sig function
}
  \item{ploidy.val}{
The estimated ploidy of the tumor sample
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
function (GD.pval, ploidy.val) 
{
    GD.status <- c(NA)
    if (ploidy.val <= 2) {
        GD.status <- ifelse(GD.pval <= 0.001, "GD", "nGD")
    }
    if (ploidy.val == 3) {
        GD.status <- ifelse(GD.pval <= 0.001, "GD", "nGD")
    }
    if (ploidy.val == 4) {
        GD.status <- ifelse(GD.pval <= 0.05, "GD", "nGD")
    }
    if (ploidy.val == 5) {
        GD.status <- ifelse(GD.pval <= 0.5, "GD", "nGD")
    }
    if (ploidy.val == 6) {
        GD.status <- ifelse(GD.pval <= 1, "GD", "nGD")
    }
    return(GD.status)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
