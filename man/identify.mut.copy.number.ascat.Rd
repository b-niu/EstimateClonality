\name{identify.mut.copy.number.ascat}
\alias{identify.mut.copy.number.ascat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify mutation local copy number
}
\description{
This function identifies local copy number for a given mutation
}
\usage{
identify.mut.copy.number.ascat(x, sub.mat.mut, sub.mat.copy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{sub.mat.mut}{
%%     ~~Describe \code{sub.mat.mut} here~~
}
  \item{sub.mat.copy}{
%%     ~~Describe \code{sub.mat.copy} here~~
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
function (x, sub.mat.mut, sub.mat.copy) 
{
    mut <- sub.mat.mut[x, , drop = FALSE]
    ww <- which(as.numeric(sub.mat.copy$Chr) == as.numeric(mut$Chr) & 
        as.numeric(sub.mat.copy$Start) <= as.numeric(mut$Start_position) & 
        as.numeric(sub.mat.copy$End) >= as.numeric(mut$Start_position))
    copy <- sub.mat.copy[ww, , drop = FALSE]
    mutation_id <- paste(mut$Patient, mut$Chr, mut$Start_position, 
        mut$Reference, sep = ":")
    ref_counts <- mut$Ref_freq
    var_counts <- mut$Variant_freq
    normal_cn <- 2
    patient <- mut$Patient
    Reference_Base <- mut$Reference
    Alternate_Base <- mut$Alternate
    if (nrow(copy) != 1) {
        minor_cn <- NA
        major_cn <- NA
        output <- data.frame(mutation_id, ref_counts, var_counts, 
            normal_cn, minor_cn, major_cn, patient, Reference_Base, 
            Alternate_Base, stringsAsFactors = FALSE)
        return(output)
    }
    minor_cn <- min(c(copy$nA, copy$nB))
    major_cn <- max(c(copy$nA, copy$nB))
    output <- data.frame(mutation_id, ref_counts, var_counts, 
        normal_cn, minor_cn, major_cn, patient, Reference_Base, 
        Alternate_Base, stringsAsFactors = FALSE)
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
