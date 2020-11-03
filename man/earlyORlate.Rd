\name{earlyORlate}
\alias{earlyORlate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clonal and temporal dissection
}
\description{
This function implementents the cancer cell fraction analysis, as well as temporal dissection of mutations
}
\usage{
earlyORlate(patient, complete.mutation.table, purity)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{patient}{
patient describes the TCGA identifier 
}
  \item{complete.mutation.table}{
This is the mutation table in matrix format. 
}
  \item{purity}{
ASCAT purity estimate [range, 0.01,1]
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
Nicholas McGranahan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# see example.script for how to use the code
}