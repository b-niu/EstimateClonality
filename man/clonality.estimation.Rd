\name{clonality.estimation}
\alias{clonality.estimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Clonality Estimation Function
}
\description{
  clonality.estimation integrates the read depth and copy number information at
  single nucleotide variant locations to determine the timing and clonality of
  each mutation.
}
\usage{
clonality.estimation(mutation.table.loc, seg.mat.loc, data.type, TCGA.barcode, ANALYSIS.DIR, sub.clonal.cut.off = 1, min.var.prop = 0.05, min.alt.reads = 5, min.depth = 30, plotting = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mutation.table.loc}{
    The location of a mutation table. 
    The table should contain the following columns:
                      c('Patient'
                      ,'Chr'
                      ,'Start_position'
                      ,'End_position'
                      ,'Reference'
                      ,'Alternate'
                      ,'Variant_freq'
                      ,'Ref_freq'
                      ,'Hugo_Symbol'
                      ,'Variant_Classification')
}
  \item{seg.mat.loc}{
The location of a segmented matrix (in RData format). This should be in the format of a list. 
}
  \item{data.type}{
Data type, e.g. TCGA_BRCA
}
  \item{TCGA.barcode}{
Sample identifier, e.g. 'TCGA-AN-A0FF'
}
  \item{ANALYSIS.DIR}{
Where do you want the folders to be saved?
}
  \item{sub.clonal.cut.off}{
What cut off should be used for the subclonal fraction, defaults to 1 (i.e. 95\% CI overlaps 1)
}
  \item{min.var.prop}{
Minimum variant proportion (defaults to 0.05)
}
  \item{min.alt.reads}{
Minimum number of alternate reads (defaults to 5)
}
  \item{min.depth}{
Minimum depth at mutated base (defaults to 30)
}
  \item{plotting}{
Do you want plotting? (defaults to TRUE)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Nothing is returned directly. However, three files (if plotting TRUE) should be returned within a folder in the specified directory. Patient.earlylate.tsv is a tab-separated table that contains a number of columns. Description of the columns is as follows:
 [1] "patient"  #Patient identifier;
 [2] "TCGA.purity"  # ASCAT tumor purity estimate;
 [3] "mutation_id" # mutation id;
 [4] "Reference_Base" # reference base; 
 [5] "Alternate_Base" # alternate base;
 [6] "ref_counts" # number of reads supporting reference base;
 [7] "var_counts" # number of reads supporting variant base;
 [8] "obs.VAF" # observed variant allele frequency;
 [9] "normal_cn" # local normal copy number at mutation site;
 [10] "minor_cn" # local minor allele copy number at mutation site;
 [11] "major_cn" # local major allele copy number at mutation site;
 [12] "mut.multi" # mutation multiplicity (number of copies of mutant allele);
[13] "absolute.ccf" # estimated cancer cell fraction;
[14] "absolute.ccf.0.05" # CI of ccf;
[15] "absolute.ccf.0.95" # CI of ccf;
[16 "prob.clonal";      
[17] "prob.subclonal";    
[18] "GD.status";        
[19] "timing" # timing based on copy number;
[20] "comb.timing" # timing based on clonality and copy number;
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicholas McGranahan
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{EstimateClonality-package},\link{earlyORlate}} 
}
\examples{
##clonality.estimation(mutation.table.loc="BRCA.mutation.table.txt"
##                     ,seg.mat.loc="tcga.brca.seg.hg19.rdata"
##                     ,data.type='TCGA_BRCA'
##                     ,TCGA.barcode="TCGA-AR-A2LR"
##                     ,ANALYSIS.DIR="example/")
}

