#' Check whether the two chromosome vectors are the same.
#'
#' @param input_chrom The vector to be checked (as A).
#' @param target_chrom The target vector (as B).
#'
#' @return A == B: "identical"; A < B: "missing"; A > B: "extra"; as well as "other".
#' @export
#'
#' @examples
#' 
#' check_chromosomal_names(input_chrom = LETTERS[1:5], target_chrom = LETTERS[1:26])
#' check_chromosomal_names(input_chrom = LETTERS[1:9], target_chrom = LETTERS[1:9])
#' 
check_chromosomal_names <- function(input_chrom, target_chrom) {
  input_chrom <- unique(as.character(input_chrom))
  target_chrom <- unique(as.character(target_chrom))
  # Check Relationships:
  all_in <- all(input_chrom %in% target_chrom)
  contain_all <- all(target_chrom %in% input_chrom)
  res <- dplyr::case_when(
    all_in & contain_all ~ "identical"
    , all_in & !contain_all ~ "missing"
    , !all_in & contain_all ~ "extra"
    , !all_in & !contain_all ~ "other"
  )
  return(res)
}
