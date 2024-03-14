print.working_cor_mat <- function(x) {
  cat("Structure Type:", x$cor_structure_type, "\n")
  cat("Dimensions:", dim(x$correlation_matrix), "\n")
  cat("Sample Type:", x$sample_type, "\n")
  cat("Working Correlation Matrix\n")
  invisible(x)
}
