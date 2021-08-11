######


# pval_thresh <- .05
# logfc_thresh <- 2
# de_vec <- 
######



## NOTES ##
# - hardcoding is.DE for now because it is already generated. Working version will take
#   p value and lfc considered significant as an argument and create is.DE column.

plotVolcano <- function(data, 
                        logfc_col, 
                        pval_col, 
                        gene_col, 
                        pval_thresh, 
                        logfc_thresh,
                        de_vec,
                        show_logfc_thresh,
                        show_pvalue_thresh,
                        thresh_color) {
  # check that columns exist
  if (!all(c(logfc_col, pval_col, gene_col) %in% colnames(data))) {
    stop("provided column names do not match dataset")
  }
  # convert pval to -log10(pval)
  data <- mutate(data,
                 log_pval = -log10(data[[pval_col]]))
  # build base of plot plot
  volcano <- ggplot(data, aes(x = .data[[logfc_col]], y = log_pval, color = de_vec)) +
    geom_point(alpha = .7)
  # if show_logfc_thresh = true add hline layer
  if (show_logfc_thresh) {
    volcano <- volcano + 
      geom_hline(yintercept = -log10(pval_thresh), col = thresh_color)
  }
  # if show_pvalue_thresh = true add vline layer
  if (show_pvalue_thresh) {
    volcano <- volcano + 
      geom_vline(xintercept = c(logfc_thresh, -logfc_thresh), col = thresh_color)
  }
  # add minimal theme
  volcanoPlot <- volcano + 
    theme_minimal()
  # display plot
  volcanoPlot
}