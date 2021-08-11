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
                        highlight_genes = NULL) {
  # check that columns exist
  if (!all(c(logfc_col, pval_col, gene_col) %in% colnames(data))) {
    stop("provided column names do not match dataset")
  }
  
  # convert pval to -log10(pval)
  data <- mutate(data,
                 log_pval = -log10(data[[pval_col]]))
  
  # build base of plot
  volcano <- ggplot(data, aes(x = .data[[logfc_col]], y = log_pval)) +
    geom_point(alpha = .6, aes(color = de_vec))
  
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
  
  # if highlight_genes isn't null, add to plot
  if (!is.null(highlight_genes)) {
    # create vector of gene ids for label
    labels <- data[data[[gene_col]] %in% highlight_genes,]
    volcano <- suppressWarnings(volcano +
        ggrepel::geom_label_repel(data = data[data[[gene_col]] %in% highlight_genes,], aes(label = labels[[gene_col]]))) 
  }
  
  # add minimal theme
  volcanoPlot <- volcano + 
    theme_minimal()
  
  # display plot
  volcanoPlot
}