######
# data <- read.table("./data/210709_KO_vs_WT.txt", sep = "\t", header = T)
log_col <- "logFC"
fdr_col <- "FDR"
gene_col <- "Genename"

######

## NOTES ##
# - hardcoding is.DE for now because it is already generated. Working version will take
#   p value and lfc considered significant as an argument and create is.DE column.

plotVolcano <- function(data, log_col, fdr_col, gene_col) {
  # check that columns exist
  if (!all(c(log_col, fdr_col, gene_col) %in% colnames(data))) {
    stop("provided column names do not match dataset")
  }
  # convert pval to -log10(pval)
  data <- mutate(data,
         log_fdr = -log10(.data[[fdr_col]]),
         is.DE_fac = as.factor(is.DE))
  # build plot
  ggplot(data, aes(x = .data[[log_col]], y = log_fdr, color = is.DE_fac)) +
    geom_point(alpha = .7)+
    geom_hline(yintercept = -log10(0.01), col = "red") +
    theme_minimal()
  #ggplotly(v, tooltip = "Genename")
}
