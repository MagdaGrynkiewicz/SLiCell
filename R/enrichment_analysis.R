#' enrichment analysis for one potential SL gene pairs
#' enrichment set - geneA altered cell lines
#' geneB sensitivity is a ranking criterion
#' @param sl_pair geneA and geneB vector
#' @param sens_matrix geneB dependency matrix in CCLE cell lines for a given knock-down experiment
#' @param alt_matrix geneA alteration matrix for CCLE cell lines
#' @param plot=FALSE option for printing enrichment plot
#' @return list of vectors of enrichment score and its p-value
#' @export

fgsea_ES_and_pval <- function(sl_pair, sens_matrix1, alt_matrix1, plot=FALSE) {

  name_gene1 <- as.character(sl_pair[1])
  name_gene2 <- as.character(sl_pair[2])

  # checking if geneA is present in alteration matrix
  if (!(name_gene1 %in% colnames(alt_matrix1))) {
    #print("geneA is not present in alteration matrix")
    return(c(NA, NA))

    # checking if geneB is present in sensitivity matrix
  } else if (!(name_gene2 %in% colnames(sens_matrix1))) {
    #print("geneB is not present in sensitivity matrix")
    return(c(NA, NA))
  }

  # list of cell lines suitable for testing - present in both matrices
  common_cell_lines <- intersect(rownames(sens_matrix1), rownames(alt_matrix1))

  # geneA alteration list
  alt_gene1 <- alt_matrix1[common_cell_lines, name_gene1]
  names(alt_gene1) <- rownames(alt_matrix1[common_cell_lines,])
  # geneB sensitivity list
  sens_gene2 = -1*unlist(sens_matrix1[common_cell_lines,name_gene2])
  names(sens_gene2) <- rownames(sens_matrix1[common_cell_lines,])

  # cell lines with gene A alteration for plotting
  alt_cls_plot <- names(alt_gene1[alt_gene1==1])
  # list of cell line names with geneA alteration
  alt_cls <- list(alt_cls_plot)
  names(alt_cls) <- c(paste0(name_gene1, "_", name_gene2))


  defaultW <- getOption("warn")
  options(warn = -1)
  res <- fgsea::fgsea(alt_cls, sens_gene2, nperm=1000)
  if (plot==TRUE) {
    ESplot <- fgsea::plotEnrichment(alt_cls_plot, na.omit(sens_gene2))
    options(warn = defaultW)
    print(ESplot)
  }


    return(c(res$ES, res$pval))


}

#' enrichment analysis for list of potential SL gene pairs
#' enrichment set - geneA altered cell lines
#' geneB sensitivity is a ranking criterion
#' @param sl_pairs list of geneA geneB pairs for verifing
#' @param sens_matrix list of geneB dependency matrices in geneA altered cell lines, default - all experimental data preloaded in the package
#' @param data_set_names names for the dependency data sets - will be visible in results
#' @param alt_matrix1 geneA alteration matrix for cell lines, default - mutation_matrix preloaded in the package
#' @return list of enrichment scores and their p-values
#' @export

verify_enrichment <- function(sl_pairs,
  sens_matrix=list(SLiCell::sens_Drive_RSA, SLiCell::sens_Drive_Ataris, SLiCell::sens_Achilles_AC, SLiCell::sens_Achilles_Demeter),
  data_set_names=c("Drive_RSA", "Drive_Ataris", "Achilles_AC", "Achilles_Demeter"),
  alt_matrix1 = SLiCell::gene_cl_mut_matrix) {

  RSA <- apply(sl_pairs, 1, SLiCell::fgsea_ES_and_pval,
    sens_matrix[[1]], alt_matrix1)
  rownames(RSA) <- c("ES", "pval")
  Ataris <- apply(sl_pairs, 1, SLiCell::fgsea_ES_and_pval,
    sens_matrix[[2]], alt_matrix1)
  rownames(Ataris) <- c("ES", "pval")
  AC <- apply(sl_pairs, 1, SLiCell::fgsea_ES_and_pval,
    sens_matrix[[3]], alt_matrix1)
  rownames(AC) <- c("ES", "pval")
  Demeter <- apply(sl_pairs, 1, SLiCell::fgsea_ES_and_pval,
    sens_matrix[[4]], alt_matrix1)
  rownames(Demeter) <- c("ES", "pval")

  results = list(sl_pairs, RSA, Ataris, AC, Demeter)
  names(results) = c("sl_pair_list", data_set_names)

  return(results)
}


#' plots 4 enrichment plots
#' @param sl_pair list of geneA geneB pairs for verifing
#' @param sens_matrix list of geneB dependency matrices in geneA altered cell lines, default - all experimental data preloaded in the package
#' @param data_set_names names for the dependency data sets - will be visible in results
#' @param alt_matrix1 geneA alteration matrix for cell lines, default - mutation_matrix preloaded in the package
#' @return list of enrichment scores and their p-values
#' @export

plot_4_enrichments <- function(sl_pair,
  sens_matrix=list(SLiCell::sens_Drive_RSA, SLiCell::sens_Drive_Ataris, SLiCell::sens_Achilles_AC, SLiCell::sens_Achilles_Demeter),
  data_set_names=c("Drive_RSA", "Drive_Ataris", "Achilles_AC", "Achilles_Demeter"),
  alt_matrix1 = SLiCell::gene_cl_mut_matrix) {

  plts <- gridExtra::grid.arrange(SLiCell::fgsea_ES_and_pval(sl_pair, sens_matrix1 = sens_matrix[[1]], alt_matrix1 = alt_matrix1, plot=T),
    SLiCell::fgsea_ES_and_pval(sl_pair, sens_matrix1 = sens_matrix[[2]], alt_matrix1 = alt_matrix1, plot=T),
    SLiCell::fgsea_ES_and_pval(sl_pair, sens_matrix1 = sens_matrix[[3]], alt_matrix1 = alt_matrix1, plot=T),
    SLiCell::fgsea_ES_and_pval(sl_pair, sens_matrix1 = sens_matrix[[4]], alt_matrix1 = alt_matrix1, plot=T),
    ncol = 2, nrow = 2)
  print(plts)

}


