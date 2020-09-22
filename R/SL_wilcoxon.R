#' perform Wilcoxon for potential geneA - geneB synthetic lethal pair
#' check if geneA altered cell lines are more dependent on geneB
#' @param sl_pair c(geneA_ENSEMBL, geneB_ENSEMBL) from a list of potential SL gene pairs
#' @param sens_matrix geneB dependency matrix in CCLE cell lines for a given knock-down experiment
#' @param alt_matrix geneA alteration matrix for CCLE cell lines
#' @param min_alt_no minimum no of cell lines which have geneA alteration, default 10
#' @return p-value of the Wilcoxon test
#' @export

run_wilcoxon_sl_pair <- function(sl_pair, sens_matrix1, alt_matrix1, plot=F, experiment_name="", min_alt_no = 10) {

  name_gene1 <- as.character(sl_pair[1])
  name_gene2 <- as.character(sl_pair[2])

  # checking if geneA is present in alteration matrix
  if (!(name_gene1 %in% colnames(alt_matrix1))) {
    print("geneA is not present in alteration matrix")
    return(NA)

  # checking if geneB is present in sensitivity matrix
  } else if (!(name_gene2 %in% colnames(sens_matrix1))) {
    print("geneB is not present in sensitivity matrix")
    return(NA)
  }

  # list of cell lines suitable for testing - present in both matrices
  common_cell_lines <- intersect(rownames(sens_matrix1), rownames(alt_matrix1))

  # geneA alteration list
  alt_gene1 <- alt_matrix1[common_cell_lines, name_gene1]
  names(alt_gene1) <- rownames(alt_matrix1[common_cell_lines,])
  # geneB sensitivity list
  sens_gene2 = sens_matrix1[common_cell_lines,name_gene2]
  names(sens_gene2) <- rownames(sens_matrix1[common_cell_lines,])

  # list of cell line names with geneA alteration
  alt_cls <- names(alt_gene1[alt_gene1==1])

  # geneB sensitivity list in geneA altered cell lines
  sens_gene2_alt_gene1_1 = sens_gene2[alt_cls]

  if (length(sens_gene2_alt_gene1_1) < min_alt_no) {
    print("number of geneA altered cell lines is lower than required min")
    return(NA)
  }

  # list of cell line names withOUT geneA alteration
  not_alt_cls <- names(alt_gene1[alt_gene1==0])

  # geneB sensitivity list in geneA NOT altered cell lines
  sens_gene2_alt_gene1_0 = sens_gene2[not_alt_cls]

  res <- wilcox.test(sens_gene2_alt_gene1_1, sens_gene2_alt_gene1_0, alternative = c("less"))
  if (plot) {
    pl <- SLiCell::wf(sl_pair, sens=sens_matrix1, alt=alt_matrix1, pval=signif(res$p.value,3), experiment=experiment_name)
    print(pl)
  }
  return(res$p.value)


}





#' verify list of potential SL gene pairs
#' performs Wilcoxon test for each pair
#' where geneA which is altered in the minimum no of cell lines for Wilcoxon test
#' and geneB is present in the chosen knock-down experiment
#' @param sl_pair_list list of geneA geneB pairs for verifing
#' @param sens_matrix list of geneB dependency matrices in geneA altered cell lines, default - all experimental data preloaded in the package
#' @param data_set_names names for the dependency data sets - will be visible in results
#' @param alt_matrix1 geneA alteration matrix for cell lines, default - mutation_matrix preloaded in the package
#' @param min_alt_no minimum no of cell lines which have geneA alteration, default 10
#' @return list of Wilcoxon rank-sum test results separately for each dependency data set
#' @export

verify_sl_pairs <- function(sl_pair_list,
  sens_matrix=list(SLiCell::sens_Drive_RSA, SLiCell::sens_Drive_Ataris, SLiCell::sens_Achilles_AC, SLiCell::sens_Achilles_Demeter),
  data_set_names=c("Drive_RSA", "Drive_Ataris", "Achilles_AC", "Achilles_Demeter"),
  alt_matrix1 = SLiCell::gene_cl_mut_matrix,
  min_alt_no = 10) {

  RSA = apply(sl_pair_list, 1, SLiCell::run_wilcoxon_sl_pair,
    sens_matrix[[1]], alt_matrix1, min_alt_no = min_alt_no)

  Ataris = apply(sl_pair_list, 1, SLiCell::run_wilcoxon_sl_pair,
    sens_matrix[[2]], alt_matrix1, min_alt_no = min_alt_no)

  AC = apply(sl_pair_list, 1, SLiCell::run_wilcoxon_sl_pair,
    sens_matrix[[3]], alt_matrix1, min_alt_no = min_alt_no)

  Demeter = apply(sl_pair_list, 1, SLiCell::run_wilcoxon_sl_pair,
    sens_matrix[[4]], alt_matrix1, min_alt_no = min_alt_no)

  results = list(sl_pair_list, RSA, Ataris, AC, Demeter)
  names(results) = c("sl_pair_list", data_set_names)

  return(results)


}
