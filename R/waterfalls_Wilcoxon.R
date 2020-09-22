#' checks primary site of a cell line
#' where geneA which is altered in the minimum no of cell lines for Wilcoxon test
#' and geneB is present in the chosen knock-down experiment
#' @param cl_id cell line id
#' @param cl_ann cell line annotations
#' @return primary site
#' @export

primary_site <- function(cl_id, cl_ann) {
  as.character(cl_ann[cl_ann[,1]==cl_id,5])
}

#' arranges single legend
#' @return grid.grab object
#' @export

grid_arrange_shared_legend <- function(..., title = "title", ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] +
                    theme(legend.position = position))$grobs 
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- unit(0.2, "npc") # percentage of screen for the legend
  gl <- lapply(plots, function(x) x +
                 theme(legend.position = "none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl), 
                                            legend,ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend, ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  
  grid.newpage()
  
  grid.draw(combined)
  grid.text(title, x=0.9, y=0.9, gp = gpar(fontsize =7))
  #vwhb
  
  grid.grab()
  

}

#' creates one waterfall plot
#' where geneA altered cell lines (for Wilcoxon test subpopulation) are colored
#' and geneB depenency is plotted
#' @param sl_pair geneA geneB pair for verifing
#' @param sens list of geneB dependency matrices in geneA altered cell lines, default - all experimental data preloaded in the package
#' @param alt geneA alteration matrix for cell lines, default - mutation_matrix preloaded in the package
#' @param experiment name of the dependency data set
#' @return waterfall plot
#' @export

wf <- function (sl_pair, sens, alt, experiment, pval=F, pallete=SLiCell::sites_colors, cl_annotations=SLiCell::cl_annotations) {
  
  name_gene1=as.character(sl_pair[1])
  name_gene2=as.character(sl_pair[2])
  
  # list of cell lines suitable for testing - present in both matrices
  common_cell_lines <- intersect(rownames(sens), rownames(alt))
  
  # geneA alteration list
  alt_gene1 <- alt[common_cell_lines, name_gene1]
  names(alt_gene1) <- rownames(alt[common_cell_lines,])
  # geneB sensitivity list
  sens_gene2 = sens[common_cell_lines,name_gene2]
  names(sens_gene2) <- rownames(sens[common_cell_lines,])
  
  
  filter_G1=as.data.frame(alt_gene1==1)
  
  sens_1=as.data.frame(sens_gene2)
  sens_2=cbind(sens_1,filter_G1)
  colnames(sens_2)=c("sens","is_alt")
  
  filter2=is.na(sens_2$sens)
  sens_2=sens_2[!filter2,]

  sens_wf=sens_2[order(sens_2$sens),]

  filter3 = sens_wf$is_alt == T
  primary_site=rep("cl_not_altered",nrow(sens_wf))
  primary_site[filter3] = sapply(rownames(sens_wf[filter3,]), SLiCell::primary_site, cl_annotations)
  
  sens_wf <- as.data.frame(cbind(cell_lines=seq(1,nrow(sens_wf),1),cl=rownames(sens_wf),sens_wf,primary_site))

  ggplot2::ggplot(data=sens_wf, ggplot2::aes(x=cell_lines, y=sens, fill = primary_site)) +
    ggplot2::geom_bar(stat="identity")+ ggplot2::scale_fill_manual(values = pallete[as.character(sens_wf$primary_site)])+
    ggplot2::ggtitle(paste0(experiment, " p_val = ", pval)) + 
    ggplot2::labs(x=paste0("cell lines, coloured when ", name_gene1, " altered"), y = paste0("dependency score ", name_gene2)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 8), axis.text.x = ggplot2::element_text(size=6),
          axis.title=ggplot2::element_text(size=7), legend.title = ggplot2::element_text(size =7), legend.key.size = ggplot2::unit(0.25, "cm"),
          legend.text = ggplot2::element_text(size = 6))
    
  #dev.off()
}

