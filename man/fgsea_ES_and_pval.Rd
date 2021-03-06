% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment_analysis.R
\name{fgsea_ES_and_pval}
\alias{fgsea_ES_and_pval}
\title{enrichment analysis for one potential SL gene pairs
enrichment set - geneA altered cell lines
geneB sensitivity is a ranking criterion}
\usage{
fgsea_ES_and_pval(sl_pair, sens_matrix1, alt_matrix1, plot = FALSE)
}
\arguments{
\item{sl_pair}{geneA and geneB vector}

\item{sens_matrix}{geneB dependency matrix in CCLE cell lines for a given knock-down experiment}

\item{alt_matrix}{geneA alteration matrix for CCLE cell lines}

\item{plot=FALSE}{option for printing enrichment plot}
}
\value{
list of vectors of enrichment score and its p-value
}
\description{
enrichment analysis for one potential SL gene pairs
enrichment set - geneA altered cell lines
geneB sensitivity is a ranking criterion
}
