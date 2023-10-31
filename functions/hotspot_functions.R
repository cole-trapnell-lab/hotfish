#' calculate spatial weights
#' functions take from graph_test
#' @k is how many neighbors in knn graph
#' @reduction_method default is UMAP
#'
calculateSpatialWeights <- function(cds, k=15, reduction_method = "UMAP"){

  lw <- monocle3:::calculateLW(cds = cds, k = k,
                                 verbose = FALSE,
                                 neighbor_graph = "knn",
                                 reduction_method = reduction_method,
                                 list(method = "nn2"))

  wc <- spdep::spweights.constants(lw, zero.policy = TRUE, adjust.n = TRUE)

  return(list(lw=lw, wc=wc))
}


#' @var is a specific label value
#' @column_name of labels of interest
#' @lw is output of calculateLW
calculateLocalG <- function(cds, lw, wc, column_name, var) {

  # absence or presence of specified variable
  # df$LG <- ifelse(df[,column_name] == var, 1,0)
  # z = df[,"LG"]
  # names(z) <- df$Cell

  cds$LG <- ifelse(pData(cds)[,column_name] == var, 1,0)
  z = pData(cds)[,"LG"]
  names(z) <- pData(cds)$Cell

  # calculate local g
  localG = spdep::localG(x=z, listw = lw, wc$n, wc$S0, zero.policy = TRUE, spChk = F)
  localG.df = localG[1:length(localG)] %>% as.data.frame()
  colnames(localG.df)[1] <- as.character(var)
  rownames(localG.df) = pData(cds)$Cell
  localG.df
}

#' select local g values of interest + merge with full df
#' only values 1 mean anything in the local G calculation
mergeLocalG <- function(localG_result, df, column) {

  localG_result$Cell = row.names(localG_result)
  wt = localG_result %>% merge(df, by=c("Cell"))
  wt = wt[wt[,column] == colnames(localG_result)[1],]
  colnames(wt)[2] <- "localG"
  wt
}

#' to turn the zscore into a pvalue
#' method : c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
getPval <- function(df, method, tail) {

  # I think this is right, but should check it...
  if (tail=="onesided"){
    df = df %>% mutate(pval = pnorm(-abs(localG)))
  } else if (tail=="twosided"){
    df = df %>% mutate(pval = 2*pnorm(-abs(localG)))
  }
  adjusted_pval = p.adjust( p = df$pval, method = method)
  df %>% mutate("adjusted_pval" = adjusted_pval)
}


#' compare col : column in cds that is your categorical label, ex: "gene_target"
#' compare: which category you want to determine hotspot, ex: "tbxta"
#' subset col: if you want to do hotspot by broad grouping, ex: "cell_type_broad"
#' subset: the specific subset, ex: "fast muscle"
calc_hotspot <- function(cds,
                         compare_col, compare,
                         subset_col = NULL, subset = NULL,
                         reduction_method="UMAP",
                         method = "bonferroni",
                         tail = "twosided") {

  if (!(is.null(subset_col) & is.null(subset))) {
    cds = cds[, tidyr::replace_na(colData(cds)[[subset_col]] == subset, F)]
  }
  pData(cds)$Cell <- pData(cds)$cell
  df = pData(cds) %>% as.data.frame()
  spatial_weights = calculateSpatialWeights(cds, reduction_method = reduction_method)

  lg <- calculateLocalG(cds = cds, lw = spatial_weights$lw,
                        wc = spatial_weights$wc, column_name = compare_col, var = compare)
  lg.df <- mergeLocalG(lg, df, compare_col) %>% getPval(method, tail)
  lg.df
}



