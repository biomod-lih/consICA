#' @title ANOVA test independent component ~ factor
#' @description ANOVA (ANalysis Of VAriance) test produced for specific 
#' independent component across each (clinical) factor as `aov(IC ~ factor)`.
#' Plot distributions of samples' weight for top 9 significant factors.
#' @param cica list compliant to `consICA()` result
#' @param Var matrix with samples' metadata. Samples in rows and factors in 
#' columns
#' @param icomp number of component to analyse
#' @param plot if plot weights distributions for top factors
#' @param mode type of plot. Can be 'violin' or 'box'
#' @param color_by_pv if TRUE plots will be colored by p-value ranges
#' @return a data.frame with 
#'     \item{factor}{name of factor}
#'     \item{p.value}{p-value for ANOVA test for factor} 
#'     \item{p.value_disp}{string for p-value printing}
#' @examples
#' data("samples_data")
#' Var <- data.frame(SummarizedExperiment::colData(samples_data))
#' cica <-  consICA(samples_data, ncomp=10, ntry=1, ncores=1, show.every=0)
#' anova <- anovaIC(cica, Var=Var, icomp = 4)
#' @export
#' @import ggplot2
anovaIC <- function(cica, Var=NULL, icomp = 1, plot = TRUE, mode = "violin",
                    color_by_pv = TRUE){
  
  if(!is.consICA(cica)) {
    message("First parameter should be compliant to `consICA()` result\n")
    return (NULL)
  }
  if(length(which(icomp == seq.int(1:ncol(cica$S)))) < 1){
    message("Component not found\n")
    return(NULL)
  }
  if(length(which(mode == c("violin", "box"))) < 1){
    message("Plot `mode` should be 'violin' or 'box'. Changed to 'violin'...\n")
    mode  <- "violin"
  }

  if (!is.null(Var)){
    pv <- double(ncol(Var))+1
    names(pv) <- names(Var)
    for (ifact in seq.int(1,ncol(Var))){
      fact <- Var[[ifact]]
      ikeep <- !is.na(fact)
      fact <- as.factor(as.character(fact[ikeep]))
      if (nlevels(fact)==1) next
      x <- cica$M[icomp,ikeep]
      res <- aov(x~fact)
      pv[ifact] <- summary(res)[[1]][1,5]
    }
    pv <- pv[getTopIdx(-log10(pv),length(pv))]
    tab <- data.frame(factor = names(pv),  
                      p.value = pv,
                      p.value_disp = sprintf("%.2e",pv))
    rownames(tab) <- NULL
    top_fact <- pv[seq.int(1,min(9, length(pv)))]
    top_fact <- names(top_fact[top_fact<0.05])
    
    if(plot & length(top_fact)>0){
      gg_d <- as.data.frame(t(cica$M))[,icomp]
      gg_df <- do.call(rbind, lapply(1:length(top_fact), function(x){
        cbind(data.frame(gg_d),top_fact[x],as.character(Var[,top_fact[x]]),
              tab[x,]$p.value, tab[x,]$p.value_disp)}))
      colnames(gg_df) <- c("weight", "factor_name", "factor",
                           "p.value","p.value_disp")
      if(color_by_pv){
        gg_df$col <- ifelse(gg_df$p.value<1e-10, "1", "2")
        gg_df$col <- ifelse(gg_df$p.value<1e-20, "3",  gg_df$col)
      } else {
        gg_df$col <-1
      }
      gg_df$factor_name_facet <- paste0(gg_df$factor_name, ", p-value = ", 
                                        gg_df$p.value_disp)
      gg_df$factor_name_facet <- factor(gg_df$factor_name_facet, 
                                        levels=unique(gg_df$factor_name_facet))
      if(mode == "violin"){
        vp <- ggplot(gg_df, aes(x = factor, y = weight)) +
          geom_violin( alpha = 0.7, aes(fill = col)) +
          geom_boxplot(width=0.08, fill="white") +
          facet_wrap(~factor_name_facet, scales = "free", ncol=3) +
          labs(title = paste0("IC", icomp, " weights for factors")) +
          theme_bw() + theme(legend.position = "none")        
        print(vp)
      }
      if(mode == "box"){
        bp <- ggplot(gg_df, aes(x = factor, y = weight)) +
          geom_boxplot( aes(fill = col)) +
          facet_wrap(~factor_name_facet, scales = "free", ncol=3) +
          labs(title = paste0("IC", icomp, " weights for factors")) +
          theme_bw() + theme(legend.position = "none")        
        print(bp)
      }
    }
    return(tab)
  }
  message("`Var` is empty. Return NULL")
  return(NULL)
}
