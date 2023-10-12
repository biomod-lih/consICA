#' @title ANOVA test for independent component across factors
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
#' # Run ANOVA for 4th independent component
#' anova <- anovaIC(cica, Var=Var, icomp = 4)
#' @export
#' @import ggplot2
anovaIC <- function(cica, Var=NULL, icomp = 1, plot = TRUE, mode = "violin",
                    color_by_pv = TRUE){
  
  if(!is.consICA(cica)) {
    message("First parameter should be compliant to `consICA()` result\n")
    return (NULL)
  }
  if(length(which(icomp == seq.int(1,ncol(cica$S)))) < 1){
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
      gg_df <- do.call(rbind, lapply(seq.int(1,length(top_fact)), function(x){
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

#' @title Similarity of two gene ontologies lists
#' @description Calculate similarity matrix of gene ontilogies (GOs) 
#' of independent components. The measure could be cosine similarity or
#' Jaccard index (see details)
#' @param GO1 list of GOs for each independent component got from `getGO()`
#' @param GO2 list of GOs for each independent component got from `getGO()`
#' @param method can be 'cosine' for non-parametric cosine similarity  or 
#' 'jaccard' for Jaccadr index. See details
#' @param fdr FDR threshold for GOs that would be used in measures. 
#' Default value is 0.01
#' @return a similarity matrix of cosine or Jaccard values, 
#' rows corresponds to independent components in `GO1`, 
#' columns to independent components in `GO2`.
#' @details 
#' Jaccard index is a measure of the similarity between two sets of data. 
#' It calculated as intersection divided by union \deqn{
#'   J(A, B) = \frac{|A \cap B|}{|A \cup B|}.
#' }
#'  Results are from 0 to 1.
#' 
#' Cosine similarity here is calculated in a non-parametric way: 
#' for two vectors of gene ontologies, the space is created as a union of GOs 
#' in both vectors.
#' Then, two rank vectors in this space created, 
#' most enriched GOs get the biggest rank and GOs from space not included in 
#' the GO vector get 0.
#' Cosine similarity is calculated between two scaled rank vectors. 
#' Such approach allows to
#' take the order of enriched GO into account.
#' Results are from -1 to 1. Zero means no similarity.
#' @author Maryna Chepeleva
#' @examples
#' \dontrun{
#' data("samples_data")
#' # Calculate ICA (run with ntry=1 for quick test, use more in real analysis)
#' cica1 <- consICA(samples_data, ncomp=5, ntry=1, show.every=0)
#' # Search enriched gene ontologies
#' cica1 <- getGO(cica1, db = "BP", ncores = 1)
#' # Calculate ICA and GOs for another dataset
#' cica2 <- consICA(samples_data[,1:100], ncomp=4, ntry=1, show.every=0) 
#' cica2 <- getGO(cica2, db = "BP", ncores = 1)

#' # Compare two lists of enriched GOs 
#' # Jaccard index
#' jc <- overlapGO(GO1 = cica1$GO$GOBP, GO2 = cica2$GO$GOBP, 
#' method = "jaccard", fdr = 0.01)
#' # Cosine similarity
#' cos_sim <- overlapGO(GO1 = cica1$GO$GOBP, GO2 = cica2$GO$GOBP, 
#' method = "cosine", fdr = 0.01)
#' }
#' @export
overlapGO <- function(GO1, GO2, method = c("cosine", "jaccard"), fdr = 0.01){
  
  method <- method[1]
  if(length(which(method == c("cosine", "jaccard"))) == 0){
    message("`method` should be 'cosine' or 'jaccard'. Return NULL.\n")
    return(NULL)
  }
  
  if(method == "jaccard"){
    jac_mat <- matrix(0, nrow = length(GO1), ncol = length(GO2))
    for(i in seq.int(1, length(GO1))){
      for(j in seq.int(1, length(GO2))){
        go_vector1 <- GO1[[i]]$pos[GO1[[i]]$pos$FDR <= fdr,]$GO.ID
        go_vector2 <- GO2[[j]]$pos[GO2[[j]]$pos$FDR <= fdr,]$GO.ID
        if(length(go_vector1) == 0 | length(go_vector2) == 0) next
        jac_mat[i,j] <- length(intersect(go_vector1, go_vector2)) / length(union(go_vector1, go_vector2))
      }
    }
    rownames(jac_mat) <- names(GO1)
    colnames(jac_mat) <- names(GO2)
    return(jac_mat)
  }
  
  if(method == "cosine"){
    cos_mat <- matrix(0, nrow = length(GO1), ncol = length(GO2))
    for(i in seq.int(1, length(GO1))){
      for(j in seq.int(1, length(GO2))){
        go_rank1 <- GO1[[i]]$pos[GO1[[i]]$pos$FDR <= fdr,]
        go_rank2 <- GO2[[j]]$pos[GO2[[j]]$pos$FDR <= fdr,]
        if(nrow(go_rank1) == 0 | nrow(go_rank2) == 0) next
        go_rank1$RANK1 <- seq.int(nrow(go_rank1),1)
        go_rank2$RANK2 <- seq.int(nrow(go_rank2),1)
        rank_df <- merge(go_rank1[,c("GO.ID","RANK1")], 
                         go_rank2[,c("GO.ID","RANK2")], by = "GO.ID",all = TRUE)     
        rank_df[is.na(rank_df)] <- 0
        rank_df$RANK1  <- (rank_df$RANK1 - min(rank_df$RANK1)) / 
          (max(rank_df$RANK1) - min(rank_df$RANK1)) * (1-0) + 0
        rank_df$RANK2  <- (rank_df$RANK2 - min(rank_df$RANK2)) / 
          (max(rank_df$RANK2) - min(rank_df$RANK2)) * (1-0) + 0
        cos_mat[i,j] <- crossprod(rank_df$RANK1, rank_df$RANK2)/
          sqrt(crossprod(rank_df$RANK1) * crossprod(rank_df$RANK2))
      }
    }
    rownames(cos_mat) <- names(GO1)
    colnames(cos_mat) <- names(GO2)
    return(cos_mat)
  }
 return(NULL)
}




