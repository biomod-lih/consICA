#' @title Save PDF report with analysis of each independent component
#' @description Save PDF report with description of each independent component 
#' (IC) consists of most affected genes, significant Go terms, survival model
#' for the component, ANOVA analysis for samples characteristics and stability
#' @param IC list compliant to `consICA()` result
#' @param Genes features list compilant to `getFeatures` output (list of 
#' dataframes `pos` for positive and `neg` for negative affecting features with
#' names of features false discovery rates columns).If NULL will generated 
#' automatically
#' @param GO list compilant to `getGO` output. If not NULL the significant GO 
#' terms will printed in report
#' @param Var matrix with samples metadata
#' @param surv dataframe with time and event values for each sample
#' @param genenames alternative gene names for printing in the report
#' @param file report filename, ends with ".pdf"
#' @param main title for each list discribes the component
#' @param show.components which compont will be shown
#' @author Petr V. Nazarov
#' @return TRUE when successfully generate report
#' @examples
#' data("samples_data")
#' cica <- consICA(samples_data, ncomp=40, ntry=10, show.every=0)
#' GOs <- NULL 
#' if(FALSE){
#' GOs <- getGO(cica, db = "BP")
#' }
#' saveReport(cica, GO=GOs, Var=samples_data$Var, surv = samples_data$Sur)
#' @export
# @import gplots
#' @importFrom pheatmap pheatmap
saveReport <- function(IC, Genes=NULL, GO=NULL, Var=NULL, surv=NULL, 
                      genenames=NULL,
                      file = sprintf("report_ICA_%d.pdf",ncol(IC$S)), 
                      main = "Component # %d (stability = %.3f)", 
                      show.components = seq.int(1,ncol(IC$S))){
  if(is.null(Genes)) Genes <- getFeatures(IC,alpha=0.01)
  if(!is.null(Var)) if (!is.data.frame(Var)) Var <- as.data.frame(Var)
  
  try(dev.off(),silent=TRUE)
  try(dev.off(),silent=TRUE)
  try(dev.off(),silent=TRUE)
  
  ncomp <- ncol(IC$S)
  pdf(file,width=8.3, height=11.7,onefile=TRUE)
  Genes <- sortFeatures(Genes)
  
  for (icomp in show.components){
    par(mfcol=c(1,1),mar=c(3,3,2,1))
    plot.new()
    title(sprintf(main,icomp,mean(IC$stab[,icomp])),cex.main=0.8)
    ## gene signature
    par(fig=c(0,0.2,0.85,1),new=TRUE,mar=c(2,2,2,1))
    plot(sort(IC$S[,icomp]),col = "#0000FF", type="l", ylab=("involvement"),
         xlab="genes",las=2,cex.axis=0.4, 
         main="Metagene\n(involvement of features)",cex.main=0.6)
    text(0,0,"negative",adj=c(0,-1),cex=0.5,col="#000088")
    text(nrow(IC$S),0,"positive",adj=c(1,1),cex=0.5,col="#880000")
    par(fig=c(0,0.24,0,0.85),new=TRUE,mar=c(2,1,0,1))
    plot.new()
    text(0,1,paste(nrow(Genes[[icomp]]$neg),"\nnegative"),font=2,adj=c(0,0),
         cex=0.5,col="#000088")
    if (nrow(Genes[[icomp]]$neg)>0){
      txt <- rownames(Genes[[icomp]]$neg)
      if (!is.null(genenames)) txt <- genenames[txt]
      for (i in seq.int(1,nrow(Genes[[icomp]]$neg))) text(0,0.99-i/80,txt[i],
                                                 col="#000088",adj=c(0,0),
                                                 cex=0.4)
    }
    rect(0.5,0,1,1,col="white",border=NA)
    text(0.5,1,paste(nrow(Genes[[icomp]]$pos),"\npositive"),font=2,adj=c(0,0),
         cex=0.5,col="#880000")
    if (nrow(Genes[[icomp]]$pos)>0){
      txt <- rownames(Genes[[icomp]]$pos)
      if (!is.null(genenames)) txt <- genenames[txt]
      for (i in seq.int(1,nrow(Genes[[icomp]]$pos))) text(0.5,0.99-i/80,txt[i],
                                              col="#880000",adj=c(0,0),cex=0.4)
    }
    
    if (!is.null(GO)){
      is_bp <- !is.null(GO$GOBP)
      is_cc <- !is.null(GO$GOCC)
      is_mf <- !is.null(GO$GOMF)
      if(is_bp){
        par(fig=c(0.2,0.6,0.7,1),new=TRUE,mar=c(2,2,2,0));plot.new()
        tab <- GO$GOBP[[icomp]]$neg[,c("Term","FDR")]
        tab <-tab[tab$FDR<0.1,]
        tab$FDR <- sprintf("%.2e",tab$FDR)
        text(0,1,sprintf("GO:BP neg : %d terms(FDR<0.1)",nrow(tab)),
             col="#000088",font=2,adj=c(0,0),cex=0.6)
        if (nrow(tab)>0) drawTable(tab[seq.int(1,min(20,nrow(tab))),],
                                   x0=0,y0=0.98,
                                   dx=c(0.8,0.2),dy=0.04,row.names=FALSE,
                                   cex=0.5,col="#000088")
        par(fig=c(0.6,1,0.7,1),new=TRUE,mar=c(2,2,2,0));plot.new()
        tab <- GO$GOBP[[icomp]]$pos[,c("Term","FDR")]
        tab <- tab[tab$FDR<0.1,]
        tab$FDR <- sprintf("%.2e",tab$FDR)
        text(0,1,sprintf("GO:BP pos : %d terms(FDR<0.1)",nrow(tab)),
             col="#880000",font=2,adj=c(0,0),cex=0.6)
        if (nrow(tab)>0) drawTable(tab[seq.int(1,min(20,nrow(tab))),],
                                   x0=0,y0=0.98,
                                   dx=c(0.8,0.2),dy=0.04,row.names=FALSE,
                                   cex=0.5,col="#880000")
      }
      if(is_cc){
        par(fig=c(0.2,0.6,0.506,0.806),new=TRUE,mar=c(2,2,2,0));plot.new()
        tab <- GO$GOCC[[icomp]]$neg[,c("Term","FDR")]
        tab <- tab[tab$FDR<0.1,]
        tab$FDR <- sprintf("%.2e",tab$FDR)
        text(0,1,sprintf("GO:CC neg : %d terms(FDR<0.1)",nrow(tab)),
             col="#000088",font=2,adj=c(0,0),cex=0.6)
        if (nrow(tab)>0) drawTable(tab[seq.int(1,min(10,nrow(tab))),],
                                   x0=0,y0=0.98,
                                   dx=c(0.8,0.2),dy=0.04,row.names=FALSE,
                                   cex=0.5,col="#000088")
        par(fig=c(0.6,1,0.506,0.806),new=TRUE,mar=c(2,2,2,0));plot.new()
        tab <- GO$GOCC[[icomp]]$pos[,c("Term","FDR")]
        tab <- tab[tab$FDR<0.1,]
        tab$FDR <- sprintf("%.2e",tab$FDR)
        text(0,1,sprintf("GO:CC pos : %d terms(FDR<0.1)",nrow(tab)),
             col="#880000",font=2,adj=c(0,0),cex=0.6)
        if (nrow(tab)>0) drawTable(tab[seq.int(1,min(10,nrow(tab))),],
                                   x0=0,y0=0.98,
                                   dx=c(0.8,0.2),dy=0.04,row.names=FALSE,
                                   cex=0.5,col="#880000")
      }      
      if(is_mf){
        par(fig=c(0.2,0.6,0.397,0.697),new=TRUE,mar=c(2,2,2,0));plot.new()
        tab <- GO$GOMF[[icomp]]$neg[,c("Term","FDR")]
        tab <- tab[tab$FDR<0.1,]
        tab$FDR <- sprintf("%.2e",tab$FDR)
        text(0,1,sprintf("GO:MF neg : %d terms(FDR<0.1)",nrow(tab)),
             col="#000088",font=2,adj=c(0,0),cex=0.6)
        if (nrow(tab)>0) drawTable(tab[seq.int(1,min(10,nrow(tab))),],
                                   x0=0,y0=0.98,
                                   dx=c(0.8,0.2),dy=0.04,row.names=FALSE,
                                   cex=0.5,col="#000088")
        par(fig=c(0.6,1,0.397,0.697),new=TRUE,mar=c(2,2,2,0));plot.new()
        tab <- GO$GOMF[[icomp]]$pos[,c("Term","FDR")]
        tab <- tab[tab$FDR<0.1,]
        tab$FDR <- sprintf("%.2e",tab$FDR)
        text(0,1,sprintf("GO:MF pos : %d terms(FDR<0.1)",nrow(tab)),
             col="#880000",font=2,adj=c(0,0),cex=0.6)
        if (nrow(tab)>0) drawTable(tab[seq.int(1,min(10,nrow(tab))),],
                                   x0=0,y0=0.98,
                                   dx=c(0.8,0.2),dy=0.04,row.names=FALSE,
                                   cex=0.5,col="#880000")
      }   
    }
    
    ix <- 1; iy <- 1
    if (!is.null(surv)){
      par(fig=c(0.2+(ix-1)*0.2,0.2+ix*0.2,0.55-iy*0.2,0.55-(iy-1)*0.2),
          new=TRUE,mar=c(4,2,2,1))
      scoreD <- c("low","high")[as.integer(IC$M[icomp,]>median(IC$M[icomp,]))+1]
      scoreD <- factor(scoreD,levels=c("low","high"))
      cox <- coxph(Surv(time = surv$time, event = surv$event) ~ IC$M[icomp,])
      pv <- summary(cox)$logtest["pvalue"]
      lhr <- log((summary(cox)$conf.int))[c(1,3,4)]
      if (pv<1e-2) {col<-c("blue","red")}else{col<-c("#666688","#886666")}
      plot(survfit(Surv(time = surv$time, event = surv$event) ~ scoreD,
                   type="kaplan-meier"),col=col,conf.int=FALSE,las=2,
           lwd=c(2,2,2),cex.axis=0.4)
      title(sprintf(
        "Cox regression:\nlogtest pv=%.1e\nLHR=%.2f (CI = %.2f, %.2f)",
        pv,lhr[1],lhr[2], lhr[3]),cex.main=0.5)
    }
    
    if (!is.null(Var)){
      #td check whetehr factor is adequate (number of levels and NAs)
      pv <- double(ncol(Var))+1
      names(pv) <- names(Var)
      for (ifact in seq.int(1,ncol(Var))){
        fact <- Var[[ifact]]
        ikeep <- !is.na(fact)
        fact <- as.factor(as.character(fact[ikeep]))
        if (nlevels(fact)==1) next
        x <- IC$M[icomp,ikeep]
        res <- aov(x~fact)
        pv[ifact] <- summary(res)[[1]][1,5]
      }
      pv <- pv[getTopIdx(-log10(pv),min(10,length(pv)))]
      tab <- data.frame(factor = names(pv), p.value = sprintf("%.2e",pv))
      par(fig=c(0.4,1,0.4,0.58),new=TRUE,mar=c(0,2,2,1))
      plot.new()
      drawTable(tab,dx=c(0.8,0.2),dy=0.08,cex=0.5,row.names=FALSE,bg="white")
      ix<-1; iy<-2
      ifact<-1
      for (ifact in seq.int(1,min(8,length(pv)))){
        factname <- names(pv)[ifact]
        fact <- Var[[factname]]
        ikeep <- !is.na(fact)
        fact <- as.factor(as.character(fact[ikeep]))
        
        fontsize <- 0.4
        if (nlevels(fact)<=5) fontsize <- 0.6
        if (nlevels(fact)>15) fontsize <- 0.3
        if (nlevels(fact)>30) fontsize <- 0.2
        
        x <- IC$M[icomp,ikeep]
        xf <- list()
        for (i in seq.int(1,nlevels(fact))) xf[[i]] <- x[fact==levels(fact)[i]]
        names(xf) <- levels(fact)
        if (ix>4) {ix<-1;iy<-iy+1}
        par(fig=abs(c(0.2+(ix-1)*0.2,0.2+ix*0.2,0.6-iy*0.2,0.6-(iy-1)*0.2)),
            new=TRUE,mar=c(4,2,2,1))
        col<-"grey"
        if (pv[ifact]<1e-5)col<-"#66AAFF"
        if (pv[ifact]<1e-10)col<-"#88FF66"
        if (pv[ifact]<1e-20)col<-"#FF8866"
        violinplot(xf, col=col,cex.axis = fontsize,colMed = "black")
        title(sprintf("%s\npv=%.1e",factname,pv[ifact]),cex.main=0.5)
        ix <- ix + 1
      }
    }
  }
  pheatmap(cor(t(IC$M))^2, main="R2 of M-matrix")
  if (!is.null(IC$stab)){
    par(mfcol=c(2,1),mar=c(4,3,3,1))
    boxplot(IC$stab,las=2,col="skyblue",main=sprintf("Stability (%d runs)",
                                                     nrow(IC$stab)),
            subj="R2 between most correlates S in multiple runs")
  }
  if (!is.null(IC$mr2)) 
    if(length(IC$mr2) > 1){
      plot(density(IC$mr2),lwd=2,col="blue",
           main="Distribution of mean R2\namong rows of M-matrix",
           xlab="Mean R2 for each single run")
    } else{
      plot(x = IC$mr2, y = 1, lwd=2,col="blue",
           main="Distribution of mean R2\namong rows of M-matrix",
           ylab="Density",
           xlab="Mean R2 for each single run")
    }
  dev.off()
  message("Report saved to ", file)
  return(TRUE)
}
