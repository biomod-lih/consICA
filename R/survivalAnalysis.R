#' @title Survival analysis based on significant IC
#' @description Cox regression (based on R package `survival`) on the weights 
#' of independent components with significant contribution in individual
#' risk model. For more see Nazarov et al. 2019 
#' In addition the function plot Kaplan-Meier diagram.
#' @param IC  list compliant to `consICA()` result
#' @param surv dataframe with time and event values for each sample. Use this
#' parameter or `time` and `event`
#' @param time survival time value for each sample 
#' @param event survival event factor for each sample (TRUE or FALSE)
#' @param fdr false discovery rate threshold for significant components 
#' involved in final model. Default value is 0.05
#' @return a list with
#'     \item{cox.model}{an object of class `coxph` representing the fit. 
#' See `coxph.object` for details}
#'     \item{hazard.score}{hazard score for significant components (fdr < `fdr` 
#'     in individual cox model)}
#' @examples
#' data("samples_data")
#' # Get deconvolution of X matrix
#' cica <-  consICA(samples_data, ncomp=10, ntry=1, show.every=0)
#' surv <- survivalAnalysis(cica, 
#'   surv = SummarizedExperiment::colData(samples_data)[,c("time", "event")]) 
#' @export
#' @importFrom survival coxph survfit Surv
survivalAnalysis <- function(IC,surv=NULL,time=NULL,event=NULL,fdr=0.05){
   
   if(!is.null(time) & !is.null(event)){
     surv <- list("time" = time, "event" = event)
   }
   if(is.null(surv)){
     message("No time-evert data!\n")
     return(NULL)
   }
   
   ZM <- IC$M
   ZM[,] <- t(scale(t(IC$M))) 
   ResR <- data.frame(matrix(nrow=nrow(IC$M), ncol=11))
   rownames(ResR) <- rownames(IC$M)
   colnames(ResR) <- c("Stab","Pval","FDR","LHR","LHRmin","LHRmax",
                       "nPos","nNeg","nPosGOBP","nNegGOBP","Pval_validation")
   if(is.null(IC$stab)){
     IC$stab <- matrix(1, nrow = 1, ncol = nrow((IC$M)))
   }
   ResR$Stab <- apply(IC$stab,2,mean)
   
   for (icomp in seq.int(1,nrow(IC$M))){
     model <- coxph(Surv(time=surv$time, event=surv$event)~ZM[icomp,])
     ResR$Pval[icomp] <- summary(model)$logtest["pvalue"]
     ResR[icomp,c("LHR","LHRmin","LHRmax")] <- 
       log((summary(model)$conf.int))[c(1,3,4)]
     model.val <- coxph(Surv(time=surv$time,event=surv$event)~ZM[icomp,])
     ResR$Pval_validation[icomp] <- summary(model.val)$logtest["pvalue"]
   }
   ResR$FDR <- p.adjust(ResR$Pval,method="fdr")
   
   d <- list() 
   # auto select 
   idx <- ResR$FDR < fdr
   d$surv <- ResR$LHR[idx]
   names(d$surv) <- rownames(ResR)[idx]
   
   R2 <- apply(IC$stab,2,mean)
   names(R2) <- colnames(IC$S)
   
   ## calculate hazard score
   ## as sum for each sample 
   d1 <- d$surv
   score1 <- lapply(names(d1), 
                    function(comp){
                      0 + d1[comp] * R2[comp] * ZM[comp,] 
                    })
   score_sumed <- matrix(unlist(score1), ncol = length(score1[[1]]))
   score_sumed <- colSums(score_sumed)
   names(score_sumed) <- names((score1[[1]])) 
   score <- list("surv" = score_sumed)
   
   i <- !is.na(surv$time)
   cox.train <- coxph(Surv(time = surv$time[i], event =surv$event[i]) ~ 
                       score$surv[i]) 
   fscore <- num2fact(score$surv,nlev=2)
   pv <- summary(cox.train)$logtest["pvalue"]
   lhr <- log((summary(cox.train)$conf.int))[c(1,3,4)]
   plot(survfit(Surv(time = surv$time[i], event = surv$event[i]) ~ 
                  fscore[i],type="kaplan-meier"),col=c("blue","red"),
        conf.int=FALSE,las=2,lwd=3,cex.axis=1)
   legend(x="bottomleft",legend=c("HS<median","HS>median"),cex=0.8,
          col=c("blue","red"),lwd=3)
   title(sprintf("logtest pv=%.1e\nLHR=%.2f (CI = %.2f, %.2f)",pv,
                 lhr[1],lhr[2], lhr[3]),cex.main=1)
   
   # avoid note: print(cox.train)
   return(list("cox.model" = cox.train,
               "hazard.score" = d$surv))
}
 
