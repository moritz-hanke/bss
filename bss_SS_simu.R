
library( parallel)
library( snow)
library( Rmpi)


library(glmnet)
library(glmnetUtils)
library(caret)
library(bestsubset)
library(tibble)
library(tidyverse)
library(mvtnorm)
library(bestridge)
library(stabs)
library(lars)

.Last <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0){
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
} 

# function to generate a block correlation structure
block_builder <- 
  function(p, rho, size){
    n.blocks <- p/size 
    out <- matrix(rep(0, p^2), ncol = p) 
    for(i in c(0:(n.blocks-1))){
      out[c(i*size+1):c(i*size+size), c(i*size+1):c(i*size+size)] <- rho
    }
    diag(out) <- 1
    return(out)
  }

N <- 1000
P <- 100
s <- 10

use_gurobi <- FALSE

#corr_type <- "block"
corr_type <- "toeplitz"
#corr_type <- "independent"

#SNR <- c(0.05, 0.09, 0.14, 0.25, 0.42, 0.71, 1.22, 2.07, 3.52, 6)
SNR <- c(0.25, 0.71, 1.22, 2.07, 6)
if(corr_type == "independent"){
  RHO <- 0
}else{
  RHO <- c(0.35, 0.7)
}

beta_position <- "spread"



Sim_n <- 10

max.k <- 15
Alpha <- seq(0.1,1,0.1)
cutoff <- 0.9


glmnet.enet <- function(x, y, q, alpha, type = c("conservative", "anticonservative"), ...) {
  if (!requireNamespace("glmnet", quietly=TRUE))
    stop("Package ", sQuote("glmnet"), " needed but not available")
  
  if (is.data.frame(x)) {
    message("Note: ", sQuote("x"),
            " is coerced to a model matrix without intercept")
    x <- model.matrix(~ . - 1, x)
  }
  
  if ("lambda" %in% names(list(...)))
    stop("It is not permitted to specify the penalty parameter ", sQuote("lambda"),
         " for lasso when used with stability selection.")
  
  ## fit model
  type <- match.arg(type)
  if (type == "conservative")
    fit <- suppressWarnings(glmnet::glmnet(x, y, pmax = q, alpha = alpha,...))
  if (type == "anticonservative")
    fit <- glmnet::glmnet(x, y, alpha = alpha, dfmax = q - 1, ...)
  
  ## which coefficients are non-zero?
  selected <- predict(fit, type = "nonzero")
  selected <- selected[[length(selected)]]
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  ## compute selection paths
  cf <- fit$beta
  sequence <- as.matrix(cf != 0)
  ## return both
  return(list(selected = ret, path = sequence))
}

bss <- function(x, y, q, type = c("conservative", "anticonservative"), ...) {
  if (!requireNamespace("bestridge", quietly=TRUE))
    stop("Package ", sQuote("bestridge"), " needed but not available")
  
  if (is.data.frame(x)) {
    message("Note: ", sQuote("x"),
            " is coerced to a model matrix without intercept")
    x <- model.matrix(~ . - 1, x)
  }
  
  if ("lambda" %in% names(list(...)))
    stop("It is not permitted to specify the penalty parameter ", sQuote("lambda"),
         " for lasso when used with stability selection.")
  
  ## fit model
  fit <- bsrr(x=X, y=Y, 
                     family="gaussian", 
                     method="sequential", 
                     tune="cv",
                     s.list=q,
                     lambda.list = 0)
  
  ## which coefficients are non-zero?
  selected <- which(fit$beta != 0)
  ## check if variables are removed again from the active set
  ## and remove these from selected
  if (any(selected < 0)) {
    idx <- which(selected < 0)
    idx <- c(idx, which(selected %in% abs(selected[idx])))
    selected <- selected[-idx]
  }
  
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  ## compute selection paths
  cf <- fit$beta
  sequence <- t(cf != 0)
  ## return both
  return(list(selected = ret, path = sequence))
}




Loop_Sim_n <- lapply(1:Sim_n, function(sim_n){
  print(sim_n)
  #Loop_Sim_n <- parLapply(cl, 1:Sim_n, function(sim_n){    
  Loop_Rho <- lapply(RHO, function(Rho){
    
    Loop_Snr <- lapply(SNR, function(snr){
      
      
      if(corr_type == "toeplitz"){
        Sigma <- Rho^toeplitz(0:(P-1))
      }else if(corr_type == "independent"){
        Sigma <- diag(P)
      }else{
        Sigma <- block_builder(P, Rho, s)
      }
      
      
      beta <- rep(0, P)
      if(beta_position == "spread"){
        beta[seq(1,P, P/s)] <- 1
      }else{
        beta[1:s] <- 1
      }
      
      
      non_zero_indices <- which(beta != 0)
      
      sigma <- sqrt( (t(beta) %*% Sigma %*% beta) / snr )
      
      
      
      n <- N #nrow(TCGA_subdata)
      seed <- round(sim_n)
      try(set.seed(seed))
      e <- rnorm(n) %*% sigma 
      
      seed <- round(sim_n)
      try(set.seed(seed+snr*1000000))
      X <- mvtnorm::rmvnorm(n = n, mean = rep(0, ncol(Sigma)), sigma = Sigma)
      X <- scale(X)
      Y <- X %*% beta + e 
      Y <- scale(Y, scale = FALSE)
      
      
      ### new BSS
      debug(
      BSS <- stabsel(x = X, y = Y,
                     fitfun = bss, cutoff = 0.9,
                     PFER = 1)
      )
      # BSS <- bsrr(x=X, y=Y, 
      #             family="gaussian", 
      #             method="sequential", 
      #             tune="cv",
      #             s.list=1:max.k,
      #             lambda.list = 0, nfolds=cv.number)
      # 
      # estimated_non_zeros_BSS <- which(BSS$beta != 0)
      # 
      # TP <- sum(estimated_non_zeros_BSS %in% non_zero_indices)
      # FP <- length(estimated_non_zeros_BSS)-TP
      # FN <- s-TP
      # F1 <- TP/(TP + 0.5*(FP+FN))
      # Precision <- TP/(TP + FP)
      # Accuracy <- TP/s
      # 
      # if(length(estimated_non_zeros_BSS) <= 1){
      #   est_det_IFIM <- 1
      #   est_det_FIM <- 1
      # }else{
      #   est_det_IFIM <- det(cov(X[,estimated_non_zeros_BSS]))
      #   est_det_FIM <- det(solve(cov(X[,estimated_non_zeros_BSS])))
      # }
      # 
      # 
      # BSS_results <- tibble(
      #   method="BSS",
      #   alpha=NA,
      #   lambda = NA,
      #   RSS = sum((Y-X %*% BSS$beta)^2),
      #   RSS_corrected = NA,
      #   est_det_IFIM = est_det_IFIM,
      #   est_det_FIM = est_det_FIM,
      #   k=length(estimated_non_zeros_BSS),
      #   TP,
      #   FP,
      #   FN,
      #   F1,
      #   Precision,
      #   Accuracy,
      #   n,
      #   p=P,
      #   s,
      #   snr,
      #   rho=Rho,
      #   sim.n = sim_n,
      #   beta_position = beta_position,
      #   corr_type = corr_type,
      #   cv.number = cv.number
      # )
      # rm(estimated_non_zeros_BSS)
      
      
      ### begin FSS
      
      stabs_FSS <- stabsel(x = X, y = Y,
                           fitfun = lars.stepwise, cutoff = 0.9,
                           PFER = 1)
      
      estimated_non_zeros_FSS <- 
        stabs_FSS$selected
      
      TP <- sum(estimated_non_zeros_FSS %in% non_zero_indices)
      FP <- length(estimated_non_zeros_FSS)-TP
      FN <- s-TP
      F1 <- TP/(TP + 0.5*(FP+FN))
      Precision <- TP/(TP + FP)
      Accuracy <- TP/s
      
      if(length(estimated_non_zeros_FSS) <= 1){
        est_det_IFIM <- 1
        est_det_FIM <- 1
      }else{
        est_det_IFIM <- det(cov(X[,estimated_non_zeros_FSS]))
        est_det_FIM <- det(solve(cov(X[,estimated_non_zeros_FSS])))
      }
      
      # FSS_beta <- 
      #   solve(t(X[,estimated_non_zeros_FSS]) %*% X[,estimated_non_zeros_FSS]) %*% 
      #   t(X[,estimated_non_zeros_FSS]) %*% Y
      
      FSS_results <- tibble(
        method="FSS",
        alpha=NA,
        lambda = NA,
        RSS = NA,
        RSS_corrected = NA,
        est_det_IFIM = est_det_IFIM,
        est_det_FIM = est_det_FIM,
        k=length(estimated_non_zeros_FSS),
        TP,
        FP,
        FN,
        F1,
        Precision,
        Accuracy,
        n,
        p=P,
        s,
        snr,
        rho=Rho,
        sim.n = sim_n,
        beta_position = beta_position,
        corr_type = corr_type,
        cutoff = cutoff
      )
      rm(estimated_non_zeros_FSS)
      
      
      
      ### begin Lasso/Enet
      
      Enet_results <- lapply(Alpha, function(i){
        
        stabs_Enet <- 
          
          stabsel(x = X, y = Y,
                  fitfun = glmnet.enet, 
                  args.fitfun=c(alpha=i), 
                  cutoff = cutoff, 
                  PFER = 1)
        estimated_non_zeros_Enet <- 
          as.numeric(stabs_Enet$selected)
        
        
        TP <- sum(estimated_non_zeros_Enet %in% non_zero_indices)
        FP <- length(estimated_non_zeros_Enet)-TP
        FN <- s-TP
        F1 <- TP/(TP + 0.5*(FP+FN))
        Precision <- TP/(TP + FP)
        Accuracy <- TP/s
        
        if(length(estimated_non_zeros_Enet) <= 1){
          est_det_IFIM <- 1
          est_det_FIM <- 1
        }else if(matrixcalc::is.singular.matrix(cov(X[,estimated_non_zeros_Enet]))){
          est_det_IFIM <- NA
          est_det_FIM <- NA
        }else{
          est_det_IFIM <- det(cov(X[,estimated_non_zeros_Enet]))
          est_det_FIM <- det(solve(cov(X[,estimated_non_zeros_Enet])))
        }
        
        tibble(
          method=paste("Enet", i, sep=" "),
          alpha=i,
          lambda = NA,
          RSS = NA,
          RSS_corrected = NA,
          est_det_IFIM = est_det_IFIM,
          est_det_FIM = est_det_FIM,
          k=TP+FP,
          TP,
          FP,
          FN,
          F1,
          Precision,
          Accuracy,
          n,
          p=P,
          s,
          snr,
          rho=Rho,
          sim.n = sim_n,
          beta_position = beta_position,
          corr_type = corr_type,
          cutoff = cutoff
        )
        
      })
      
      Enet_results <- do.call(rbind, Enet_results)
      
      
      out <- rbind(#BSS_results,
        FSS_results,
        Enet_results)
      
      
      out
      
    })
    Loop_Snr <- do.call(rbind, Loop_Snr)
    Loop_Snr
    
    
  })
  
  Loop_Rho <- do.call(rbind, Loop_Rho)
  Loop_Rho
  
})
Loop_Sim_n <- do.call(rbind, Loop_Sim_n)

Loop_Sim_n$method <- 
  factor(Loop_Sim_n$method, 
         levels = c(paste("Enet", seq(0.1, 1, 0.1)), "FSS"), 
         labels = c(paste("Enet", seq(0.1, 1, 0.1)), "FSS"))


ggplot(Loop_Sim_n,
       aes(x = as.factor(snr), y = F1, fill=method)) +
  facet_wrap( ~rho) + 
  geom_boxplot() +
  scale_fill_manual(values=c(
    colorRampPalette(c("#FF99CC", "#B266FF"))(9),
    "#FF3333",
    "#0080FF"
  )) + 
  ylim(0,1)


