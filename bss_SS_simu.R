
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



Sim_n <- 100

max.k <- 15
Alpha <- seq(0.1,1,0.1)

B <- 200
cutoff <- 0.9


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
      
      q <- floor(sqrt(0.8*P))
      ### new BSS
      
     
      
      Loop_BSS_B <- lapply(1:B, function(b){
        
        set.seed(b)
        n_sample <- sample(1:N, N/2)
        
        BSS_b <- bsrr(x=X[n_sample,], y=Y[n_sample],
                      family="gaussian",
                      method="sequential",
                      tune="gic",
                      s.list=1:min(max.k, q),
                      lambda.list = 0)
        
        BSS_b$beta.all[[1]] != 0
      })
      
      prob_selction_BSS <- Reduce("+", Loop_BSS_B)/B
      estimated_non_zeros_BSS <- 
        which(apply(prob_selction_BSS, 1, function(i){any(i >= cutoff)}))
      
      TP <- sum(estimated_non_zeros_BSS %in% non_zero_indices)
      FP <- length(estimated_non_zeros_BSS)-TP
      FN <- s-TP
      F1 <- TP/(TP + 0.5*(FP+FN))
      Precision <- TP/(TP + FP)
      Accuracy <- TP/s
      
      
      
      BSS_results <- tibble(
        method="BSS",
        alpha=NA,
        lambda = NA,
        RSS = NA,
        RSS_corrected = NA,
        est_det_IFIM = NA,
        est_det_FIM = NA,
        k=length(estimated_non_zeros_BSS),
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
      rm(estimated_non_zeros_BSS)
      
      
      # FSS
      Loop_FSS_B <- lapply(1:B, function(b){
        
        set.seed(b)
        n_sample <- sample(1:N, N/2)
        
        FSS_b <- 
          lars::lars(x=X[n_sample,], y=Y[n_sample], 
                     max.steps = min(max.k q), type = "stepwise",
                     use.Gram=FALSE)
        
        maxk_selected <- unlist(FSS_b$actions)
        
        selection_matrix <- matrix(rep(0, P*max.k), ncol=max.k)
        
        for(i in seq_along(maxk_selected)){
          selection_matrix[maxk_selected[1:i],i] <- 1
        }
        
        selection_matrix
        
      })
      
      prob_selction_FSS <- Reduce("+", Loop_FSS_B)/B
      estimated_non_zeros_FSS <- 
        which(apply(prob_selction_FSS, 1, function(i){any(i >= cutoff)}))
      
      TP <- sum(estimated_non_zeros_FSS %in% non_zero_indices)
      FP <- length(estimated_non_zeros_FSS)-TP
      FN <- s-TP
      F1 <- TP/(TP + 0.5*(FP+FN))
      Precision <- TP/(TP + FP)
      Accuracy <- TP/s
      
    
      
      FSS_results <- tibble(
        method="FSS",
        alpha=NA,
        lambda = NA,
        RSS = NA,
        RSS_corrected = NA,
        est_det_IFIM = NA,
        est_det_FIM = NA,
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
        
        Enet <- glmnet(x=X, y=Y, alpha = i, 
               nlambda = 1000, 
               intercept = F)
        
        
        
        n_non_zeros <- apply(Enet$beta, 2, function(x){
          sum(x != 0)
        })
        
        
        lambda_indices <- which(n_non_zeros <= q)
        
        Enet_lambdas <- Enet$lambda[lambda_indices]
        
        Loop_Enet_B <- lapply(1:B, function(b){
          
          set.seed(b)
          n_sample <- sample(1:N, N/2)
          
          Enet_b <- 
            glmnet(x=X, y=Y, alpha = i,
                   lambda = Enet_lambdas, 
                   intercept = F)
          
          Enet_b$beta != 0
          
          
        })
        
        prob_selction_Enet <- Reduce("+", Loop_Enet_B)/B
        estimated_non_zeros_Enet <- 
          which(apply(prob_selction_Enet, 1, function(i){any(i >= cutoff)}))
        
        TP <- sum(estimated_non_zeros_Enet %in% non_zero_indices)
        FP <- length(estimated_non_zeros_Enet)-TP
        FN <- s-TP
        F1 <- TP/(TP + 0.5*(FP+FN))
        Precision <- TP/(TP + FP)
        Accuracy <- TP/s
        
        
        
        tibble(
          method=paste("Enet", i, sep=" "),
          alpha=i,
          lambda = NA,
          RSS = NA,
          RSS_corrected = NA,
          est_det_IFIM = NA,
          est_det_FIM = NA,
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
      
      
      out <- rbind(
        BSS_results,
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
         levels = c(paste("Enet", seq(0.1, 1, 0.1)), "FSS", "BSS"), 
         labels = c(paste("Enet", seq(0.1, 1, 0.1)), "FSS", "BSS"))


ggplot(Loop_Sim_n,
       aes(x = as.factor(snr), y = F1, fill=method)) +
  facet_wrap( ~rho) + 
  geom_boxplot() +
  scale_fill_manual(values=c(
    colorRampPalette(c("#FF99CC", "#B266FF"))(9),
    "#FF3333",
    "#0080FF",
    "#00CC00"
  )) + 
  ylim(0,1)

saveRDS(Loop_Sim_n, paste("~/Documents/BIPS/bss/data/SS_", P,"_",N,"_", 
                          beta_position,
                          "_",
                          corr_type,
                          ".RDS",
                          sep="")))
