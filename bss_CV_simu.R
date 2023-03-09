#' Simulate a low/high dimensional block setting and use Best subset selection (BSS)
#' to find the best subset based on different time limits and subset sizes. 
#' 
#' @param N number of observations
#' @param P number of variables
#' @param s number of non-zero coefficients
#' @param SNR signal-to-noise-ratio
#' @param RHO correlation strength between variables 
#' @param Sim_n number of simulations
#' @param max.k maximal subset size
#' @param Time.limits vector of different time limits for each subset size
#' @param mc number of workers for parallel computation


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
cv.number <- 10
# 
# mc<-105 #Optional: number of cores can be determined automatically
# 
# cl <- makeCluster(mc, type="MPI")
# 
# 
# clusterExport(cl, c("N",
#                     "P",
#                     "s",
#                     "SNR",
#                     "RHO",
#                     "Sim_n",  
#                     "max.k",
#                     "beta_position",
#                     "corr_type",
#                     "block_builder",
#                     "use_gurobi",
#                     "cv.number"))
# 
# clusterEvalQ(cl, {
#   library(glmnet)
#   library(bestsubset)
#   library(tibble)
#   library(tidyverse)
#   library(mvtnorm)
#   library(bestridge)
# })


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
      BSS <- bsrr(x=X, y=Y, 
                  family="gaussian", 
                  method="sequential", 
                  tune="cv",
                  s.list=1:max.k,
                  lambda.list = 0, nfolds=cv.number)
      
      estimated_non_zeros_BSS <- which(BSS$beta != 0)
      
      TP <- sum(estimated_non_zeros_BSS %in% non_zero_indices)
      FP <- length(estimated_non_zeros_BSS)-TP
      FN <- s-TP
      F1 <- TP/(TP + 0.5*(FP+FN))
      Precision <- TP/(TP + FP)
      Accuracy <- TP/s
      
      if(length(estimated_non_zeros_BSS) <= 1){
        est_det_IFIM <- 1
        est_det_FIM <- 1
      }else{
        est_det_IFIM <- det(cov(X[,estimated_non_zeros_BSS]))
        est_det_FIM <- det(solve(cov(X[,estimated_non_zeros_BSS])))
      }
      
      
      BSS_results <- tibble(
        method="BSS",
        alpha=NA,
        lambda = NA,
        RSS = sum((Y-X %*% BSS$beta)^2),
        RSS_corrected = NA,
        est_det_IFIM = est_det_IFIM,
        est_det_FIM = est_det_FIM,
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
        cv.number = cv.number
      )
      rm(estimated_non_zeros_BSS)
      
      
      ### begin FSS
      #### new
      data_YX <- cbind(Y,X)
      colnames(data_YX) <- c("Y", paste("Var", 1:ncol(X), sep=""))
      train.control <- trainControl(method = "cv", number = cv.number)
      step.model <- train(y = data_YX[,1], x = data_YX[,-1], 
                          method = "leapForward",
                          trace = T, 
                          trControl = train.control,
                          tuneGrid = expand.grid(nvmax = seq(1, max.k, 1))
                          
      )
      
      step.model$bestTune
      
      
      FSS <- regsubsets(x = X, y = Y, nvmax = max.k, 
                        method = "forward", matrix=T)
      estimated_non_zeros_FSS <- 
        which(
          summary(FSS)$which[as.numeric(step.model$bestTune), -1] == T
          )
      
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
      
      FSS_beta <- 
        solve(t(X[,estimated_non_zeros_FSS]) %*% X[,estimated_non_zeros_FSS]) %*% 
        t(X[,estimated_non_zeros_FSS]) %*% Y
      
      FSS_results <- tibble(
        method="FSS",
        alpha=NA,
        lambda = NA,
        RSS = sum((Y-X[,estimated_non_zeros_FSS] %*% FSS_beta)^2),
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
        cv.number = cv.number
      )
      rm(estimated_non_zeros_FSS)
      
      
      
      ### begin Lasso/Enet
      #### new
      cv_Enet <- cva.glmnet(x=X, y=Y, alpha = Alpha, nfolds = cv.number)
      
      best_Enets <- lapply(seq_along(Alpha), function(i){
        glmnet(x=X, y=Y, alpha = Alpha[i], 
               lambda = cv_Enet$modlist[[i]]$lambda.1se, 
               intercept = F)
      })
      
      Enet_results <- lapply(seq_along(Alpha), function(i){
        estimated_non_zeros_Enet <- 
          which(best_Enets[[i]]$beta != 0)
        
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
          method=paste("Enet", Alpha[i], sep=" "),
          alpha=Alpha[i],
          lambda = best_Enets[[i]]$lambda,
          RSS = sum((Y-X %*% best_Enets[[i]]$beta)^2),
          RSS_corrected = 
            sum((Y-X %*% ((1+(1-Alpha[i])*best_Enets[[i]]$lambda)*
                            best_Enets[[i]]$beta))^2),
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
          cv.number = cv.number
        )
        
        
      })
      
      Enet_results <- do.call(rbind, Enet_results)
      
      
      out <- rbind(BSS_results,
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

# save results
saveRDS(Loop_Sim_n, 
        paste("/home/hanke/Programme/RSS_all_methods/RSS_",P,"_",N,"_", 
              beta_position,
              "_",
              corr_type,
              ".RDS",
              sep=""))

stopCluster(cl)
mpi.quit()


ggplot(Loop_Sim_n_n1000_p100_toep_spread,
       aes(x = as.factor(snr), y = F1, fill=method)) +
  facet_wrap( ~rho) + 
  geom_boxplot() +
  scale_fill_manual(values=c(
    colorRampPalette(c("#FF99CC", "#B266FF"))(9),
    "#FF3333",
    "#0080FF",
    "#00CC00"
  )) 



library(stabs)
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

glmnet.fss <- function(x, y, q, alpha, type = c("conservative", "anticonservative"), ...) {
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
      
      
      ### new BSS
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





### GUROBI-Zeugs
library( parallel)
library( snow)
library( Rmpi)


library(glmnet)
library(bestsubset)
library(mvtnorm)
library(tibble)
library(tidyverse)


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

CORR_TYPE <- "toeplitz"
DIM <- c("low", "high")
BETA_POS <- c("spread")

s <- 10

SNR <- c(0.05, 0.09, 0.14, 0.25, 0.42, 0.71, 1.22, 2.07, 3.52, 6)

Sim_n <- 100

ALPHA <- c(0.1, 
           0.2, 0.3, 0.4, 
           0.5, 
           0.6, 0.7, 0.8, 
           0.9, 1)

nLambda <- 1000
max.k <- 15 

RHO <- c(0.35, 0.7)


mc<-105 #Optional: number of cores can be determined automatically
cl <- makeCluster(mc, type="MPI")


clusterExport(cl, c("s",
                    "RHO",
                    "SNR", 
                    "CORR_TYPE",
                    "DIM",
                    "BETA_POS",
                    "Sim_n", 
                    "ALPHA",  
                    "max.k", 
                    "nLambda"))

clusterEvalQ(cl, {
  library(glmnet)
  library(bestsubset)
  library(tibble)
  library(tidyverse)
  library(mvtnorm)
})



Loop_Dim <- lapply(DIM, function(Dim){
  Loop_Beta_pos <- lapply(BETA_POS, function(Beta_pos){
    Loop_Corr_type <- lapply(CORR_TYPE, function(Corr_type){
      Loop_Rho <- lapply(RHO, function(rho){
        Loop_Snr <- lapply(SNR, function(snr){
          Loop_Sim_n <- parLapply(cl, 1:Sim_n, function(sim_n){

            if(Dim == "low"){
              p <- 100
              n <- 1000
            }else if(Dim == "high"){
              p <- 1000
              n <- 100
            }

            if(Beta_pos == "adjacent"){
              non_zero_indices <- 1:s
            }else if(Beta_pos == "spread"){
              non_zero_indices <- seq(1,p, p/s)
            }

            if(Corr_type == "block"){
              # Block-design
              Sigma <- diag(rep(1, p))
              for(i in 1:(p/s)){
                j <- s*(i-1)
                Sigma[(j+1):(j+s), (j+1):(j+s)] <- rho
                }
              diag(Sigma) <- 1
            }else if(Corr_type == "toeplitz"){
              Sigma <- rho^toeplitz(0:(p-1))
            }else if(Corr_type == "independent"){
              Sigma <- diag(p)
            }


              beta <- rep(0, p)
              beta[non_zero_indices] <- 1
              
              
              sigma <- sqrt( (t(beta) %*% Sigma %*% beta) / snr )
              
              
              seed <- round(sim_n+snr*1000000)
              try(set.seed(seed))
              e <- rnorm(n) %*% sigma 
              
              seed <- round(sim_n+rho*1000)
              X <- mvtnorm::rmvnorm(n = n, mean = rep(0, ncol(Sigma)), sigma = Sigma)
              Y <- X %*% beta + e 
              Y <- scale(Y, scale = FALSE)
              
              ### begin BSS
              BSS <- bs(x=X, y=Y, intercept=FALSE, time.limit=180, k=1:max.k, verbose=F)
              
              BSS_results <- 
                lapply(1:max.k, function(x){
                  RSS <- sum((Y - X %*% BSS$beta[,x])^2)
                  
                  TP <- sum(which(abs(BSS$beta[,x]) > 0.00001) %in% non_zero_indices)
                  FP <- x-TP
                  FN <- s-TP
                  F1 <- TP/(TP + 0.5*(FP+FN)) 
                  Precision <- TP/(TP + FP)
                  Accuracy <- TP/s
                  
                  tibble(
                    # method = paste("Enet ", alpha, sep=""),
                    method="BSS",
                    
                    # alpha=alpha,
                    alpha=NA,
                    
                    k=x,
                    RSS = RSS,
                    TP,
                    FP,
                    FN,
                    F1,
                    Precision,
                    Accuracy,
                    n,
                    p,
                    s,
                    snr,
                    corr_type = Corr_type,
                    rho,
                    beta_pos = Beta_pos,
                    beta_switch = NA,

                    status = BSS$status[x],
                    # status = NA,
                    
                    sim.n = sim_n
                  )
                  
                })
               
              BSS_results <- 
                do.call(rbind, BSS_results)
              
              # end BSS
              
              
              # begin FSS
              FSS <- fs(x=X, y=Y, maxsteps= max.k, intercept=FALSE, verbose=F)
              
              
              FSS_results <- 
                lapply(1:max.k, function(x){
                  RSS <- sum((Y - X %*% FSS$beta[,x+1])^2)
                  TP <- sum(which(abs(FSS$beta[,x+1]) > 0.00001) %in% non_zero_indices)
                  FP <- x-TP
                  FN <- s-TP
                  F1 <- TP/(TP + 0.5*(FP+FN)) 
                  Precision <- TP/(TP + FP)
                  Accuracy <- TP/s
                  
                  tibble(
                    # method = paste("Enet ", alpha, sep=""),
                    method="FSS",
                    
                    # alpha=alpha,
                    alpha=NA,
                    
                    k=x,
                    RSS = RSS,
                    TP,
                    FP,
                    FN,
                    F1,
                    Precision,
                    Accuracy,
                    n,
                    p,
                    s,
                    snr,
                    corr_type = Corr_type,
                    rho,
                    beta_pos = Beta_pos,
                    beta_switch = NA,
                    
                    # status = BSS$status[x],
                    status = NA,
                    
                    sim.n = sim_n
                  )
                  
                })
              
              FSS_results <- 
                do.call(rbind, FSS_results)
              
              # end FSS
              
              
              
              ### begin Enet
              
              Enet_results <- lapply(ALPHA, function(alpha){

                fit_enet <- glmnet(X, Y, alpha = alpha, nlambda = nLambda)
                

                non_zero_betas <- 
                  apply(fit_enet$beta, 2, function(x){as.numeric(x != 0)})

                beta_switches <- c(sapply(1:(ncol(non_zero_betas)-1), function(i){
                  non_zero_diffs <- non_zero_betas[,i+1]-non_zero_betas[,i+1]
                  sum(non_zero_diffs < 0)
                }), 0)
                

                

                # enet_betas <- as.numeric(apply(fit_enet$beta, 2, function(x){
                #   sum(x != 0) 
                # }))
                
                # bigger_subset <- sapply(2:length(enet_betas), function(i){
                #   enet_betas[i] > enet_betas[i-1]
                # })
                
                # unique_enet_subsets <- c(seq_along(enet_betas)[-1])[bigger_subset]
                


                enet_alpha <- lapply(1:ncol(fit_enet$beta), function(x){
                  
                  RSS <- sum((Y - X %*% fit_enet$beta[,x])^2)

                  TP <- sum(which(abs(fit_enet$beta[,x]) > 0.000001) %in% non_zero_indices)
                  FP <- sum(abs(fit_enet$beta[,x]) > 0.000001) - TP
                  FN <- s - TP
                  F1 <- TP/(TP + 0.5*(FP+FN)) 
                  Precision <- TP/(TP + FP)
                  Accuracy <- TP/s
                  
                  tibble(
                    method = paste("Enet ", alpha, sep=""),
                    # method="BSS",
                    
                    alpha=alpha,
                    # alpha=NA,
                    
                    k=sum(abs(fit_enet$beta[,x]) > 0.0000001),
                    RSS = RSS,
                    TP,
                    FP,
                    FN,
                    F1,
                    Precision,
                    Accuracy,
                    n,
                    p,
                    s,
                    snr,
                    corr_type = Corr_type,
                    rho,
                    beta_pos = Beta_pos,
                    beta_switch = beta_switches[x],
                    
                    # status = BSS$status[x],
                    status = NA,
                    
                    sim.n = sim_n
                  )
                  
                  
                })
                
                enet_alpha <- do.call(rbind, enet_alpha)
              })


              Enet_results <- do.call(rbind, Enet_results)
              ### end Enet
              
              rbind(BSS_results, FSS_results, Enet_results)
             
              
              
              # write_csv(BSS$status[1], 
              #           path ="~/test_gurobi.txt", 
              #           append = TRUE) 
              
              
              
              
          })
          Loop_Sim_n <- do.call(rbind, Loop_Sim_n)
          Loop_Sim_n 
          
        })
        Loop_Snr <- do.call(rbind, Loop_Snr)
        Loop_Snr
        

      })
      Loop_Rho <- do.call(rbind, Loop_Rho)
        Loop_Rho

        saveRDS(Loop_Rho,
                paste("/home/hanke/Programme/RSS_all_methods/SimuSynthetic_RSS_",
                      Corr_type, "_",
                      Dim, "_",
                      Beta_pos,
                      ".RDS",
                      sep=""))
    })
    
  })
  
})



stopCluster(cl)
mpi.quit()
