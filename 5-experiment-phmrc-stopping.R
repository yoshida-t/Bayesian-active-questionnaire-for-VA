######################################################
#
# Experiment comparing stopping rules
#    using PHMRC data
#
######################################################


source("functions.R")

## ---------------------------------------------------------  ##
## Prepare full PHMRC data
## ---------------------------------------------------------  ##
library(openVA)
PHMRC <- read.csv(getPHMRC_url("adult"))

# turn into binary data
binarydata <- ConvertData.phmrc(input = PHMRC, input.test = PHMRC, cause = "gs_text34")
causes <- as.character(unique(PHMRC$gs_text34))

tmp <- binarydata$output[, -c(1, 2)]
X0 <- matrix(0, dim(tmp)[1], dim(tmp)[2])
X0[tmp == "Y"] <- 1
X0[tmp == "."] <- NA
head(X0)
Y0 <- match(as.character(binarydata$output[,2]), causes)

sympnames <- colnames(binarydata$output)[-c(1, 2)]
# load better symptom names
symplist <- read.csv("data/PHMRC_symptoms.csv")
sympnames_new <- symplist[match(sympnames, symplist[,1]), 3]

colnames(X0) <- sympnames_new 

library(tidyverse)
library(caret)

C = length(unique(Y0))
n = dim(X0)[1]
p = dim(X0)[2]

set.seed(0)
folds = createFolds(Y0, k = 10, returnTrain = TRUE)

N2 = 100 # the number of samples drawn from posterior distribution 
p1st_cand = c(0.6, 0.7, 0.8, 0.9, 0.95)
d_cand = c(NA, 0.75, 0.5, 0.25, 0)
st2 = p


for (i in c(1:10)) {
  trainid = folds[[i]] # get the i-th fold
  smlted = simulateXY(n, p, C, listXY=list(X0, Y0), trainid=trainid)
  alpha=smlted$alpha; ac=smlted$ac; bc=smlted$bc; Pi=smlted$Pi; theta=smlted$theta; Y0=smlted$Y; X0=smlted$X; nc=smlted$nc; X0_train=smlted$X_train; Y0_train=smlted$Y_train; X0_test=smlted$X_test; Y0_test=smlted$Y_test
  samp = post_conjugate(X0_train, Y0_train, alpha, ac, bc, N2, C)
  
  save(smlted,file = paste0("output/fold", i, "_smlted.RData"))
  save(samp,  file = paste0("output/fold", i, "_samp.RData"))
  
  Pihat = colMeans(samp$Pi)
  thetahat = rowMeans(samp$theta, dims=2)
  
  ite_rand = 1
  method = list("rand", "PWKL")
  
  # use point estimate
  PPmeans_mean_2 = experiment_stopping_mat(X0_test, 
                                           Y0_test, 
                                           N=1, Pihat, thetahat, 
                                           st2, ite_rand, method, 
                                           probs= c(0.025, 0.975), 
                                           p1st_cand, 
                                           d_cand, 
                                           garbagefirst=FALSE,  use_posterior_draw=FALSE)
  
  tab_acc_len2_2 = experiment_stopping_tab(Y0_test, C, p1st_cand, d_cand, PPmeans_mean_2)
  file_name = paste0("output/fold",i, "_PPmeans_mean_list.RData")
  tab_name = paste0("output/fold", i, "_tab_acc_len2_2.RData")
  save(PPmeans_mean_2, file = file_name)
  save(tab_acc_len2_2, file = tab_name)
  
  
  # use posteriors with r_prop = 0.7
  PPmeans_post_0.7_2 = experiment_stopping_mat(X0_test, 
                                               Y0_test, 
                                               N=1, Pihat, thetahat, 
                                               st2, ite_rand, method, 
                                               probs = c(0.025, 0.975), 
                                               p1st_cand, 
                                               d_cand, 0.7, 
                                               samp, garbagefirst=FALSE, use_posterior_draw=TRUE)
  
  tab_acc_len3_0.7_2 = experiment_stopping_tab(Y0_test, C, p1st_cand, d_cand, PPmeans_post_0.7_2)
  
  file_name = paste0("output/fold",i, "_PPmeans_post_0.7_list.RData")
  tab_name = paste0("output/fold", i, "tab_acc_len3_0.7_2.RData")
  save(PPmeans_post_0.7_2, file = file_name)
  save(tab_acc_len3_0.7_2, file = tab_name)
  
  
  # use posterior with r_prop = 0.5
  PPmeans_post_0.5_2 = experiment_stopping_mat(X0_test, 
                                               Y0_test, 
                                               N=1, Pihat, thetahat, 
                                               st2, ite_rand, method, 
                                               probs = c(0.025, 0.975), 
                                               p1st_cand, 
                                               d_cand, 0.5, 
                                               samp, garbagefirst=FALSE, use_posterior_draw=TRUE)
  tab_acc_len3_0.5_2 = experiment_stopping_tab(Y0_test, C, p1st_cand, d_cand, PPmeans_post_0.5_2)
  
  file_name = paste0("output/fold",i, "_PPmeans_post_0.5_list.RData")
  tab_name = paste0("output/fold", i, "tab_acc_len3_0.5_2.RData")
  save(PPmeans_post_0.5_2, file = file_name)
  save(tab_acc_len3_0.5_2, file = tab_name)
}




