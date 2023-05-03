######################################################
#
# Experiment comparing accuracy over iterations (order penalty)
#
######################################################
source("functions.R")

ppmean_jump21  <- NULL
ppmean_jump22  <- NULL

### 
set.seed(1)
lambda_list = c(2,10)
h_list = c(2,10)

for(case in 1:2){
  
  if(case == 1){
    ### load simulated data (misspecified)
    load('input_data/misspecified.RData')
    smlted = smlted_misspecified
    file = 'fig/sim_nonstop_misspecified.pdf'
    
  }else{
    ### load simulated data (correct, different theta)
    load('input_data/correct.RData')
    smlted = smlted_correct
    file = 'fig/sim_nonstop_correct.pdf'
  }
  
  
  ##################################################################################################################
  ##################################################################################################################
  ##################################################################################################################
  
  ### load parameters
  alpha=smlted$alpha; ac=smlted$ac; bc=smlted$bc; Pi=smlted$Pi; theta=smlted$theta; Y=smlted$Y; X=smlted$X; nc=smlted$nc; X_train=smlted$X_train; Y_train=smlted$Y_train; X_test0=smlted$X_test; Y_test0=smlted$Y_test
  
  ### subset testing data
  test_size = 200
  testid2  = sample(1:(n-train_size), size=test_size, replace=FALSE)
  X_test = X_test0[testid2, ]
  Y_test = Y_test0[testid2]
  
  
  ### set experiments parameters. When ite_rand=1, "rand" method gives serial order.
  method=list("rand","serial" ,"PWKL", "PWKL_penalized", "PWKL_post", "PWKL_penalized_post")
  probs = c(0.025, 0.975)
  
  ### obtain posterior samples
  samp = post_conjugate(X_train, Y_train, alpha, ac, bc, N2, C)
  Pihat = colMeans(samp$Pi)
  thetahat = rowMeans(samp$theta, dims=2)
  
  ### run experiments
  h = h_list[1]
  ppmean_jump21[[case]] = experiment_orderpenalty(X_test, Y_test, Pihat, thetahat, st, h, lambda_list, method, samp)
  
  h = h_list[2]
  ppmean_jump22[[case]] = experiment_orderpenalty(X_test, Y_test, Pihat, thetahat, st, h, lambda_list, method, samp)
  
}

output <- NULL
ind = 0:st2
for(case in 1:2){
  sim <- c( "Misspecified", "Correct")[case]
  tmp1 <- data.frame(itr = ind, acc = ppmean_jump21[[case]]$pcv_rand, type = "Random Order", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[1]))
  tmp2 <- data.frame(itr = ind, acc = ppmean_jump21[[case]]$pcv_serial, type = "Fixed Order", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[1]))
  tmp3 <- data.frame(itr = ind, acc = ppmean_jump21[[case]]$pcv_PWKL, type = "Active Order: Point Estimate, lambda = 0", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[1]))
  tmp4 <- data.frame(itr = ind, acc = ppmean_jump21[[case]]$pcv_Penalized_PWKL_mat[,1], type = "Active Order: Point Estimate, lambda = 2", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[1]))
  tmp5 <- data.frame(itr = ind, acc = ppmean_jump21[[case]]$pcv_Penalized_PWKL_mat[,2], type = "Active Order: Point Estimate, lambda = 10", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[1]))
  tmp6 <- data.frame(itr = ind, acc = ppmean_jump21[[case]]$pcv_PWKL_post, type = "Active Order: Posterior Predictive, lambda = 0", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[1]))
  tmp7 <- data.frame(itr = ind, acc = ppmean_jump21[[case]]$pcv_Penalized_PWKL_post_mat[,1], type = "Active Order: Posterior Predictive, lambda = 2", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[1]))
  tmp8 <- data.frame(itr = ind, acc = ppmean_jump21[[case]]$pcv_Penalized_PWKL_post_mat[,2], type = "Active Order: Posterior Predictive, lambda = 10", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[1]))
  tmp <- rbind(tmp1, tmp2, tmp3, tmp6, tmp4, tmp5, tmp7, tmp8)
  output <- rbind(output, tmp)
}
for(case in 1:2){
  sim <- c( "Misspecified", "Correct")[case]
  tmp1 <- data.frame(itr = ind, acc = ppmean_jump22[[case]]$pcv_rand, type = "Random Order", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[2]))
  tmp2 <- data.frame(itr = ind, acc = ppmean_jump22[[case]]$pcv_serial, type = "Fixed Order", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[2]))
  tmp3 <- data.frame(itr = ind, acc = ppmean_jump22[[case]]$pcv_PWKL, type = "Active Order: Point Estimate, lambda = 0", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[2]))
  tmp4 <- data.frame(itr = ind, acc = ppmean_jump22[[case]]$pcv_Penalized_PWKL_mat[,1], type = "Active Order: Point Estimate, lambda = 2", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[2]))
  tmp5 <- data.frame(itr = ind, acc = ppmean_jump22[[case]]$pcv_Penalized_PWKL_mat[,2], type = "Active Order: Point Estimate, lambda = 10", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[2]))
  tmp6 <- data.frame(itr = ind, acc = ppmean_jump22[[case]]$pcv_PWKL_post, type = "Active Order: Posterior Predictive, lambda = 0", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[2]))
  tmp7 <- data.frame(itr = ind, acc = ppmean_jump22[[case]]$pcv_Penalized_PWKL_post_mat[,1], type = "Active Order: Posterior Predictive, lambda = 2", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[2]))
  tmp8 <- data.frame(itr = ind, acc = ppmean_jump22[[case]]$pcv_Penalized_PWKL_post_mat[,2], type = "Active Order: Posterior Predictive, lambda = 10", sim = sim, n = paste0("n_train = ", train_size), h = paste0("h = ", h_list[2]))
  tmp <- rbind(tmp1, tmp2, tmp3, tmp6, tmp4, tmp5, tmp7, tmp8)
  output <- rbind(output, tmp)
}
output$n <- train_size
output$h <- factor(output$h, levels = c("h = 2", "h = 10"))
output$type <- factor(output$type, levels = c("Random Order",
                                              "Fixed Order", 
                                              "Active Order: Point Estimate, lambda = 0",
                                              "Active Order: Posterior Predictive, lambda = 0",
                                              "Active Order: Point Estimate, lambda = 2",
                                              "Active Order: Point Estimate, lambda = 10",
                                              "Active Order: Posterior Predictive, lambda = 2",
                                              "Active Order: Posterior Predictive, lambda = 10"
))
output$fixed <- grepl("Fixed", output$type)
output$lambda <- "no penalization"
output$lambda[grepl("lambda = 2", output$type)] <- "lambda = 2"
output$lambda[grepl("lambda = 10", output$type)] <- "lambda = 10"
output$lambda <- factor(output$lambda, c("no penalization", "lambda = 2", "lambda = 10"))
output$case <- "Low noise"
output$case[output$h == "h = 2"] <- "High noise"
output$case <- factor(output$case, levels = c("Low noise", "High noise"))
output$method <- "Random Order"
output$method[grepl("Fixed", output$type)] <- "Fixed Order"
output$method[grepl("Point Estimate", output$type)] <- "Active Order: Point Estimate"
output$method[grepl("Posterior Predictive", output$type)] <- "Active Order: Posterior Predictive"
output$method <- factor(output$method, c("Fixed Order", "Active Order: Point Estimate", "Active Order: Posterior Predictive"))

fixed_only <- subset(output, method == "Fixed Order")[, c("itr", "acc", "sim", "case")]

library(ggplot2)
g <- ggplot(subset(output, type != "Random Order" & method != "Fixed Order")) + 
  geom_line(aes(x = itr, y = acc, color = lambda, group = type), linewidth = 1, alpha = 0.9) + 
  geom_line(data = fixed_only, aes(x = itr, y = acc), linewidth = 1, color = "#f03b20") + 
  facet_grid(method~ sim + case) + 
  scale_color_manual("Penalization", values = c("#74a9cf", "#2b8cbe", "#045a8d")) +
  guides(color=guide_legend(nrow=1,byrow=F)) + 
  theme_bw() + 
  theme(legend.position = "bottom") +
  xlab("Number of Questioned Asked") + 
  ylab("Probability of Correct Classification")
g
ggsave(g, file = "fig/3-experiment-penalty.pdf", width = 6, height = 6)
