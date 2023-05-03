######################################################
#
# Experiment comparing accuracy over iterations
#    using PHMRC data
#
######################################################

source("functions.R")

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

n0 = dim(X0)[1]
p = dim(X0)[2]
C = max(Y0)
alpha <- rep(1, C)
ac <- bc <- matrix(1, C, p)
Pi <- rep(1/C, C)

set.seed(1)
rand <- sample(1:n0)

K <- 10
out1 <- out2 <- NULL

for(k in 1:K){
       X_train <- X0[rand %% K != k - 1, ]
       X_test <- X0[rand %% K == k - 1, ]
       Y_train <- Y0[rand %% K != k - 1]
       Y_test <- Y0[rand %% K == k - 1]
       N2 <- dim(X_test)[1]
       
       ### obtain posterior samples
       samp = post_conjugate(X_train, Y_train, alpha, ac, bc, N2, C)
       Pihat = colMeans(samp$Pi)
       thetahat = rowMeans(samp$theta, dims=2)

       # Point estimate
       ppmean_post21 = active_q_selection(X_test, Y_test, N = 1, Pihat, thetahat, p, 1, method = list("rand", "PWKL") , probs=c(0.025, 0.975))      
       # Pscore
       ppmean_post31 = simu1cp_post(samp, Pi = Pihat, theta = thetahat, X_test, Y_test, N = 1, p, method = list("rand", "PWKL") , probs=c(0.025, 0.975))
       # save results
       out1[[k]] <- ppmean_post21
       out2[[k]] <- ppmean_post31
       save(ppmean_post21, ppmean_post31, file = paste0("output/PHMRC-CV-", k, ".RData"))
}

# get average curves
output <- NULL
K <- 10
ind <- 0:168
for(k in 1:K){
       load(paste0("output/PHMRC-CV-", k, ".RData"))
       tmp1 <- data.frame(itr = ind, acc = ppmean_post21$pcv_rand, type = "Fixed Order", fold = k)
       tmp2 <- data.frame(itr = ind, acc = ppmean_post21$pcv_PWKL, type = "Active Order: Point Estimate", fold = k)
       tmp3 <- data.frame(itr = ind, acc = ppmean_post31$pcv_PWKL, type = "Active Order: Posterior Predictive", fold = k)
       tmp <- rbind(tmp1, tmp2, tmp3)
       output <- rbind(output, tmp)
}
output <- aggregate(acc~type + itr, data = output, FUN = mean)
output$type <- factor(output$type, levels = c("Fixed Order", 
                             "Active Order: Point Estimate", 
                             "Active Order: Posterior Predictive"))
library(ggplot2)
g <- ggplot(output) + 
       geom_line(aes(x = itr, y = acc, color = type), linewidth = 1.2) + 
       scale_color_manual("", values = c("#000000",  "#377eb8", "#4daf4a")) +
       theme_bw() + 
       theme(legend.position = c(0.7, 0.2)) +
       xlab("Number of Questioned Asked") + 
       ylab("Probability of Correct Classification")
g
ggsave(g, file = "fig/5-experiment-phmrc-nonstop.pdf", width = 5, height = 4)




# # plot
# pdf(file, width = 8, height = 6, onefile = F, paper = "special")
# par(mfrow=c(1, 1))
# ind = 0:p
# plot(NULL, xlim=range(ind), ylim=c(0,1),  xlab="number of questions asked", ylab="pattern recovery rate", main=bquote(p==~.(p)~","~C==.(C)~','~n[train]==.(length(Y_train))))
# lines(ind,ppmean_post21$pcv_rand, col="black", lwd = 2)
# lines(ind,ppmean_post21$pcv_PWKL, col="darkorange", lwd = 2)
# lines(ind,ppmean_post31$pcv_PWKL, col="dodgerblue", lwd = 2)
# legend("topright", 
#        legend=c("fixed order",
#                 "active order: point estimate", 
#                 "active order: posterior predictive"), 
#        col=c("black", "darkorange", "dodgerblue"), 
#        lty=c(1, 1, 1), cex=1)
# dev.off()

