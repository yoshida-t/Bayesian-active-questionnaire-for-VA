######################################################
#
# Experiment comparing accuracy over iterations
#
######################################################
source("functions.R")

ppmean_post11 <- ppmean_post21 <- ppmean_post31 <- NULL
ppmean_post12 <- ppmean_post22 <- ppmean_post32 <- NULL

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
	set.seed(1)
	test_size = 200
	testid2  = sample(1:(n-train_size), size=test_size, replace=FALSE)
	X_test = X_test0[testid2, ]
	Y_test = Y_test0[testid2]


	### set experiments parameters. When ite_rand=1, "rand" method gives serial order.
	ite_rand = 1
	probs=c(0.025, 0.975)
	method=list("rand", "PWKL") 


	### question selection with two different train_size
	train_size1 = 200
	train_size2 = 1000

	# train_size 100
	set.seed(1)
	subsetid1 = sample(1:train_size, size=train_size1, replace=FALSE)
	X_train_sub1 = X_train[subsetid1,]
	Y_train_sub1 = Y_train[subsetid1]

	### obtain posterior samples
	samp = post_conjugate(X_train_sub1, Y_train_sub1, alpha, ac, bc, N2, C)
	Pihat = colMeans(samp$Pi)
	thetahat = rowMeans(samp$theta, dims=2)

	ppmean_post11[[case]] = active_q_selection(X_test, Y_test, N, Pi, theta, st2, ite_rand, method, probs )            # True Pi and theta
	ppmean_post21[[case]] = active_q_selection(X_test, Y_test, N, Pihat, thetahat, st2, ite_rand, method, probs )      # Point estimate
	ppmean_post31[[case]] = simu1cp_post(samp, Pi, theta, X_test, Y_test, N, st2, method, probs)                       # Pscore
	# ppmean_post41 = simu1cp_post(samp, Pi, theta, X_test, Y_test, N, st2, method, probs,use_each_draw=TRUE)  # Probabilistic

	# train_size 300
	set.seed(1)
	subsetid2 = sample(1:train_size, size=train_size2, replace=FALSE)
	X_train_sub2 = X_train[subsetid2,]
	Y_train_sub2 = Y_train[subsetid2]

	### obtain posterior samples
	samp = post_conjugate(X_train_sub2, Y_train_sub2, alpha, ac, bc, N2, C)
	Pihat = colMeans(samp$Pi)
	thetahat = rowMeans(samp$theta, dims=2)

	ppmean_post12[[case]] = active_q_selection(X_test, Y_test, N, Pi, theta, st2, ite_rand, method, probs, garbagefirst=FALSE)
	ppmean_post22[[case]] = active_q_selection(X_test, Y_test, N, Pihat, thetahat, st2, ite_rand, method, probs, garbagefirst=FALSE)
	ppmean_post32[[case]] = simu1cp_post(samp, Pi, theta, X_test, Y_test, N, st2, method, probs, garbagefirst=FALSE)
	# ppmean_post42 = simu1cp_post(samp, Pi, theta, X_test, Y_test, N, st2, method, probs, garbagefirst=FALSE, use_each_draw=TRUE)
}

output <- NULL
ind = 0:st2
for(case in 1:2){
	sim <- c( "Misspecified", "Correct")[case]
	tmp1 <- data.frame(itr = ind, acc = ppmean_post11[[case]]$pcv_rand, type = "Fixed Order", sim = sim, n = paste0("n_train = ", train_size1))
	tmp2 <- data.frame(itr = ind, acc = ppmean_post11[[case]]$pcv_PWKL, type = "Active Order: Oracle", sim = sim, n = paste0("n_train = ", train_size1))

	# no oracle for misspecified model
	if(case == 1) tmp2 <- NULL

	tmp3 <- data.frame(itr = ind, acc = ppmean_post21[[case]]$pcv_PWKL, type = "Active Order: Point Estimate", sim = sim, n = paste0("n_train = ", train_size1))
	tmp4 <- data.frame(itr = ind, acc = ppmean_post31[[case]]$pcv_PWKL, type = "Active Order: Posterior Predictive", sim = sim, n = paste0("n_train = ", train_size1))
	tmp <- rbind(tmp1, tmp2, tmp3, tmp4)
	output <- rbind(output, tmp)
}
for(case in 1:2){
	sim <- c( "Misspecified", "Correct")[case]
	tmp1 <- data.frame(itr = ind, acc = ppmean_post12[[case]]$pcv_rand, type = "Fixed Order", sim = sim, n = paste0("n_train = ", train_size2))
	tmp2 <- data.frame(itr = ind, acc = ppmean_post12[[case]]$pcv_PWKL, type = "Active Order: Oracle", sim = sim, n = paste0("n_train = ", train_size2))

	# no oracle for misspecified model
	if(case == 1) tmp2 <- NULL

	tmp3 <- data.frame(itr = ind, acc = ppmean_post22[[case]]$pcv_PWKL, type = "Active Order: Point Estimate", sim = sim, n = paste0("n_train = ", train_size2))
	tmp4 <- data.frame(itr = ind, acc = ppmean_post32[[case]]$pcv_PWKL, type = "Active Order: Posterior Predictive", sim = sim, n = paste0("n_train = ", train_size2))
	tmp <- rbind(tmp1, tmp2, tmp3, tmp4)
	output <- rbind(output, tmp)
}
output$n <- factor(output$n, levels = c("n_train = 200", "n_train = 1000"))
output$type <- factor(output$type, levels = c("Fixed Order", 
				 "Active Order: Oracle",
				 "Active Order: Point Estimate", 
				 "Active Order: Posterior Predictive"))
output$oracle <- grepl("Oracle", output$type)

library(ggplot2)
g <- ggplot(output) + 
	geom_line(aes(x = itr, y = acc, color = type, linetype = oracle)) + 
	facet_grid(sim ~ n) + 
	scale_color_manual("", values = c("#000000", "#e41a1c", "#377eb8", "#4daf4a")) +
	scale_linetype_manual(values = c(1, 2)) + 
	guides(linetype = "none", color=guide_legend(nrow=2,byrow=F)) + 
	theme_bw() + 
	theme(legend.position = "bottom") +
	xlab("Number of Questioned Asked") + 
	ylab("Probability of Correct Classification")
ggsave(g, file = "fig/1-experiment-nonstop.pdf", width = 5, height = 4)


# # plot
# pdf(file, width = 12, height = 6, onefile = F, paper = "special")
# par(mfrow=c(1, 2))
# ind = 0:st2
# plot(NULL, xlim=range(ind), ylim=c(0,1),  xlab="Number of questions asked", ylab="Average probability of correct classification", main=bquote(n[train]==.(train_size1)))
# lines(ind,ppmean_post21$pcv_rand, col="black", lwd = 2)
# lines(ind,ppmean_post11$pcv_PWKL, col="red", lty = 2, lwd = 2)
# lines(ind,ppmean_post21$pcv_PWKL, col="darkorange", lwd = 2)
# lines(ind,ppmean_post31$pcv_PWKL, col="dodgerblue", lwd = 2)
# plot(NULL, xlim=range(ind), ylim=c(0,1),  xlab="Number of questions asked", ylab="Average probability of correct classification", main=bquote(n[train]==.(train_size2)))
# lines(ind,ppmean_post22$pcv_rand, col="black", lwd = 2)
# lines(ind,ppmean_post12$pcv_PWKL, col="red", lty = 2, lwd = 2)
# lines(ind,ppmean_post22$pcv_PWKL, col="darkorange", lwd = 2)
# lines(ind,ppmean_post32$pcv_PWKL, col="dodgerblue", lwd = 2)
# legend("bottomright", 
# 		legend=c("fixed order", 
# 				 "active order: oracle",
# 				 "active order: point estimate", 
# 				 "active order: point estimate"), 
# 		col=c("black", "red", "darkorange", "dodgerblue"), 
# 		lty=c(1, 2, 1, 1), cex=1)
# dev.off()

