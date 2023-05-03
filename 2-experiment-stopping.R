######################################################
#
# Experiment comparing different stopping rules
#
######################################################

source("functions.R")

### load simulated data (misspecified)
load('input_data/misspecified.RData')
smlted = smlted_misspecified
model = 'sim_misspecified'

### load simulated data (correct, different theta)
load('input_data/correct.RData')
smlted = smlted_correct
model = 'sim_correct'


### load parameters
alpha=smlted$alpha; ac=smlted$ac; bc=smlted$bc; Pi=smlted$Pi; theta=smlted$theta; Y=smlted$Y; X=smlted$X; nc=smlted$nc; X_train=smlted$X_train; Y_train=smlted$Y_train; X_test0=smlted$X_test; Y_test0=smlted$Y_test

### set seed
set.seed(1)

### subset testing data
test_size = 200
testid2  = sample(1:(n-train_size), size=test_size, replace=FALSE)
X_test = X_test0[testid2, ]
Y_test = Y_test0[testid2]


### set experiments parameters. When ite_rand=1, "rand" method gives serial order.
ite_rand = 1
probs=c(0.025, 0.975)
method=list("rand", "PWKL") 

### obtain posterior samples
samp = post_conjugate(X_train, Y_train, alpha, ac, bc, N2, C)
Pihat = colMeans(samp$Pi)
thetahat = rowMeans(samp$theta, dims=2)

### set criteria
p1st_cand = c(0.8, 0.9, 0.95)
d_cand = c(0.5, 0.25, 0)
p2nd = p2nd_formula(p1st, C, 0.75)
p_bounds = c(p1st, p2nd)


type_names = c("mean", "posterior (r=0.7)")
output_p1d = list(p1st = c(0.8, 0.9, 0.95), d=c(0.5, 0.25, 0))

### run stopping algorithms
# point estimate
PPmeans_mean = experiment_stopping_mat(X_test, Y_test, N, Pihat, thetahat, st2, ite_rand, method, probs, p1st_cand, d_cand, garbagefirst=FALSE)
table_stop_mean = experiment_stopping_tab(Y_test, C, p1st_cand, d_cand, PPmeans_mean)
# PScore (r=0.7)
PPmeans_post_0.7 = experiment_stopping_mat(X_test, Y_test, N, Pihat, thetahat, st2, ite_rand, method, probs, p1st_cand, d_cand, 0.7, samp, garbagefirst=FALSE, use_posterior_draw=TRUE)
table_stop_post0.7 = experiment_stopping_tab(Y_test, C, p1st_cand, d_cand, PPmeans_post_0.7)
# PScore (r=0.5)
PPmeans_post_0.5 = experiment_stopping_mat(X_test, Y_test, N, Pihat, thetahat, st2, ite_rand, method, probs, p1st_cand, d_cand, 0.5, samp, garbagefirst=FALSE, use_posterior_draw=TRUE)
table_stop_post0.5 = experiment_stopping_tab(Y_test, C, p1st_cand, d_cand, PPmeans_post_0.5)

### save output
PPmeans_list = list(PPmeans_mean, PPmeans_post_0.7, PPmeans_post_0.5)
save(n, train_size, st, N, p, C, N2, st2, trainid, X_train, X_test, Y_train, X_test, PPmeans_list, file = paste('output/', model, '_stopping.RData', sep=''))
# save tables
latex(table_stop_mean, file=paste('table/table_',model,'_mean.tex', sep=''), label=paste('tab;table_stop_',model,'_mean', sep=''), title="", digits=3, here=F, caption="Classification Accuracy for Latent Class and Test Length (with $\\hat{\\btheta}$ and $\\hat{\\bpi}$).")
latex(table_stop_post0.7, file=paste('table/table_',model,'_post0.7.tex', sep=''), label=paste('tab;table_stop_',model,'_post0.7', sep=''), title="", digits=3, here=F, caption="Classification Accuracy for Latent Class and Test Length by using posterior draw ($r=0.7$)")
latex(table_stop_post0.5, file=paste('table/table_',model,'_post0.5.tex', sep=''), label=paste('tab;table_stop_',model,'_post0.5', sep=''), title="", digits=3, here=F, caption="Classification Accuracy for Latent Class and Test Length by using posterior draw ($r=0.5$)")
# save plots
pdf(paste('fig/', model, '_scatter.pdf', sep=''), width = 6, height = 10, onefile = F, paper = "special")
stopping_scatter(p, PPmeans_list, quant_cand, mains, output_p1d)
dev.off()
