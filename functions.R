library(LaplacesDemon)
library(Rcpp)
library(Hmisc)
library(ggjoy)
library(dplyr)
library(gridExtra)
library(ggplot2)
sourceCpp('cpp/PWKL.cpp')

simulateXY = function(n, p, C, n_o=n, listXY=NA, trainid=NULL, generate_option=0){
  alpha = rep(1, C)
  ac = array(1, dim=c(C,p))
  bc = array(1, dim=c(C,p))
  
  Pi = as.numeric(rdirichlet(1, alpha))
  Y = sample(1:C, n, replace=TRUE, prob=Pi)
  X = array(NA, dim=c(n,p))
  if(generate_option == 'Misspecified (K=2)' | generate_option == 'Misspecified (K=3)'){
    if(generate_option == 'Misspecified (K=2)'){
      K = 2
    }else if(generate_option == 'Misspecified (K=3)'){
      K = 3
    }
    theta = rep(list(NA),K)
    id_tmp = sample(n, n)
    for(z in 1:K){
      theta[[z]] = array(rbeta(C*p,1,1), dim=c(C,p))
      id_k = id_tmp[((z-1)*(n/K)+1):((n/K)*z)]
      for(m in id_k){
        X[m,] = rbern(p, theta[[z]][Y[m],])
      }
    }
  }else{
    theta = array(rbeta(C*p,1,1), dim=c(C,p))# *array(rep(Pi,p),dim=c(C,p))
    for(m in 1:n){
      X[m,] = rbern(p, theta[Y[m],])
    }
  }
  if(generate_option == 'Different theta'){
    theta = array(NA, dim=c(C,p))
    for(c in 1:C){
      if(c <= 3){
        theta[c,] = array(rbeta(C*p, 0.5,0.5), dim=c(1,p))
      }else if(c >= 4 & c <= 6){
        theta[c,] = array(rbeta(C*p, 2, 2), dim=c(1,p))
      }else if(c >= 7 & c <= 10){
        theta[c,] = array(rbeta(C*p, 1, 3), dim=c(1,p))
      }
    }
    for(m in 1:n){
      X[m,] = rbern(p, theta[Y[m],])
    }
  }else{
    theta = array(rbeta(C*p,1,1), dim=c(C,p))# *array(rep(Pi,p),dim=c(C,p))
    for(m in 1:n){
      X[m,] = rbern(p, theta[Y[m],])
    }
  }
  if(!is.na(listXY[1])){
    X = listXY[[1]]
    Y = listXY[[2]]
  }
  nc = rep(NA, C)
  for(c in 1:C){
    nc[c] = sum(Y==c)
  }
  
  n_u = n - n_o
  Y_o = Y[1:n_o]
  
  smltd = list(alpha=alpha, ac=ac, bc=bc, Pi=Pi, theta=theta, Y=Y, X=X, nc=nc, n_u=n_u, Y_o=Y_o)
  if(!is.null(trainid)){
    if(trainid[1]=='half'){
      trainid = 1:floor(n/2)
    }
    smltd$X_train = X[trainid,]
    smltd$Y_train = Y[trainid]
    smltd$X_test = X[-trainid,]
    smltd$Y_test = Y[-trainid]
  }
  smltd
 }
post_conjugate = function(X, Y, alpha, ac, bc, N2, C){
  p = dim(X)[2]
  nc = rep(NA, C)
  X_NA2zero = X
  X_NA2zero[is.na(X_NA2zero)] = 0
  for(c in 1:C){
    nc[c] = sum(Y==c)
  }
  addmat = array(NA, dim=c(C,p))
  for(c in 1:C){
    addmat[c,] = colSums(matrix(X_NA2zero[Y==c,],ncol=p))
  }
  acn = ac+addmat
  bcn = bc+matrix(rep(nc,p),ncol=p, byrow=F)-addmat
  samp = NULL
  samp$Pi = rdirichlet(N2, alpha+nc)
  samp$theta = array(NA, dim=c(C, p, N2))
  for(s in 1:N2){
    samp$theta[, , s] = matrix(rbeta(C*p, acn, bcn), byrow=F, nrow=C)
  }
  samp
}

PWKL = function(X, Pi, theta, st, i){
  C = length(Pi)
  p = dim(theta)[2]
  st2 = min(st,p)
  unasked = 1:p
  asked = rep(NA, st2)
  lpYit_pre = log(Pi)
  for(ite in 1:st2){
    Score = PWKL_Score_cpp(theta, unasked, lpYit_pre)
    j = na.omit(unasked[Score == max(Score[!is.na(Score)])])[1]
    asked[ite] = j
    unasked = unasked[! unasked %in% j]
    if(is.na(X[i,j])){
      lpYit_pre_ = lpYit_pre
    }else{
      lpYit_pre_ = lpYit_pre + log( ifelse(rep(X[i,j]==1,C),theta[,j],1-theta[,j]) )
    }
    lpYit_pre = lpYit_pre_ - log(sum(exp(lpYit_pre_)))
  }
  asked
}
postp_Y = function(X, Y, Pi, theta, asked=1:p, i){
  C = length(Pi)
  st = length(asked)
  Post = NULL
  Post$lpost_Y = array(NA, dim=c(st+1,C))
  Post$pp_Y = array(NA, dim=c(st+1,C))
  Post$classified = rep(NA, st+1)
  Post$correct_classification = rep(NA, st+1)
  lpY = log(Pi)
  Post$lpost_Y[1, ] = lpY
  Post$pp_Y[1, ] = Pi
  Post$classified[1] = which.max(lpY)
  Post$correct_classification[1] = ifelse(which.max(lpY) == Y[i], 1, 0)
  for(ite in 1:st){
    j = asked[ite]
    if(!is.na(X[i,j]==1)){
      lpY = lpY + log( ifelse(rep(X[i,j]==1,C),theta[,j],1-theta[,j]) )
    }
    Post$lpost_Y[ite+1, ] = lpY
    pY = exp(lpY)
    Post$pp_Y[ite+1, ] = pY/sum(pY)
    classified = which.max(lpY)
    Post$classified[ite+1] = classified
    Post$correct_classification[ite+1] = ifelse(classified == Y[i], 1, 0)
  }
  Post$pp_Yi = Post$pp_Y[,Y[i]]
  Post
}
# extract mean and interval estimate
mprob = function(pp_Yis, probs=c(0.025, 0.975)){
  ppmean_method = cbind(mean=apply(pp_Yis, 2, mean), lower=apply(pp_Yis, 2, quantile, probs=probs[1], na.rm=TRUE), upper=apply(pp_Yis, 2, quantile, probs=probs[2], na.rm=TRUE))
  ppmean_method
}
# Active question selection from (X, Y)
active_q_selection = function(X, Y, N, Pi, theta, st, ite_rand, method=list("rand","SHE", "PWKL"), probs, garbagefirst=FALSE, true_Pi=NULL, true_theta=NULL){
  if(is.null(true_Pi)){
    true_Pi = Pi
  }
  if(is.null(true_theta)){
    true_theta = theta
  }
  n = dim(X)[1]
  p = dim(X)[2]
  pp_Yis_rand = array(NA, dim=c(n,st+1,N*ite_rand))
  pp_Yis_randm = array(NA, dim=c(n,st+1))
  pp_Yis_SHE = array(NA, dim=c(n,st+1,N))
  pp_Yis_PWKL = array(NA, dim=c(n,st+1,N))
  correct_classification_rand = array(NA, dim=c(n,st+1,N))
  correct_classification_SHE = array(NA, dim=c(n,st+1,N))
  correct_classification_PWKL = array(NA, dim=c(n,st+1,N))
  ppmean = NULL
  if("rand" %in% method){
    for(l in 1:N){
      for(ir in 1:ite_rand){
        if(garbagefirst){
          asked_rand = c(1:10,sample(11:p))
        }else if(ite_rand == 1){
          asked_rand = 1:p
        }else{
          asked_rand = sample(1:p)
        }
        for(i in 1:n){
          res_rand = postp_Y(X, Y, true_Pi, true_theta, asked_rand, i)
          pp_Yis_randm[i, ] = res_rand$pp_Yi[1:(st+1)]
          correct_classification_rand[i, ,l] = res_rand$correct_classification
        }
        pp_Yis_rand[, ,(l-1)*ite_rand+ir] = pp_Yis_randm
      }
    }
    ppmean$rand = mprob(pp_Yis_rand, probs)
    ppmean$pcv_rand = colSums(correct_classification_rand)/n
  }
  if("SHE" %in% method){
    for(l in 1:N){
      for(i in 1:n){
        asked_SHE = SHE_cpp(X, Pi, true_theta, st, i)
        res_SHE = postp_Y(X, Y, true_Pi, true_theta, asked_SHE, i)
        pp_Yis_SHE[i, ,l] = res_SHE$pp_Yi[1:(st+1)]
        correct_classification_SHE[i, ,l] = res_SHE$correct_classification
      }
    }
    ppmean$SHE = mprob(pp_Yis_SHE, probs)
    ppmean$pcv_SHE = colSums(correct_classification_SHE)/n
    if("rand" %in% method){
      diff_SHE = rep(pp_Yis_SHE,ite_rand)-pp_Yis_rand
      ppmean$diff_SHE = mprob(diff_SHE, probs)
      ppmean$prop_SHE =  apply(apply(diff_SHE, 2, function(di) di>0), 2, sum)/(n*N*ite_rand)
    }
  }
  if("PWKL" %in% method){
    for(l in 1:N){
      for(i in 1:n){
        asked_PWKL = PWKL(X, Pi, true_theta, st, i)
        res_PWKL = postp_Y(X, Y, true_Pi, true_theta, asked_PWKL, i)
        pp_Yis_PWKL[i, ,l] = res_PWKL$pp_Yi[1:(st+1)]
        correct_classification_PWKL[i, ,l] = res_PWKL$correct_classification
      }
    }
    ppmean$PWKL = mprob(pp_Yis_PWKL, probs)
    ppmean$pcv_PWKL = colSums(correct_classification_PWKL)/n
    if("rand" %in% method){
      diff_PWKL = rep(pp_Yis_PWKL,ite_rand)-pp_Yis_rand
      ppmean$diff_PWKL = mprob(diff_PWKL, probs=probs)
      ppmean$prop_PWKL =  apply(apply(diff_PWKL, 2, function(di) di>0), 2, sum)/(n*N*ite_rand)
    }
  }
  ppmean
}

## Experiments
simu1cp = function(Pi, theta, n, N, st2, method=list("rand","SHE", "PWKL"), probs=c(0.025, 0.975), garbagefirst=FALSE){
  ite_rand = 1 # iteration for random order
  C = length(Pi)
  p = dim(theta)[2]
  for(l in 1:N){
    Y = sample(1:C, n, replace=TRUE, prob=Pi)
    X = array(NA, dim=c(n,p))
    for(m in 1:n){
      X[m,] = rbern(p, theta[Y[m],])
    }
  }
  ppmean = active_q_selection(X, Y, N, Pi, theta, st2, ite_rand, method, probs, garbagefirst=FALSE)
  ppmean
}
scenario_theta0 = function(mp, mC, scen=1){
  garbagefirst = FALSE
  if(scen == 1){ # theta ~ Unif(0,1)
    theta0 = array(runif(mC*mp), dim=c(mC,mp))
  }else if(scen == 2){
    theta0 = array(runif(mC*mp), dim=c(mC,mp))
    for(c in (1:10)){
      theta0[c,] = rbeta(mp, 1, 99)
    }
  }else if(scen == 3){
    theta0 = array(runif(mC*mp), dim=c(mC,mp))
    for(j in (1:10)){
      theta0[,j] = rbeta(mC, 1, 99)
    }
  }else if(scen == 4){
    theta0 = array(runif(mC*mp), dim=c(mC,mp))
    for(j in (1:50)){
      theta0[,j] = rbeta(mC, 1, 99)
    }
  }else if(scen == 5){
    theta0 = array(runif(mC*mp), dim=c(mC,mp))
    for(j in (1:80)){
      theta0[,j] = rbeta(mC, 1, 99)
    }
  }else if(scen == 6){
    theta0 = array(runif(mC*mp), dim=c(mC,mp))
    for(j in (1:20)){
      aj = runif(1,0,0.95)
      theta0[,j] = runif(mC, aj, aj+0.05)
    }
  }else if(scen == 7){
    theta0 = array(runif(mC*mp, 0.2, 0.4), dim=c(mC,mp))
    for(j in (1:20)){
      aj = runif(1,0,0.95)
      theta0[,j] = runif(mC, aj, aj+0.05)
    }
  }else if(scen == 8){
    theta0 = array(rbeta(mC*mp, 5, 10), dim=c(mC,mp))
    for(j in (1:20)){
      aj = runif(1,0,0.95)
      theta0[,j] = runif(mC, aj, aj+0.05)
    }
  }
  else if(scen == 9){
    theta0 = array(runif(mC*mp), dim=c(mC,mp))
    for(j in (1:20)){
      for(c in sample(mC,20)){
        aj = runif(1,0,0.95)
        theta0[c,j] = runif(1, aj, aj+0.05)
      }
    }
  }else if(scen == 10){
    theta0 = array(rbeta(mC*mp, 5, 10), dim=c(mC,mp))
    for(j in (1:20)){
      for(c in sample(mC,20)){
        aj = runif(1,0,0.95)
        theta0[c,j] = runif(1, aj, aj+0.05)
      }
    }
  }
  list(theta0, garbagefirst)
}
simu1 = function(n, N, st, pp, CC, scen=1, method=list("rand","SHE", "PWKL"), probs=c(0.025, 0.975)){
  ppmean_pc = NULL
  mp = max(pp)
  mC = max(CC)
  st2 = min(st,mp)
  ppmean_pc_mint = array(NA, dim=c(st2+1, 3, length(pp), length(CC)))
  if("rand" %in% method){
    ppmean_pc$rand = ppmean_pc_mint
  }
  if("SHE" %in% method){
    ppmean_pc$SHE = ppmean_pc_mint
    if("rand" %in% method){
      ppmean_pc$diff_SHE = ppmean_pc_mint
      ppmean_pc$prop_SHE = array(NA, dim=c(st2+1, length(pp), length(CC)))
    }
  }
  if("PWKL" %in% method){
    ppmean_pc$PWKL = ppmean_pc_mint
    if("rand" %in% method){
      ppmean_pc$diff_PWKL = ppmean_pc_mint
      ppmean_pc$prop_PWKL = array(NA, dim=c(st2+1, length(pp), length(CC)))
    }
  }
  Pi0 = as.numeric(rdirichlet(1, rep(0.5,mC)))
  scn = scenario_theta0(mp, mC, scen)
  theta0 = scn[[1]]
  garbagefirst = scn[[2]]
  for(m in 1:length(pp)){
    for(q in 1:length(CC)){
      p = pp[m]
      C = CC[q]
      Pi = Pi0[1:C]
      theta = theta0[1:C,1:p]
      ind2 = 1:(st2+1)
      ppmean = simu1cp(Pi, theta, n, N, st2, method, probs, garbagefirst)
      if("rand" %in% method){
        ppmean_pc$rand[ind2, , m, q] = ppmean$rand
      }
      if("SHE" %in% method){
        ppmean_pc$SHE[ind2, , m, q] = ppmean$SHE
        if("rand" %in% method){
          ppmean_pc$diff_SHE[ind2, , m, q] = ppmean$diff_SHE
          ppmean_pc$prop_SHE[ind2, m, q] = ppmean$prop_SHE
        }
      }
      if("PWKL" %in% method){
        ppmean_pc$PWKL[ind2, , m, q] = ppmean$PWKL
        if("rand" %in% method){
          ppmean_pc$diff_PWKL[ind2, , m, q] = ppmean$diff_PWKL
          ppmean_pc$prop_PWKL[ind2, m, q] = ppmean$prop_PWKL
        }
      }
    }
  }
  ppmean_pc
}
simu2 = function(n, N, st, pp, CC, SCEN, method=list("rand","SHE", "PWKL"), probs=c(0.025, 0.975)){
  ppmean_pcs = NULL
  if("rand" %in% method){
    ppmean_pcs$rand = array(NA, c(st+1, 3, length(pp), length(CC), SCEN))
    for(scen in 1:SCEN){
      ppmean_pcs$rand[, , , ,scen] = simu1(n, N, st, pp, CC, scen, method="rand", probs)$rand[,,1,1]
    }
  }
  if("SHE" %in% method){
    ppmean_pcs$SHE = array(NA, c(st+1, 3, length(pp), length(CC), SCEN))
    for(scen in 1:SCEN){
      ppmean_pcs$SHE[, , , ,scen] = simu1(n, N, st, pp, CC, scen, method="SHE", probs)$SHE[,,1,1]
    }
  }
  if("PWKL" %in% method){
    ppmean_pcs$PWKL = array(NA, c(st+1, 3, length(pp), length(CC), SCEN))
    for(scen in 1:SCEN){
      ppmean_pcs$PWKL[, , , ,scen] = simu1(n, N, st, pp, CC, scen, method="PWKL", probs)$PWKL[,,1,1]
    }
  }
  ppmean_pcs
}

## PScore
PWKL_post = function(X, samp, st, i){
  C = dim(samp$Pi)[2]
  p = dim(samp$theta)[2]
  N2 = dim(samp$Pi)[1]
  st2 = min(st,p)
  unasked = 1:p
  asked = rep(NA, st2)
  lpYit_pre_ = log(samp$Pi)
  lpYit_pre = lpYit_pre_ - matrix(rep(log(rowSums(exp(lpYit_pre_))),C),nrow=N2)
  for(ite in 1:st2){
    Score = array(NA, dim=c(N2,length(unasked)))
    for(s in 1:N2){
      Score[s,] = PWKL_Score_cpp(samp$theta[,,s], unasked, lpYit_pre[s,])
    }
    mean_Score = colMeans(Score)
    j = na.omit(unasked[mean_Score == max(mean_Score[!is.na(mean_Score)])])[1]
    asked[ite] = j
    unasked = unasked[! unasked %in% j]
    if(!is.na(X[i,j])){
      lpYit_pre_ = lpYit_pre + t(matrix(log( ifelse(rep(X[i,j]==1,C*N2),samp$theta[,j,],1-samp$theta[,j,]) ), byrow=F, nrow=C))
    }
    lpYit_pre = lpYit_pre_ - matrix(rep(log(rowSums(exp(lpYit_pre_))),C),nrow=N2)
  }
  asked
}

postp_Y_post = function(X, Y, samp, asked=1:p, i){ # postp_Y using posterior samples
  C = dim(samp$Pi)[2]
  st = length(asked)
  p = dim(samp$theta)[2]
  N2 = dim(samp$Pi)[1]
  Post = NULL
  Post$lpost_Y = array(NA, dim=c(st+1,C))
  Post$pp_Y = array(NA, dim=c(st+1,C))
  Post$classified = rep(NA, st+1)
  Post$correct_classification = rep(NA, st+1)
  pY = colMeans(samp$Pi)
  lpY = log(pY)
  lpY_mat = array(NA, dim=c(st+1,C,N2)); lpY_mat[1,,] = t(samp$Pi)
  Post$lpost_Y[1, ] = lpY
  Post$pp_Y[1, ] = pY/sum(pY)
  Post$classified[1] = which.max(lpY)
  Post$correct_classification[1] = ifelse(which.max(lpY) == Y[i], 1, 0)
  for(ite in 1:st){
    j = asked[ite]
    for(s in 1:N2){
      if(!is.na(X[i,j]==1)){
        lpY_mat[ite+1,,s] = lpY_mat[ite,,s] + log( ifelse(rep(X[i,j]==1,C),samp$theta[,j,s],1-samp$theta[,j,s]) )
      }else{
        lpY_mat[ite+1,,s] = lpY_mat[ite,,s]
      }
    }
    lpY = log(rowMeans(exp(lpY_mat[ite+1,,])))
    Post$lpost_Y[ite+1, ] = lpY
    pY = exp(lpY)
    Post$pp_Y[ite+1, ] = pY/sum(pY)
    classified = which.max(lpY)
    Post$classified[ite+1] = classified
    Post$correct_classification[ite+1] = ifelse(classified == Y[i], 1, 0)
  }
  Post$pp_Yi = Post$pp_Y[,Y[i]]
  Post
}

# using posterior samples and choose next question probabilistically
PWKL_post_probabilistic = function(X, samp, st, i){
  C = dim(samp$Pi)[2]
  p = dim(samp$theta)[2]
  N2 = dim(samp$Pi)[1]
  st2 = min(st,p)
  unasked = 1:p
  asked = rep(NA, st2)
  lpYit_pre = log(samp$Pi)
  for(ite in 1:st2){
    next_item_cand = rep(NA, N2)
    for(s in 1:N2){
      Score = PWKL_Score_cpp(samp$theta[,,s], unasked, lpYit_pre[s,])
      next_item_cand[s] = na.omit(unasked[Score == max(Score[!is.na(Score)])])[1]
    }
    j = as.numeric(sample(names(table(next_item_cand)),1,prob=table(next_item_cand)/N2, replace=TRUE))
    asked[ite] = j
    unasked = unasked[! unasked %in% j]
    if(!is.na(X[i,j])){
      lpYit_pre_ = lpYit_pre + t(matrix(log( ifelse(rep(X[i,j]==1,C*N2),samp$theta[,j,],1-samp$theta[,j,]) ), byrow=F, nrow=C))
    }
    lpYit_pre = lpYit_pre_ - matrix(rep(log(rowSums(exp(lpYit_pre_))),C),nrow=N2)
  }
  asked
}

simu1cp_post = function(samp, Pi, theta, X, Y, N, st2, method=list("rand","SHE", "PWKL"), probs=c(0.025, 0.975), garbagefirst=FALSE, use_each_draw=FALSE){
  ite_rand = 1 
  n = dim(X)[1]
  C = dim(samp$Pi)[2]
  p = dim(X)[2]
  ppmean_post = NULL
  ppmean_rand_samp = array(NA, dim=c(1,st2+1, N2))
  ppmean_post$rand = active_q_selection(X, Y, N=1, NA, NA, st2, ite_rand, method="rand", probs, garbagefirst=FALSE, true_Pi=Pi, true_theta=theta)$rand
  ppmean_post$pcv_rand = active_q_selection(X, Y, N=1, NA, NA, st2, ite_rand, method="rand", probs, garbagefirst=FALSE, true_Pi=Pi, true_theta=theta)$pcv_rand
  pp_Yis_PWKL = array(NA, dim=c(n,st2+1,N))
  correct_classification_PWKL = array(NA, dim=c(n,st2+1,N))
  for(l in 1:N){
    for(i in 1:n){
      if(use_each_draw){
        asked_PWKL = PWKL_post_probabilistic(X, samp, st2, i)
      }else{
        asked_PWKL = PWKL_post(X, samp, st2, i)
      }
      res_PWKL = postp_Y_post(X, Y, samp, asked_PWKL, i)
      pp_Yis_PWKL[i,,l] = res_PWKL$pp_Yi[1:(st2+1)]
      correct_classification_PWKL[i, ,l] = res_PWKL$correct_classification
    }
  }
  ppmean_post$PWKL = mprob(pp_Yis_PWKL, probs)
  ppmean_post$pcv_PWKL = colSums(correct_classification_PWKL)/n
  ppmean_post
}


###
### stopping
active_q_selection_stopping = function(X, Y, N, Pi, theta, st, ite_rand, method=list("rand","SHE", "PWKL"), probs, p_bounds, garbagefirst=FALSE, true_Pi=NULL, true_theta=NULL){
  st2 = min(st, p)
  if(is.null(true_Pi)){
    true_Pi = Pi
  }
  if(is.null(true_theta)){
    true_theta = theta
  }
  n = dim(X)[1]
  p = dim(X)[2]
  pp_Yis_rand = array(NA, dim=c(n,st2+1,N*ite_rand))
  pp_Yis_randm = array(NA, dim=c(n,st2+1))
  pp_Yis_SHE = array(NA, dim=c(n,st2+1,N))
  ppmean = NULL
  ppmean$test_length_PWKL = rep(NA,n)
  ppmean$assigned_PWKL = rep(NA,n)
  p1st = p_bounds[1]
  p2nd = p_bounds[2]
  if(is.na(p_bounds[2])){
    p2nd = 1
  }
  if("PWKL" %in% method){
    for(i in 1:n){
      asked_PWKL = PWKL(X, Pi, theta, st2, i)
      res_PWKL = postp_Y(X, Y, true_Pi, true_theta, asked_PWKL, i)
      test_length = 0
      largestPPLS = sort(res_PWKL$pp_Y[test_length+1,], decreasing=T)[1:2]
      while(largestPPLS[1]<p1st | largestPPLS[2]>p2nd){
        test_length = test_length + 1
        largestPPLS = sort(res_PWKL$pp_Y[test_length+1,], decreasing=T)[1:2]
        if(test_length >= p){
          break
        }
      }
      ppmean$test_length_PWKL[i] = test_length
      ppmean$assigned_PWKL[i] = which.max(res_PWKL$pp_Y[test_length+1,])
    }
  }
  ppmean
}
p2nd_formula = function(p1st, C, d){
  p2nd = (1-p1st)*((C-2)*d+1)/(C-1)
  p2nd
}
# Use each draw
active_q_selection_stopping_posterior_draw = function(X, Y, N, st, ite_rand, method=list("rand","SHE", "PWKL"), probs, p_bounds, r_prop, samp, garbagefirst=FALSE, true_Pi=NULL, true_theta=NULL){
  st2 = min(st, p)
  n = dim(X)[1]
  p = dim(X)[2]
  N2 = dim(samp$Pi)[1]
  pp_Yis_rand = array(NA, dim=c(n,st2+1,N*ite_rand))
  pp_Yis_randm = array(NA, dim=c(n,st2+1))
  pp_Yis_SHE = array(NA, dim=c(n,st2+1,N))
  ppmean = NULL
  ppmean$test_length_PWKL = rep(NA,n)
  ppmean$assigned_PWKL = rep(NA,n)
  p1st = p_bounds[1]
  p2nd = p_bounds[2]
  if(is.na(p_bounds[2])){
    p2nd = 1
  }
  TF = rep(NA, N2)
  PPLS_tab = array(NA, dim=c(st2+1,C,N2))
  if("PWKL" %in% method){
    for(i in 1:n){
      for(s in 1:N2){
        if(is.null(true_Pi)){
          true_Pi = samp$Pi[s,]
        }
        if(is.null(true_theta)){
          true_theta = samp$theta[,,s]
        }
        asked_PWKL = PWKL(X, samp$Pi[s,], samp$theta[,,s], st2, i)
        res_PWKL = postp_Y(X, Y, true_Pi, true_theta, asked_PWKL, i)
        PPLS_tab[,,s] = res_PWKL$pp_Y
      }
      test_length = 0
      c_star = which.max(rowMeans(PPLS_tab[test_length+1,,]))
      c_s = apply(PPLS_tab[test_length+1,-c_star,], 2, which.max)
      for(s in 1:N2){
        TF[s] = PPLS_tab[test_length+1,c_star,s]>p1st && PPLS_tab[test_length+1,c_s[s],s]<p2nd
      }
      while(sum(TF)/N2<r_prop){
        test_length = test_length + 1
        c_star = which.max(rowMeans(PPLS_tab[test_length+1,,]))
        c_s = apply(PPLS_tab[test_length+1,-c_star,], 2, which.max)
        for(s in 1:N2){
          TF[s] = PPLS_tab[test_length+1,c_star,s]>p1st && PPLS_tab[test_length+1,c_s[s],s]<p2nd
        }
        if(test_length >= p){
          break
        }
      }
      ppmean$test_length_PWKL[i] = test_length
      ppmean$assigned_PWKL[i] = which.max(res_PWKL$pp_Y[test_length+1,])
    }
  }
  ppmean
}
experiment_stopping_mat = function(X, Y, N, Pi, theta, st2, ite_rand, method, probs, p1st_cand, d_cand, r_prop, samp, garbagefirst=FALSE, use_posterior_draw=FALSE){
  n = dim(X)[1]
  p1st_cand_len = length(p1st_cand)
  d_cand_len = length(d_cand)
  PPmeans = NULL
  PPmeans$results = matrix(list(NA), ncol=p1st_cand_len, nrow=d_cand_len)
  for(p1i in 1:p1st_cand_len){
    for(p2i in 1:d_cand_len){
      p1st = p1st_cand[p1i]
      d = d_cand[p2i]
      p2nd = p2nd_formula(p1st, C, d)
      p2nd_disp = paste(round(p2nd,4), "$","(d=","$", d, ")")
      if(is.na(d)){
        p2nd = p1st
        p2nd_disp = '-'
      }
      p_bounds = c(p1st, p2nd)
      if(use_posterior_draw){
        ppmean = active_q_selection_stopping_posterior_draw(X, Y, N, st2, ite_rand, method, probs, p_bounds, r_prop, samp, garbagefirst=FALSE)
      }else{
        ppmean = active_q_selection_stopping(X, Y, N, Pi, theta, st2, ite_rand, method, probs, p_bounds, garbagefirst=FALSE)
      }
      PPmeans$results[p1i,p2i]=list(ppmean)
    }
  }
  rownames(PPmeans$results) = paste("p1st = ", p1st_cand, sep="")
  colnames(PPmeans$results) = paste("d = ", d_cand, sep="")
  PPmeans$n_test = length(Y)
  PPmeans$p = dim(theta)[2]
  PPmeans$C = dim(theta)[1]
  PPmeans$Y_true = Y
  PPmeans$X = X
  PPmeans$N = N
  PPmeans$p1st_cand = p1st_cand
  PPmeans$d_cand = d_cand
  if(use_posterior_draw){
    PPmeans$samp = samp
  }else{
    PPmeans$Pi = Pi
    PPmeans$theta = theta
  }
  PPmeans
}
experiment_stopping_tab = function(Y, C, p1st_cand, d_cand, PPmeans){
  n = length(Y)
  p1st_cand_len = length(p1st_cand)
  d_cand_len = length(d_cand)
  tab_acc_len = array(NA, dim=c(p1st_cand_len*d_cand_len, 7))
  for(p1i in 1:p1st_cand_len){
    for(p2i in 1:d_cand_len){
      p1st = p1st_cand[p1i]
      d = d_cand[p2i]
      p2nd = p2nd_formula(p1st, C, d)
      p2nd_disp = paste(round(p2nd,4), "$","(d=","$", d, ")")
      if(is.na(d)){
        p2nd = p1st
        p2nd_disp = '-'
      }
      ppmean = PPmeans$results[[p1i, p2i]]
      test_length_PWKL = ppmean$test_length_PWKL
      row_ind = d_cand_len*(p1i-1)+p2i
      tab_acc_len[row_ind, 1] = p1st
      tab_acc_len[row_ind, 2] = p2nd_disp
      tab_acc_len[row_ind, 3] = format(round(sum(ppmean$assigned_PWKL==Y)/n, 2), nsmall = 2)
      tab_acc_len[row_ind, 4] = as.character(format(round(mean(test_length_PWKL), 2), nsmall=3))
      tab_acc_len[row_ind, 5] = as.character(format(round(quantile(test_length_PWKL, 0.05), 3), nsmall=2))
      tab_acc_len[row_ind, 6] = as.character(format(round(quantile(test_length_PWKL, 0.5), 3), nsmall=2))
      tab_acc_len[row_ind, 7] = as.character(format(round(quantile(test_length_PWKL, 0.95), 3), nsmall=2))
    }
  }
  tab_acc_len = data.frame(tab_acc_len)
  colnames(tab_acc_len) = c("$p_{1{\\rm st}}$", "$p_{2{\\rm nd}}$", "Accuracy", "mean", "5th", "median", "95th")
  tab_acc_len
}
experiment_stopping = function(X, Y, N, Pi, theta, st2, ite_rand, method, probs, p1st_cand, d_cand, r_prop, garbagefirst=FALSE, samp=FALSE){
  n = dim(X)[1]
  p1st_cand_len = length(p1st_cand)
  d_cand_len = length(d_cand)
  test_length_PWKL = array(NA, dim=c(p1st_cand_len, d_cand_len, n))
  accuracy_PWKL = array(NA, dim=c(p1st_cand_len, d_cand_len))
  for(p1i in 1:p1st_cand_len){
    for(p2i in 1:d_cand_len){
      p1st = p1st_cand[p1i]
      d = d_cand[p2i]
      p2nd = p2nd_formula(p1st, C, d)
      p2nd_disp = paste(round(p2nd,4), "$","(d=","$", d, ")")
      if(is.na(d)){
        p2nd = p1st
        p2nd_disp = '-'
      }
      p_bounds = c(p1st, p2nd)
      if(isFALSE(samp)){
        ppmean = active_q_selection_stopping(X, Y, N, Pi, theta, st2, ite_rand, method, probs, p_bounds, garbagefirst=FALSE)
      }else{
        ppmean = active_q_selection_stopping_posterior_draw(X, Y, N, st2, ite_rand, method, probs, p_bounds, r_prop, samp)
      }
      accuracy_PWKL[p1i,p2i] = sum(ppmean$assigned_PWKL==Y)/n
      test_length_PWKL[p1i,p2i,] = ppmean$test_length_PWKL
    }
  }
  list(accuracy = accuracy_PWKL, test_length = test_length_PWKL)
}
experiment_stopping_ridgeplot = function(method, p1st_cand, d_cand, accu_length){
  p1st_cand_len = length(p1st_cand)
  d_cand_len = length(d_cand)
  test_length_PWKL = accu_length$test_length
  accuracy_PWKL = accu_length$accuracy
  plots = list()
  for(p1i in 1:p1st_cand_len){
    neach = dim(test_length_PWKL)[3]
    tl = data.frame(array(NA, dim=c(neach*d_cand_len,3)))
    p1st = p1st_cand[p1i]
    p2nds = rep(NA,d_cand_len)
    for(p2i in 1:d_cand_len){
      d = d_cand[p2i]
      p2nd = p2nd_formula(p1st, C, d)
      p2nd_disp = round(p2nd, 4)
      if(is.na(d)){
        p2nd = p1st
        p2nd_disp = '--'
      }
      p2nds[p2i] = p2nd_disp
      tl[((p2i-1)*neach+1):(p2i*neach),1] = p2nd_disp
      tl[((p2i-1)*neach+1):(p2i*neach),2] = accuracy_PWKL[p1i,p2i]
      tl[((p2i-1)*neach+1):(p2i*neach),3] = test_length_PWKL[p1i,p2i,]
    }
    colnames(tl) = c('p2nd', 'accuracy','length')
    tl$p2nd = factor(tl$p2nd, levels=names(table(tl$p2nd))[c(2:d_cand_len,1)])
    plots[[p1i]] = ggplot(tl, aes(x=length, y=p2nd, fill=accuracy)) +
      geom_joy(scale=0.9, rel_min_height=0, panel_scaling=TRUE) +
      scale_y_discrete(expand = c(0.001, 0)) +
      ylab(bquote(p['2nd'])) + 
      xlab('') + 
      xlim(c(0, p)) +
      ggtitle(bquote(p['1st']==.(p1st))) + 
      labs(fill='accuracy') +
      scale_fill_gradient(limits = range(0.5, 1)) +
      theme_minimal() + 
      theme(plot.title = element_text(hjust = 0.5), axis.ticks.y = element_blank(), panel.grid.major.y=element_blank(), axis.text.y = element_text(size=10, vjust=-1.5), axis.line.y = element_blank(), axis.title.y = element_text(size=13, margin = margin(r = 10)))
  }
  do.call("grid.arrange", c(plots,ncol=1))
}
stopping_scatter = function(p, PPmeans_list, quant_cand, mains, output_p1d=NULL){
  ll = length(PPmeans_list)
  lquant_cand = length(quant_cand)
  xlim1_ = array(NA, dim=c(ll, lquant_cand))
  ylim1_ = array(NA, dim=c(ll, lquant_cand))
  if(is.null(output_p1d)){
    p1st_cand = PPmeans_list[[1]]$p1st_cand
    d_cand = PPmeans_list[[1]]$d_cand
  }else{
    p1st_cand = output_p1d$p1st
    d_cand = output_p1d$d
  }
  n = length(PPmeans_list[[1]]$results[[1,1]][[1]])
  lpp_points = length(p1st_cand) * length(d_cand)
  length_mat = array(NA, dim=c(n, lpp_points, ll))
  accu_mat = array(NA, dim=c(lpp_points, ll))
  trueY = PPmeans_list[[1]]$Y_true
  ind_p1st = which(PPmeans_list[[1]]$p1st_cand %in% p1st_cand)
  ind_d = which(PPmeans_list[[1]]$d_cand %in% d_cand)
  for(l in 1:ll){
    for(ppl in 1:lpp_points){
      PPmeans_array = array(PPmeans_list[[l]]$results[ind_p1st, ind_d])
      length_mat[, ppl, l] = PPmeans_array[[ppl]]$test_length_PWKL
      accu_mat[ppl, l] = sum(PPmeans_array[[ppl]]$assigned_PWKL == trueY)/n
    }
  }
  xlim1_ = apply(length_mat, c(2,3), quantile, prob=quant_cand)
  xlim1 = range(xlim1_)
  ylim1 = range(accu_mat)
  xlim1[1]=xlim1[1]-1; xlim1[2]=xlim1[2]+1
  ylim1[1]=ylim1[1]-0.05; ylim1[2]=min(ylim1[2]+0.05, 1.01)
  par(mfrow=c(2,ceiling(lquant_cand/2)))
  for(qq in 1:lquant_cand){
      plot(NULL, xlim=xlim1, ylim=ylim1, xlab="number of questions asked", ylab="mean probability of correct classification", main=mains[qq])
    for(l in 1:ll){
      quant = quant_cand[qq]
      length_vec = apply(length_mat[,,l], 2, quantile, prob=quant)
      accu_vec = accu_mat[,l]
      points(length_vec, accu_vec, col=1+l)
    }
    legend("bottomright", legend=c("posterior mean","posterior samples (r=0.7)", "posterior samples (r=0.5)"), col=seq(2,4), pch=1)
  }
}
createResultsTable = function(PPmeans_list, file=filename, type_names = c("mean", "0.7", "0.5")){
  n = PPmeans_list[[1]]$n_test
  C = PPmeans_list[[1]]$C
  p1st_cand = PPmeans_list[[1]]$p1st_cand
  d_cand = PPmeans_list[[1]]$d_cand
  p1st_cand_len = length(p1st_cand)
  d_cand_len = length(d_cand)
  type_len = length(PPmeans_list)
  p2nd_disps = rep(NA, p1st_cand_len*d_cand_len)
  for(p1i in 1:p1st_cand_len){
    p1st = p1st_cand[p1i]
    for(p2i in 1:d_cand_len){
      d = d_cand[p2i]
      p2nd = p2nd_formula(p1st, C, d)
      p2nd_disp = paste(round(p2nd,4), "$","(d=","$", d, ")")
      if(is.na(d)){
        p2nd = p1st
        p2nd_disp = '-'
      }
      p2nd_disps[(p1i-1)*d_cand_len+p2i] = p2nd_disp
    }
  }
  results_tab = array(NA, dim=c(type_len*p1st_cand_len*d_cand_len*n, 7))
  results_tab[, 1] = rep(type_names, each=p1st_cand_len*d_cand_len*n)
  results_tab[, 2] = rep(p1st_cand, type_len, each=d_cand_len*n)
  results_tab[, 3] = rep(p2nd_disps, type_len, each=n)
  results_tab[, 4] = rep(1:n, type_len*p1st_cand_len*d_cand_len)
  for(typei in 1:type_len){
    PPmeans = PPmeans_list[[typei]]
    for(p1i in 1:p1st_cand_len){
      p1st = p1st_cand[p1i]
      for(p2i in 1:d_cand_len){
        ppmean = PPmeans$results[[p1i, p2i]]
        row_ind_mat=expand.grid(1:n, 1:d_cand_len, 1:p1st_cand_len, 1:type_len)
        row_ind = ((typei-1)*p1st_cand_len*d_cand_len*n + (p1i-1)*d_cand_len*n + (p2i-1)*n + 1):((typei-1)*p1st_cand_len*d_cand_len*n + (p1i-1)*d_cand_len*n + p2i*n)
        results_tab[row_ind, 5] = ppmean$test_length_PWKL
        results_tab[row_ind, 6] = PPmeans$Y_true
        results_tab[row_ind, 7] = ppmean$assigned_PWKL
      }
    }
  }
  results_tab = data.frame(results_tab)
  colnames(results_tab) = c("Type", "$p_{1{\\rm st}}$", "$p_{2{\\rm nd}}$", "ID", "Question length", "True CoD", "Predicted CoD")
  if(!is.null(file)){
    write.csv(x = results_tab, file = file, row.names = FALSE)
  }
  results_tab
}

## Jumping
PWKL_orderpenalty = function(X, Pi, theta, st, i, h, lambda, samp=NULL, use_posterior_draw=FALSE, Score_penalizing=FALSE){
  C = length(Pi)
  p = dim(theta)[2]
  st2 = min(st,p)
  unasked = 1:p
  asked = rep(NA, st2)
  if(!use_posterior_draw){
    lpYit_pre = log(Pi)
    lpYit_pre_ = lpYit_pre
  }else{
    if(is.null(samp)){
      stop('samp is NULL')
    }
    lpYit_pre_post = log(samp$Pi)
    lpYit_pre_post_ = lpYit_pre_post
  }
  X_noisei = X[i,]
  for(ite in 1:st2){
    X_tmp = X
    X_tmp[i,] = X_noisei
    if(!use_posterior_draw){
      Score = PWKL_Score_cpp(theta, unasked, lpYit_pre)
    }else{
      if(is.null(samp)){
        stop('samp is NULL')
      }
      Score_m = array(NA, dim=c(N2,length(unasked)))
      for(s in 1:N2){
        Score_m[s,] = PWKL_Score_cpp(samp$theta[,,s], unasked, lpYit_pre_post[s,])
      }
      Score = colMeans(Score_m)
    }
    Score_penalized = Score
    if(Score_penalizing & ite>1){
      Score_penalized = Score - abs(unasked-asked[ite-1])/p * lambda * max(Score)
    }
    j = na.omit(unasked[Score_penalized == max(Score_penalized[!is.na(Score_penalized)])])[1]
    if(ite>1){
      q = abs(asked[ite-1]-j)/p * (1/h)
      if(rbinom(1,1,q)==1){
        X_noisei[j] = ifelse(X[i,j]==1, 0, 1)
      }else{
        X_noisei[j] = ifelse(X[i,j]==1, 1, 0)
      }
    }
    asked[ite] = j
    unasked = unasked[! unasked %in% j]
    if(!use_posterior_draw){
      if(!is.na(X[i,j])){
        lpYit_pre_ = lpYit_pre + log( ifelse(rep(X_noisei[j]==1,C),theta[,j],1-theta[,j]) )
      }
      if(is.na(lpYit_pre_[1])){
        # break
      }
      lpYit_pre = lpYit_pre_ - log(sum(exp(lpYit_pre_)))
    }else{
      if(!is.na(X[i,j])){
        lpYit_pre_post_ = lpYit_pre_post + t(matrix(log( ifelse(rep(X[i,j]==1,C*N2),samp$theta[,j,],1-samp$theta[,j,]) ), byrow=F, nrow=C))
      }
      if(is.na(lpYit_pre_post_[1,1])){
        # break
      }
      lpYit_pre_post = lpYit_pre_post_ - matrix(rep(log(rowSums(exp(lpYit_pre_post_))),C),nrow=N2)
    }
  }
  list(asked=asked, X_noisei=X_noisei)
}

noiserate = function(X_noise_mat, X){
  n = dim(X)[1]
  p = dim(X)[2]
  sum((X-X_noise_mat)[!is.na(X)]^2)/(n*p)
}

X_noise_for_given_order = function(X, st2, h, asked){
  n = dim(X)[1]
  p = dim(X)[2]
  X_noise = X
  qq = abs(diff(asked))/p * (1/h)
  for(i in 1:n){
    for(ite in 2:st2){
      j = asked[ite]
      if(rbinom(1,1,qq[ite-1])==1){
        X_noise[i,j] = ifelse(X[i,j]==1, 0, 1)
      }else{
        X_noise[i,j] = ifelse(X[i,j]==1, 1, 0)
      }
    }
  }
  X_noise
}
pp_Yism_gen = function(X_noise, Y, Pi, theta, asked, st2){
  n = dim(X_noise)[1]
  res = NULL
  res$pp_Yism = array(NA, dim=c(n,st2+1))
  res$correct_classification_m = array(NA, dim=c(n,st2+1))
  for(i in 1:n){
    res_rand = postp_Y(X_noise, Y, Pi, theta, asked, i)
    res$pp_Yism[i, ] = res_rand$pp_Yi[1:(st2+1)]
    res$correct_classification_m[i,] = res_rand$correct_classification
  }
  res
}
X_noise_ppYism_gen_PWKL = function(X, Y, Pi, theta, st2, h, lambda, samp=NULL, use_posterior_draw=FALSE, Score_penalizing=FALSE){
  n = dim(X)[1]
  p = dim(X)[2]
  st2 = min(st, p)
  X_noise_PWKL_mat = array(NA, dim=c(n,p))
  PWKL_out = NULL
  PWKL_out$pp_Yis_PWKL = array(NA, dim=c(n,st2+1))
  PWKL_out$correct_classification_m = array(NA, dim=c(n,st2+1))
  X_noise_PWKL = X
  for(i in 1:n){
    PWKL_tmp = PWKL_orderpenalty(X, Pi, theta, st2, i, h, lambda, samp, use_posterior_draw, Score_penalizing)
    asked_PWKL = PWKL_tmp$asked
    X_noise_PWKL[i,] = PWKL_tmp$X_noisei
    if(!use_posterior_draw){
      res_PWKL = postp_Y(X_noise_PWKL, Y, Pi, theta, asked_PWKL, i)
    }else{
      res_PWKL = postp_Y_post(X_noise_PWKL, Y, samp, asked_PWKL, i)
    }
    PWKL_out$pp_Yis_PWKL[i,] = res_PWKL$pp_Yi[1:(st2+1)]
    PWKL_out$correct_classification_m[i,] = res_PWKL$correct_classification
  }
  X_noise_PWKL_mat = X_noise_PWKL
  PWKL_out$X_noise_PWKL_mat = X_noise_PWKL_mat
  PWKL_out
}

experiment_orderpenalty = function(X, Y, Pi, theta, st, h, lambda_list, method, samp=NULL){
  n = dim(X)[1]
  p = dim(X)[2]
  st2 = min(st, p)
  
  ppmean = NULL
  ppmean$pwrong = NULL
  ir = 1
  
  # random order
  if("rand" %in% method){
    asked_rand = sample(1:p)
    X_noise_rand = X_noise_for_given_order(X, st2, h, asked_rand)
    res_rand= pp_Yism_gen(X_noise_rand, Y, Pi, theta, asked_rand, st2)
    ppmean$rand = mprob(res_rand$pp_Yism, probs)
    ppmean$pwrong = cbind(ppmean$pwrong, noiserate(X_noise_rand, X))
    colnames(ppmean$pwrong)[length(ppmean$pwrong)] = 'rand'
    ppmean$pcv_rand = colSums(res_rand$correct_classification_m)/n
  }
  
  # serial order
  if("serial" %in% method){
    asked_serial = 1:p
    X_noise_serial = X_noise_for_given_order(X, st2, h, asked_serial)
    res_serial = pp_Yism_gen(X_noise_serial, Y, Pi, theta, asked_serial, st2)
    ppmean$serial = mprob(res_serial$pp_Yism, probs)
    ppmean$pwrong = cbind(ppmean$pwrong, noiserate(X_noise_serial, X))
    colnames(ppmean$pwrong)[length(ppmean$pwrong)] = 'serial'
    ppmean$pcv_serial = colSums(res_serial$correct_classification_m)/n
  }
  
  # PWKL
  if("PWKL" %in% method){
    PWKL_tmp2 = X_noise_ppYism_gen_PWKL(X, Y, Pi, theta, st2, h, lambda, use_posterior_draw=FALSE, Score_penalizing=FALSE)
    ppmean$PWKL = mprob(PWKL_tmp2$pp_Yis_PWKL, probs)
    ppmean$pwrong = cbind(ppmean$pwrong, noiserate(PWKL_tmp2$X_noise_PWKL_mat, X))
    colnames(ppmean$pwrong)[length(ppmean$pwrong)] = 'PWKL'
    ppmean$pcv_PWKL = colSums(PWKL_tmp2$correct_classification_m)/n
  }
  
  # PWKL (penalized for jumping)
  if("PWKL_penalized" %in% method){
    len_lambda_list = length(lambda_list)
    ppmean_PWKL_penalized_mat = array(NA, dim=c(st2+1, 3, len_lambda_list))
    ppmean$pcv_Penalized_PWKL_mat = array(NA, dim=c(st2+1, len_lambda_list))
    for(ll in 1:len_lambda_list){
      lambda = lambda_list[ll]
      PWKL_tmp3 = X_noise_ppYism_gen_PWKL(X, Y, Pi, theta, st2, h, lambda, use_posterior_draw=FALSE, Score_penalizing=TRUE)
      ppmean_PWKL_penalized_mat[, , ll] = mprob(PWKL_tmp3$pp_Yis_PWKL, probs)
      ppmean$pwrong = cbind(ppmean$pwrong, noiserate(PWKL_tmp3$X_noise_PWKL_mat, X))
      colnames(ppmean$pwrong)[length(ppmean$pwrong)] = paste0('Penalized PWKL', ' (lambda=', lambda, ')')
      ppmean$pcv_Penalized_PWKL_mat[,ll] = colSums(PWKL_tmp3$correct_classification_m)/n
    }
  }
  
  # PWKL (with posterior samples)
  if("PWKL_post" %in% method){
    PWKL_tmp4 = X_noise_ppYism_gen_PWKL(X, Y, Pi, theta, st2, h, lambda, samp, use_posterior_draw=TRUE, Score_penalizing=FALSE)
    ppmean$PWKL_post = mprob(PWKL_tmp4$pp_Yis_PWKL, probs)
    ppmean$pwrong = cbind(ppmean$pwrong, noiserate(PWKL_tmp4$X_noise_PWKL_mat, X))
    colnames(ppmean$pwrong)[length(ppmean$pwrong)] = 'PWKL'
    ppmean$pcv_PWKL_post = colSums(PWKL_tmp4$correct_classification_m)/n
  }
  
  # PWKL (penalized for jumping, with posterior samples)
  if("PWKL_penalized_post" %in% method){
    len_lambda_list = length(lambda_list)
    ppmean_PWKL_penalized_post_mat = array(NA, dim=c(st2+1, 3, len_lambda_list))
    ppmean$pcv_Penalized_PWKL_post_mat = array(NA, dim=c(st2+1, len_lambda_list))
    for(ll in 1:len_lambda_list){
      lambda = lambda_list[ll]
      PWKL_tmp5 = X_noise_ppYism_gen_PWKL(X, Y, Pi, theta, st2, h, lambda, samp, use_posterior_draw=TRUE, Score_penalizing=TRUE)
      ppmean_PWKL_penalized_post_mat[, , ll] = mprob(PWKL_tmp5$pp_Yis_PWKL, probs)
      ppmean$pwrong = cbind(ppmean$pwrong, noiserate(PWKL_tmp5$X_noise_PWKL_mat, X))
      colnames(ppmean$pwrong)[length(ppmean$pwrong)] = paste0('Penalized PWKL', ' (lambda=', lambda, ')')
      ppmean$pcv_Penalized_PWKL_post_mat[,ll] = colSums(PWKL_tmp5$correct_classification_m)/n
    }
  }
  ppmean
}
