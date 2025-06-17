
source(".../functions/scCOSMiX_functions.R")

library(mvtnorm)



cells_per_patient <- 200
proportion_CD4 <- 2466 / 4032

n_patients <- 5 # 10, 15
B <- 500

n <- n_patients * cells_per_patient

offset1 <- rep(log(5e5), n)
offset2 <- rep(log(5e5), n)

ourdf <- data.frame(patientID = rep(paste("patient", 1:n_patients, sep=""), each=cells_per_patient),
                    celltype = "CD8",
                    y1 = 0,
                    y2 = 0)
ourdf[sample(1:n, size=floor(proportion_CD4*n)), "celltype"] <- "CD4"

ourdf$celltype <- factor(ourdf$celltype, levels = c("CD4", "CD8"))
ourdf$patientID <- factor(ourdf$patientID)

ourdf$interaction <- interaction(ourdf$patientID, ourdf$celltype, sep = "_")


## Define fixed effects design matrix X
X <- model.matrix(~ celltype, data = ourdf)


## Define random effects design matrices Z_theta
Z_mu1 <- model.matrix(~ interaction -1, data = ourdf)
Z_mu2 <- model.matrix(~ interaction -1, data = ourdf)
Z_sig1 <- model.matrix(~ 1, data = ourdf)
Z_sig2 <- model.matrix(~ 1, data = ourdf)
Z_rho <- model.matrix(~ interaction -1, data = ourdf)
Z_p1 <- model.matrix(~ patientID -1, data = ourdf)
Z_p2 <- model.matrix(~ patientID -1, data = ourdf)



beta1 <- c(-7.43, -1.26)
beta2 <- c(-8.26, -0.65)
alpha1 <- c(0.51, 0.65)
alpha2 <- c(0.85, 0.26)
tau <- c(0.35, -0.34)
kappa1 <- c(logit(0.001), 0) # c(-2.44, 0.63)
kappa2 <- c(logit(0.001), 0) # c(-2.09, 0.76)

paramvec <- c(beta1, beta2, alpha1, alpha2, tau, kappa1, kappa2)
names(paramvec) <- c("beta01","beta11","beta02","beta12","alpha01","alpha11",
                     "alpha02","alpha12","tau0","tau1","p01","p11","p02","p12")

paramvec[c("p11", "p12")] <- c(sum(paramvec[c("p01", "p11")]), sum(paramvec[c("p02", "p12")]))
paramvec[c("p01", "p11", "p02", "p12")] <- sigmoid(paramvec[c("p01", "p11", "p02", "p12")])



eq1 <- y1 ~ celltype + s(patientID, bs="re", by=celltype)
eq2 <- y2 ~ celltype + s(patientID, bs="re", by=celltype)
eq3 <- ~ celltype 
eq4 <- ~ celltype 
eq5 <- ~ celltype + s(patientID, bs="re", by=celltype)
eq6 <- ~ celltype + s(patientID, bs="re")
eq7 <- ~ celltype + s(patientID, bs="re")
eqlist <- list(eq1, eq2, eq3, eq4, eq5, eq6, eq7)




retmat <- data.frame(coverage=rep(0, length(paramvec)), MSE=0, MBE=0, CI_width=0)
rownames(retmat) <- names(paramvec)

tempmat <- data.frame(CIlow = rep(0, length(paramvec)), CIup = 0)
rownames(tempmat) <- rownames(retmat)


for (b in 1:B){
  
  print(paste("b =", b))
  
  ## Define random effect coeffs
  u_mu1 <- c(rnorm(n=n_patients, 0, sd=0.25), rnorm(n=n_patients, 0, sd=0.25))
  u_mu2 <- c(rnorm(n=n_patients, 0, sd=0.25), rnorm(n=n_patients, 0, sd=0.25))
  u_sig1 <- 0
  u_sig2 <- 0
  u_rho <- c(rnorm(n=n_patients, 0, sd=0.05), rnorm(n=n_patients, 0, sd=0.05))
  u_p1 <- rnorm(n=n_patients, 0, sd=0.125)
  u_p2 <- rnorm(n=n_patients, 0, sd=0.125)
  
  
  ## Define mu1, mu2, sig1, sig2, rho, p1, p2
  mu1 <- exp(c(X%*%beta1 + Z_mu1%*%u_mu1 + offset1))
  mu2 <- exp(c(X%*%beta2 + Z_mu2%*%u_mu2 + offset2))
  sig1 <- exp(c(X%*%alpha1 + Z_sig1%*%u_sig1))
  sig2 <- exp(c(X%*%alpha2 + Z_sig2%*%u_sig2))
  rho <- tanh(c(X%*%tau + Z_rho%*%u_rho))
  p1 <- sigmoid(c(X%*%kappa1 + Z_p1%*%u_p1))
  p2 <- sigmoid(c(X%*%kappa2 + Z_p2%*%u_p2))
  
  
  ## generate y1, y2
  for (i in 1:n){
    
    bivgauss <- rmvnorm(n = 1, mean = c(0,0), sigma = rbind(c(1, rho[i]), c(rho[i], 1)))[1,]
    
    prob_bivgauss <- pnorm(bivgauss)
    
    ourdf[i, "y1"] <- qnbinom(p=prob_bivgauss[1], mu = mu1[i], size=1/sig1[i]) * rbinom(n=1, size=1, prob=1-p1[i])
    ourdf[i, "y2"] <- qnbinom(p=prob_bivgauss[2], mu = mu2[i], size=1/sig2[i]) * rbinom(n=1, size=1, prob=1-p2[i])
    
  }
  
  
  result_out <- sccosmix(formula = eqlist, data=ourdf, offset1=offset1, offset2=offset2)
  
  result_summ <- summary.custom(result_out)
  rownames(result_summ) <- rownames(retmat)

  zinf1_idcs <- result_out$VC$X1.d2 + result_out$VC$X2.d2 + result_out$VC$X3.d2 + result_out$VC$X4.d2 + result_out$VC$X5.d2 + 1:2
  zinf2_idcs <- result_out$VC$X1.d2 + result_out$VC$X2.d2 + result_out$VC$X3.d2 + result_out$VC$X4.d2 + result_out$VC$X5.d2 + result_out$VC$X6.d2 + 1:2
  
  ests <- result_out$coefficients
  varcovmat <- result_out$Vb
  
  result_summ[c("p01", "p11"), "coefficients"] <- c(rbind(c(1,0), c(1,1))%*%ests[zinf1_idcs])
  result_summ[c("p02", "p12"), "coefficients"] <- c(rbind(c(1,0), c(1,1))%*%ests[zinf2_idcs])
  
  result_summ[c("p01", "p11"), "std.errors"] <- sqrt(diag(rbind(c(1,0), c(1,1))%*%result_out$Vb[zinf1_idcs, zinf1_idcs]%*%t(rbind(c(1,0), c(1,1)))))
  result_summ[c("p01", "p11"), "std.errors"] <- sqrt(diag(rbind(c(1,0), c(1,1))%*%result_out$Vb[zinf2_idcs, zinf2_idcs]%*%t(rbind(c(1,0), c(1,1)))))
  
  tempmat$CIlow <- result_summ$coefficients - 1.96 * result_summ$std.errors
  tempmat$CIup <- result_summ$coefficients + 1.96 * result_summ$std.errors  
  
  tempmat[c("p01", "p11", "p02", "p12"), ] <- sigmoid(as.matrix(tempmat[c("p01", "p11", "p02", "p12"), ]))
  result_summ[c("p01", "p11", "p02", "p12"), "coefficients"] <- sigmoid(as.matrix(result_summ[c("p01", "p11", "p02", "p12"), "coefficients"]))

  retmat$coverage <- 1*((tempmat$CIlow < result_summ$coefficients) & (tempmat$CIup > result_summ$coefficients)) / B
  retmat$MSE <- retmat$MSE + (paramvec - result_summ$coefficients)^2 / B
  retmat$MBE <- retmat$MBE + abs(paramvec - result_summ$coefficients) / B
  retmat$CI_width <- retmat$CI_width + (tempmat$CIup - tempmat$CIlow) / B


}
 




retmat[c("beta01","beta11","beta02","beta12","alpha01","alpha11","alpha02","alpha12","tau0","tau1"), ]

 




