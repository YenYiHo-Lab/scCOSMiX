
source(".../functions/scCOSMiX_functions.R")

library(mvtnorm)



cells_per_patient <- 250
proportion_pretreatment <- 2215 / 5018

tau1s_use <- seq(0, 0.4, length=10)
npatients_use <- c(5, 10, 15)
B <- 500


beta1 <- c(-7.26, 0.45)
beta2 <- c(-7.10, 0.20)
alpha1 <- c(-0.46, -0.19)
alpha2 <- c(-0.44, -0.04)
tau <- c(0.06, -0.20)
kappa1 <- c(-4.35, 1.68)
kappa2 <- c(-5.59, 1.04)


eq1 <- y1 ~ timepoint + s(patientID, bs="re", by=timepoint)
eq2 <- y2 ~ timepoint + s(patientID, bs="re", by=timepoint)
eq3 <- ~ timepoint 
eq4 <- ~ timepoint 
eq5 <- ~ timepoint + s(patientID, bs="re", by=timepoint)
eq6 <- ~ timepoint + s(patientID, bs="re")
eq7 <- ~ timepoint + s(patientID, bs="re")
eqlist <- list(eq1, eq2, eq3, eq4, eq5, eq6, eq7)




retmat <- expand.grid(tau1 = tau1s_use, npatients = npatients_use)
retmat$rejections <- 0

for (pat_idx in 1:length(npatients_use)){
  
  n_patients <- npatients_use[pat_idx]
  n <- n_patients * cells_per_patient

  offset1 <- rep(log(3e4), n)
  offset2 <- rep(log(3e4), n)
  
  ourdf <- data.frame(patientID = rep(paste("patient", 1:n_patients, sep=""), each=cells_per_patient),
                      timepoint = "post",
                      y1 = 0,
                      y2 = 0)
  ourdf[sample(1:n, size=floor(proportion_pretreatment*n)), "timepoint"] <- "pre"
  
  ourdf$timepoint <- factor(ourdf$timepoint, levels = c("pre", "post"))
  ourdf$patientID <- factor(ourdf$patientID)
  
  ourdf$interaction <- interaction(ourdf$patientID, ourdf$timepoint, sep = "_")
  
  
  
  ## Define fixed effects design matrix X
  X <- model.matrix(~ timepoint, data = ourdf)
  
  
  ## Define random effects design matrices Z_theta
  Z_mu1 <- model.matrix(~ interaction -1, data = ourdf)
  Z_mu2 <- model.matrix(~ interaction -1, data = ourdf)
  Z_sig1 <- model.matrix(~ 1, data = ourdf)
  Z_sig2 <- model.matrix(~ 1, data = ourdf)
  Z_rho <- model.matrix(~ interaction -1, data = ourdf)
  Z_p1 <- model.matrix(~ patientID -1, data = ourdf)
  Z_p2 <- model.matrix(~ patientID -1, data = ourdf)

  for (tau1_idx in 1:length(tau1s_use)){

    tau <- c(0, tau1s_use[tau1_idx])

    for (b in 1:B){
      
      print(paste("npatients =", n_patients, "tau1 =", tau[2], "b =", b))
      
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
      
      retmat[(retmat$tau1 == tau1s_use[tau1_idx]) & (retmat$npatients == n_patients), "rejections"] <-
        retmat[(retmat$tau1 == tau1s_use[tau1_idx]) & (retmat$npatients == n_patients), "rejections"] + 1*(summary(result_out)$tableP5[2,4] < 0.05) / B
      
    }
    
    
  }
  
  
}





retmat


