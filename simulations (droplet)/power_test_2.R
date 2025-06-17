
source(".../functions/CSCORE_functions.R")
source(".../functions/sctransformrho_functions.R")
source(".../functions/CoCoA_functions.R")
source(".../functions/CNM_functions.R")
library(GJRM)
library(scDECO)


library(mvtnorm)
library(pracma)


methods_use <-  c("CSCORE","sctransform-rho","CoCoA","CNM.full","GJRM","scDECO","scCOSMiX")



cells_per_patient <- 250
proportion_pretreatment <- 2215 / 5018

tau1s_use <- seq(0, 0.4, length=2) # 10
n_patients <- 10
n <- n_patients * cells_per_patient

offset1 <- rep(log(3e4), n)
offset2 <- rep(log(3e4), n)

zinfs_use <- c(0, 0.15, 0.30)
B <- 1 # 500

ourdf <- data.frame(patientID = rep(paste("patient", 1:n_patients, sep=""), each=cells_per_patient),
                    timepoint = "post",
                    y1 = 0,
                    y2 = 0)
ourdf[sample(1:n, size=floor(proportion_pretreatment*n)), "timepoint"] <- "pre"

ourdf$timepoint <- factor(ourdf$timepoint, levels = c("pre", "post"))
ourdf$patientID <- factor(ourdf$patientID)

ourdf$interaction <- interaction(ourdf$patientID, ourdf$timepoint, sep = "_")

ourdf$logn_counts <- offset1


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
eqlist <- list(eq1, eq2, eq3, eq4, eq5)



retmat <- expand.grid(tau1 = tau1s_use, method = methods_use, zinf = zinfs_use)
retmat$rejections <- 0

dataset_list <- list()

for (zinf_idx in 1:length(zinfs_use)){

  zinf_amnt <- zinfs_use[zinf_idx]
  
  kappa1 <- c(logit(zinf_amnt+0.001), 0)
  kappa2 <- c(logit(zinf_amnt+0.001), 0)
  
  for (tau1_idx in 1:length(tau1s_use)){

    tau <- c(0, tau1s_use[tau1_idx])

    for (b in 1:B){
      
      print(paste("zinf =", zinf_amnt, "tau1 =", tau[2], "b =", b))
      
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
      
      dataset_list[[paste("zinfidx_", zinf_idx, "_tauidx_", tau1_idx, "_bidx_", b, sep="")]] <- ourdf
      
      
      for (method_idx in 1:(length(methods_use)-1)){
        
        if (methods_use[method_idx] == "CSCORE"){
          
          what <- csCORE_pval(
                ourdf = ourdf,
                genecols = colnames(ourdf)[3:4],
                metacols = c("timepoint", "patientID", "logn_counts"),
                groupcol = "timepoint",
                logcountcol = "logn_counts",
                B = 100
          )

          retmat[(retmat$tau1 == tau[2]) & (retmat$method == "CSCORE") & (retmat$zinf == zinf_amnt), "rejections"] <-
            retmat[(retmat$tau1 == tau[2]) & (retmat$method == "CSCORE") & (retmat$zinf == zinf_amnt), "rejections"] + 1*(what < 0.05) / B

        }

        if (methods_use[method_idx] == "sctransform-rho"){

          what <- scrho_pval(
            ourdf = ourdf,
            genecols = colnames(ourdf)[3:4],
            metacols = c("timepoint", "patientID", "logn_counts"),
            groupcol = "timepoint",
            logcountcol = "logn_counts",
            B = 100
          )

          retmat[(retmat$tau1 == tau[2]) & (retmat$method == "sctransform-rho") & (retmat$zinf == zinf_amnt), "rejections"] <-
            retmat[(retmat$tau1 == tau[2]) & (retmat$method == "sctransform-rho") & (retmat$zinf == zinf_amnt), "rejections"] + 1*(what < 0.05) / B


        }
        
        
        if (methods_use[method_idx] == "CoCoA"){

          tempdf <- ourdf

          tempdf$timepoint_bern <- 1*(tempdf$timepoint == "post")
          
          tempdf$y1.norm <- gamlss::gamlss(formula = y1 ~ timepoint + offset(logn_counts), 
                                           sigma.formula = ~ timepoint,
                                           data = tempdf, family = gamlss.dist::NBI, trace=FALSE)$residuals
          
          tempdf$y2.norm <- gamlss::gamlss(formula = y2 ~ timepoint + offset(logn_counts), 
                                           sigma.formula = ~ timepoint, 
                                           data = tempdf, family = gamlss.dist::NBI, trace=FALSE)$residuals
          
          tempdf$Tpar <- NA
          tempdf$ZTpar <- NA
          
          cocoa_df = data.frame(X = tempdf$y1.norm, Y = tempdf$y2.norm, Z=tempdf$timepoint_bern)
          
          
          what <- get_params_mle(dat_xyz = cocoa_df)
          
          retmat[(retmat$tau1 == tau[2]) & (retmat$method == "CoCoA") & (retmat$zinf == zinf_amnt), "rejections"] <- 
            retmat[(retmat$tau1 == tau[2]) & (retmat$method == "CoCoA") & (retmat$zinf == zinf_amnt), "rejections"] + 1*(what$out2$p[nrow(what$out2)] < 0.05) / B


        }
        
        
        if (methods_use[method_idx] == "CNM.full"){
          
          tempdf$rownum <- 1:nrow(tempdf)
          
          ourdat_long <- tempdf %>%
            pivot_longer(cols = c("y1.norm", "y2.norm"),
                         names_to = "Margin",
                         values_to = "Value") %>%
            mutate(Margin = ifelse(Margin == "y1.norm", "Margin 1", "Margin 2"))
          
          zcor_matrix <- model.matrix(~ timepoint_bern, data = tempdf)
          
          cnm_res <- summary(geese(Value ~ Margin,
                                   id = rownum,
                                   data = ourdat_long,
                                   family = gaussian,
                                   corstr = "exchangeable",
                                   sformula = ~ Margin,
                                   zcor = zcor_matrix,
                                   cor.link = "fisherz"))$correlation
          
          retmat[(retmat$tau1 == tau[2]) & (retmat$method == "CNM.full") & (retmat$zinf == zinf_amnt), "rejections"] <- 
            retmat[(retmat$tau1 == tau[2]) & (retmat$method == "CNM.full") & (retmat$zinf == zinf_amnt), "rejections"] + 1*(cnm_res$p[2] < 0.05) / B
          
        }
        
        
        if (methods_use[method_idx] == "GJRM"){


          what <- gjrm(formula = eqlist, data=ourdf, margins = c("NBI", "NBI"), model = "B", offset1 = offset1, offset2 = offset2)

          retmat[(retmat$tau1 == tau[2]) & (retmat$method == "GJRM") & (retmat$zinf == zinf_amnt), "rejections"] <-
            retmat[(retmat$tau1 == tau[2]) & (retmat$method == "GJRM") & (retmat$zinf == zinf_amnt), "rejections"] + 1*(summary(what)$tableP5[2,4] < 0.05) / B

        }
        
        if (methods_use[method_idx] == "scDECO"){

          mcmc.out <- scdeco.cop(y=ourdf[, c("y1", "y2")], x=1*(ourdf$timepoint == "post"), marginals=c("ZINB", "ZINB"), 
                                 n.mcmc=100, burn=0, thin=1, offset1 = offset1, offset2 = offset2) # update offset to 6000, burn to 1000 for real test
          
          pval <- 2*(min(mean(mcmc.out[,"tau1"] < 0), mean(mcmc.out[, "tau1"] > 0)))
          
          retmat[(retmat$tau1 == tau[2]) & (retmat$method == "scDECO") & (retmat$zinf == zinf_amnt), "rejections"] <- 
            retmat[(retmat$tau1 == tau[2]) & (retmat$method == "scDECO") & (retmat$zinf == zinf_amnt), "rejections"] + 1*(pval < 0.05) / B
          
        }
        

      }

      
    }
    
    
  }
  
  
}




library(stringr)

source(".../functions/scCOSMiX_functions.R")

eqlist <- list(eq1, eq2, eq3, eq4, eq5, eq6, eq7)

for (k in 1:length(dataset_list)){
  
  tempdf <- dataset_list[[k]]
  
  dataset_name <- names(dataset_list)[k]
  
  print(dataset_name)
  
  idcs <- as.integer(str_extract_all(dataset_name, "\\d+")[[1]])
  
  what <- sccosmix(eqlist, data = tempdf, offset1 = offset1, offset2 = offset2)

  retmat[(retmat$tau1 == tau1s_use[idcs[2]]) & (retmat$method == "scCOSMiX") & (retmat$zinf == zinfs_use[idcs[1]]), "rejections"] <-
    retmat[(retmat$tau1 == tau1s_use[idcs[2]]) & (retmat$method == "scCOSMiX") & (retmat$zinf == zinfs_use[idcs[1]]), "rejections"] + 1*(summary(what)$tableP5[2,4] < 0.05) / B
  
  
}


retmat







