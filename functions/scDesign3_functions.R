

library(scDesign3)


fit_marginal <- function (data, predictor = "gene", mu_formula, sigma_formula, 
          family_use, n_cores, usebam = FALSE, edf_flexible = FALSE, 
          parallelization = "mcmapply", BPPARAM = NULL, trace = FALSE, 
          simplify = FALSE, filter_cells = FALSE) 
{
  count_mat <- data$count_mat
  dat_cov <- data$dat
  filtered_gene <- data$filtered_gene
  feature_names <- colnames(count_mat)
  matches <- regexpr("k\\s*=\\s*([0-9]+)", mu_formula, perl = TRUE)
  extracted_value <- regmatches(mu_formula, matches)
  extracted_K <- as.numeric(sub("k\\s*=\\s*", "", extracted_value))
  if (identical(extracted_K, numeric(0))) {
    extracted_K <- 0
  }
  num <- 100
  if (dim(count_mat)[2] > num & extracted_K >= 200 & edf_flexible == 
      TRUE) {
    edf_fitting <- TRUE
    edf_gini_genes <- sample(seq_len(dim(count_mat)[2]), 
                             num)
    edf_gini_count_mat <- count_mat[, edf_gini_genes]
    edf_gini_feature_names <- feature_names[edf_gini_genes]
    edf_flexible_genes <- seq_len(dim(count_mat)[2])[-edf_gini_genes]
    edf_flexible_count_mat <- count_mat[, -edf_gini_genes]
    edf_flexible_feature_names <- feature_names[-edf_gini_genes]
  }
  else {
    edf_fitting <- FALSE
  }
  if (length(family_use) == 1) {
    if (edf_fitting == TRUE) {
      edf_gini_family_use <- rep(family_use, length(edf_gini_feature_names))
      edf_flexible_family_use <- rep(family_use, length(edf_flexible_feature_names))
    }
    family_use <- rep(family_use, length(feature_names))
  }
  if (length(family_use) != length(feature_names)) {
    stop("The family_use must be either a single string or a vector with the same length as all features!")
  }
  fit_model_func <- function(gene, family_gene, dat_use, mu_formula, 
                             sigma_formula, predictor, count_mat, edf = NULL) {
    if (!is.null(edf)) {
      mu_formula_ex <- sub("(k\\s*=).*", "\\1", mu_formula)
      mu_formula = paste0(mu_formula_ex, round(edf[[gene]]), 
                          ")")
    }
    mgcv_formula <- stats::formula(paste0(predictor, "~", 
                                          mu_formula))
    mu_mgcvform <- grepl("s\\(", mu_formula) | grepl("te\\(", 
                                                     mu_formula)
    usebam <- usebam & mu_mgcvform
    if (usebam) {
      fitfunc = mgcv::bam
    }
    else {
      fitfunc = mgcv::gam
    }
    if (mu_mgcvform) {
      terms <- attr(stats::terms(mgcv_formula), "term.labels")
      terms_smooth <- terms[which(grepl("s\\(", terms))]
      if (usebam) {
        terms_smooth_update <- sapply(terms_smooth, 
                                      function(x) {
                                        paste0("ba(~", x, ", method = 'fREML', gc.level = 0, discrete = TRUE)")
                                      })
        if (length(terms_smooth) == length(terms)) {
          mu_formula <- stats::formula(paste0(predictor, 
                                              "~", paste0(terms_smooth_update, collapse = "+")))
        }
        else {
          terms_linear <- terms[which(!grepl("s\\(", 
                                             terms))]
          terms_update <- c(terms_linear, terms_smooth_update)
          mu_formula <- stats::formula(paste0(predictor, 
                                              "~", paste0(terms_update, collapse = "+")))
        }
      }
      else {
        terms_smooth_update <- sapply(terms_smooth, 
                                      function(x) {
                                        paste0("ga(~", x, ", method = 'REML')")
                                      })
        if (length(terms_smooth) == length(terms)) {
          mu_formula <- stats::formula(paste0(predictor, 
                                              "~", paste0(terms_smooth_update, collapse = "+")))
        }
        else {
          terms_linear <- terms[which(!grepl("s\\(", 
                                             terms))]
          terms_update <- c(terms_linear, terms_smooth_update)
          mu_formula <- stats::formula(paste0(predictor, 
                                              "~", paste0(terms_update, collapse = "+")))
        }
      }
    }
    else {
      mu_formula <- stats::formula(paste0(predictor, "~", 
                                          mu_formula))
    }
    sigma_mgcvform <- grepl("s\\(", sigma_formula) | grepl("te\\(", 
                                                           sigma_formula)
    if (sigma_mgcvform) {
      temp_sigma_formula <- stats::formula(paste0(predictor, 
                                                  "~", sigma_formula))
      terms <- attr(stats::terms(temp_sigma_formula), 
                    "term.labels")
      terms_smooth <- terms[which(grepl("s\\(", terms))]
      if (usebam) {
        terms_smooth_update <- sapply(terms_smooth, 
                                      function(x) {
                                        paste0("ba(~", x, ", method = 'fREML', gc.level = 0, discrete = TRUE)")
                                      })
        if (length(terms_smooth) == length(terms)) {
          sigma_formula <- stats::formula(paste0("~", 
                                                 paste0(terms_smooth_update, collapse = "+")))
        }
        else {
          terms_linear <- terms[which(!grepl("s\\(", 
                                             terms))]
          terms_update <- c(terms_linear, terms_smooth_update)
          sigma_formula <- stats::formula(paste0("~", 
                                                 paste0(terms_update, collapse = "+")))
        }
      }
      else {
        terms_smooth_update <- sapply(terms_smooth, 
                                      function(x) {
                                        paste0("ga(~", x, ", method = 'REML')")
                                      })
        if (length(terms_smooth) == length(terms)) {
          sigma_formula <- stats::formula(paste0("~", 
                                                 paste0(terms_smooth_update, collapse = "+")))
        }
        else {
          terms_linear <- terms[which(!grepl("s\\(", 
                                             terms))]
          terms_update <- c(terms_linear, terms_smooth_update)
          sigma_formula <- stats::formula(paste0("~", 
                                                 paste0(terms_update, collapse = "+")))
        }
      }
    }
    else {
      sigma_formula <- stats::formula(paste0("~", sigma_formula))
    }
    dat_use$gene <- count_mat[, gene]
    add_log <- function(function_name, type, message) {
      new_l <- logs
      new_log <- list(function_name = function_name, type = type, 
                      message = message)
      new_l[[length(new_l) + 1]] <- new_log
      logs <<- new_l
    }
    logs <- list()
    if (!is.null(filtered_gene) & gene %in% filtered_gene) {
      add_log("fit_marginal", "warning", paste0(gene, 
                                                "is expressed in too few cells."))
      return(list(fit = NA, warning = logs, time = c(NA, 
                                                     NA)))
    }
    if (filter_cells) {
      all_covariates <- all.vars(mgcv_formula)[-1]
      dat_cova <- dat_use[, all_covariates]
      check_factor <- all(sapply(dat_cova, is.factor))
      if (length(all_covariates) > 0 & check_factor) {
        remove_idx_list <- lapply(all_covariates, function(x) {
          curr_x <- tapply(dat_use$gene, dat_use[, x], 
                           sum)
          zero_group <- which(curr_x == 0)
          if (length(zero_group) == 0) {
            return(list(idx = NA, changeFormula = FALSE))
          }
          else {
            type <- names(curr_x)[zero_group]
            if (length(type) == length(unique(dat_use[, 
                                                      x])) - 1) {
              return(list(idx = NA, changeFormula = TRUE))
            }
            return(list(idx = which(dat_use[, x] %in% 
                                      type), changeFormula = FALSE))
          }
        })
        names(remove_idx_list) <- all_covariates
        remove_idx <- lapply(remove_idx_list, function(x) x$idx)
        remove_cell <- unlist(remove_idx)
        if (all(is.na(remove_cell))) {
          remove_cell <- NA
        }
        else {
          remove_cell <- unique(stats::na.omit(remove_cell))
        }
        if (length(remove_cell) > 0 && !any(is.na(remove_cell))) {
          dat_use <- dat_use[-remove_cell, ]
        }
        changeFormula <- sapply(remove_idx_list, function(x) x$changeFormula)
        if (length(which(changeFormula)) > 0) {
          changeVars <- names(which(changeFormula))
          formulaUpdate <- paste0(changeVars, collapse = "-")
          mgcv_formula <- stats::update.formula(mgcv_formula, 
                                                stats::as.formula(paste0("~.-", formulaUpdate)))
          mu_formula <- stats::update.formula(mu_formula, 
                                              stats::as.formula(paste0("~.-", formulaUpdate)))
          sigmaVars <- which(changeVars %in% as.character(sigma_formula))
          if (length(sigmaVars) > 0) {
            formulaUpdate <- paste0(changeVars[sigmaVars], 
                                    collapse = "-")
          }
          sigma_formula = stats::update.formula(sigma_formula, 
                                                stats::as.formula(paste0("~.-", formulaUpdate)))
        }
      }
      else {
        remove_cell <- NA
      }
    }
    else {
      remove_cell <- NA
    }
    time_list <- c(NA, NA)
    if (family_gene == "binomial") {
      mgcv.fit <- withCallingHandlers(tryCatch({
        start.time <- Sys.time()
        res <- fitfunc(formula = mgcv_formula, data = dat_use, 
                       family = "binomial", discrete = usebam)
        end.time <- Sys.time()
        time <- as.numeric(end.time - start.time)
        time_list[1] <- time
        res
      }, error = function(e) {
        add_log("gam", "error", toString(e))
        NULL
      }), warning = function(w) {
        add_log("gam", "warning", toString(w))
      })
      if (sigma_formula != "~1") {
        gamlss.fit <- withCallingHandlers(tryCatch({
          start.time = Sys.time()
          res <- gamlss::gamlss(formula = mu_formula, 
                                data = dat_use, family = gamlss.dist::BI, 
                                control = gamlss::gamlss.control(trace = FALSE, 
                                                                 c.crit = 0.1))
          end.time = Sys.time()
          time = as.numeric(end.time - start.time)
          time_list[2] <- time
          res
        }, error = function(e) {
          add_log("gamlss", "error", toString(e))
          NULL
        }), warning = function(w) {
          add_log("gamlss", "warning", toString(w))
        })
      }
      else {
        gamlss.fit <- NULL
      }
    }
    else if (family_gene == "poisson") {
      mgcv.fit <- withCallingHandlers(tryCatch({
        start.time <- Sys.time()
        res <- fitfunc(formula = mgcv_formula, data = dat_use, 
                       family = "poisson", discrete = usebam)
        end.time <- Sys.time()
        time <- as.numeric(end.time - start.time)
        time_list[1] <- time
        res
      }, error = function(e) {
        add_log("gam", "error", toString(e))
        NULL
      }), warning = function(w) {
        add_log("gam", "warning", toString(w))
      })
      if (sigma_formula != "~1") {
        gamlss.fit <- withCallingHandlers(tryCatch({
          start.time = Sys.time()
          res <- gamlss::gamlss(formula = mu_formula, 
                                data = dat_use, family = gamlss.dist::PO, 
                                control = gamlss::gamlss.control(trace = FALSE, 
                                                                 c.crit = 0.1))
          end.time = Sys.time()
          time = as.numeric(end.time - start.time)
          time_list[2] <- time
          res
        }, error = function(e) {
          add_log("gamlss", "error", toString(e))
          NULL
        }), warning = function(w) {
          add_log("gamlss", "warning", toString(w))
        })
      }
      else {
        gamlss.fit <- NULL
      }
    }
    else if (family_gene == "gaussian") {
      mgcv.fit <- withCallingHandlers(tryCatch({
        start.time <- Sys.time()
        res <- fitfunc(formula = mgcv_formula, data = dat_use, 
                       family = "gaussian", discrete = usebam)
        end.time <- Sys.time()
        time <- as.numeric(end.time - start.time)
        time_list[1] <- time
        res
      }, error = function(e) {
        add_log("gam", "error", toString(e))
        NULL
      }), warning = function(w) {
        add_log("gam", "warning", toString(w))
      })
      if (sigma_formula != "~1") {
        gamlss.fit <- withCallingHandlers(tryCatch({
          start.time = Sys.time()
          res <- gamlss::gamlss(formula = mu_formula, 
                                sigma.formula = sigma_formula, data = dat_use, 
                                family = gamlss.dist::NO, control = gamlss::gamlss.control(trace = FALSE, 
                                                                                           c.crit = 0.1))
          end.time = Sys.time()
          time = as.numeric(end.time - start.time)
          time_list[2] <- time
          res
        }, error = function(e) {
          add_log("gamlss", "error", toString(e))
          NULL
        }), warning = function(w) {
          add_log("gamlss", "warning", toString(w))
        })
      }
      else {
        gamlss.fit <- NULL
      }
    }
    else if (family_gene == "nb") {
      mgcv.fit <- withCallingHandlers(tryCatch({
        start.time <- Sys.time()
        res <- fitfunc(formula = mgcv_formula, data = dat_use, 
                       family = "nb", discrete = usebam)
        end.time <- Sys.time()
        time <- as.numeric(end.time - start.time)
        time_list[1] <- time
        res
      }, error = function(e) {
        add_log("gam", "error", toString(e))
        NULL
      }), warning = function(w) {
        add_log("gam", "warning", toString(w))
      })
      if (sigma_formula != "~1") {
        gamlss.fit <- withCallingHandlers(tryCatch({
          start.time = Sys.time()
          res <- gamlss::gamlss(formula = mu_formula, 
                                sigma.formula = sigma_formula, data = dat_use, 
                                family = gamlss.dist::NBI, control = gamlss::gamlss.control(trace = FALSE, 
                                                                                            c.crit = 0.1))
          end.time = Sys.time()
          time = as.numeric(end.time - start.time)
          time_list[2] <- time
          res
        }, error = function(e) {
          add_log("gamlss", "error", toString(e))
          NULL
        }), warning = function(w) {
          add_log("gamlss", "warning", toString(w))
        })
      }
      else {
        gamlss.fit <- NULL
      }
    }
    else if (family_gene == "zip") {
      mgcv.fit <- withCallingHandlers(tryCatch({
        start.time <- Sys.time()
        res <- fitfunc(formula = mgcv_formula, data = dat_use, 
                       family = "poisson", discrete = usebam)
        end.time <- Sys.time()
        time <- as.numeric(end.time - start.time)
        time_list[1] <- time
        res
      }, error = function(e) {
        add_log("gam", "error", toString(e))
        NULL
      }), warning = function(w) {
        add_log("gam", "warning", toString(w))
      })
      gamlss.fit <- withCallingHandlers(tryCatch({
        start.time = Sys.time()
        res <- gamlss::gamlss(formula = mu_formula, 
                              sigma.formula = mu_formula, data = dat_use, 
                              family = gamlss.dist::ZIP, control = gamlss::gamlss.control(trace = FALSE, 
                                                                                          c.crit = 0.1))
        end.time = Sys.time()
        time = as.numeric(end.time - start.time)
        time_list[2] <- time
        res
      }, error = function(e) {
        add_log("gamlss", "error", toString(e))
        NULL
      }), warning = function(w) {
        add_log("gamlss", "warning", toString(w))
      })
    }
    else if (family_gene == "zinb") {
      mgcv.fit <- withCallingHandlers(tryCatch({
        start.time <- Sys.time()
        res <- fitfunc(formula = mgcv_formula, data = dat_use, 
                       family = "nb", discrete = usebam)
        end.time <- Sys.time()
        time <- as.numeric(end.time - start.time)
        time_list[1] <- time
        res
      }, error = function(e) {
        add_log("gam", "error", toString(e))
        NULL
      }), warning = function(w) {
        add_log("gam", "warning", toString(w))
      })
      gamlss.fit <- withCallingHandlers(tryCatch({
        start.time = Sys.time()
        res <- gamlss::gamlss(formula = as.formula(paste(predictor, 
                                                         " ~ Group + Group:Patient + offset(logn_counts)", 
                                                         sep = "")), sigma.formula = ~Group, nu.formula = ~Group + 
                                Group:Patient, data = dat_use, family = gamlss.dist::ZINBI, 
                              control = gamlss::gamlss.control(trace = FALSE, 
                                                               c.crit = 0.1))
        end.time = Sys.time()
        time = as.numeric(end.time - start.time)
        time_list[2] <- time
        res
      }, error = function(e) {
        add_log("gamlss", "error", toString(e))
        NULL
      }), warning = function(w) {
        add_log("gamlss", "warning", toString(w))
      })
    }
    else {
      stop("The regression distribution must be one of gaussian, poisson, nb, zip or zinb!")
    }
    if (!"gamlss" %in% class(gamlss.fit)) {
      if (sigma_formula != "~1") {
        message(paste0(gene, " uses mgcv::gam due to gamlss's error!"))
        if (!is.null(gamlss.fit)) {
          if (is.null(warn)) {
            warn = gamlss.fit
          }
          else {
            warn = c(warn, gamlss.fit)
          }
        }
      }
      fit <- mgcv.fit
    }
    else {
      mean_vec <- stats::predict(gamlss.fit, type = "response", 
                                 what = "mu", data = dat_use)
      theta_vec <- stats::predict(gamlss.fit, type = "response", 
                                  what = "sigma", data = dat_use)
      if_infinite <- (sum(is.infinite(mean_vec + theta_vec)) > 
                        0)
      if_overmax <- (max(mean_vec, na.rm = TRUE) > 10 * 
                       max(dat_use$gene, na.rm = TRUE))
      if (family_gene %in% c("nb", "zinb")) {
        if_overdisp <- (max(theta_vec, na.rm = TRUE) > 
                          1000)
      }
      else {
        if_overdisp <- FALSE
      }
      if (if_infinite | if_overmax | if_overdisp) {
        add_log("fit_marginal", "warning", paste0(gene, 
                                                  " gamlss returns abnormal fitting values!"))
        fit <- mgcv.fit
      }
      else if (stats::AIC(mgcv.fit) - stats::AIC(gamlss.fit) < 
               -Inf) {
        message(paste0(gene, "'s gamlss AIC is not signifincantly smaller than gam!"))
        fit <- mgcv.fit
      }
      else {
        fit <- gamlss.fit
      }
    }
    if (simplify) {
      fit <- simplify_fit(fit)
    }
    if (trace) {
      return(list(fit = fit, warning = logs, time = time_list, 
                  removed_cell = remove_cell))
    }
    return(list(fit = fit, removed_cell = remove_cell))
  }
  paraFunc <- parallel::mcmapply
  if (.Platform$OS.type == "windows") {
    BPPARAM <- BiocParallel::SnowParam()
    parallelization <- "bpmapply"
  }
  if (parallelization == "bpmapply") {
    paraFunc <- BiocParallel::bpmapply
  }
  if (parallelization == "pbmcmapply") {
    paraFunc <- pbmcapply::pbmcmapply
  }
  if (edf_fitting == FALSE) {
    if (parallelization == "bpmapply") {
      if (class(BPPARAM)[1] != "SerialParam") {
        BPPARAM$workers <- n_cores
      }
      model_fit <- suppressMessages(paraFunc(fit_model_func, 
                                             gene = feature_names, family_gene = family_use, 
                                             MoreArgs = list(dat_use = dat_cov, mu_formula = mu_formula, 
                                                             sigma_formula = sigma_formula, predictor = predictor, 
                                                             count_mat = count_mat), SIMPLIFY = FALSE, 
                                             BPPARAM = BPPARAM))
    }
    else {
      model_fit <- suppressMessages(paraFunc(fit_model_func, 
                                             gene = feature_names, family_gene = family_use, 
                                             mc.cores = n_cores, MoreArgs = list(dat_use = dat_cov, 
                                                                                 mu_formula = mu_formula, sigma_formula = sigma_formula, 
                                                                                 predictor = predictor, count_mat = count_mat), 
                                             SIMPLIFY = FALSE))
    }
  }
  else {
    if (parallelization == "bpmapply") {
      if (class(BPPARAM)[1] != "SerialParam") {
        BPPARAM$workers <- n_cores
      }
      model_fit_edf_gini <- suppressMessages(paraFunc(fit_model_func, 
                                                      gene = edf_gini_feature_names, family_gene = edf_gini_family_use, 
                                                      MoreArgs = list(dat_use = dat_cov, mu_formula = mu_formula, 
                                                                      sigma_formula = sigma_formula, predictor = predictor, 
                                                                      count_mat = edf_gini_count_mat), SIMPLIFY = FALSE, 
                                                      BPPARAM = BPPARAM))
    }
    else {
      model_fit_edf_gini <- suppressMessages(paraFunc(fit_model_func, 
                                                      gene = edf_gini_feature_names, family_gene = edf_gini_family_use, 
                                                      mc.cores = n_cores, MoreArgs = list(dat_use = dat_cov, 
                                                                                          mu_formula = mu_formula, sigma_formula = sigma_formula, 
                                                                                          predictor = predictor, count_mat = edf_gini_count_mat), 
                                                      SIMPLIFY = FALSE))
    }
    edf <- rep(NA, length(model_fit_edf_gini))
    for (i in 1:length(model_fit_edf_gini)) {
      res_ind <- model_fit_edf_gini[i]
      if (lengths(res_ind) == 2) {
        res_ind <- res_ind[[names(res_ind)]]
        edf[i] <- sum(res_ind$fit$edf)
      }
    }
    edf_gini_count_gini <- apply(log(edf_gini_count_mat + 
                                       1), MARGIN = 2, FUN = gini)
    edf_gini_df <- data.frame(edf = edf, gini = edf_gini_count_gini)
    lm_edf_gini <- stats::lm(edf ~ gini, data = edf_gini_df)
    edf_flexible_count_gini <- apply(log(edf_flexible_count_mat + 
                                           1), MARGIN = 2, FUN = gini)
    edf_flexible_df <- data.frame(gini = edf_flexible_count_gini)
    edf_flexible_predicted <- stats::predict(lm_edf_gini, 
                                             edf_flexible_df, se.fit = TRUE, interval = "confidence", 
                                             level = 0.95)
    edf_flexible_predicted_upr <- edf_flexible_predicted$fit[, 
                                                             3]
    if (parallelization == "bpmapply") {
      if (class(BPPARAM)[1] != "SerialParam") {
        BPPARAM$workers <- n_cores
      }
      model_fit_edf_flexible <- suppressMessages(paraFunc(fit_model_func, 
                                                          gene = edf_flexible_feature_names, family_gene = edf_flexible_family_use, 
                                                          MoreArgs = list(dat_use = dat_cov, mu_formula = mu_formula, 
                                                                          sigma_formula = sigma_formula, predictor = predictor, 
                                                                          count_mat = edf_flexible_count_mat, edf = edf_flexible_predicted_upr), 
                                                          SIMPLIFY = FALSE, BPPARAM = BPPARAM))
    }
    else {
      model_fit_edf_flexible <- suppressMessages(paraFunc(fit_model_func, 
                                                          gene = edf_flexible_feature_names, family_gene = edf_flexible_family_use, 
                                                          mc.cores = n_cores, MoreArgs = list(dat_use = dat_cov, 
                                                                                              mu_formula = mu_formula, sigma_formula = sigma_formula, 
                                                                                              predictor = predictor, count_mat = edf_flexible_count_mat, 
                                                                                              edf = edf_flexible_predicted_upr), SIMPLIFY = FALSE))
    }
    model_fit <- vector(mode = "list", length = length(feature_names))
    names(model_fit) <- feature_names
    for (index in names(model_fit_edf_gini)) {
      model_fit[[index]] <- model_fit_edf_gini[[index]]
    }
    for (index in names(model_fit_edf_flexible)) {
      model_fit[[index]] <- model_fit_edf_flexible[[index]]
    }
  }
  return(model_fit)
}



tmpfun <- get("fit_marginal", envir = asNamespace("scDesign3"))

environment(fit_marginal) <- environment(tmpfun)

assignInNamespace("fit_marginal", fit_marginal, ns = "scDesign3")


