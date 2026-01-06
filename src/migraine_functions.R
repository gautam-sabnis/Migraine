prep_fboli_data <- function(data){

    fboli <- data |> dplyr::select(c('NetworkFilename', 'MouseID', 'Tx', 'Sex', 0:104 |> as.character()))
    colnames(fboli)[5:109] <- paste0('minuteMax', colnames(fboli)[5:109]) |> as.character()
    fboli$MouseID <- as.factor(fboli$MouseID)
    fboli$Tx <- as.factor(fboli$Tx)
    Txs <- unique(fboli$Tx)
    mIDs <- unique(fboli$MouseID)

    df_equal <- data.frame()
    df_unequal <- data.frame()
    pb <- txtProgressBar(min = 0, max = length(Txs), style = 3)

    for (t in seq(Txs)){

  	    df <- fboli[fboli$Tx %in% Txs[t],]
  	    tmp <- list()
  	    invisible(suppressMessages({(lapply(seq(nrow(df)), function(x) {
            tmp[[x]] <<- reshape2::melt(df[df$NetworkFilename == unique(df$NetworkFilename)[x], which(names(df) %in% c('MouseID','Sex','Tx', sapply(seq(0,104), function(x) paste0("minuteMax",x))))]);
			tmpp <- tmp[[x]]$value;
            while (any(diff(tmpp) < 0)){
                ind <- which(diff(tmpp) < 0)
				for (i in ind){
					tmpp[i + 1] <- tmpp[i]
				}
			}
			tmp[[x]]$value <<- tmpp;}))}))
  		df_melt <- do.call(rbind,tmp)
  		df_melt$value <- as.numeric(df_melt$value)
  		df_melt$variable <- as.integer(df_melt$variable)
  		df_melt$Sex <- as.factor(rep(df$Sex,each = 105)) 
        df_melt$Tx <- as.factor(rep(df$Tx,each = 105))
  		df_equal <- rbind(df_equal, df_melt)

  		tmp <- list()
  		invisible(lapply(seq(length(unique(df_melt$MouseID))), function(x) {
            tmpp <<- df_melt[df_melt$MouseID == unique(df_melt$MouseID)[x],]
            tmp[[x]] <<- tmpp[!duplicated(tmpp$value),]}))
        df_melt2 <- do.call(rbind,tmp)
  		df_melt2$variable <- as.integer(df_melt2$variable)
  		df_unequal <- rbind(df_unequal, df_melt2)

  		setTxtProgressBar(pb, t)
  		close(pb)
        }
    df_equal <- nlme::groupedData(value ~ variable|MouseID, data = df_equal, order.groups=FALSE)
    df_unequal <- nlme::groupedData(value ~ variable|MouseID, data = df_unequal, order.groups=FALSE)
    df_equal <- df_equal[complete.cases(df_equal),]
    df_unequal <- df_unequal[complete.cases(df_unequal),]
    df_unequal <- df_unequal[!duplicated(df_unequal[c("MouseID", "variable")]),]
    names(df_unequal)[names(df_unequal) == "time"] <- "variable"
    df_unequal <- df_unequal |> filter(!variable == 105)
    #df_unequal <- df_unequal |> filter(!variable == last(variable))
    df_unequal <- bind_rows(df_unequal, df_equal |> group_by(MouseID) |> filter(variable == last(variable)))
    
    p <- ggplot(df_unequal, aes(x = variable, y = value, fill = MouseID)) + geom_point(size = 2, alpha = 0.3) + geom_line(alpha = 0.3) + theme_minimal(base_size = 14) + theme(legend.position = "none") + facet_wrap(. ~ Tx, scales = "free") + labs(x = "Time (in minutes)", y = "Count")

    return(list(df_equal = df_equal, df_unequal = df_unequal, plot = p))

}


fit_bayes_growth_model <- function(data, trt, type){
    df_trt <- data[data$Tx %in% trt, ]
    df_trt$MouseID <- droplevels(df_trt$MouseID)
    df_trt$Sex <- as.factor(as.numeric(df_trt$Sex))
    if (type == "Gompertz"){
    	fit <- brm(bf(value ~ theta1*exp(-exp(-theta2*(variable - theta3))), theta1 ~ (1|MouseID), theta2 ~ (1|MouseID), theta3 ~ (1|MouseID), nl = TRUE), prior = c(prior(normal(0,5), class = "b", nlpar = "theta1", lb = 0), prior(normal(0,1), class = "b", nlpar = "theta2", lb = 0), prior(normal(0,7), class = "b", nlpar = "theta3", lb = 0)), data = df_trt, family = gaussian(), iter = 2000, cores = 4, backend = "cmdstanr", threads = threading(2), silent = 2, control = list(max_treedepth = 10), , sample_prior = "yes")
    } else if (type == "Logistic"){
    	fit <- brm(bf(value ~ theta1/(1 + exp(-theta2*(variable - theta3))), theta1 ~ (1|MouseID), theta2 ~ (1|MouseID), theta3 ~ (1|MouseID), nl = TRUE), prior = c(prior(normal(0,5), class = "b", nlpar = "theta1", lb = 0), prior(normal(0,1), class = "b", nlpar = "theta2",lb = 0), prior(normal(0,7), class = "b", nlpar = "theta3", lb = 0)), data = df_trt, family = gaussian(), iter = 2000, cores = 4, backend = "cmdstanr", threads = threading(2), silent = 2, control = list(max_treedepth = 10), , sample_prior = "yes")
    } else {
    	fit <- brm(bf(value ~ theta1*(1 + psi*exp(-theta2*(variable - theta3)))^(-1/psi), theta1 ~ (1|MouseID), theta2 ~ (1|MouseID), theta3 ~ (1|MouseID), psi ~ 1, nl = TRUE), prior = c(prior(normal(0,4), class = "b", nlpar = "theta1"), prior(normal(0,4), class = "b", nlpar = "theta2"), prior(normal(0,4), class = "b", nlpar = "theta3"), prior(normal(0,4), class = "b", nlpar = "psi")), data = df_trt, family = gaussian(), iter = 2000, cores = 4, backend = "cmdstanr", threads = threading(2), refresh = 0, silent = 2)
    }

    return(fit)
}

prep_vignette <- function(modelfit, data){

    fit <- modelfit
    Txs <- unique(data$Tx)

    #Population level
    p <- conditional_effects(fit) |> plot() #fixed effect
    p1 <- p[[1]] + theme_bw(base_size = 14) + labs(x = "Time (in mins)", y = "Fecal Boli Count")

    #Animal Level
    conditions <- data |> dplyr::select(MouseID, Tx) |> distinct() |> as.data.frame() 
    rownames(conditions) <- paste0(conditions$MouseID, "_", conditions$Tx)
    conditions$MouseID <- droplevels(conditions$MouseID)


    p <- plot(conditional_effects(fit, conditions = conditions[sample(10),], re_formula = NULL, method = "predict"), points = T, ncol = 5)
    
    p2 <- p[[1]] + theme_bw(base_size = 14) + labs(x = "Time (in mins)", y = "Fecal Boli Count") + scale_x_continuous(limits = c(0, 60)) + theme(strip.text.x = element_text(size = 7)) 

    return(list(p1 = p1, p2 = p2))

}

extract_thetas <- function(modelfit, data){

    fit <- modelfit
    fe <- fixef(fit)
    re <- ranef(fit)
    rnames <- rownames(re[[1]])
    re <- re[[1]]

    df_fe <- matrix(fixef(fit)[c(1,2,3),1], nrow = 1, ncol = 3) %x% rep(1, length(mIDs)) # %x% is the kronecker product
    df_re <- matrix(re[,1,1:3], nrow = length(mIDs), ncol = 3) #animal level

    dftheta <- df_fe + df_re
    colnames(dftheta) <- c("theta1", "theta2", "theta3")
    dftheta <- as_tibble(dftheta)
    dftheta <- dftheta |> mutate(MouseID = rnames) |> left_join(data |> distinct(MouseID, .keep_all = TRUE) |> dplyr::select(MouseID, Tx, Sex)) #add metadata

    dftheta <- dftheta |> dplyr::select(MouseID, theta1, theta2, theta3)
    return(list(data = dftheta))

}

#naive version of sparse pca and sparse lda for dimensionality reduction and visualization
plot_sparse_dimRed <- function(data, features, remove_timepoint, remove_tx, sex){

    dt <- data |> dplyr::select(all_of(c(meta, features)))
    dt <- dt |> filter(!Timepoint %in% remove_timepoint)
    dt <- dt |> filter(!Tx %in% remove_tx)
    dt <- dt |> filter(Sex %in% sex)

    pca_recipe <- recipe(~ ., data = dt |> dplyr::select(all_of(features))) |> step_nzv(all_predictors()) |> step_normalize(all_predictors()) |> step_pca(all_predictors(), threshold = 0.99) |> prep()
    zero_var_features <- pca_recipe |> tidy(number = 1) |> pull(terms)
    
    #corr_filter <- recipe(~ ., data = dt |> dplyr::select(all_of(c(setdiff(features, zero_var_features))))) |> step_corr(all_predictors(), threshold = 0.9) |> prep()
    #high_corr_features <- corr_filter |> tidy(number = 1) |> pull(terms)

    remove_features <- c(zero_var_features)
    X <- dt |> dplyr::select(setdiff(features, remove_features))   

    pca.results <- spca(X, scale = TRUE, ncomp = ncol(X)-1)
    dfpc <- pca.results$X%*%pca.results$loadings[[1]] |> bind_cols(dt |> dplyr::select(MouseID, Tx, Timepoint, Sex)) |> mutate(Tx = factor(Tx))

    p1 <- dfpc |> ggplot(aes(x = PC1, y = PC2)) + geom_point(size = 3, aes(color = Tx, shape = Sex)) + stat_ellipse(geom = "polygon", aes(fill = Tx), alpha = 0.3) + ggsci::scale_color_nejm() + ggsci::scale_fill_nejm() + labs(x = paste0("PC1 (", round(pca.results$prop_expl_var$X[1], 2)*100, "%)"), y = paste0("PC2 (", round(pca.results$prop_expl_var$X[2], 2)*100, "%)")) + theme_minimal(base_size = 14) + facet_grid(. ~ Timepoint)

    Y <- dt |> pull(Tx) |> as.factor()
    res <- splsda(X, Y, multilevel = dt$mID)
    dff <- res$X%*%res$loadings[[1]] |> bind_cols(dt |> dplyr::select(mID, Tx, Timepoint, Sex)) |> mutate(Tx = factor(Tx))
    
    p2 <- dff |> ggplot(aes(x = comp1, y = comp2)) + geom_point(size = 3, aes(color = Tx, shape = Sex)) + stat_ellipse(geom = "polygon", aes(fill = Tx), alpha = 0.3, level = 0.99) + ggsci::scale_color_nejm() + ggsci::scale_fill_nejm() + labs(x = 'Component 1', y = 'Component 2') + theme_minimal(base_size = 14) + facet_grid(. ~ Timepoint)

    return(list(p1 = p1, p2 = p2, pca.results = pca.results, lda.results = dff))


}

prep_data_univariate <- function(data, sex, remove_tx, remove_timepoint, features){

    dt <- data |> dplyr::select(all_of(c(meta, features)))
    dt <- dt |> filter(!Timepoint %in% remove_timepoint)
    dt <- dt |> filter(!Tx %in% remove_tx)
    dt <- dt |> filter(Sex %in% sex)
    
    return(list(data = dt))

}

glmmLasso_fun <- function(data, features, lambdas){

    dt <- data |> dplyr::select(all_of(c(meta, features))) |> mutate_at('mID', factor) 

    pca_recipe <- recipe(~ ., data = dt |> dplyr::select(all_of(features))) |> step_nzv(all_predictors()) |> step_normalize(all_predictors()) |> step_pca(all_predictors(), threshold = 0.99) |> prep()
    zero_var_features <- pca_recipe |> tidy(number = 1) |> pull(terms)

    remove_features <- c(zero_var_features)
    features <- setdiff(features, remove_features)
    dt <- dt |> dplyr::select(!all_of(remove_features))

    ff <- function(lambda){
        fit <- glmmLasso(fix = (paste(paste("Tx ~ ", paste(features, collapse = " + "), sep = ""))) |> as.formula(), rnd = list(mID=~1), data = dt, family = cumulative(), lambda = lambda, control = list(standardize = TRUE))
        return(fit$bic) 
    }

    res_list <- furrr::future_map(lambdas, ff, .progress = TRUE)
    res <- data.frame(lambda = lambdas, bic = unlist(res_list))
    p <- res |> ggplot(aes(x = lambda, y = bic)) + geom_line() + geom_point() + theme_classic(base_size = 12) + labs(x = latex2exp::TeX("$\\lambda$"), y = "BIC") + scale_x_continuous(breaks = lambdas)
    return(list(res = res, plot = p))

}

fit_glmmLasso <- function(data, features, lambda){
    
    dt <- data |> dplyr::select(all_of(c(meta, features)))

    fit <- glmmLasso(fix = (paste(paste("Tx ~ ", paste(features, collapse = " + "), sep = ""))) |> as.formula(), rnd = list(mID=~1), data = dt, family = cumulative(), lambda = lambda)
    return(list(fit = fit, data = dt, features = features, lambda = lambda))

}

plot_fit <- function(fit, data, features){

    dt <- data |> dplyr::select(all_of(c(meta, features)))

    theta_hat <- fit$coefficients[-c(1:3)]
    tau <- fit$coefficients[c(1:3)]
    df_theta <- data.frame(Features = names(theta_hat), theta_hat = theta_hat |> as.numeric())
    df_theta$Feature <- factor(df_theta$Feature, levels = names(fit$features))
    p1 <- ggplot(df_theta, aes(x = Feature, y = -theta_hat)) + geom_bar(stat = "identity", color = TRUE) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(y = latex2exp::TeX("Coefficients, $\\theta$"))

    X <- dt |> dplyr::select(all_of(features)) |> as.matrix()
    theta_hat <- as.numeric(theta_hat)
    df_univariate <- dt |> dplyr::select(Tx) |> mutate(score = as.numeric(X%*%theta_hat)) |> mutate(Tx = factor(Tx, ordered = TRUE))
    p2 <- ggplot(df_univariate, aes(x = -score, color = Tx)) + geom_density(lwd = 1.1) + ggsci::scale_color_jama(name = "Seizure") + labs(x = "Univariate Migraine Scale", y = "Density") + scale_y_continuous(expand = c(0,0)) + geom_vline(xintercept = tau[1], linetype = "dashed") + geom_vline(xintercept = tau[2], linetype = "dashed") + geom_vline(xintercept = tau[3], linetype = "dashed") + theme_classic(base_size = 16) + guides(color=guide_legend(title="Tx"))

    p3 <- ggplot(df_univariate |> filter(!Tx %in% c(2.5, 5)), aes(x = -score, color = Tx)) + geom_density(lwd = 1.1) + ggsci::scale_color_jama(name = "Seizure") + labs(x = "Univariate Migraine Scale", y = "Density") + theme_classic(base_size = 16) + guides(color=guide_legend(title="Tx"))

    return(list(p1 = p1, p2 = p2, p3 = p3, fit = fit))

}

