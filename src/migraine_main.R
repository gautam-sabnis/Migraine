libs <- c('dplyr', 'ggplot2', 'readr', 'brms', 'mixOmics', 'tidymodels', 'glmmLasso', 'future')
sapply(libs, require, character.only = TRUE)

setwd("/Users/sabnig/Documents/WIP/Misc/Migraine")
source("Scripts/migraine_functions.R")
source("Scripts/utils.R")

data <- read_csv("Data/final_nextflow_dataset.csv")
data <- data |> rename(Sex = 'sex') |> mutate(mID = MouseID) |> mutate(MouseID = paste0(MouseID, '_', Timepoint)) |> relocate(mID, .before = 'MouseID') |> rename_with(~ gsub(" ", "_", .)) 
data <- data |> mutate(MouseID = as.character(MouseID))
mIDs <- unique(data$MouseID)
Txs <- unique(data$Tx)

dfboli <- prep_fboli_data(data = data)
dfboli$plot


#Identify mice with 0 fecal boli counts through out the experiment. For such animals, only the asymptote is defined, but slope and displacement are not defined. We may need to filter these animals out. 
zerofb_mIDs <- dfboli$df_unequal |> filter(variable == last(variable)) |> filter(value == 0) |> pull(MouseID)

fit0 <- fit_bayes_growth_model(data = dfboli$df_unequal, trt = Txs[1:4], type = "Gompertz")
fit1 <- fit_bayes_growth_model(data = dfboli$df_unequal, trt = Txs[1:4], type = "Logistic")
loo(fit0, fit1)

dftheta <- extract_thetas(modelfit = fit1, data = dfboli$df_unequal)
dftheta <- dftheta$data

data <- data |> dplyr::select(!c(0:104 |> as.character())) |> left_join(dftheta, by = "MouseID") |> rename(Asymptote = 'theta1', Slope = 'theta2', Displacement = 'theta3')
gait <- colnames(data)[12:95]
engineered <- colnames(data)[96:112]
jabs <- colnames(data)[113:272]
fboli <- colnames(data)[273:275]
features <- c(gait, engineered, jabs, fboli)

#Naive Sparse PCA and SPLS-DA
res <- plot_sparse_dimRed(data = data, features = c(gait, engineered, jabs, fboli), remove_timepoint = setdiff(unique(data$Timepoint), 'D09'), remove_tx = c(2.5, 5.0), sex = c('M'))
res$p1
res$p2

#Univariate migraine scale
dt <- prep_data_univariate(data = data, features = c(gait, engineered, jabs, fboli), remove_timepoint = c('D00', 'D21'), remove_tx = c(2.5, 5.0), sex = c('M'))$data

plan(multisession, workers = availableWorkers() |> length())
lambdas <- seq(0, 3, by = 1)
tune_lasso <- glmmLasso_fun(data = dt, features = features, lambdas = lambdas)


#Boosted models 
require(mboost)

na_features <- which(apply(data, 2, function(x) sum(is.na(x))) > 0) |> names()
zero_var_features <- which(apply(data[, features], 2, function(x) sd(x)) < 1e-16) |> names()
remove_features <- c(na_features, zero_var_features)
features <- setdiff(features, remove_features) 
data <- data |> dplyr::select(!all_of(remove_features)) |> mutate_at('mID', factor)
data$Tx <- factor(data$Tx, ordered = TRUE)

data <- data.frame(data)
data$Tx <- factor(data$Tx, order = TRUE)
mod <- mboost((paste(paste("Tx ~ brandom(mID) + ", paste(features, collapse = " + "), sep = ""))) |> as.formula(), data = data, family = PropOdds(), baselearner = 'btree', control = boost_control(mstop = 5000)) 

data2 <- data |> filter(Tx %in% c(0, 10) & Timepoint == "D09")
#data2 <- data |> filter(Tx %in% c(0, 10) & Timepoint %in% setdiff(unique(data$Timepoint), c("D00", "D21"))) 
data2$Tx <- droplevels(data2$Tx)
data2$mID <- droplevels(data2$mID)

data2 <- data2 |> filter(Sex == 'M')
data2$Tx <- droplevels(data2$Tx)
data2$mID <- droplevels(data2$mID)

na_features <- which(apply(data2, 2, function(x) sum(is.na(x))) > 0) |> names()
zero_var_features <- which(apply(data2[, features], 2, function(x) sd(x)) < 1e-16) |> names()
remove_features <- c(na_features, zero_var_features)
features <- setdiff(features, remove_features) 
data2 <- data2 |> dplyr::select(!all_of(remove_features))

data2 <- data2 |> mutate(across(all_of(features), ~as.vector(scale(., center = TRUE, scale = TRUE))))

mod <- mboost((paste(paste("Tx ~ ", paste(features, collapse = " + "), sep = ""))) |> as.formula(), data = data2, family = Binomial(), baselearner = 'bols', control = boost_control(mstop = 5000)) 

coef <- lapply(coef(mod), function(x) x[[2]]) |> unlist() 
names(coef) <- gsub('bols|[()&]', '', names(coef))
coef <- coef[-1]
dftmp_coef <- data.frame(Feature = names(coef), Coefficient = coef) 
df_coef <- data.frame(Feature = features)
df_coef <- df_coef |> left_join(dftmp_coef, by = 'Feature')
df_coef$Coefficient[df_coef$Coefficient %in% NA] <- 0
df_coef2 <- df_coef |> filter(abs(Coefficient) > 0.05) 
df_coef2$Feature <- factor(df_coef2$Feature, levels = df_coef2$Feature)
df_coef2 <- df_coef2 |> mutate(Coefficient = round(Coefficient, 2))

ggplot(df_coef2, aes(x = Feature, y = Coefficient)) + geom_bar(stat = "identity", color = 'black') + theme_minimal(base_size = 14) + labs(x = "Feature", y = latex2exp::TeX("Coefficient,  $\\theta$")) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) #+ ggbreak::scale_y_break(c(-7, -93), scales = 1) + ggbreak::scale_y_break(c(7, 11), scales = 1)
ggsave("Results/Feature_coefficients_D09_M.pdf", width = 6.9, height = 5.8)
ggsave("Results/Feature_coefficients_D09.pdf", width = 8.9, height = 5.6)

X <- data2 |> dplyr::select(all_of(features)) |> as.matrix()
col_fun2 <- colorRamp2(c(1,2,3,4), c("#7570B3", '#ffcb69', "#D95F02", "#1B9E77"), transparency = 0)
my_features <- data.frame(Type = ifelse(features %in% gait, 1, ifelse(features %in% engineered, 2, ifelse(features %in% jabs, 3, 4))), Feature = features)
my_features <- my_features |> filter(Feature %in% df_coef$Feature) |> dplyr::select(!Feature)
colnames(my_features) <- ""
ht_cluster <- Heatmap(t(as.matrix(my_features)),cluster_rows = FALSE, cluster_columns = FALSE,col = col_fun2, show_heatmap_legend = FALSE,column_names_gp = gpar(fontsize = 10),row_names_gp = gpar(fontsize = 10),row_names_side = "left")


dftmp <- data.frame(score = X %*% df_coef$Coefficient) |> bind_cols(data2 |> dplyr::select(mID, Tx, Timepoint))
ggplot(dftmp, aes(x = score, color = Tx)) + geom_density(lwd = 1.1) + ggsci::scale_color_jama(name = "Tx") + labs(x = "Composite Score", y = "Density") + theme_classic(base_size = 14) + guides(color=guide_legend(title="Tx")) + expand_limits(y = 0)
ggsave("Results/Composite_score_D09_M.pdf", width = 5.02, height = 4.51)

cv10f <- cv(model.weights(mod), type = "kfold") #type = "bootstrap"
tmp <- cvrisk(mod, folds = cv10f, papply = lapply)

mstop <- mstop(tmp)
mIDs <- unique(data2$mID)

ff <- function(x){
    datat <- data2 |> filter(!mID == mIDs[x]) 
    datat$mID <- droplevels(datat$mID)
    feature_means <- apply(datat[, features], 2, mean)
    feature_sds <- apply(datat[, features], 2, sd)
    datat <- datat |> mutate(across(all_of(features), ~as.vector(scale(., center = TRUE, scale = TRUE))))
    mod <- mboost((paste(paste("Tx ~ ", paste(features, collapse = " + "), sep = ""))) |> as.formula(), data = datat, family = Binomial(), baselearner = 'bols', control = boost_control(mstop = 5000)) 
    datan <- data2 |> filter(mID == mIDs[x]) 
    datan[features] <- sweep(datan[features], 2, feature_means, FUN = '-')
    datan[features] <- sweep(datan[features], 2, feature_sds, FUN = '/')
    yhat <- predict(mod, datan, type = 'class') #type = 'class'
    return(list(mod = mod, yhat = yhat))
}

require(future)
plan(multisession, workers = availableWorkers() |> length())
res_list <- furrr::future_map(mIDs, ff, .progress = TRUE, .options=furrr::furrr_options(seed = TRUE))
dfres <- data.frame(pred = lapply(res_list, function(x) x[['yhat']]) |> unlist())
dfres <- dfres |> bind_cols(data2 |> dplyr::select(mID, Tx))

tab <- dfres |> dplyr::select(pred, Tx) |> table() |> as_tibble() 
plot_confusion_matrix(tab, target_col = "Tx", prediction_col = "pred", counts_col = 'n', place_x_axis_above = FALSE, palette = "Oranges") + theme_minimal(base_size = 14) + labs(x = 'True Genotype', y = 'Predicted Genotype')
dev.print(pdf, "Results/conf_mat_D09.pdf", width = 5.02, height = 4.51)


ev <- evaluate(data = dfres |> mutate_at(c('Tx', 'pred'), as.character), target_col = "Tx", prediction_cols = "pred", type = "binomial", id_method = "mean")

#An alternative: glmnet
data2 <- data |> filter(Tx %in% c(0, 10) & Timepoint == "D09")
#data2 <- data |> filter(Tx %in% c(0, 10) & Timepoint %in% setdiff(unique(data$Timepoint), c("D00", "D21"))) 
data2$Tx <- droplevels(data2$Tx)
data2$mID <- droplevels(data2$mID)

data2 <- data2 |> filter(Sex == 'M')
data2$Tx <- droplevels(data2$Tx)
data2$mID <- droplevels(data2$mID)

na_features <- which(apply(data2, 2, function(x) sum(is.na(x))) > 0) |> names()
zero_var_features <- which(apply(data2[, features], 2, function(x) sd(x)) < 1e-16) |> names()
remove_features <- c(na_features, zero_var_features)
features <- setdiff(features, remove_features) 
data2 <- data2 |> dplyr::select(!all_of(remove_features))

#data2 <- data2 |> mutate(across(all_of(features), ~as.vector(scale(., center = TRUE, scale = TRUE))))

X <- data2 |> dplyr::select(all_of(features)) |> mutate(across(all_of(features), ~as.vector(scale(., center = TRUE, scale = TRUE)))) |> as.matrix()
y <- data2$Tx 
mIDs <- unique(data2$mID)

fit <- glmnet::cv.glmnet(x = X, y = y, family = "binomial", nfolds = 10, intercept = FALSE)
df_coef <- coef(fit, s = "lambda.min") |> as.matrix() |> as.data.frame() |> rownames_to_column(var = "Feature") |> rename(Coefficient = 's0') |> mutate(Coefficient = round(Coefficient, 2)) |> filter(abs(Coefficient) > 0.05) 
df_coef$Feature <- factor(df_coef$Feature, levels = df_coef$Feature)

ggplot(df_coef, aes(x = Feature, y = Coefficient)) + geom_bar(stat = "identity", color = 'black') + theme_minimal(base_size = 14) + labs(x = "Feature", y = latex2exp::TeX("Coefficient,  $\\theta$")) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("Results/Feature_coefficients_D09.pdf", width = 9.6, height = 5.4)

ff <- function(x){
    datat <- data2 |> filter(!mID == mIDs[x]) 
    datat$mID <- droplevels(datat$mID)
    feature_means <- apply(datat[, features], 2, mean)
    feature_sds <- apply(datat[, features], 2, sd)
    datat <- datat |> mutate(across(all_of(features), ~as.vector(scale(., center = TRUE, scale = TRUE))))
    Xtr <- datat |> dplyr::select(all_of(features)) |> as.matrix()
    ytr <- datat$Tx 
    fit <- glmnet::cv.glmnet(x = Xtr, y = ytr, family = "binomial", nfolds = 10, intercept = FALSE)
    datan <- data2 |> filter(mID == mIDs[x]) 
    datan[features] <- sweep(datan[features], 2, feature_means, FUN = '-')
    datan[features] <- sweep(datan[features], 2, feature_sds, FUN = '/')
    Xte <- datan |> dplyr::select(all_of(features)) |> as.matrix()
    yhat <- predict(fit, newx = Xte, s = "lambda.min", type = 'class') #type = 'class'
    return(list(fit = fit, yhat = yhat))
}

require(future)
plan(multisession, workers = availableWorkers() |> length())
res_list <- furrr::future_map(mIDs, ff, .progress = TRUE, .options=furrr::furrr_options(seed = TRUE))
dfres <- data.frame(pred = lapply(res_list, function(x) x[['yhat']]) |> unlist())
dfres <- dfres |> bind_cols(data2 |> dplyr::select(mID, Tx))

tab <- dfres |> dplyr::select(pred, Tx) |> table() |> as_tibble() 
plot_confusion_matrix(tab, target_col = "Tx", prediction_col = "pred", counts_col = 'n', place_x_axis_above = FALSE, palette = "Oranges") + theme_minimal(base_size = 14) + labs(x = 'True Genotype', y = 'Predicted Genotype')
dev.print(pdf, "Results/conf_mat_D09.pdf", width = 5.02, height = 4.51)


ev <- evaluate(data = dfres |> mutate_at(c('Tx', 'pred'), as.character), target_col = "Tx", prediction_cols = "pred", type = "binomial", id_method = "mean")