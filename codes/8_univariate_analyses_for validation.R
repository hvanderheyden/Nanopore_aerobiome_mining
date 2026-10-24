#######################################################
############ Diversity ################################
#######################################################

library("phyloseq")
library("pheatmap")
library("ggplot2")
library("tidyverse")
library("vegan")
library("ggpubr")

Runs_PRF <- readRDS("R_objects/RunsR_R.rds");Runs_PRF

################################################################################
################################################################################

Runs_PRF_filtered <- prune_samples(sample_sums(Runs_PRF) > 0, Runs_PRF)

# Set seed for reproducible results
set.seed(1978)
# Calculate bacterial alpha diversity metrics (Shannon, Chao1)
Runs_PRF_a_div <- estimate_richness(Runs_PRF_filtered, measures = c("Shannon","Observed", "Simpson"))

head(Runs_PRF_a_div)

# Extract sample data from phyloseq object to use for analysis
Runs_metadata <- as(sample_data(Runs_PRF_filtered), "data.frame")

# Create new metadata columns for Chao1 and Shannon
Runs_metadata$Shannon <- Runs_PRF_a_div$Shannon
Runs_metadata$Observed <- Runs_PRF_a_div$Observed
Runs_metadata$Simpson <- Runs_PRF_a_div$Simpson

head(Runs_metadata)


###################### Shannon diversity #########################

Runs_metadata %>%
  group_by(BBCH_grouped) %>%
  summarise(
    mean = mean(Shannon),
    median = median(Shannon),
    q1 = quantile(Shannon, 0.25),
    q3 = quantile(Shannon, 0.75),
    min = min(Shannon),
    max = max(Shannon),
    .groups = "drop"
  )

# Test for normality and perform Wilcox test
shapiro.test(Runs_metadata$Shannon)
Runs_wilcox_results_Shannon <- pairwise.wilcox.test(Runs_metadata$Shannon, 
                                                    Runs_metadata$BBCH_grouped, 
                                                    p.adjust.method = "BH")

# Convert results to a data frame
Runs_pvalues_shannon <- Runs_wilcox_results_Shannon$p.value
Runs_sig_data_Shannon <- data.frame(
  group1 = rep(rownames(Runs_pvalues_shannon), times = ncol(Runs_pvalues_shannon)),
  group2 = rep(colnames(Runs_pvalues_shannon), each = nrow(Runs_pvalues_shannon)),
  p_value = as.vector(Runs_pvalues_shannon)
)

Runs_sig_data_Shannon <- na.omit(Runs_sig_data_Shannon)  # Remove NA values

# Add significance labels
Runs_sig_data_Shannon$label <- ifelse(is.na(Runs_sig_data_Shannon$p_value), "ns",
                                      ifelse(Runs_sig_data_Shannon$p_value < 0.001, "***",
                                             ifelse(Runs_sig_data_Shannon$p_value < 0.01, "**",
                                                    ifelse(Runs_sig_data_Shannon$p_value < 0.05, "*", "ns"))))

# Add a y.position for plotting the significance
Runs_sig_data_Shannon$y.position <- c(6, 6.5, 7)  # Customize based on your data

# Plot Shannon index by BBCH for Ascodiomycota

Runs_shannon <- ggplot(Runs_metadata, aes(x = BBCH_grouped, y = Shannon)) +
  geom_boxplot(aes(fill = BBCH_grouped), alpha = 0.6, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), size = 0.2) +
  labs(x = "BBCH stage", y = "Shannon diversity") +
  scale_fill_manual(values = c("cyan4", "grey", "aquamarine2", "#90353BFF")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  stat_pvalue_manual(
    data = Runs_sig_data_Shannon,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = Runs_sig_data_Shannon$y.position,
    tip.length = 0.01,
    size = 2.5
  ); Runs_shannon


###################### Observed Diversity ########################

as.data.frame (Runs_metadata %>%
  group_by(BBCH_grouped) %>%
  summarise(
    mean = signif(mean(Observed),5),
    median = median(Observed),
    q1 = quantile(Observed, 0.25),
    q3 = quantile(Observed, 0.75),
    min = min(Observed),
    max = max(Observed),
    .groups = "drop"
  ))

shapiro.test(Runs_metadata$Observed)

Runs_wilcox_results_Observed <- pairwise.wilcox.test(Runs_metadata$Observed, 
                                                     Runs_metadata$BBCH_grouped, 
                                                     p.adjust.method = "BH")

# Convert results to a data frame
Runs_pvalues_Observed <- Runs_wilcox_results_Observed$p.value
Runs_sig_data_Observed <- data.frame(
  group1 = rep(rownames(Runs_pvalues_Observed), times = ncol(Runs_pvalues_Observed)),
  group2 = rep(colnames(Runs_pvalues_Observed), each = nrow(Runs_pvalues_Observed)),
  p_value = as.vector(Runs_pvalues_Observed)
)

Runs_sig_data_Observed <- na.omit(Runs_sig_data_Observed)  # Remove NA values

# Add significance labels
Runs_sig_data_Observed$label <- ifelse(is.na(Runs_sig_data_Observed$p_value), "ns",
                                       ifelse(Runs_sig_data_Observed$p_value < 0.001, "***",
                                              ifelse(Runs_sig_data_Observed$p_value < 0.01, "**",
                                                     ifelse(Runs_sig_data_Observed$p_value < 0.05, "*", "ns"))))

# Add a y.position for plotting the significance
Runs_sig_data_Observed$y.position <- c(1150, 1250, 1350)  # Customize based on your data

# Plot Shannon index by BBCH for Ascodiomycota

# Plot Observed diversity by BBCH 
Runs_observed <- ggplot(Runs_metadata, aes(x = BBCH_grouped, y = Observed)) +
  geom_boxplot(aes(fill = BBCH_grouped), alpha = 0.6, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), size = 0.2) +
  labs(x = "BBCH stage", y = "Observed richness") +
  scale_fill_manual(values = c("cyan4", "grey", "aquamarine2", "#90353BFF")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  ylim(0, 1350)+
  stat_pvalue_manual(
    data = Runs_sig_data_Observed,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = Runs_sig_data_Observed$y.position,
    tip.length = 0.01,
    size = 2.5
  ); Runs_observed

##################################################################
###################### PCoA ######################################
##################################################################

Runs_PRF_filtered_otu_matrix <- t(otu_table(Runs_PRF_filtered))
dim(Runs_PRF_filtered_otu_matrix)

Runs_PRF_filtered_dist_matrix <- vegdist(Runs_PRF_filtered_otu_matrix, method = "hellinger")

# Perform PCoA
Runs_pcoa_result <- capscale(Runs_PRF_filtered_dist_matrix ~ 1)
Runs_eigenvalues <- Runs_pcoa_result$CA$eig
Runs_prop_explained <- Runs_eigenvalues / sum(Runs_eigenvalues);Runs_prop_explained

# Extract PCoA scores
Runs_pcoa_scores <- scores(Runs_pcoa_result)$sites

row.names(Runs_metadata) == row.names(scores(Runs_pcoa_result)$sites)

Runs_metadata$PCoA1 <- scores(Runs_pcoa_result)$sites[,1]      
Runs_metadata$PCoA2 <- scores(Runs_pcoa_result)$sites[,2]


# Plot  beta diversity by sample time (Figure 1B)
Runs_PCoA <- Runs_metadata %>% 
  ggplot(aes(PCoA1, PCoA2, shape = BBCH_grouped, 
             fill = BBCH_grouped, 
             color = BBCH_grouped))+
  geom_point(size=2)+
  stat_ellipse(geom = "polygon", alpha = 1/6, linetype = 2, size = 0.3, 
               aes(fill = BBCH_grouped, color = BBCH_grouped)) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  theme_bw() + 
  xlab("PCo1 (32.8%)")+
  ylab("PCo2 (10.9%)")+
  theme(legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 9),
        legend.position = "bottom",
        axis.title = element_text(face = "bold", size = 9, color = "black"),
        axis.text = element_text(size = 9)) +
  scale_color_manual(values = c("black","black", "black", "black")) + 
  scale_fill_manual(values = c("cyan4",'grey',"aquamarine2", "#90353BFF")) + 
  guides(size = FALSE);Runs_PCoA


####################################################################################################################################

diversity <- ggarrange(Runs_shannon, 
                       Runs_observed,
                       Runs_PCoA,
                       ncol = 3,
                       # nrow = 3,
                       labels = c("A", "B", "E"),
                       hjust = -1,
                       legend="bottom",
                       common.legend = TRUE,
                       align="v");diversity


################################################################################



################################################################################
################ Permanova #####################################################

set.seed(1432)
permanova_Runs <- adonis2(Runs_PRF_filtered_dist_matrix ~ (BBCH_grouped + DOY + Site)^2 + 
                            Mean_temp_5 + Wind_direction + Relative_humidity + Cum_rain_5 + 
                            Wind_speed + Solar_radiation,
                          data = Runs_metadata, strata=Runs_metadata$YEAR,
                          permutations = 9999); permanova_Runs




################################################################################
########################## Linear mixed effect model ###########################
################################################################################

# Main function to run LMM scenarios and compare models
run_lmm_scenarios_v2 <- function(
    metadata, 
    dataset_name = "Dataset", 
    response_var = "Observed",
    log_transform = TRUE,
    bbch_var = "BBCH_grouped", 
    time_var = "DOY",          
    meteo_vars = c("Mean_temp_5", "Wind_direction", 
                   "Relative_humidity", "Cum_rain_5", 
                   "Wind_speed", "Solar_radiation"),
    by_year = TRUE,
    include_year = FALSE,
    save_results = FALSE,
    model_type = "lmer",          
    glmer_family = NULL           
) {
  # Load required packages
  required_packages <- c("AICcmodavg", "ggeffects", "ggplot2", "car", "dplyr", "lme4", "MuMIn")
  sapply(required_packages, require, character.only = TRUE)
  
  # Step 1: Generate model formulas for different scenarios
  get_model_formulas <- function(include_year) {
    year_term <- if(include_year) " + YEAR" else ""
    resp <- if (log_transform) paste0("log(", response_var, " + 1)") else response_var
    list_formulas <- list(
      Null = as.formula(paste0(resp, " ~ 1", year_term, " + (1|Site)")),
      Tempo = as.formula(paste0(resp, " ~ ",  year_term, time_var, " + (1|Site)")),
      growth = as.formula(paste0(resp, " ~ ",  year_term, bbch_var," + (1|Site)")),
      growth_tempo = as.formula(paste0(resp, " ~ ", bbch_var, " + ", time_var, year_term, " + (1|Site)"))
    )
    meteo_grid <- unlist(lapply(1:length(meteo_vars), function(n) combn(meteo_vars, n, simplify=FALSE)), recursive=FALSE)
    for (i in seq_along(meteo_grid)) {
      v <- meteo_grid[[i]]
      name <- paste("Meteo", paste(sort(v), collapse="_"), sep="_")
      meteo_terms <- paste(v, collapse = " + ")
      fml_txt <- paste0(resp, " ~ ", meteo_terms, year_term, " + (1|Site)")
      list_formulas[[name]] <- as.formula(fml_txt)
    }
    return(list_formulas)
  }
  
  # Step 2: Fit a model (lmer or glmer) to the data
  fit_model <- function(formula, data) {
    if(model_type == "lmer") {
      lme4::lmer(formula, data = data, REML = FALSE)
    } else if(model_type == "glmer") {
      if (is.null(glmer_family)) stop("Please specify 'glmer_family' for glmer. Example: Gamma(link=\"log\").")
      lme4::glmer(formula, data = data, family = glmer_family)
    }
  }
  
  # Step 3: Fit all models and compare using AICc, R2, RMSE, etc.
  fit_and_compare_models <- function(data, include_year) {
    formulas <- get_model_formulas(include_year)
    models <- list()
    stats <- list()
    for (name in names(formulas)) {
      fml <- formulas[[name]]
      mod <- tryCatch(fit_model(fml, data), error = function(e) NULL)
      if (!is.null(mod)) {
        models[[name]] <- mod
        aicc <- tryCatch(AICcmodavg::AICc(mod), error=function(e) NA_real_)
        r2 <- tryCatch(MuMIn::r.squaredGLMM(mod), error=function(e) c(NA_real_, NA_real_))
        loglik <- tryCatch(logLik(mod), error=function(e) NA_real_)
        pred <- tryCatch(predict(mod), error=function(e) NA_real_)
        obs <- tryCatch(model.response(model.frame(mod)), error=function(e) NA_real_)
        rmse <- if (!any(is.na(pred)) && !any(is.na(obs))) sqrt(mean((obs - pred)^2, na.rm=TRUE)) else NA_real_
        var_exp <- r2[1] * 100
        stats[[name]] <- data.frame(
          Model = name,
          AICc = aicc,
          R2_marginal = r2[1],
          R2_conditional = r2[2],
          Variance_explained_pct = var_exp,
          LogLikelihood = as.numeric(loglik),
          RMSE = rmse,
          row.names = NULL
        )
      }
    }
    stats_tab <- do.call(rbind, stats)
    meteo_idx <- grep("^Meteo_", stats_tab$Model)
    get_best_model <- function(mod_name) {
      if(is.na(mod_name)) return(NA)
      idx <- which(stats_tab$Model == mod_name)
      if(length(idx) == 1) return(stats_tab$Model[idx]) else return(NA)
    }
    best_meteo <- NA
    if(length(meteo_idx)) {
      min_aicc <- which.min(stats_tab$AICc[meteo_idx])
      if(length(min_aicc) == 1) best_meteo <- stats_tab$Model[meteo_idx][min_aicc]
    }
    best <- list(
      Null = get_best_model("Null"),
      Temporel = get_best_model("Tempo"),
      Growth = get_best_model("growth"),
      Growth_temporal = get_best_model("growth_tempo"),
      Meteo = best_meteo
    )
    return(list(aic_tab=stats_tab, models=models, best=best))
  }
  
  # Step 4: Predict model output along a given axis (variable)
  predict_on_axis <- function(model, axis_name, data_orig = NULL, group_with = NULL, year_val = NULL) {
    require(ggeffects)
    mod <- model
    if(!is.null(data_orig) && axis_name %in% colnames(data_orig)) {
      xseq <- seq(
        min(data_orig[[axis_name]], na.rm=TRUE),
        max(data_orig[[axis_name]], na.rm=TRUE),
        length.out=200
      )
      pred <- tryCatch(
        as.data.frame(ggeffects::ggpredict(mod, terms = setNames(list(xseq), axis_name))),
        error = function(e) NULL
      )
    } else {
      pred <- tryCatch(
        as.data.frame(ggeffects::ggpredict(mod, terms = paste0(axis_name, " [all]"))),
        error = function(e) NULL
      )
    }
    if(!is.null(pred) & !is.null(year_val)) pred$YEAR <- year_val
    pred
  }
  
  # Step 5: Loop over years (or globally), fit models, store results and predictions
  results <- list()
  pred_combined_temp <- NULL
  years <- if(by_year) sort(unique(metadata$YEAR)) else "global"
  
  obs_var <- if (log_transform) paste0("log_", response_var) else response_var
  response_label <- if (log_transform) paste0("log(", response_var, " + 1)") else response_var
  
  if (log_transform && !(obs_var %in% colnames(metadata))) {
    metadata[[obs_var]] <- log(metadata[[response_var]] + 1)
  }
  
  for (yr in years) {
    data_y <- if(by_year) dplyr::filter(metadata, YEAR==yr) else metadata
    if (log_transform && !(obs_var %in% colnames(data_y))) {
      data_y[[obs_var]] <- log(data_y[[response_var]] + 1)
    }
    comp <- fit_and_compare_models(data_y, include_year = (!by_year && include_year))
    scenario_preds <- list()
    obs_y <- data_y
    
    # Step 6: For each scenario, get predictions from best model
    for (sc in names(comp$best)) {
      modname <- comp$best[[sc]]
      if (!is.null(modname) && !is.na(modname) && (modname %in% names(comp$models))) {
        bestmod <- comp$models[[modname]]
        axis_x <- switch(sc,
                         Null             = NULL,
                         Temporel         = time_var,
                         Spatio_temporel  = time_var,
                         Meteo            = { 
                           fixeds <- attr(terms(bestmod), "term.labels")
                           meteo_found <- intersect(fixeds, meteo_vars)
                           if(length(meteo_found)) meteo_found[1] else NULL
                         }
        )
        pred <- if (!is.null(axis_x)) predict_on_axis(bestmod, axis_x, obs_y, group_with=NULL, year_val=yr) else NULL
        scenario_preds[[sc]] <- pred
        if (!(sc %in% c("Null", "Meteo")) && !is.null(pred)) {
          pred_combined_temp <- rbind(pred_combined_temp, pred)
        }
      }
    }
    results[[as.character(yr)]] <- list(models=comp$models, aic=comp$aic_tab, preds=scenario_preds, best=comp$best, observed=obs_y)
    if (save_results && !is.null(comp$aic_tab)) write.csv(comp$aic_tab, paste0(dataset_name, "_AIC_table_", yr, ".csv"), row.names = FALSE)
  }
  
  # Step 7: Collect selection stats for best models
  selection_stats <- lapply(results, function(res) {
    x <- res$aic
    do.call("rbind", lapply(res$best, function(modname){
      if(is.null(modname) || is.na(modname)) return(NULL)
      x[x$Model==modname,]
    }))
  })
  names(selection_stats) <- names(results)
  
  # Step 8: Plot observed vs predicted values if available
  if(!is.null(pred_combined_temp)) {
    obs_all <- do.call(rbind, lapply(results, function(res) { res$observed }))
    
    y_obs <- if (response_var == "Observed" && any(grepl("^log\\(", names(obs_all)))) {
      exp(obs_all[[paste0("log(Observed + 1)")]]) - 1
    } else if (exists("log_transform") && log_transform) {
      exp(obs_all[[paste0("log_", response_var)]]) - 1
    } else {
      obs_all[[response_var]]
    }
    
    obs_all$y_obs_bt <- y_obs
    
    # Plot
    p <- ggplot() + 
      geom_point(data=obs_all, aes_string(x=time_var, y="y_obs_bt", color="factor(YEAR)"), alpha=0.5) +
      geom_line(data=pred_combined_temp, aes(x=x, y=predicted, color=factor(YEAR)), size=1) +
      geom_ribbon(
        data = pred_combined_temp,
        aes(x = x, ymin = conf.low, ymax = conf.high, fill = factor(YEAR), group = YEAR),
        alpha = 0.30, inherit.aes = FALSE
      ) +
      scale_fill_manual(values = c("#354823", "#273046")) + 
      scale_colour_manual(values = c("#354823", "#273046")) + 
      labs(x = time_var, y = response_var, color = "Year", fill = "Year") +
      theme_classic()
    print(p)
  }
  
  # Step 9: Return results, selection stats, predictions, and plot
  out <- list(
    results_by_block = results,
    selection_stats  = selection_stats,
    scenario_preds   = pred_combined_temp,
    plot             = p
  )
  invisible(out)
}

################################################################################

library("phyloseq")
Runs_PRF <- readRDS("R_objects/RunsR_R.rds");Runs_PRF

# Set seed for reproducible results
set.seed(1978)
# Calculate bacterial alpha diversity metrics (Shannon, Chao1)
Runs_PRF_a_div <- estimate_richness(Runs_PRF, measures = c("Shannon","Observed", "Simpson"))

head(Runs_PRF_a_div)

# Extract sample data from phyloseq object to use for analysis
Runs_PRF_metadata <- as(sample_data(Runs_PRF), "data.frame")

# Create new metadata columns for Chao1 and Shannon
Runs_PRF_metadata$Shannon <- Runs_PRF_a_div$Shannon
Runs_PRF_metadata$Observed <- Runs_PRF_a_div$Observed
Runs_PRF_metadata$Simpson <- Runs_PRF_a_div$Simpson

head(Runs_PRF_metadata)

Runs_PRF_metadata$Week <- as.numeric(Runs_PRF_metadata$Week)

################################################################################
# Fitting linear mixed effect models 

res_lmm_obs <- run_lmm_scenarios_v2(
  Runs_PRF_metadata,
  response_var = "Observed",
  log_transform = TRUE ,
  bbch_var = "BBCH_grouped",
  time_var = "DOY",
  meteo_vars = c("Mean_temp_5", "Wind_direction", 
                 "Relative_humidity", "Cum_rain_5", 
                 "Wind_speed", "Solar_radiation"), 
  by_year = TRUE,
  include_year = FALSE, 
  save_results = TRUE,
  model_type = "lmer"
)

res_lmm_obs$plot
res_lmm_obs$selection_stats

res_lmm_shan <- run_lmm_scenarios_v2(
  Runs_PRF_metadata,
  response_var = "Shannon",
  log_transform = TRUE ,
  bbch_var = "BBCH_grouped",
  time_var = "DOY",
  meteo_vars = c("Mean_temp_5", "Wind_direction", 
                 "Relative_humidity", "Cum_rain_5", 
                 "Wind_speed", "Solar_radiation"), 
  by_year = TRUE,
  include_year = FALSE, 
  save_results = TRUE,
  model_type = "lmer"
)

res_lmm_shan$plot
res_lmm_shan$selection_stats


library(ggpubr)
lme_graph <- ggarrange(res_lmm_shan$plot, 
          res_lmm_obs$plot,
          ncol = 2,
          nrow = 1,
          labels = c("C", "D"),
          legend="bottom",
          common.legend = TRUE,
          align="v",
          heights = c(0.7, 0.7));lme_graph


################# final graph ##################

ggarrange(diversity, 
          lme_graph,
          ncol = 1,
          nrow = 2)


################################################################################

ggsave(file="Figures/Fig5_Diversity.pdf", width=8, 
       height=7, units="in", dpi=900)


