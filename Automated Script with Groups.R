# Import Necessary Packages
library(netmeta)
library(readxl)
library(dplyr)
library(stringr)
library(meta)
library(grid)
library(NMA)
library(openxlsx)
library(metafor)
library(tidyr)

# =========================
# Helper function
# =========================
extract_mean_sd <- function(x) {
  x <- str_replace_all(x, ",", ".")
  nums <- str_extract_all(x, "[0-9.]+")[[1]]
  
  if (length(nums) == 1) {
    c(mean = as.numeric(nums[1]), sd = NA)
  } else {
    c(mean = as.numeric(nums[1]), sd = as.numeric(nums[2]))
  }
}

# =========================
# Outcome groups & figures
# =========================

body_comp <- c("BMI", "BMI Z Score", "Weight", "Body Fat", "FM", "FFM")
cardio    <- c("DBP", "SBP", "VO2")
metabolic <- c("Cholestrol", "LDL", "HDL", "TG", "Glucose", "Insulin", "HOMA-IR")

fig_map <- list(
  body = list(forest = 1, funnel = 2, network = 3, transitivity = 4, loo = 5, baujat = 6, inf = 7),
  cardio = list(forest = 8, funnel = 9, network = 10, transitivity = 11, loo = 12, baujat = 13, inf = 14),
  metabolic = list(forest = 15, funnel = 16, network = 17, transitivity = 18, loo = 19, baujat = 20, inf = 21)
)


# =========================
# Main analysis function
# =========================
run_analysis <- function(file_path, results) {
  
  outcome <- tools::file_path_sans_ext(basename(file_path))
  message("Processing: ", outcome)
  
  # -------------------------
  # Read data
  # -------------------------
  data_raw <- read_excel(file_path)
  
  data <- data_raw %>%
    rowwise() %>%
    mutate(
      mean1 = extract_mean_sd(`Treatment 1 Outcome`)[1],
      sd1   = extract_mean_sd(`Treatment 1 Outcome`)[2],
      mean2 = extract_mean_sd(`Treatment 2 Outcome`)[1],
      sd2   = extract_mean_sd(`Treatment 2 Outcome`)[2]
    ) %>%
    ungroup() %>%
    rename(
      Study  = `STUDY ID`,
      n1     = `Treatment 1 sample size`,
      n2     = `Treatment 2 Sample size`,
      Treat1 = `Treatment 1`,
      Treat2 = `Treatment 2`
    ) %>%
    filter(!is.na(sd1), !is.na(sd2),
           !grepl("hiit", Treat2, ignore.case = TRUE))
  
  data$Study <- make.unique(data$Study)
  
  # Effect size
  data$TE   <- data$mean1 - data$mean2
  data$seTE <- sqrt((data$sd1^2 / data$n1) + (data$sd2^2 / data$n2))
  
  # -------------------------
  # Transitivity
  # -------------------------
  hf_arm <- data %>%
    rename(
      study = Study,
      trt1  = Treat1,
      trt2  = Treat2,
      d1    = mean1,
      d2    = mean2,
      hiit  = `HIIT Treatment Type`,
      age   = `Mean Age`
    ) %>%
    pivot_longer(
      cols = c(trt1, trt2, n1, n2, d1, d2),
      names_to = c(".value", "arm"),
      names_pattern = "([a-z]+)([12])"
    )
  
  hf_nma <- setup(
    study   = study,
    trt     = trt,
    d       = d,
    n       = n,
    z       = c(age, hiit),
    measure = "OR",
    ref     = "Control (No Exercise)",
    data    = hf_arm
  )
  
  # -------------------------
  # Network meta-analysis
  # -------------------------
  nma_result <- netmeta(
    TE = TE,
    seTE = seTE,
    treat1 = Treat1,
    treat2 = Treat2,
    studlab = Study,
    sm = "SMD",
    reference.group = "Control (No Exercise)",
    data = data
  )
  
  nettable(
    nma_result,
    digits = 2,
    path = paste0("Network Tables/", outcome, " nettable.xlsx"),
    overwrite = TRUE
  )
  
  # -------------------------
  # Pairwise meta-analysis
  # -------------------------
  data_meta <- data.frame(
    study = data$Study,
    treatment = data$Treat1,
    control = data$Treat2,
    mean_int = data$mean1,
    sd_int = data$sd1,
    n_int = data$n1,
    mean_ctrl = data$mean2,
    sd_ctrl = data$sd2,
    n_ctrl = data$n2
  )
  
  data_meta$subgroup <- ifelse(grepl("no exercise", tolower(data_meta$control)), "Non-Exercise Control",
                               ifelse(grepl("mict", tolower(data_meta$control)), "MICT",
                                      ifelse(grepl("low intensity", tolower(data_meta$control)), "Low Intensity",
                                             ifelse(grepl("nutrition|l-cit|diet|chocolate", tolower(data_meta$control)), "Nutrition",
                                                    ifelse(grepl("supra|sprint|sit", tolower(data_meta$control)), "Sprint Variants",
                                                           ifelse(grepl("miit|sbe|intermitten|football", tolower(data_meta$control)), "Mix Exercise",
                                                                  "Other"))))))
  
  meta_result <- metacont(
    n.e = n_int,
    mean.e = mean_int,
    sd.e = sd_int,
    n.c = n_ctrl,
    mean.c = mean_ctrl,
    sd.c = sd_ctrl,
    studlab = study,
    data = data_meta,
    sm = "SMD",
    method.tau = "REML",
    subgroup = subgroup,
    random = TRUE
  )

  
  # -------------------------
  # Leave-one-out (Influence analysis)
  # -------------------------
  loo_result <- metainf(meta_result)
  loo_sum <- summary(loo_result)
  
  loo_df <- data.frame(
    Study = loo_sum$studlab,
    TE    = loo_sum$TE,
    seTE = loo_sum$seTE,
    lower = loo_sum$lower,
    upper = loo_sum$upper,
    pval = loo_sum$pval,
    tau2 = loo_sum$tau2,
    tau = loo_sum$tau,
    I2 = loo_sum$I2
  )
  
  write.xlsx(
    loo_df,
    file = paste0("Leave-One-Out/", outcome, " leave1out.xlsx"),
    overwrite = TRUE
  )
  
  # -------------------------
  # Baujat plot
  # -------------------------
  baujat_obj <- meta_result
  
  # -------------------------
  # Influence Diagnostic Plots
  # -------------------------
  rma_obj <- rma(
    yi  = meta_result$TE,
    sei = meta_result$seTE,
    method = "REML"
  )
  
  inf_obj <- influence(rma_obj)
  
  # Simpan semua
  return(list(
    outcome = outcome,
    meta   = meta_result,
    nma    = nma_result,
    transitivity = hf_nma,
    loo    = loo_result,
    baujat = meta_result,
    inf    = inf_obj
  ))
}

# Forest plot gabungan
plot_forest_group <- function(outcomes, results, fig_no, title) {
  
  pdf(
    file = paste0("Figures/Figure_", fig_no, "_Forest.pdf"),
    width = 9,
    height = 3 * length(outcomes)  # otomatis tinggi
  )
  
  par(mfrow = c(length(outcomes), 1), mar = c(4, 4, 2, 2))
  
  for (o in outcomes) {
    meta::forest(
      results[[o]]$meta,
      main = o,
      xlab = "Standardized Mean Difference",
      fontsize = 7
    )
  }
  
  dev.off()
}


# Funnel plot gabungan
plot_funnel_group <- function(outcomes, results, fig_no, title) {
  
  png(
    filename = paste0("Figures/Figure_", fig_no, "_Funnel.png"),
    width = 4000, height = 5000, res = 300
  )
  
  par(mfrow = c(ceiling(length(outcomes)/2), 2))
  
  for (o in outcomes) {
    funnel(
      results[[o]]$meta,
      main = o
    )
  }
  dev.off()
}

# Network graph gabungan
plot_network_group <- function(outcomes, results, fig_no, title) {
  
  png(
    filename = paste0("Figures/Figure_", fig_no, "_Network.png"),
    width = 4000, height = 4000, res = 300
  )
  
  par(mfrow = c(ceiling(length(outcomes)/2), 2))
  
  for (o in outcomes) {
    netgraph(
      results[[o]]$nma,
      main = o,
      plastic = FALSE
    )
  }
  dev.off()
}

plot_transitivity_group <- function(outcomes, results, fig_no, title) {
  
  png(
    filename = paste0("Figures/Figure_", fig_no, "_Transitivity.png"),
    width = 4000, height = 5000, res = 300
  )
  
  par(mfrow = c(length(outcomes), 1), mar = c(4, 4, 2, 2))
  
  for (o in outcomes) {
    transitivity(results[[o]]$transitivity,age)
    title(o, line = 1, cex.main = 1)
  }
  
  mtext(title, outer = TRUE, line = -1.5, cex = 1.6, font = 2)
  dev.off()
}

plot_loo_group <- function(outcomes, results, fig_no, title) {
  
  pdf(
    file = paste0("Figures/Figure_", fig_no, "_LeaveOneOut.pdf"),
    width = 9,
    height = 3 * length(outcomes)  # otomatis tinggi
  )
  
  par(mfrow = c(length(outcomes), 1), mar = c(4, 4, 2, 2))
  
  for (o in outcomes) {
    meta::forest(
      results[[o]]$loo,
      xlab     = "SMD after omitting one study",
      main     = o,
      fontsize = 7,
    )}
  dev.off()
}


plot_baujat_group <- function(outcomes, results, fig_no, title) {
  
  png(
    filename = paste0("Figures/Figure_", fig_no, "_Baujat.png"),
    width = 4000, height = 4000, res = 300
  )
  
  par(mfrow = c(ceiling(length(outcomes)/2), 2))
  
  for (o in outcomes) {
    baujat(results[[o]]$baujat)
    title(o, line = 1, cex.main = 1)
  }
  
  mtext(title, outer = TRUE, line = -1.5, cex = 1.6, font = 2)
  dev.off()
}

plot_influence_group <- function(outcomes, results, fig_no, title) {
  
  png(
    filename = paste0("Figures/Figure_", fig_no, "_Influence.png"),
    width = 4000, height = 5000, res = 300
  )
  
  par(mfrow = c(length(outcomes), 1))
  
  for (o in outcomes) {
    plot(results[[o]]$inf)
  }
  
  mtext(title, outer = TRUE, line = -1.5, cex = 1.6, font = 2)
  dev.off()
}


# =========================
# RUN FOR ALL DATASETS
# =========================
results <- list()

files <- list.files("Datasets", pattern = "\\.xlsx$", full.names = TRUE)

for (f in files) {
  res <- run_analysis(f)
  results[[res$outcome]] <- res
}

## Body Composition
plot_forest_group(body_comp, results,
                  fig_map$body$forest,
                  "Forest Plot of Body Composition Outcomes")

plot_funnel_group(body_comp, results,
                  fig_map$body$funnel,
                  "Funnel Plot of Body Composition Outcomes")

plot_network_group(body_comp, results,
                   fig_map$body$network,
                   "Network Graph of Body Composition Outcomes")

plot_transitivity_group(body_comp, results, fig_map$body$transitivity,
                        "")

plot_loo_group(body_comp, results, fig_map$body$loo,
               "Leave-one-out Sensitivity Analysis of Body Composition Outcomes")

plot_baujat_group(body_comp, results, fig_map$body$baujat,
                  "Baujat plots of Body Composition Outcomes")

plot_influence_group(body_comp, results, fig_map$body$inf,
                     "")


## Cardiorespiratory
plot_forest_group(cardio, results,
                  fig_map$cardio$forest,
                  "Forest Plot of Cardiorespiratory Outcomes")

plot_funnel_group(cardio, results,
                  fig_map$cardio$funnel,
                  "Funnel Plot of Cardiorespiratory Outcomes")

plot_network_group(cardio, results,
                   fig_map$cardio$network,
                   "Network Graph of Cardiorespiratory Outcomes")

plot_transitivity_group(cardio, results, fig_map$cardio$transitivity,
                        "")

plot_loo_group(cardio, results, fig_map$cardio$loo,
               "Leave-one-out Sensitivity Analysis of Cardiorespiratory Outcomes")

plot_baujat_group(cardio, results, fig_map$cardio$baujat,
                  "Baujat plots of Cardiorespiratory Outcomes")

plot_influence_group(cardio, results, fig_map$cardio$inf,
                     "")


## Metabolic
plot_forest_group(metabolic, results,
                  fig_map$metabolic$forest,
                  "Forest Plot of Metabolic Outcomes")

plot_funnel_group(metabolic, results,
                  fig_map$metabolic$funnel,
                  "Funnel Plot of Metabolic Outcomes")

plot_network_group(metabolic, results,
                   fig_map$metabolic$network,
                   "Network Graph of Metabolic Outcomes")

plot_transitivity_group(metabolic, results, fig_map$metabolic$transitivity,
                        "")

plot_loo_group(metabolic, results, fig_map$metabolic$loo,
               "Leave-one-out Sensitivity Analysis of Metabolic Outcomes")

plot_baujat_group(metabolic, results, fig_map$metabolic$baujat,
                  "Baujat plots of Metabolic Outcomes")

plot_influence_group(metabolic, results, fig_map$metabolic$inf,
                     "")
