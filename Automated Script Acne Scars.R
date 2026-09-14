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
library(gemtc)
library(rjags)

dir.create("SUCRA")
dir.create("SUCRA/sucra plot") 

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
# Main analysis function
# =========================
run_analysis <- function(file_path) {
  
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
    filter(!is.na(sd1), !is.na(sd2))
    
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
      d2    = mean2
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
    measure = "OR",
    ref     = "Subcision only",
    data    = hf_arm
  )
  
  pdf(paste0("TRANSITIVITY GRAPHS/", outcome, ".pdf"), 15, 20)
  transitivity(hf_nma, age)
  dev.off()
  
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
    reference.group = "Subcision only",
    data = data
  )
  
  nettable(
    nma_result,
    digits = 2,
    path = paste0("Network Tables/", outcome, " nettable.xlsx"),
    overwrite = TRUE
  )
  
  pdf(paste0("Network Graphs/net plot ", outcome, ".pdf"), 15, 20)
  netgraph(nma_result, plastic = FALSE)
  dev.off()
  
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
  
  data_meta$subgroup <- ifelse(grepl("subcision|mn|mn only|subcision+mn|mnrf", tolower(data_meta$control)), "Surgical/Mechanical Procedures",
                               ifelse(grepl("co2 laser|fractional co2 laser|co2 laser only|co2+mn|co2+prp|co2+laser", tolower(data_meta$control)), "Laser Therapy",
                                      ifelse(grepl("prp|subcision+prp|mn+prp|co2+prp|subcision+ha", tolower(data_meta$control)), "Biologic Based Treatments",
                                             ifelse(grepl("tca", tolower(data_meta$control)), "Chemical Therapy",
                                                    ifelse(grepl("carboxytherapy", tolower(data_meta$control)), "Gas Therapy",
                                                           "Other")))))
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
  
  pdf(paste0("Forest Plots/forest plot ", outcome, ".pdf"), 15, 20)
  forest(meta_result,
         xlab = paste("Standardized Mean Difference (", outcome, ")", sep = ""),
         subgroup = TRUE,
         fontsize = 7)
  dev.off()
  
  pdf(paste0("Funnel Plots/funnel plot ", outcome, ".pdf"), 15, 20)
  funnel(meta_result, main = paste("Funnel Plot -", outcome))
  dev.off()
  
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
  
  pdf(paste0("Leave-One-Out/leave1out plot ", outcome, ".pdf"), 15, 20)
  forest(
    loo_result,
    xlab = paste("SMD after omitting one study"),
    fontsize = 7
  )
  dev.off()
  
  # -------------------------
  # Baujat plot
  # -------------------------
  pdf(paste0("Baujat Plots/baujat plot ", outcome, ".pdf"), 10, 8)
  baujat(meta_result)
  dev.off()
  
  # -------------------------
  # Influence Diagnostic Plots
  # -------------------------
  rma_obj <- rma(
    yi  = meta_result$TE,
    sei = meta_result$seTE,
    method = "REML"
  )
  
  inf <- influence(rma_obj)
  
  png(
    filename = paste0("Influence Diagnostic/influence_plot_", outcome, ".png"),
    width = 4000,
    height = 3000,
    res = 300
  )
  plot(inf)
  dev.off()
  
  # -------------------------
  # SUCRA (Bayesian - gemtc)
  # -------------------------
  
  library(gemtc)
  
  # Convert ke arm-based untuk gemtc
  gemtc_data <- data %>%
    select(Study, Treat1, Treat2, mean1, mean2, sd1, sd2, n1, n2) %>%
    rename(
      study = Study,
      trt1 = Treat1,
      trt2 = Treat2,
      y1 = mean1,
      y2 = mean2,
      sd1 = sd1,
      sd2 = sd2,
      n1 = n1,
      n2 = n2
    ) %>%
    pivot_longer(
      cols = c(trt1, trt2, y1, y2, sd1, sd2, n1, n2),
      names_to = c(".value", "arm"),
      names_pattern = "([a-z]+)([12])"
    )
  
  # Build network
  network <- mtc.network(data.ab = gemtc_data)
  
  # Model (continuous → normal)
  model <- mtc.model(
    network,
    linearModel = "random",
    likelihood = "normal",
    link = "identity",
    n.chain = 4
  )
  
  # Run MCMC
  mcmc <- mtc.run(
    model,
    n.adapt = 5000,
    n.iter = 20000,   # bisa dinaikkan kalau final
    thin = 10
  )
  
  # Ranking probability
  rp <- rank.probability(mcmc)
  
  # SUCRA
  sucra_res <- sucra(rp)
  
  # Save SUCRA
  sucra_df <- data.frame(
    Treatment = names(sucra_res),
    SUCRA = as.numeric(sucra_res)
  )
  
  write.xlsx(
    sucra_df,
    file = paste0("SUCRA/", outcome, " sucra_gemtc.xlsx"),
    overwrite = TRUE
  )
  
  # Plot SUCRA
  pdf(paste0("SUCRA/sucra plot", outcome, ".pdf"), 10, 8)
  plot(sucra_res)
  dev.off()
}

# =========================
# RUN FOR ALL DATASETS
# =========================
files <- list.files("Datasets", pattern = "\\.xlsx$", full.names = TRUE)
lapply(files, run_analysis)
