# =========================
# SUCRA ANALYSIS ONLY
# =========================

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(openxlsx)
library(gemtc)
library(rjags)

dir.create("SUCRA", showWarnings = FALSE)
dir.create("SUCRA/sucra plot", showWarnings = FALSE)

# =========================
# Helper function (SAMA)
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

clean_treatment <- function(x) {
  x %>%
    str_replace_all("\\+", "_plus_") %>%   # + → _plus_
    str_replace_all("\\s+", "_") %>%       # spasi → _
    str_replace_all("[^A-Za-z0-9_]", "")   # hapus karakter lain
}

# =========================
# SUCRA FUNCTION
# =========================
run_sucra <- function(file_path) {
  
  outcome <- tools::file_path_sans_ext(basename(file_path))
  message("Processing SUCRA: ", outcome)
  
  # -------------------------
  # Read & prepare data
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
    mutate(
      Treat1 = clean_treatment(Treat1),
      Treat2 = clean_treatment(Treat2)
    ) %>%
    filter(!is.na(sd1), !is.na(sd2))
  
  data$Study <- make.unique(data$Study)
  
  # -------------------------
  # Convert ke arm-based
  # -------------------------
  arm1 <- data %>%
    transmute(
      study = Study,
      treatment = Treat1,
      mean = mean1,
      std.dev = sd1,
      sampleSize = n1
    )
  
  arm2 <- data %>%
    transmute(
      study = Study,
      treatment = Treat2,
      mean = mean2,
      std.dev = sd2,
      sampleSize = n2
    )
  
  gemtc_data <- bind_rows(arm1, arm2)
  
  # -------------------------
  # Build network
  # -------------------------
  network <- mtc.network(data.ab = gemtc_data)
  
  # -------------------------
  # Model (continuous data!)
  # -------------------------
  model <- mtc.model(
    network,
    linearModel = "random",
    likelihood = "normal",
    link = "identity",
    n.chain = 4
  )
  
  # -------------------------
  # Run MCMC
  # -------------------------
  mcmc <- mtc.run(
    model,
    n.adapt = 5000,
    n.iter = 20000,
    thin = 10
  )
  
  # -------------------------
  # Ranking & SUCRA
  # -------------------------
  rp <- rank.probability(mcmc)
  sucra_res <- sucra(rp)
  
  # -------------------------
  # Save result
  # -------------------------
  sucra_df <- data.frame(
    Treatment = names(sucra_res),
    SUCRA = as.numeric(sucra_res)
  )
  
  write.xlsx(
    sucra_df,
    file = paste0("SUCRA/", outcome, " sucra.xlsx"),
    overwrite = TRUE
  )
  
  # -------------------------
  # Plot
  # -------------------------
  pdf(paste0("SUCRA/sucra plot/", outcome, ".pdf"), 10, 8)
  plot(sucra_res)
  dev.off()
}

# =========================
# RUN ALL FILES
# =========================
files <- list.files("Datasets", pattern = "\\.xlsx$", full.names = TRUE)
lapply(files, run_sucra)