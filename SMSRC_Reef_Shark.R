
#Necessary to install all these packages and load them to R
library(ellipse)
library(Matrix)
library(lme4)
library(emmeans)
library(tidyverse)
library(lmerTest)
library(reshape2)
library(viridis)
library(ggplot2)
library(nicheROVER)
library(SIBER)
library(dplyr)
library(simmr)
library(emmeans)
library(nlme)
library(rjags)
library(tidyr)
library(RColorBrewer)
library(readxl)
library(ggridges) 
library(gridExtra) 
library(coda)
library(MixSIAR)


# set your working directory (where your data file is)
setwd("C:/Users/F00076/Desktop/Megafauna/Reef Shark Mark publication")

#upload the dataset
data <- read_excel("Stable_isotopes_MM.xlsx")


#Shark Trophic position####

# Constants
delta_N_trophic_fractionation <- 2.3   # from Hussey et al. 2010
TP_baseline <- 2                      # trophic position of primary consumers

# 1. Calculate δ15N baseline from C. sordidus at Rowleys
delta15N_b <- data %>%
  filter(Species == "C. sordidus", Reef == "Rowleys") %>%
  summarise(mean_baseline_d15N = mean(d15N, na.rm = TRUE)) %>%
  pull(mean_baseline_d15N)

# 2. Filter sharks and calculate individual TP
tp_sharks <- data %>%
  filter(Species == "C. amblyrhynchos") %>%
  mutate(
    TP = TP_baseline + (d15N - delta15N_b) / delta_N_trophic_fractionation
  )

# 3. Summary statistics
tp_summary <- tp_sharks %>%
  summarise(
    mean_TP = mean(TP, na.rm = TRUE),
    sd_TP   = sd(TP, na.rm = TRUE),
    min_TP  = min(TP, na.rm = TRUE),
    max_TP  = max(TP, na.rm = TRUE),
    n       = n() )

# 4. Print outputs
   
# Individual values
tp_sharks %>%
  select("Sample ID", TP) %>%
  print()

print(tp_summary)   # Summary statistics

write.csv(tp_sharks, "TP_sharks.csv", row.names = FALSE)



# Basic bar plot: Trophic Position per individual shark
# Ensure it's a factor (to get individual bars)
tp_sharks$`Sample ID` <- factor(tp_sharks$`Sample ID`, levels = unique(tp_sharks$`Sample ID`))

# Plot
# Define correct order, skipping sh8
shark_order <- c("sh1", "sh2", "sh3", "sh4", "sh5", "sh6", "sh7", "sh9", "sh10")

# Apply ordered factor to Sample ID
tp_sharks$`Sample ID` <- factor(tp_sharks$`Sample ID`, levels = shark_order)

# Plot
Shark_TP_plot<-ggplot(tp_sharks, aes(x = `Sample ID`, y = TP)) +
  geom_col(fill = "lightblue") +
  theme_minimal() +
  labs(
    title = "Trophic Position of C. amblyrhynchos Individuals at the Rowley Shoals",
    x = "Sample ID",
    y = "Trophic Position (TP)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 14)
  )

Shark_TP_plot

tiff("Shark_TP_plot.tiff", units="cm", width=20, height=10, res=300)
Shark_TP_plot
dev.off()

# %CSD Shark Diet ####
# 1) Settings
predator <- "C. amblyrhynchos"

# The 9 sources used in the original figure (order matters for plotting):
preys <- c("C. microrhinos", 
           "C. sordidus", 
           "L. kasmira",
           "Z. scopas",
           "M. grandoculis",
           "P. vaiuli",
           "S. vulpinus",
           "F. flavissimus",
           "L. bohar",
           "L. decussatus", 
           "L. gibbous")

# TEF functions (Caut 2009; Cardona 2012)
tef_c13 <- function(x) -0.213 * x - 2.848
tef_n15 <- function(x) -0.261 * x + 4.895

# Element-specific TEF SDs (roomier than 0.5/0.5; adjust if you know the exact values used)
TEF_SD_C <- 0.8
TEF_SD_N <- 1.0

# MCMC controls (more robust than defaults)
mcmc_ctrl <- list(n.chain = 4, burn = 10000, iter = 50000, thin = 40)


# 2) Load & harmonize data
df <- readxl::read_excel("Stable_isotopes_MM.xlsx") %>%
  rename(
    Species = tidyselect::any_of("Species"),
    Reef    = tidyselect::any_of("Reef"),
    d13C    = tidyselect::any_of(c("d13C","delta13C","C13")),
    d15N    = tidyselect::any_of(c("d15N","delta15N","N15"))) 


 # {OPTIONAL} Pool the three Lutjanus spp. *before* choosing the source list

 #   %>% mutate(
 #   Species = ifelse(Species %in% c("L. bohar","L. decussatus","L. gibbous"),
 #                    "L. bo./dec./gib.", Species))


# 3) Mixtures (predator individuals @ Rowleys)
mix_data <- df %>%
  filter(Species == predator, Reef == "Rowleys") %>%
  select(d13C, d15N) %>%
  drop_na()

stopifnot(nrow(mix_data) > 0)
mix <- as.matrix(mix_data)


# 4) Source means & SDs (Rowleys only) for the 9 sources
source_data <- df %>%
  filter(Reef == "Rowleys", Species %in% preys) %>%
  group_by(Species) %>%
  summarise(
    mean_d13C = mean(d13C, na.rm = TRUE),
    sd_d13C   = sd(d13C,   na.rm = TRUE),
    mean_d15N = mean(d15N, na.rm = TRUE),
    sd_d15N   = sd(d15N,   na.rm = TRUE),
    n         = sum(!is.na(d13C) & !is.na(d15N)),
    .groups = "drop"
  ) %>%
  filter(n > 0) %>%
  # enforce plotting/model order = 'preys'
  arrange(match(Species, preys))

# Guard against zero/NA SDs
source_data <- source_data %>%
  mutate(
    sd_d13C = ifelse(is.na(sd_d13C) | sd_d13C == 0, 1e-6, sd_d13C),
    sd_d15N = ifelse(is.na(sd_d15N) | sd_d15N == 0, 1e-6, sd_d15N)
  )

s_names <- source_data$Species
s_means <- as.matrix(source_data %>% select(mean_d13C, mean_d15N))
s_sds   <- as.matrix(source_data %>% select(sd_d13C,    sd_d15N))


# 5) TEFs (Caut/Cardona) on *source means* with element-specific SDs
c_means <- cbind(
  tef_c13(source_data$mean_d13C),
  tef_n15(source_data$mean_d15N)
)
c_sds <- matrix(c(TEF_SD_C, TEF_SD_N), nrow = nrow(c_means), ncol = 2, byrow = TRUE)


# 6) simmr model
simmr_in <- simmr_load(
  mixtures         = mix,
  source_names     = s_names,
  source_means     = s_means,
  source_sds       = s_sds,
  correction_means = c_means,
  correction_sds   = c_sds
)

simmr_out <- simmr_mcmc(simmr_in, mcmc_control = mcmc_ctrl)

# Convergence diagnostics (should be ~1)
summary(simmr_out, type = "diagnostics") 


# 7) Summaries (Mean, SD, 95% CrI) in %
# Extract posterior p's (stack all chains)
post_list <- simmr_out$output
post_mat  <- do.call(rbind, lapply(post_list, function(ch)
  ch$BUGSoutput$sims.list$p))

colnames(post_mat) <- simmr_out$input$source_names

posterior_long <- as.data.frame(post_mat) %>%
  tidyr::pivot_longer(everything(), names_to = "Source", values_to = "Prop") %>%
  mutate(
    Source = factor(Source, levels = preys),
    Prop   = Prop * 100
  )

summary_table <- posterior_long %>%
  group_by(Source) %>%
  summarise(
    Mean = mean(Prop),
    SD   = sd(Prop),
    L95  = quantile(Prop, 0.025),
    U95  = quantile(Prop, 0.975),
    .groups = "drop"
  ) %>%
  arrange(match(Source, preys))

print(summary_table, n = nrow(summary_table))


summary_table <- summary_table %>%
  rename(Species = Source)

Shark_diet<-summary_table




# Reorder factor levels by Mean (descending)
summary_table <- summary_table %>%
  mutate(Species = factor(Species, levels = Species[order(-Mean)]))

# Plot


Shark_diet_plot<-ggplot(summary_table, aes(x = Species, y = Mean)) +
  geom_col(width = 0.7, fill = "lightblue") +
  geom_errorbar(aes(ymin = pmax(0, Mean - SD),
                    ymax = pmin(100, Mean + SD)),
                width = 0.15) +
  labs(
    x = "Prey species",
    y = "Mean % contribution to diet",
    title = bquote("Estimated diet of " * italic(.(predator)) * " at Rowley Shoals")
  ) +
  scale_x_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank())

Shark_diet_plot

tiff("Shark_diet_plot.tiff", units="cm", width=20, height=10, res=300)
Shark_diet_plot
dev.off()






#Teleost Trophic position####
species_keep <- c("C. microrhinos", "C. sordidus", "L. kasmira",
                  "Z. scopas", "M. grandoculis", "P. vaiuli",
                  "S. vulpinus", "F. flavissimus", "L. bohar",
                  "L. decussatus", "L. gibbous")

Delta_d15N <- 3.4   # trophic enrichment factor, Vander Zanden & Rasmussen (2001)
TP0 <- 2            # baseline = primary consumer TP

# Reef-specific baselines: mean δ15N of C. sordidus per reef
reef_baselines <- data %>%
  filter(Reef %in% c("Rowleys", "Scott"),
         Species == "C. sordidus") %>%
  summarise(
    d15N_baseline = mean(d15N, na.rm = TRUE),
    .by = Reef)

# Calculate TP for each individual of selected species
tp_individuals <- data %>%
  filter(Reef %in% c("Rowleys", "Scott"),
         Species %in% species_keep) %>%
  left_join(reef_baselines, by = "Reef") %>%
  mutate(
    TP = ifelse(
      !is.na(d15N) & !is.na(d15N_baseline),
      TP0 + (d15N - d15N_baseline) / Delta_d15N,
      NA_real_ )
  )

tp_individuals %>%
  select(Species, Reef, d15N, d15N_baseline, TP)   # inspect results

# 3) Summarize TP per Species × Reef
tp_summary <- tp_individuals %>%
  group_by(Species, Reef) %>%
  summarise(
    mean_TP = mean(TP, na.rm = TRUE),
    sd_TP   = sd(TP, na.rm = TRUE),
    n       = sum(!is.na(TP)),
    .groups = "drop"
  )

tp_summary   # dataframe with mean TP per species × reef

TP_teleosts<-tp_summary

write.csv(TP_teleosts, "TP_teleosts.csv", row.names = FALSE)

#View(TP_teleosts)

#%ΔTP Estimation, Scott relative to Rowleys####
# Collapse so we have one mean TP per species × reef
tp_summary2 <- tp_summary %>%
  group_by(Species, Reef) %>%
  summarise(mean_TP = mean(mean_TP, na.rm = TRUE), .groups = "drop")

# Pivot to wide format (Rowleys vs Scott)
tp_wide <- tp_summary2 %>%
  pivot_wider(names_from = Reef, values_from = mean_TP)

# Calculate %ΔTP (Scott relative to Rowleys)
tp_diff <- tp_wide %>%
  mutate(
    perc_diff_TP = 100 * ((Scott - Rowleys) / Rowleys)
  )


tp_diff

write.csv(tp_diff, "TP_percentage_difference_teleosts", row.names = FALSE)



#View(tp_diff)

#PLOT Mean trophic position Rowleys and shoals####
#colors
reef_cols <- c(
  "Scott"   = '#F8766D',
  "Rowleys" = "#0033C9")



plot_3_df <- TP_teleosts %>%
  inner_join(tp_diff, by = "Species") 


# Vector of species to drop
dropping <- c("C. microrhinos", 
           "C. sordidus", 
           "Z. scopas",
           "P. vaiuli",
           "S. vulpinus",
           "F. flavissimus")

# Drop them from dataframe df
df_filtered <- plot_3_df %>%
  filter(!Species %in% dropping)

df_filtered


#Compute p-values from t-tests on log10(TP) between reefs, per species

pvals_df <- tp_individuals %>%
  filter(Species %in% unique(df_filtered$Species),
         Reef %in% c("Rowleys","Scott"),
         is.finite(TP)) %>%
  mutate(log_TP = log10(TP)) %>%
  group_by(Species) %>%
  summarise(
    # require both reefs and at least 2 obs per reef
    p_value = {
      sub <- cur_data()
      ok <- all(c("Rowleys","Scott") %in% sub$Reef) &&
        min(table(sub$Reef)) >= 2
      if (ok) t.test(log_TP ~ Reef, data = sub, var.equal = TRUE)$p.value else NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(p_label = paste0("p = ", signif(p_value, 3)))

# 2) Build label positions and the +% labels from your plotting dataframe
labels_df <- df_filtered %>%
  group_by(Species) %>%
  summarise(
    perc_diff_TP = unique(perc_diff_TP),
    y_pos = max(mean_TP + sd_TP) + 0.3,   # a bit above error bars
    .groups = "drop"
  ) %>%
  left_join(pvals_df, by = "Species") %>%
  mutate(perc_label = paste0("+", round(perc_diff_TP, 1), "%"))

# 3) Plot with %ΔTP and p-values stacked above each species group
df_filtered$Species <- factor(df_filtered$Species,
                              levels = c("L. kasmira",
                                         "M. grandoculis",
                                         "L. decussatus",
                                         "L. gibbous",
                                         "L. bohar"))

labels_df$Species <- factor(labels_df$Species,
                            levels = levels(df_filtered$Species))


#plot
TP_teleosts_plot<-ggplot(df_filtered, aes(x = Species, y = mean_TP, fill = Reef)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7, color = "black", alpha = 0.9) +
  geom_errorbar(aes(ymin = mean_TP - sd_TP, ymax = mean_TP + sd_TP, color = Reef),
                position = position_dodge(width = 0.8),
                width = 0.25, linewidth = 0.7) +
  scale_fill_manual(values = reef_cols) +
  scale_color_manual(values = reef_cols, guide = "none") +
  geom_text(data = labels_df,
            aes(x = Species, y = y_pos, label = paste0("+", round(perc_diff_TP, 1), "%")),
            inherit.aes = FALSE, vjust = 0, fontface = "bold", size = 3.5) +
  geom_text(data = labels_df,
            aes(x = Species, y = y_pos + 0.25, label = p_label),
            inherit.aes = FALSE, vjust = 0, size = 3.2) +
  labs(
    title = "Mean trophic position (± SD) by species and reef",
    x = "Species",
    y = "Mean trophic position (TP)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "italic"),
    axis.title  = element_text(size = 12, face = "bold"),
    plot.title  = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "top",
    legend.title    = element_blank()
  )

TP_teleosts_plot

tiff("TP_differeces_plot.tiff", units="cm", width=20, height=14, res=300)
TP_teleosts_plot
dev.off()



#Species lenghts means per reef####
length_means <- data %>%
  mutate(length_mm = as.numeric(length_mm)) %>%       # convert chr → numeric
  filter(Reef %in% c("Rowleys", "Scott")) %>%
  filter(!is.na(length_mm)) %>%                      # drop missing values
  group_by(Species, Reef) %>%
  summarise(
    mean_length_mm = mean(length_mm, na.rm = TRUE),
    sd_length_mm   = sd(length_mm, na.rm = TRUE),    # standard deviation
    n = n(),
    .groups = "drop"
  )

length_means   # dataframe of mean lengths + SD

write.csv(length_means, "length_means.csv", row.names = FALSE)


#Species weight means per reef####
weight_means <- data %>%
  mutate(
    weight_kg_num = readr::parse_number(gsub(",", ".", weight_kg))  # clean & convert
  ) %>%
  filter(Reef %in% c("Rowleys", "Scott")) %>%
  filter(!is.na(weight_kg_num)) %>%
  group_by(Species, Reef) %>%
  summarise(
    mean_weight_kg = mean(weight_kg_num, na.rm = TRUE),
    sd_weight_kg   = sd(weight_kg_num, na.rm = TRUE),   # add standard deviation
    n = dplyr::n(),
    .groups = "drop"
  )

weight_means   # dataframe with mean, sd, and sample size

write.csv(weight_means, "weight_means.csv", row.names = FALSE)


#T-tests (tweo way) for lenght/weight####
species_keep <- c("C. microrhinos","C. sordidus","L. kasmira","Z. scopas",
                  "M. grandoculis","P. vaiuli","S. vulpinus","F. flavissimus",
                  "L. bohar","L. decussatus","L. gibbous")

dat <- data %>%
  filter(Reef %in% c("Rowleys","Scott"),
         Species %in% species_keep) %>%
  mutate(
    length_mm = as.numeric(length_mm),
    weight_kg_num = parse_number(gsub(",", ".", weight_kg)),
    Reef = factor(Reef, levels = c("Rowleys","Scott")),
    Species = factor(Species)
  )


data$length_mm <- as.numeric(data$length_mm)
data$weight_kg_num <- readr::parse_number(gsub(",", ".", data$weight_kg))

# Unique species list
species_list <- unique(data$Species)

# T-test length
ttest_length <- lapply(species_list, function(sp) {
  subdat <- subset(data, Species == sp & Reef %in% c("Rowleys","Scott") & !is.na(length_mm))
  if(length(unique(subdat$Reef)) == 2) {
    res <- t.test(length_mm ~ Reef, data = subdat)
    c(Species = sp,
      p.value = res$p.value,
      mean_Rowleys = mean(subdat$length_mm[subdat$Reef=="Rowleys"], na.rm=TRUE),
      mean_Scott   = mean(subdat$length_mm[subdat$Reef=="Scott"], na.rm=TRUE))
  }
})
ttest_length <- do.call(rbind, ttest_length)
ttest_length

# T-test weight
ttest_weight <- lapply(species_list, function(sp) {
  subdat <- subset(data, Species == sp & Reef %in% c("Rowleys","Scott") & !is.na(weight_kg_num))
  if(length(unique(subdat$Reef)) == 2) {
    res <- t.test(weight_kg_num ~ Reef, data = subdat)
    c(Species = sp,
      p.value = res$p.value,
      mean_Rowleys = mean(subdat$weight_kg_num[subdat$Reef=="Rowleys"], na.rm=TRUE),
      mean_Scott   = mean(subdat$weight_kg_num[subdat$Reef=="Scott"], na.rm=TRUE))
  }
})
ttest_weight <- do.call(rbind, ttest_weight)
ttest_weight



#Linear Mixed-Effects Model (LMM)####
#to test the overall Reef effect (Rowleys vs Scott) across species (random intercept for Species), 
#for length and weight, and per-species Welch t-tests for Reef differences 
#with Holm-adjusted p-values (multiple testing correction).


# Be robust to length column name: take length_mm if present, else Lenth_mm
length_col <- if ("length_mm" %in% names(data)) {
  "length_mm"
} else if ("Lenth_mm" %in% names(data)) {
  "Lenth_mm"
} else {
  stop("No length column named 'length_mm' or 'Lenth_mm' found.")
}

dat <- data %>%
  filter(Reef %in% c("Rowleys","Scott"),
         Species %in% species_keep) %>%
  mutate(
    # Length to numeric (works whether it's 'length_mm' or 'Lenth_mm')
    !!length_col := as.numeric(.data[[length_col]]),
    # Weight: clean character like "1,23 kg" → 1.23
    weight_kg_num = parse_number(gsub(",", ".", weight_kg)),
    Reef   = factor(Reef, levels = c("Rowleys","Scott")),
    Species = factor(Species)
  )

# For convenience create plain vectors
dat$length_num <- dat[[length_col]]  # numeric length column used below

##  LMMs: overall reef effect across species 
# LENGTH
fit_length <- lmer(length_num ~ Reef + (1 | Species), data = dat, REML = TRUE)
cat("\n===== LMM: LENGTH =====\n")
print(summary(fit_length))   # Fixed effect for Reef tests Rowleys vs Scott overall
print(anova(fit_length))     # F-test for Reef

# WEIGHT
fit_weight <- lmer(weight_kg_num ~ Reef + (1 | Species),
                   data = dat %>% filter(!is.na(weight_kg_num)),
                   REML = TRUE)
cat("\n===== LMM: WEIGHT =====\n")
print(summary(fit_weight))
print(anova(fit_weight))






species_list <- levels(dat$Species)

# Helper to safely build a row (as named vector) from a t-test
make_t_row <- function(sp, res, mean_R, mean_S, n_R, n_S) {
  c(Species = sp,
    mean_Rowleys = mean_R,
    mean_Scott   = mean_S,
    n_Rowleys    = n_R,
    n_Scott      = n_S,
    estimate_diff = unname(res$estimate[[1]] - res$estimate[[2]]), # may be NA if names vary
    t_stat       = unname(res$statistic),
    df           = unname(res$parameter),
    p_value      = unname(res$p.value))
}

## LENGTH t-tests
tt_rows_len <- lapply(species_list, function(sp) {
  subdat <- subset(dat, Species == sp & !is.na(length_num))
  if (length(unique(subdat$Reef)) < 2) return(NULL)
  
  mean_R <- mean(subdat$length_num[subdat$Reef=="Rowleys"], na.rm = TRUE)
  mean_S <- mean(subdat$length_num[subdat$Reef=="Scott"],   na.rm = TRUE)
  n_R    <- sum(subdat$Reef=="Rowleys" & !is.na(subdat$length_num))
  n_S    <- sum(subdat$Reef=="Scott"   & !is.na(subdat$length_num))
  
  res <- t.test(length_num ~ Reef, data = subdat)  # Welch by default
  make_t_row(sp, res, mean_R, mean_S, n_R, n_S)
})
ttest_length <- if (length(Filter(Negate(is.null), tt_rows_len)) > 0) {
  as.data.frame(do.call(rbind, tt_rows_len), stringsAsFactors = FALSE)
} else data.frame()

# Convert numeric columns back to numeric (they come in as characters)
num_cols_len <- c("mean_Rowleys","mean_Scott","n_Rowleys","n_Scott",
                  "estimate_diff","t_stat","df","p_value")
if (nrow(ttest_length) > 0) ttest_length[num_cols_len] <- lapply(ttest_length[num_cols_len], as.numeric)

# Holm correction
if (nrow(ttest_length) > 0) {
  ttest_length$p_adj_holm <- p.adjust(ttest_length$p_value, method = "holm")
}

cat("\n===== Per-species Welch t-tests: LENGTH =====\n")
print(ttest_length)

## ---- WEIGHT t-tests
tt_rows_wt <- lapply(species_list, function(sp) {
  subdat <- subset(dat, Species == sp & !is.na(weight_kg_num))
  if (length(unique(subdat$Reef)) < 2) return(NULL)
  
  mean_R <- mean(subdat$weight_kg_num[subdat$Reef=="Rowleys"], na.rm = TRUE)
  mean_S <- mean(subdat$weight_kg_num[subdat$Reef=="Scott"],   na.rm = TRUE)
  n_R    <- sum(subdat$Reef=="Rowleys" & !is.na(subdat$weight_kg_num))
  n_S    <- sum(subdat$Reef=="Scott"   & !is.na(subdat$weight_kg_num))
  
  res <- t.test(weight_kg_num ~ Reef, data = subdat)  # Welch by default
  make_t_row(sp, res, mean_R, mean_S, n_R, n_S)
})
ttest_weight <- if (length(Filter(Negate(is.null), tt_rows_wt)) > 0) {
  as.data.frame(do.call(rbind, tt_rows_wt), stringsAsFactors = FALSE)
} else data.frame()

num_cols_wt <- c("mean_Rowleys","mean_Scott","n_Rowleys","n_Scott",
                 "estimate_diff","t_stat","df","p_value")
if (nrow(ttest_weight) > 0) ttest_weight[num_cols_wt] <- lapply(ttest_weight[num_cols_wt], as.numeric)

if (nrow(ttest_weight) > 0) {
  ttest_weight$p_adj_holm <- p.adjust(ttest_weight$p_value, method = "holm")
}

cat("\n===== Per-species Welch t-tests: WEIGHT =====\n")
print(ttest_weight)

## (Optional) Save tables
# write.csv(ttest_length, "ttest_length_by_species.csv", row.names = FALSE)
# write.csv(ttest_weight, "ttest_weight_by_species.csv", row.names = FALSE)

#PLOT of models####
# Optional reef colors
reef_cols <- c(
  "Scott"   = '#F8766D',
  "Rowleys" = "#0033C9")

# LENGTH summary per species
len_summary <- dat %>%
  filter(!is.na(length_num)) %>%
  group_by(Species, Reef) %>%
  summarise(
    mean = mean(length_num, na.rm = TRUE),
    sd   = sd(length_num, na.rm = TRUE),
    n    = n(),
    se   = sd / sqrt(n),   # standard error
    .groups = "drop"
  )

p_len <- ggplot(len_summary, aes(x = Reef, y = mean, color = Reef)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
  facet_wrap(~ Species, scales = "free_y") +
  scale_color_manual(values = reef_cols) +
  labs(
    title = "Mean LENGTH (±SE) per species by Reef",
    x = "Reef",
    y = "Mean length (mm)"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_len


# WEIGHT summary per species
wt_summary <- dat %>%
  filter(!is.na(weight_kg_num)) %>%
  group_by(Species, Reef) %>%
  summarise(
    mean = mean(weight_kg_num, na.rm = TRUE),
    sd   = sd(weight_kg_num, na.rm = TRUE),
    n    = n(),
    se   = sd / sqrt(n),
    .groups = "drop"
  )

# Plot
p_wt <- ggplot(wt_summary, aes(x = Reef, y = mean, color = Reef)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
  facet_wrap(~ Species, scales = "free_y") +
  scale_color_manual(values = reef_cols) +
  labs(
    title = "Mean WEIGHT (±SE) per species by Reef",
    x = "Reef",
    y = "Mean weight (kg)"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_wt


#%ΔW Niche width  change per species and reef####
  meso_spp <- c("L. bohar", "L. decussatus", "L. gibbous",
                "L. kasmira", "M. grandoculis", "P. vaiuli",
                "S. vulpinus", "F. flavissimus")

# Compute niche width (δ13C range) and height (δ15N range) per Species × Reef 
niche_table <- data %>%
  filter(Reef %in% c("Rowleys", "Scott"),
         Species %in% meso_spp) %>%
  filter(!is.na(d13C), !is.na(d15N)) %>%        # keep complete isotope pairs
  group_by(Species, Reef) %>%
  summarise(
    n          = n(),
    niche_width_d13C  = diff(range(d13C)),      # δ13C max - min
    niche_height_d15N = diff(range(d15N)),      # δ15N max - min
    .groups = "drop"
  ) %>%
  arrange(Species, Reef)

# Print per-reef niche metrics
niche_table

# % change in niche width (Scott vs Rowleys) for each species
niche_wide <- niche_table %>%
  select(Species, Reef, niche_width_d13C) %>%
  pivot_wider(names_from = Reef, values_from = niche_width_d13C)

niche_change <- niche_wide %>%
  mutate(
    perc_change_width_Scott_vs_Rowleys =
      ifelse(!is.na(Rowleys) & Rowleys != 0,
             100 * ((Scott - Rowleys) / Rowleys),
             NA_real_)  # avoid divide-by-zero or missing Rowleys
  ) %>%
  arrange(Species)

# Print % change table
niche_change

# (Optional) merge counts back in for context 
n_counts <- niche_table %>%
  select(Species, Reef, n) %>%
  pivot_wider(names_from = Reef, values_from = n, names_prefix = "n_")

niche_change_with_n <- niche_change %>%
  left_join(n_counts, by = "Species")

# View with sample sizes per reef
niche_change_with_n


#PLOTS - Linear Regressions####
TP_scott <- TP_teleosts %>%
  filter(Reef == "Scott")

TP_scott              #dataframe with trophic positions
tp_diff                  #dataframe with trophic positions % changes
niche_change_with_n      #dataframe with niche % changes
Shark_diet               #dataframe with shark diet %


combined <- TP_scott %>%
  inner_join(tp_diff, by = "Species") %>%
  inner_join(niche_change_with_n, by = "Species") %>%
  inner_join(Shark_diet, by = "Species")




#keeping species to match the  species used in the publication draft
mesopredators_to_keep <- c( "L. kasmira",
           "M. grandoculis",
           "L. bohar",
           "L. decussatus", 
           "L. gibbous")


combined_new <- combined %>%
  filter(Species %in% mesopredators_to_keep)


#Tidy data for each model
  d1 <- combined_new %>%
    select(Species, mean_TP, perc_change_width_Scott_vs_Rowleys) %>%
    rename(x = mean_TP, y = perc_change_width_Scott_vs_Rowleys) %>%
    filter(is.finite(x), is.finite(y))
  
  d2 <- combined_new %>%
    select(Species, perc_diff_TP, perc_change_width_Scott_vs_Rowleys) %>%
    rename(x = perc_diff_TP, y = perc_change_width_Scott_vs_Rowleys) %>%
    filter(is.finite(x), is.finite(y))
  
  d3 <- combined_new %>%    # assuming 'Mean' is % diet contribution
    select(Species, Mean, perc_diff_TP) %>%
    rename(x = Mean, y = perc_diff_TP) %>%
    filter(is.finite(x), is.finite(y))
  
  #Helper to compute model labels 
  model_label <- function(dat){
    fit <- lm(y ~ x, data = dat)
    s   <- summary(fit)
    slope <- unname(coef(fit)[2])
    r2    <- s$r.squared
    pval  <- coef(s)[2, "Pr(>|t|)"]
    list(fit = fit,
         lab = sprintf("slope = %.3f\nR² = %.3f\np = %.3g", slope, r2, pval))
  }
  
  m1 <- model_label(d1)
  m2 <- model_label(d2)
  m3 <- model_label(d3)
  
#Journal-friendly theme 
  theme_pub <- theme_classic(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", size = 12, hjust = 0),
      axis.title   = element_text(size = 11),
      axis.text    = element_text(size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
  
  pt_sz <- 2.6
  line_w <- 0.9
  
  # Build the three panels
  L1 <- ggplot(d1, aes(x = x, y = y)) +
    geom_point(size = pt_sz) +
    geom_smooth(method = "lm", se = TRUE, linewidth = line_w) +
    labs(
      title = "A) % niche-width change vs trophic position",
      x = "Mean trophic position (TP)",
      y = "% change in niche width (Scott vs Rowleys)"
    ) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1, label = m1$lab) +
    theme_pub
  
  L2 <- ggplot(d2, aes(x = x, y = y)) +
    geom_point(size = pt_sz) +
    geom_smooth(method = "lm", se = TRUE, linewidth = line_w) +
    labs(
      title = "B) % niche-width change vs ΔTP (%)",
      x = "ΔTP (%) (Scott relative to Rowleys)",
      y = "% change in niche width (Scott vs Rowleys)"
    ) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 4.5, label = m2$lab) +
    theme_pub
  
  L3 <- ggplot(d3, aes(x = x, y = y)) +
    geom_point(size = pt_sz) +
    geom_smooth(method = "lm", se = TRUE, linewidth = line_w) +
    labs(
      title = "C) ΔTP (%) vs prey contribution (%)",
      x = "Mean diet contribution (%CSD)",
      y = "ΔTP (%) (Scott relative to Rowleys)"
    ) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1, label = m3$lab) +
    theme_pub
  
  
# Arrange vertically without patchwork 
 grid.arrange(L1, L2, L3, ncol = 1, heights = c(1, 1, 1))
 
 
 tiff("Linear_regressions_plot.tiff", units="cm", width=12, height=29, res=300)
 grid.arrange(L1, L2, L3, ncol = 1, heights = c(1, 1, 1))
 dev.off()
 
 
 
  
#Ranges of C13 values for mesopredators####

# 1. Define your custom palette
 reef_cols <- c(
   "Scott"   = '#F8766D',
   "Rowleys" = "#0033C9")

# 2. Species to drop
drop_spp <- c("C. amblyrhynchos","H. atra","C. sordidus",
              "P. vaiuli","S. vulpinus","C. margaritifer",
              "C. microrhinos","Z. scopas")

# 3. Filter to just the two reefs of interest and drop unwanted species
data2 <- data %>% 
  filter(
    Reef %in% c("Scott", "Rowleys"),
    !Species %in% drop_spp
  ) %>%
  mutate(Species = factor(Species, levels = sort(unique(Species))))

# 4. Prepare summary for plot 3
summary_df <- data2 %>%
  group_by(Species, Reef) %>%
  summarise(
    min_d13C  = min(d13C, na.rm = TRUE),
    max_d13C  = max(d13C, na.rm = TRUE),
    mean_d13C = mean(d13C, na.rm = TRUE),
    .groups    = "drop"
  )

# 5a. Classic boxplot: δ13C by Species, coloured by Reef
p1 <- ggplot(data2, aes(x = Species, y = d13C, fill = Reef)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 1) +
  scale_fill_manual(values = reef_cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "δ¹³C distribution by Species & Reef",
    x     = "Species",
    y     = expression(delta^{13}*C)
  )

# 5b. Violin + jitter: show full density + individual points
p2 <- ggplot(data2, aes(x = Species, y = d13C, fill = Reef)) +
  geom_violin(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(
    aes(colour = Reef),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 1, alpha = 0.6
  ) +
  scale_fill_manual(values = reef_cols) +
  scale_colour_manual(values = reef_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Violin plot of δ¹³C by Species & Reef",
    x     = "Species",
    y     = expression(delta^{13}*C)
  )

# 5c. Range + mean plot: error bars from min→max with mean point
p3 <- ggplot(summary_df, aes(x = Species, group = Reef)) +
  geom_errorbar(
    aes(ymin = min_d13C, ymax = max_d13C, colour = Reef),
    position = position_dodge(width = 0.5), width = 0.2
  ) +
  geom_point(
    aes(y = mean_d13C, colour = Reef),
    position = position_dodge(width = 0.5), size = 2
  ) +
  scale_colour_manual(values = reef_cols) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Range (min–max) & mean δ¹³C by Species & Reef",
    x     = "Species",
    y     = expression(delta^{13}*C)
  )

# 5d. Ridge density plot: overlapping δ13C densities per species
p4 <- ggplot(data2, aes(x = d13C, y = Species, fill = Reef)) +
  geom_density_ridges(alpha = 0.6, scale = 1) +
  scale_fill_manual(values = reef_cols) +
  theme_minimal() +
  labs(
    title = "Density ridges of δ¹³C by Species & Reef",
    x     = expression(delta^{13}*C),
    y     = "Species"
  )

# 6. Arrange all four plots in a 2×2 grid
grid.arrange(p1, p2, p3, p4, ncol = 2)



tiff("Ranges_C13.tiff", units="cm", width=26, height=22, res=300)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()



#We can pick whichever plot we waant for the publication


#Ranges of N15 values for mesopredators####

# 1. Define your custom palette
reef_cols <- c(
  "Scott"   = '#F8766D',
  "Rowleys" = "#0033C9")


# 2. Species to drop
drop_spp <- c("C. amblyrhynchos","H. atra","C. sordidus",
              "P. vaiuli","S. vulpinus","C. margaritifer",
              "C. microrhinos","Z. scopas")

# 3. Filter to just the two reefs of interest and drop unwanted species
data2 <- data %>% 
  filter(
    Reef %in% c("Scott", "Rowleys"),
    !Species %in% drop_spp
  ) %>%
  mutate(Species = factor(Species, levels = sort(unique(Species))))

# 4. Prepare summary for plot 3
summary_df <- data2 %>%
  group_by(Species, Reef) %>%
  summarise(
    min_d15N  = min(d15N, na.rm = TRUE),
    max_d15N  = max(d15N, na.rm = TRUE),
    mean_d15N = mean(d15N, na.rm = TRUE),
    .groups    = "drop"
  )

# 5a. Classic boxplot: δ13C by Species, coloured by Reef
N1 <- ggplot(data2, aes(x = Species, y = d15N, fill = Reef)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 1) +
  scale_fill_manual(values = reef_cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "δ¹⁵N distribution by Species & Reef",
    x     = "Species",
    y     = expression(delta^{15}*N)
  )

# 5b. Violin + jitter: show full density + individual points
N2 <- ggplot(data2, aes(x = Species, y = d15N, fill = Reef)) +
  geom_violin(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(
    aes(colour = Reef),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 1, alpha = 0.6
  ) +
  scale_fill_manual(values = reef_cols) +
  scale_colour_manual(values = reef_cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Violin plot of δ¹⁵N by Species & Reef",
    x     = "Species",
    y     = expression(delta^{15}*N)
  )

# 5c. Range + mean plot: error bars from min→max with mean point
N3 <- ggplot(summary_df, aes(x = Species, group = Reef)) +
  geom_errorbar(
    aes(ymin = min_d15N, ymax = max_d15N, colour = Reef),
    position = position_dodge(width = 0.5), width = 0.2
  ) +
  geom_point(
    aes(y = mean_d15N, colour = Reef),
    position = position_dodge(width = 0.5), size = 2
  ) +
  scale_colour_manual(values = reef_cols) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Range (min–max) & mean δ¹⁵N by Species & Reef",
    x     = "Species",
    y     = expression(delta^{15}*N)
  )

# 5d. Ridge density plot: overlapping δ¹⁵N densities per species
N4 <- ggplot(data2, aes(x = d15N, y = Species, fill = Reef)) +
  geom_density_ridges(alpha = 0.6, scale = 1) +
  scale_fill_manual(values = reef_cols) +
  theme_minimal() +
  labs(
    title = "Density ridges of δ¹⁵N by Species & Reef",
    x     = expression(delta^{15}*N),
    y     = "Species"
  )

# 6. Arrange all four plots in a 2×2 grid
grid.arrange(N1, N2, N3, N4, ncol = 2)



tiff("Ranges_N13.tiff", units="cm", width=26, height=22, res=300)
grid.arrange(N1, N2, N3, N4, ncol = 2)
dev.off()

#Standard ellipses####

# Prepare data: filter missing values, exclude H. atra and S. vulpinus, and create identifiers
data2 <- data %>%
  filter(!is.na(d13C), !is.na(d15N),
         !Species %in% c("C. amblyrhynchos", "F. flavissimus", "H. atra","C. sordidus",
                         "P. vaiuli","S. vulpinus","C. margaritifer","C. microrhinos","Z. scopas")) %>% #excluded species
  mutate(
    Species = factor(Species),
    Reef    = factor(Reef),
    group   = paste(Species, Reef, sep = "_")
  )

# Compute convex hull for each Species × reef group
hull_data <- data2 %>%
  group_by(group) %>%
  slice(chull(d13C, d15N)) %>%
  ungroup()

# Isotopic niche plots with colored hulls and points
ellipse_plot<-ggplot(data2, aes(x = d13C, y = d15N)) +
  # Colored hull polygons
  geom_polygon(
    data = hull_data,
    aes(group = group, fill = group),
    color = "black",
    size = 0.5,
    alpha = 0.4
  ) +
  # Points filled with the same group color
  geom_point(
    aes(fill = group),
    shape = 21,
    color = "black",
    size = 2,
    stroke = 0.5) +
  #xlim(-17.5,-8.5)+ylim(8.5,12.5)
  # Facet: reefs in rows, species in columns
  facet_grid(Reef ~ Species, switch = "y") +
  coord_equal() +
  # Publication-ready theme
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", size = 0.5),
    strip.background = element_rect(fill = "grey90", color = "black", size = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.ticks = element_line(color = "black"),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  ) +
  labs(
    title = "Isotopic Niche Areas per Species and Reef",
    x = expression(delta^{13}*"C (‰)"),
    y = expression(delta^{15}*"N (‰)")
  )


ellipse_plot

tiff("Ellipse_plot.tiff", units="cm", width=25, height=12, res=300)
ellipse_plot
dev.off()



#95% Standard Ellipses ####

# === Load Data ===
df <- read_excel("Stable_isotopes_MM.xlsx", sheet = "Sheet1")


suppressPackageStartupMessages({
  library(tidyverse)
})

#columns (edit if needed) 
iso_c_col   <- "d13C"      # carbon
iso_n_col   <- "d15N"      # nitrogen
group_col   <- "Reef"      # grouping to overlap within each species
species_col <- "Species"   # species column

# your species 
preys <- c( "L. bohar", "L. decussatus",    "L. gibbous",     "L. kasmira",
            "M. grandoculis")

# ---- your requested colors 
reef_cols <- c(
  "Scott"   = "#F8766D",
  "Rowleys" = "#0033C9"
)

# prep & sanity 
stopifnot(all(c(iso_c_col, iso_n_col, group_col, species_col) %in% names(data)))

df <- data %>%
  dplyr::rename(x = !!iso_c_col,
                y = !!iso_n_col,
                grp = !!group_col,
                spp = !!species_col) %>%
  drop_na(x, y, grp, spp) %>%
  filter(spp %in% preys) %>%
  mutate(
    grp = factor(grp),
    spp = factor(spp, levels = preys)  # keep facet order as given
  ) %>%
  group_by(spp, grp) %>%
  filter(n() >= 3) %>%                # need >=3 per group for covariance
  ungroup()

# Optional: keep only species with at least two groups present
# df <- df %>%
#   group_by(spp) %>% filter(n_distinct(grp) >= 2) %>% ungroup()

# plot: nicer styling, your palette, outlines only 
p_all <- ggplot(df, aes(x = x, y = y)) +
  # points with subtle black outline
  geom_point(aes(fill = grp),
             shape = 21, size = 2.6, alpha = 0.85,
             color = "black", stroke = 0.35) +
  
  # 95% ellipse FILL (translucent)
  stat_ellipse(
    aes(fill = grp),
    level = 0.95, type = "norm",
    geom = "polygon", alpha = 0.30, color = NA
  ) +
  # 95% ellipse OUTLINE
  stat_ellipse(
    aes(color = grp),
    level = 0.95, type = "norm",
    geom = "path", linewidth = 1.05
  ) +
  
  facet_wrap(~ spp, ncol = 3, scales = "fixed") +
  coord_equal() +
  scale_color_manual(values = reef_cols, name = group_col) +
  scale_fill_manual(values = reef_cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  labs(
    title = "δ¹³C–δ¹⁵N: 95% Standard Ellipses by Species",
    subtitle = "Filled ellipses: 95% ML; points are individuals",
    x = expression(delta^{13}*C~("\u2030")),
    y = expression(delta^{15}*N~("\u2030"))
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(linewidth = 0.25, linetype = "dotted"),
    panel.background  = element_rect(fill = "white", color = NA)) # frame
    
  
p_all


tiff("Niche_overlap__plot.tiff", units="cm", width=25, height=12, res=300)
p_all
dev.off()



#Niche overalps % estimations####
preys <- c("L. bohar","L. decussatus","L. gibbous",
           "L. kasmira","M. grandoculis")

nsamp  <- 1000   # posterior draws for parameters per group
nprob  <- 1000   # Monte Carlo reps for overlap probability
alpha  <- 0.95   # 95% niche region

# Restrict to the two reefs and isotopes needed
colnames(data)

df2 <- subset(data,
              Species %in% preys & Reef %in% c("Rowleys","Scott"),
              select = c(Species, Reef, d13C, d15N))

for (sp in preys) {
  dat <- subset(df2, Species == sp)
  
  # Need both reefs present
  if (length(unique(dat$Reef)) < 2L) {
    cat("\n", sp, ": skipped (needs both Rowleys and Scott)\n", sep = "")
    next
  }
  
  # Build posteriors for each reef (2D: d13C, d15N)
  posts <- tapply(1:nrow(dat), dat$Reef, function(ii) {
    niw.post(nsamples = nsamp, X = as.matrix(dat[ii, c("d13C","d15N")]))
  })
  
  # Overlap draws (directional A→B)
  over.stat <- overlap(posts, nreps = nsamp, nprob = nprob, alpha = alpha)
  
  # Mean overlap (%) for the single alpha provided
  over.mean <- apply(over.stat, c(1, 2), mean) * 100
  
  cat("\n=== ", sp, " — Overlap % (A in B, α = ", alpha*100, "%) ===\n", sep = "")
  print(round(over.mean, 1))
  cat("Rows = A (from), Cols = B (to). Directional & asymmetric.\n")
}



results_list <- list()  # empty list to collect

for (sp in preys) {
  dat <- subset(df2, Species == sp)
  
  # Need both reefs present
  if (length(unique(dat$Reef)) < 2L) {
    message(sp, " skipped (needs both Rowleys and Scott)")
    next
  }
  
  # Build posteriors
  posts <- tapply(1:nrow(dat), dat$Reef, function(ii) {
    niw.post(nsamples = nsamp, X = as.matrix(dat[ii, c("d13C","d15N")]))
  })
  
  # Overlap draws
  over.stat <- overlap(posts, nreps = nsamp, nprob = nprob, alpha = alpha)
  
  # Mean overlap (%)
  over.mean <- apply(over.stat, c(1, 2), mean) * 100
  
  # Save as tidy table
  res_df <- as.data.frame(over.mean)
  res_df$From <- rownames(over.mean)
  res_df$Species <- sp
  
  results_list[[sp]] <- res_df
}

# Bind all species into one data frame
results_all <- do.call(rbind, results_list)

# Export to CSV
write.csv(results_all, "overlap_results.csv", row.names = FALSE)



