# -----------------------------
# Setup
# -----------------------------
library(tidyverse)
library(dagitty)
library(performance)


# -----------------------------
# Input Root and Path
# -----------------------------
root <- "/Users/miftahfaizah/Library/CloudStorage/OneDrive-UniversityofLeeds/PHD JOURNEY/YBMAP/CN_DATASET/CN_analysis/YBMAP_P3b_review"
sm_path    <- file.path(root, "2_Preprocessing/ND_26/bids_data_phenotypes/phone_social_media.tsv")
erp_path   <- file.path(root, "2_Preprocessing/ND_26/preprocessing_output/p3b_erp_summary_autoreject_nd_26_allch.csv")
parts_path <- file.path(root, "2_Preprocessing/ND_26/bids_data_phenotypes/participants.tsv")    # age, sex
socio_path <- file.path(root, "2_Preprocessing/ND_26/bids_data_phenotypes/socioeconomic.tsv")   # kss, imd_decile, self_ses, family_ses
sdq_path <- file.path(root, "2_Preprocessing/ND_26/bids_data_phenotypes/sdq.tsv")
stai_path <- file.path(root, "2_Preprocessing/ND_26/bids_data_phenotypes/stai.tsv")
swemwbs_path <- file.path(root, "2_Preprocessing/ND_26/bids_data_phenotypes/swemwbs.tsv")
rcads_path <- file.path(root, "2_Preprocessing/ND_26/bids_data_phenotypes/rcads.tsv")
qc_path <- file.path(root, "2_Preprocessing/ND_26/manual_channel_qc_nd26.csv")
behaviour_path <- file.path(root, "2_Preprocessing/ND_26/behaviour_p3b_nd26.csv")


# -----------------------------
# Load and Clean Data
# -----------------------------
# 1. sm data
sm <- read_tsv(sm_path, na=c("", "NA", "N/A", "n/a"), show_col_types = FALSE)
midpoints <- c(0.25, 0.75, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 8)   # index = code 1..9
n_school <- 190 # UK school calendar weighting: 190 school days,
n_nonschool <- 365 - n_school   # 175 non-school days
sm <- sm %>% 
  mutate(
    weekday_hrs      = midpoints[as.integer(as.numeric(sm_weekdays_ord))],
    weekend_hrs      = midpoints[as.integer(as.numeric(sm_weekend_holiday_ord))],
    sm_composite_hrs = (n_school * weekday_hrs + n_nonschool * weekend_hrs) / 365  
  )


# 2. eeg data
erp <- read_csv(erp_path, show_col_types = FALSE)
erp <- erp %>% 
  rename(participant_id = "subject") %>%
  filter(channels %in% c("O1", "O2")) %>% 
  filter(n_target_clean >= 20 & n_nontarget_clean >= 80)

roi <- erp %>% 
  group_by(participant_id) %>%
  summarise(p3b_roi = mean(target_minus_non_target_uv, na.rm=TRUE), .groups = "drop")
            

exclude <- read_csv(qc_path, show_col_types = FALSE)
exclude <- exclude %>% 
  rename(participant_id = "subject") %>% 
  filter(exclude_subject == "TRUE")

# 3. Other phenotype data
participants  <- read_tsv(parts_path, na = c("", "NA", "N/A", "n/a", "Not found"), show_col_types = FALSE) %>% select(participant_id, age, sex)
sdq           <- read_tsv(sdq_path, na = c("", "NA", "N/A", "n/a", "Not found"), show_col_types = FALSE) 
stai          <- read_tsv(stai_path, na = c("", "NA", "N/A", "n/a", "Not found"), show_col_types = FALSE) %>% select(participant_id, stai_s_total)
swemwbs       <- read_tsv(swemwbs_path, na = c("", "NA", "N/A", "n/a", "Not found"), show_col_types = FALSE) %>% select(participant_id, swemwbs_total)
socioeconomic <- read_tsv(socio_path, na = c("", "NA", "N/A", "n/a", "Not found"), show_col_types = FALSE) 
rcads         <- read_tsv(rcads_path, na = c("", "NA", "N/A", "n/a", "Not found"), show_col_types = FALSE) %>% select(participant_id, rcads25_anxiety_raw, rcads25_depression_raw)


# 4. behaviour_p3b_data
behaviour <- read_csv(behaviour_path, na = c("","NA","n/a","N/A"), show_col_types = FALSE) %>%
  rename(participant_id = subject) %>%
  select(participant_id, acc_overall, acc_target, acc_nontarget,
         rt_target_mean, rt_overall_mean, invalid_rate)

# 5. Maindataset
Data <- roi %>% 
  anti_join(exclude, by= "participant_id") %>%
  inner_join(sm, by= "participant_id") %>%
  left_join(participants, by= "participant_id") %>%
  left_join(sdq, by= "participant_id") %>%
  left_join(stai, by= "participant_id") %>%
  left_join(swemwbs, by= "participant_id") %>%
  left_join(socioeconomic, by= "participant_id") %>%
  left_join(rcads, by= "participant_id") %>% 
  left_join(behaviour, by = "participant_id")


Data <- Data %>% 
  filter(!is.na(p3b_roi), !is.na(sm_composite_hrs))

names(Data)
glimpse(Data)

num_vars <- c("age", "imd_decile", "self_ses", "family_ses", "kss",
              "sdq_internalizing", "sdq_externalizing", "sdq_total_difficulties", 
              "stai_s_total", "swemwbs_total", "rcads25_anxiety_raw", "rcads25_depression_raw",
              "acc_overall","acc_target","acc_nontarget","rt_target_mean","rt_overall_mean","invalid_rate")

Data <- Data %>% mutate(across(any_of(num_vars), as.numeric))
summary(Data$sm_composite_hrs)
hist(Data$sm_composite_hrs)

# -----------------------------
# Analyis
# -----------------------------
# 1. main correlation
Primary_correlation <- cor.test(Data$sm_composite_hrs, Data$p3b_roi, method = "spearman")
Primary_correlation

# 2. regression and covariate adjusted (exploratory)
model_unadjusted <- lm(p3b_roi ~ sm_composite_hrs, data = Data)
summary(model_unadjusted)
check_model(model_unadjusted)
cooks <- cooks.distance(model_unadjusted)
which(cooks > 4/nobs(model_unadjusted)) 

# --- Sensitivity: drop the Cook's-D influential point ---
infl <- which(cooks > 4/nobs(model_unadjusted))   # subject flagged (largest P3b)
Data_sens <- Data[-infl, ]

cor.test(Data_sens$sm_composite_hrs, Data_sens$p3b_roi, method = "spearman")  # rho robust?
summary(lm(p3b_roi ~ sm_composite_hrs, data = Data_sens))                     # β robust?

# robust regression on FULL data
summary(MASS::rlm(p3b_roi ~ sm_composite_hrs, data = Data))


# Adjusted model (exploratory)
  # imd_decile = confounder (SES -> SM use and -> P3b)
  # age, kss   = precision covariates (predict P3b, not SM); not confounders
  # Behaviour & mental-health variables are NOT adjusted (mediator/collider) - reported separately.
adjusted_variables <- c("p3b_roi", "sm_composite_hrs", "age", "kss", "imd_decile")
Data_complete <- Data %>% tidyr::drop_na(all_of(adjusted_variables))
nrow(Data_complete)                          
model_adjusted <- lm(p3b_roi ~ sm_composite_hrs + age +kss + imd_decile, data = Data_complete)
summary(model_adjusted) 
performance::check_collinearity(model_adjusted)
performance::check_model(model_adjusted)

# 3. behaviour x sm, behaviour x p3b, and behaviour as covariate 
cor.test(Data$sm_composite_hrs, Data$acc_target,     method = "spearman")
cor.test(Data$sm_composite_hrs, Data$rt_target_mean, method = "spearman")
cor.test(Data$p3b_roi, Data$acc_target,    method = "spearman")
cor.test(Data$p3b_roi, Data$acc_nontarget, method = "spearman")
summary(lm(p3b_roi ~ sm_composite_hrs + acc_target, data = Data))
summary(lm(p3b_roi ~ sm_composite_hrs + acc_nontarget, data = Data))


# 4. exploratory mental health as predictor and moderator
mh_unadjusted <- lm(p3b_roi ~ sdq_total_difficulties, data = Data)
summary(mh_unadjusted)                    

# -----------------------------
# PLOT
# -----------------------------
dir.create("figures_pilot", showWarnings = FALSE)

p1 <- ggplot(Data, aes(x = sm_composite_hrs, y = p3b_roi)) +
  geom_point(alpha = 0.6, color = "#1f77b4") +
  geom_smooth(method = MASS::rlm, se = FALSE, color = "#d62728", linewidth = 0.9) +   # robust (Huber) fit, outlier-resistant
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 4,
           label = sprintf("rho = %.2f, p = %.3f, n = %d",
                           Primary_correlation$estimate, Primary_correlation$p.value, nrow(Data))) +
  labs(x = "Social-media use (hours/day, time-weighted composite)",
       y = expression("P3b difference wave (target - non-target, " * mu * "V)"),
       title = "Higher social-media use, smaller occipital P3b") +
  theme_minimal(base_size = 13)
ggsave("figures_pilot/fig_scatter_composite_p3b.png", p1, width = 6.4, height = 4.6, dpi = 300)


Data$influential <- FALSE
Data$influential[infl] <- TRUE

p2 <- ggplot(Data, aes(sm_composite_hrs, p3b_roi)) +
  geom_point(aes(color = influential), alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "#d62728") +                 # all data
  geom_smooth(data = subset(Data, !influential), method = "lm",
              se = FALSE, color = "grey40", linetype = "dashed") +              # without outlier
  scale_color_manual(values = c(`FALSE` = "#1f77b4", `TRUE` = "red"),
                     labels = c("retained", "influential (Cook's D)"), name = NULL) +
  labs(x = "Social-media use (hours/day)", y = expression("P3b difference wave (" * mu * "V)"),
       title = "Sensitivity: trend with vs without influential point") +
  theme_minimal(base_size = 13)
ggsave("figures_pilot/fig_sensitivity_influential.png", p2, width = 6.4, height = 4.6, dpi = 300)


p3 <- ggplot(Data, aes(sm_composite_hrs)) +
  geom_histogram(binwidth = 1, boundary = 0, fill = "#4C78A8", color = "white") +
  labs(x = "Social-media use (hours/day)", y = "Number of participants",
       title = "Distribution of time-weighted composite social media") +
  theme_minimal()
ggsave("figures_pilot/fig_hist_sm_composite.png", p3, width = 6.0, height = 4.0, dpi = 300)

diss <- bind_rows(
  data.frame(var="SM composite", rho=cor(Data$sm_composite_hrs, Data$p3b_roi, method="spearman", use="complete.obs")),
  data.frame(var="SDQ total",    rho=cor(Data$sdq_total_difficulties, Data$p3b_roi, method="spearman", use="complete.obs")),
  data.frame(var="STAI anxiety", rho=cor(Data$stai_s_total, Data$p3b_roi, method="spearman", use="complete.obs")),
  data.frame(var="SWEMWBS",      rho=cor(Data$swemwbs_total, Data$p3b_roi, method="spearman", use="complete.obs"))
)

p4 <- ggplot(diss, aes(reorder(var, rho), rho, fill = var=="SM composite")) +
  geom_col() + coord_flip() +
  scale_fill_manual(values=c("grey70","#d62728"), guide="none") +
  labs(x=NULL, y="Spearman rho with P3b", title="Specificity: SM vs mental health") + theme_minimal()
ggsave("figures_pilot/fig_dissociation_sm_mh.png", p4, width = 6.0, height = 4.0, dpi = 300)

sessionInfo() 
