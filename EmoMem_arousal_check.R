library(readxl)
library(readr)
library(dplyr)
library(stringr)
library(lme4)
library(ggplot2)
library(patchwork)

trial_data <- read_csv("//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/memory_variables/Malo_Mem_cc280_trialwise.csv")
head(trial_data)

txt_data <- read.table("//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/EM_CC520168_140507_studyres1.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(txt_data) <- c("Col1", "Stimulus_File", "Object_File", "Col4", "Col5", "Col6")

# Extract stimulus ID from the filename using regex
stim_labels <- gsub("\\D", "", txt_data$Stimulus_File)
trial_data$stim_id <- rep(stim_labels, length(unique(trial_data$CCID)))

descriptor_data <- read_excel("//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/IAPS_ratings/Libkuman_2007_Ratings.xls")
descriptor_data$stim_id <- as.character(as.integer(descriptor_data$`1APS#`))

df_merged <- merge(
  trial_data,
  descriptor_data[, c("stim_id", "Libkuman_val_M", "Libkuman_arousal_M")],
  by = "stim_id",
  all.x = TRUE
)


df_merged$Age_z <- scale(df_merged$Age)
df_merged$Libkuman_val_M_z <- scale(df_merged$Libkuman_val_M)
df_merged$Libkuman_arousal_M_z <- scale(df_merged$Libkuman_arousal_M)



model_all <- glmer(Correct_Recall ~ Libkuman_val_M_z * Age_z + (1 | CCID) + (1 | stim_id), data = df_merged, family = binomial)
summary(model_all)


model_int <- glmer(Correct_Recall ~ Libkuman_val_M_z * Age_z * Libkuman_arousal_M_z + (1 | CCID) + (1 | stim_id), data = df_merged, family = binomial)
summary(model_int)


df_low_arousal <- df_merged %>%
  group_by(Valence) %>%
  mutate(arousal_median = median(Libkuman_arousal_M, na.rm = TRUE)) %>%
  filter(Libkuman_arousal_M <= arousal_median) %>%
  ungroup()

model_low <- glmer(Correct_Recall ~ Libkuman_val_M_z * Age_z + (1 | CCID) + (1 | stim_id), data = df_low_arousal, family = binomial)
summary(model_low)

df_high_arousal <- df_merged %>%
  group_by(Valence) %>%
  mutate(arousal_median = median(Libkuman_arousal_M, na.rm = TRUE)) %>%
  filter(Libkuman_arousal_M >= arousal_median) %>%
  ungroup()

model_high <- glmer(Correct_Recall ~ Libkuman_val_M_z * Age_z + (1 | CCID) + (1 | stim_id), data = df_high_arousal, family = binomial)
summary(model_high)








