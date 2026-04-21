library(readxl)
library(ggplot2)
library(dplyr)
library(broom)

ratings <- read.delim("//cbsu/data/Group/Camcan/sandbox/cc07/Project3_IAPS_RSA/Gruhn-BRM-2008a/Gruhn-BRM-2008/Gruhn-Scheibe_2008_PictureData.txt")


# plot RDMs
o_RDM <- outer(ratings$val_o, ratings$val_o, "-")
y_RDM <- outer(ratings$val_y, ratings$val_y, "-")

o_RDM <- as.data.frame(as.table(o_RDM))
y_RDM <- as.data.frame(as.table(y_RDM))

colnames(o_RDM) <- c("Row", "Col", "Value")
colnames(y_RDM) <- c("Row", "Col", "Value")



ggplot(y_RDM, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "yellow",
    high = "red",
    midpoint = 0
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) +
  labs(
    title = "Heatmap of vRDM, young group",
    fill = "Value"
  )

ggplot(o_RDM, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "yellow",
    high = "red",
    midpoint = 0
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) +
  labs(
    title = "Heatmap of vRDM, old group",
    fill = "Value"
  )


#get indices of neutral, negative, and positive ratings
neg_df <- data.frame(index = which(ratings$val_y >= 1 & ratings$val_y <= 3.4), valence_y = ratings$val_y[which(ratings$val_y >= 1 & ratings$val_y <= 3.4)], valence_o = ratings$val_o[which(ratings$val_y >= 1 & ratings$val_y <= 3.4)])
neu_df <- data.frame(index = which(ratings$val_y >= 3.5 & ratings$val_y <= 6.4), valence_y = ratings$val_y[which(ratings$val_y >= 3.5 & ratings$val_y <= 6.4)], valence_o = ratings$val_o[which(ratings$val_y >= 3.5 & ratings$val_y <= 6.4)])
pos_df <- data.frame(index = which(ratings$val_y >= 6.5 & ratings$val_y <= 9), valence_y = ratings$val_y[which(ratings$val_y >= 6.5 & ratings$val_y <= 9)], valence_o = ratings$val_o[which(ratings$val_y >= 6.5 & ratings$val_y <= 9)])



#I checked, none of these have missing values here
dmneg_o <- outer(neg_df$valence_o, neg_df$valence_o, "-")
dmneu_o <- outer(neu_df$valence_o, neu_df$valence_o, "-")
dmpos_o <- outer(pos_df$valence_o, pos_df$valence_o, "-")
dmnegneu_o <- outer(neg_df$valence_o, neu_df$valence_o, "-")
dmnegpos_o <- outer(neu_df$valence_o, pos_df$valence_o, "-")
dmneupos_o <- outer(neu_df$valence_o, pos_df$valence_o, "-")
old_negdv <- log(abs(dmneg_o[upper.tri(dmneg_o)]) + 1)
old_neudv <- log(abs(dmneu_o[upper.tri(dmneu_o)]) + 1)
old_posdv <- log(abs(dmpos_o[upper.tri(dmpos_o)]) + 1)
old_negneudv <- log(abs(dmnegneu_o[upper.tri(dmnegneu_o)]) + 1)
old_negposdv <- log(abs(dmnegpos_o[upper.tri(dmnegpos_o)]) + 1)
old_neuposdv <- log(abs(dmneupos_o[upper.tri(dmneupos_o)]) + 1)


dmneg_y <- outer(neg_df$valence_y, neg_df$valence_y, "-")
dmneu_y <- outer(neu_df$valence_y, neu_df$valence_y, "-")
dmpos_y <- outer(pos_df$valence_y, pos_df$valence_y, "-")
dmnegneu_y <- outer(neg_df$valence_y, neu_df$valence_y, "-")
dmnegpos_y <- outer(neu_df$valence_y, pos_df$valence_y, "-")
dmneupos_y <- outer(neu_df$valence_y, pos_df$valence_y, "-")
young_negdv <- log(abs(dmneg_y[upper.tri(dmneg_y)]) + 1)
young_neudv <- log(abs(dmneu_y[upper.tri(dmneu_y)]) + 1)
young_posdv <- log(abs(dmpos_y[upper.tri(dmpos_y)]) + 1)
young_negneudv <- log(abs(dmnegneu_y[upper.tri(dmnegneu_y)]) + 1)
young_negposdv <- log(abs(dmnegpos_y[upper.tri(dmnegpos_y)]) + 1)
young_neuposdv <- log(abs(dmneupos_y[upper.tri(dmneupos_y)]) + 1)


##
MNeg <- lm(old_negdv ~ 0 + young_negdv)
summary(MNeg)
confint(MNeg, level=0.95)

MNeu <- lm(old_neudv ~ 0 + young_neudv)
summary(MNeu)
confint(MNeu, level=0.95)

MPos <- lm(old_posdv ~ 0 + young_posdv)
summary(MPos)
confint(MPos, level=0.95)

models <- list(
  Neg = MNeg,
  Neu = MNeu,
  Pos = MPos
)

results <- lapply(names(models), function(name) {
  model <- models[[name]]
  
  tidy(model, conf.int = TRUE) %>%
    mutate(
      Model = name,
      df = df.residual(model)
    )
}) %>%
  bind_rows()

results %>%
  mutate(
    t_value_1 = (estimate - 1) / std.error,
    p_value_1 = 2 * pt(abs(t_value_1), df = df, lower.tail = FALSE)
  ) %>%
  mutate(
    estimate = round(estimate, 2),
    std.error = round(std.error, 2),
    conf.low = round(conf.low, 2),
    conf.high = round(conf.high, 2),
    p_value_1 = ifelse(p_value_1 < .001, "< .001", round(p_value_1, 3))
  )




## check across valences: 
MNegneu <- lm(old_negneudv ~ 0 + young_negneudv)
summary(MNegneu)

MNegpos <- lm(old_negposdv ~ 0 + young_negposdv)
summary(MNegpos)

MNeupos <- lm(old_neuposdv ~ 0 + young_neuposdv)
summary(MNeupos)


## check without log:
onegdv <- abs(dmneg_o[upper.tri(dmneg_o)])
oneudv <- abs(dmneu_o[upper.tri(dmneu_o)])
oposdv <- abs(dmpos_o[upper.tri(dmpos_o)])
ynegdv <- abs(dmneg_y[upper.tri(dmneg_y)])
yneudv <- abs(dmneu_y[upper.tri(dmneu_y)])
yposdv <- abs(dmpos_y[upper.tri(dmpos_y)])

MNegR <- lm(onegdv ~ 0 + ynegdv)
summary(MNegR)
confint(MNegR, level=0.95)

MNeuR <- lm(oneudv ~ 0 + yneudv)
summary(MNeuR)
confint(MNeuR, level=0.95)

MPosR <- lm(oposdv ~ 0 + yposdv)
summary(MPosR)
confint(MPosR, level=0.95)


## check directly for the ratings:
MNeg2 <- lm(log(valence_o) ~ 0 + log(valence_y), neg_df)
summary(MNeg2)

MNeu2 <- lm(log(valence_o) ~ 0 + log(valence_y), neu_df)
summary(MNeu2)

MPos2 <- lm(log(valence_o) ~ 0 + log(valence_y), pos_df)
summary(MPos2)

