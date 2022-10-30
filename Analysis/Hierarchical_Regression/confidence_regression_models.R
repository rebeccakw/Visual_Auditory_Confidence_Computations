setwd("/Users/rebeccawest/Dropbox/Documents/MATLAB/confidence-master_2/new_data")
confidencedata <- read.csv(file = 'confidencedata_formatted.csv')

### FUNCTION #####
overlap_function <- function(mu1,sigma1,mu2,sigma2) {
  c <- (0*sigma1^2 - sigma2*(0*sigma2 + sigma1*sqrt((0 - 0)^2 + 
                                                      2*(sigma1^2 - sigma2^2)*log(sigma1/sigma2))))/(sigma1^2 - sigma2^2)
  lowerbound = mu1 - c;
  upperbound = mu1 + c;
  bounds <- c(upperbound, lowerbound)
  return(bounds)
}

### GET BOUNDARIES #####
taskA_visual_bounds <- 0
taskA_auditory_bounds <- 2700
taskB_visual_bounds <- overlap_function(0, 12, 0, 3)
taskB_auditory_bounds <- overlap_function(2700, 500, 2700, 125)

### CALCULATE EVIDENCE FOR SINGLE BOUNDARY #####
confidencedata$truecat_sd[confidencedata$stim_type == "Visual" & 
                            confidencedata$task_type == "Single Boundary"] <- 5
confidencedata$truecat_mean[confidencedata$stim_type == "Visual" & 
                              confidencedata$task_type == "Single Boundary" &
                              confidencedata$stim_orientation < taskA_visual_bounds] <- -4
confidencedata$truecat_mean[confidencedata$stim_type == "Visual" & 
                              confidencedata$task_type == "Single Boundary" &
                              confidencedata$stim_orientation > taskA_visual_bounds] <- 4

confidencedata$truecat_sd[confidencedata$stim_type == "Auditory" & 
                            confidencedata$task_type == "Single Boundary"] <- 475
confidencedata$truecat_mean[confidencedata$stim_type == "Auditory" & 
                              confidencedata$task_type == "Single Boundary" &
                              confidencedata$stim_orientation < taskA_auditory_bounds] <- 2300
confidencedata$truecat_mean[confidencedata$stim_type == "Auditory" & 
                              confidencedata$task_type == "Single Boundary" &
                              confidencedata$stim_orientation > taskA_auditory_bounds] <- 3100

### CALCULATE EVIDENCE FOR DUAL BOUNDARY #####
confidencedata$truecat_mean[confidencedata$stim_type == "Visual" & 
                              confidencedata$task_type == "Dual Boundary"] <- 0
confidencedata$truecat_sd[confidencedata$stim_type == "Visual" & 
                            confidencedata$task_type == "Dual Boundary" &
                            abs(confidencedata$stim_orientation) < taskB_visual_bounds[2]] <- 3
confidencedata$truecat_sd[confidencedata$stim_type == "Visual" & 
                            confidencedata$task_type == "Dual Boundary" &
                            abs(confidencedata$stim_orientation) > taskB_visual_bounds[2]] <- 12


confidencedata$truecat_mean[confidencedata$stim_type == "Auditory" & 
                              confidencedata$task_type == "Dual Boundary"] <- 2700
confidencedata$truecat_sd[confidencedata$stim_type == "Auditory" & 
                            confidencedata$task_type == "Dual Boundary" &
                            confidencedata$stim_orientation > taskB_auditory_bounds[1] & 
                            confidencedata$stim_orientation < taskB_auditory_bounds[2]] <- 125
confidencedata$truecat_sd[confidencedata$stim_type == "Auditory" & 
                            confidencedata$task_type == "Dual Boundary" &
                            confidencedata$stim_orientation < taskB_auditory_bounds[1] | 
                            confidencedata$stim_orientation > taskB_auditory_bounds[2]] <- 500


### CALCULATE EVIDENCE FOR SINGLE BOUNDARY VISUAL #####
vis_single_cat1_index <- confidencedata$stim_type == "Visual" & 
  confidencedata$task_type == "Single Boundary" &
  confidencedata$stim_orientation < taskA_visual_bounds 
vis_single_cat2_index <- confidencedata$stim_type == "Visual" & 
  confidencedata$task_type == "Single Boundary" &
  confidencedata$stim_orientation > taskA_visual_bounds 

confidencedata$evidence[vis_single_cat1_index] <- dnorm(confidencedata$stim_orientation[vis_single_cat1_index], -4, 5)/
  (dnorm(confidencedata$stim_orientation[vis_single_cat1_index], -4, 5) + 
     dnorm(confidencedata$stim_orientation[vis_single_cat1_index], 4, 5))

confidencedata$evidence[vis_single_cat2_index] <- dnorm(confidencedata$stim_orientation[vis_single_cat2_index], 4, 5)/
  (dnorm(confidencedata$stim_orientation[vis_single_cat2_index], -4, 5) + 
     dnorm(confidencedata$stim_orientation[vis_single_cat2_index], 4, 5))

### CALCULATE EVIDENCE FOR SINGLE BOUNDARY AUDITORY #####
aud_single_cat1_index <- confidencedata$stim_type == "Auditory" & 
  confidencedata$task_type == "Single Boundary" &
  confidencedata$stim_orientation < taskA_auditory_bounds 
aud_single_cat2_index <- confidencedata$stim_type == "Auditory" & 
  confidencedata$task_type == "Single Boundary" &
  confidencedata$stim_orientation > taskA_auditory_bounds 

confidencedata$evidence[aud_single_cat1_index] <- dnorm(confidencedata$stim_orientation[aud_single_cat1_index], 2300, 475)/
  (dnorm(confidencedata$stim_orientation[aud_single_cat1_index], 2300, 475) + 
     dnorm(confidencedata$stim_orientation[aud_single_cat1_index], 3100, 475))

confidencedata$evidence[aud_single_cat2_index] <- dnorm(confidencedata$stim_orientation[aud_single_cat2_index], 3100, 475)/
  (dnorm(confidencedata$stim_orientation[aud_single_cat2_index], 2300, 475) + 
     dnorm(confidencedata$stim_orientation[aud_single_cat2_index], 3100, 475))


### CALCULATE EVIDENCE FOR DUAL BOUNDARY VISUAL #####
vis_dual_cat1_index <- confidencedata$stim_type == "Visual" & 
  confidencedata$task_type == "Dual Boundary" &
  abs(confidencedata$stim_orientation) < taskB_visual_bounds[2]
vis_dual_cat2_index <- confidencedata$stim_type == "Visual" & 
  confidencedata$task_type == "Dual Boundary" &
  abs(confidencedata$stim_orientation) > taskB_visual_bounds[2]

confidencedata$evidence[vis_dual_cat1_index] <- dnorm(confidencedata$stim_orientation[vis_dual_cat1_index], 0, 3)/
  (dnorm(confidencedata$stim_orientation[vis_dual_cat1_index], 0, 3) + 
     dnorm(confidencedata$stim_orientation[vis_dual_cat1_index], 0, 12))

confidencedata$evidence[vis_dual_cat2_index] <- dnorm(confidencedata$stim_orientation[vis_dual_cat2_index], 0, 12)/
  (dnorm(confidencedata$stim_orientation[vis_dual_cat2_index], 0, 3) + 
     dnorm(confidencedata$stim_orientation[vis_dual_cat2_index], 0, 12))


### CALCULATE EVIDENCE FOR DUAL BOUNDARY AUDITORY #####
aud_dual_cat1_index <- confidencedata$stim_type == "Auditory" & 
  confidencedata$task_type == "Dual Boundary" &
  confidencedata$stim_orientation > taskB_auditory_bounds[1] & 
  confidencedata$stim_orientation < taskB_auditory_bounds[2]
aud_dual_cat2_index <- confidencedata$stim_type == "Auditory" & 
  confidencedata$task_type == "Dual Boundary" &
  confidencedata$stim_orientation < taskB_auditory_bounds[1] | 
  confidencedata$stim_orientation > taskB_auditory_bounds[2]

confidencedata$evidence[aud_dual_cat1_index] <- dnorm(confidencedata$stim_orientation[aud_dual_cat1_index], 2700, 125)/
  (dnorm(confidencedata$stim_orientation[aud_dual_cat1_index], 2700, 125) + 
     dnorm(confidencedata$stim_orientation[aud_dual_cat1_index], 2700, 500))

confidencedata$evidence[aud_dual_cat2_index] <- dnorm(confidencedata$stim_orientation[aud_dual_cat2_index], 2700, 500)/
  (dnorm(confidencedata$stim_orientation[aud_dual_cat2_index], 2700, 125) + 
     dnorm(confidencedata$stim_orientation[aud_dual_cat2_index], 2700, 500))

confidencedata$subject_name <- as.factor(confidencedata$subject_name)
bin_size = 36
bin_numbers = 160

## MODEL SINGLE BOUDARY VISUAL ###
confidencedata_vis_single <- confidencedata %>%
  filter(task_type == "Single Boundary", stim_type == "Visual")

confidencedata_vis_single <- confidencedata_vis_single %>%
  mutate(m_stim_rel = scale(stim_reliability),
         m_evidence = scale(evidence), 
         int = m_stim_rel*m_evidence)

model_vis_single = lmer(resp_confidence ~ 1 + m_stim_rel + m_evidence + int + (1 + m_stim_rel + m_evidence + int | subject_name), 
                        data = confidencedata_vis_single)

confidencedata_vis_single <- data.frame(confidencedata_vis_single, 
                                        pred_group = predict(model_vis_single, re.form = NA),
                                        pred_individ = predict(model_vis_single))

confidencedata_vis_single_binned <- confidencedata_vis_single %>%
  group_by(stim_reliability_level) %>%
  arrange(evidence, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(order = rep(1:bin_numbers, each = bin_size)) %>%
  group_by(stim_reliability_level, order) %>%
  mutate(bin = mean(evidence),
            mean_category = mean(resp_category - 1),
            mean_confidence = mean(resp_confidence),
            mean_resp = mean(resp_buttonid),
            no = n(),
            std_error = sd(resp_buttonid)/sqrt(8))


## MODEL SINGLE BOUDARY AUDITORY ###
confidencedata_aud_single <- confidencedata %>%
  filter(task_type == "Single Boundary", stim_type == "Auditory")

confidencedata_aud_single <- confidencedata_aud_single %>%
  mutate(m_stim_rel = scale(stim_reliability),
         m_evidence = scale(evidence), 
         int = m_stim_rel*m_evidence)

model_aud_single = lmer(resp_confidence ~ 1 + m_stim_rel + m_evidence + int + (1 + m_stim_rel + m_evidence + int | subject_name), 
                        data = confidencedata_aud_single)

confidencedata_aud_single <- data.frame(confidencedata_aud_single, 
                                        pred_group = predict(model_aud_single, re.form = NA),
                                        pred_individ = predict(model_aud_single))

confidencedata_aud_single_binned <- confidencedata_aud_single %>%
  group_by(stim_reliability_level) %>%
  arrange(evidence, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(order = rep(1:bin_numbers, each = bin_size)) %>%
  group_by(stim_reliability_level, order) %>%
  mutate(bin = mean(evidence),
         mean_category = mean(resp_category - 1),
         mean_confidence = mean(resp_confidence),
         mean_resp = mean(resp_buttonid),
         no = n(),
         std_error = sd(resp_buttonid)/sqrt(8))


## MODEL DUAL BOUDARY VISUAL ###
confidencedata_vis_dual <- confidencedata %>%
  filter(task_type == "Dual Boundary", stim_type == "Visual")

confidencedata_vis_dual <- confidencedata_vis_dual %>%
  mutate(m_stim_rel = scale(stim_reliability),
         m_evidence = scale(evidence), 
         int = m_stim_rel*m_evidence)

model_vis_dual = lmer(resp_confidence ~ 1 + m_stim_rel + m_evidence + int + (1 + m_stim_rel + m_evidence + int | subject_name), 
                        data = confidencedata_vis_dual)

confidencedata_vis_dual <- data.frame(confidencedata_vis_dual, 
                                        pred_group = predict(model_vis_dual, re.form = NA),
                                        pred_individ = predict(model_vis_dual))

confidencedata_vis_dual_binned <- confidencedata_vis_dual %>%
  group_by(stim_reliability_level) %>%
  arrange(evidence, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(order = rep(1:bin_numbers, each = bin_size)) %>%
  group_by(stim_reliability_level, order) %>%
  mutate(bin = mean(evidence),
         mean_category = mean(resp_category - 1),
         mean_confidence = mean(resp_confidence),
         mean_resp = mean(resp_buttonid),
         no = n(),
         std_error = sd(resp_buttonid)/sqrt(8))


## MODEL DUAL BOUDARY AUDITORY ###
confidencedata_aud_dual <- confidencedata %>%
  filter(task_type == "Dual Boundary", stim_type == "Auditory")

confidencedata_aud_dual <- confidencedata_aud_dual %>%
  mutate(m_stim_rel = scale(stim_reliability),
         m_evidence = scale(evidence), 
         int = m_stim_rel*m_evidence)

model_aud_dual = lmer(resp_confidence ~ 1 + m_stim_rel + m_evidence + int + (1 + m_stim_rel + m_evidence + int | subject_name), 
                      data = confidencedata_aud_dual)

confidencedata_aud_dual <- data.frame(confidencedata_aud_dual, 
                                      pred_group = predict(model_aud_dual, re.form = NA),
                                      pred_individ = predict(model_aud_dual))

confidencedata_aud_dual_binned <- confidencedata_aud_dual %>%
  group_by(stim_reliability_level) %>%
  arrange(evidence, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(order = rep(1:bin_numbers, each = bin_size)) %>%
  group_by(stim_reliability_level, order) %>%
  mutate(bin = mean(evidence),
         mean_category = mean(resp_category - 1),
         mean_confidence = mean(resp_confidence),
         mean_resp = mean(resp_buttonid),
         no = n(),
         std_error = sd(resp_buttonid)/sqrt(8))


all_binned_data <- rbind(confidencedata_vis_single_binned,
      confidencedata_aud_single_binned,
      confidencedata_vis_dual_binned,
      confidencedata_aud_dual_binned)

all_binned_data$stim_type <- factor(all_binned_data$stim_type, levels = c("Visual", "Auditory"))
all_binned_data$task_type <- factor(all_binned_data$task_type, levels = c("Single Boundary", "Dual Boundary"))

confidence_plots <- ggplot(data = all_binned_data) +
  geom_point(aes(x = bin, y = mean_confidence, col =as.factor(stim_reliability_level)), size = 1.5,
             alpha = .8) + 
  geom_line(size = 2, aes(x = evidence, y = pred_group, col = as.factor(stim_reliability_level))) + 
  viridis::scale_color_viridis(discrete = TRUE, begin = 0.95, end = 0) +
  facet_grid(~task_type + stim_type, scales = 'free') +
  xlab("Evidence for Given Category") + ylab("Confidence") +
  labs(colour = "Stimulus Reliability") +
  theme(axis.text=element_text(size=15), axis.title=element_text(size = 18), strip.text =element_text(size= 0.5),
        legend.position = "none") + ylim(0.9, 4)
