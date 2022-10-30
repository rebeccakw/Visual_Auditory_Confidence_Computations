setwd("/Users/s4323621/Dropbox/Documents/MATLAB/confidence-master_2/new_data")
category_confidencedata <- read.csv(file = 'confidencedata_formatted.csv')

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

### CALCULATE EVIDENCE FOR SINGLE BOUNDARY VISUAL #####
vis_single_cat1_index <- category_confidencedata$stim_type == "Visual" & 
  category_confidencedata$task_type == "Single Boundary" &
  category_confidencedata$stim_orientation < taskA_visual_bounds 
vis_single_cat2_index <- category_confidencedata$stim_type == "Visual" & 
  category_confidencedata$task_type == "Single Boundary" &
  category_confidencedata$stim_orientation > taskA_visual_bounds 

category_confidencedata$evidence[vis_single_cat1_index] <- dnorm(category_confidencedata$stim_orientation[vis_single_cat1_index], 4, 5)/
  (dnorm(category_confidencedata$stim_orientation[vis_single_cat1_index], -4, 5) + 
     dnorm(category_confidencedata$stim_orientation[vis_single_cat1_index], 4, 5))

category_confidencedata$evidence[vis_single_cat2_index] <- dnorm(category_confidencedata$stim_orientation[vis_single_cat2_index], 4, 5)/
  (dnorm(category_confidencedata$stim_orientation[vis_single_cat2_index], -4, 5) + 
     dnorm(category_confidencedata$stim_orientation[vis_single_cat2_index], 4, 5))

### CALCULATE EVIDENCE FOR SINGLE BOUNDARY AUDITORY #####
aud_single_cat1_index <- category_confidencedata$stim_type == "Auditory" & 
  category_confidencedata$task_type == "Single Boundary" &
  category_confidencedata$stim_orientation < taskA_auditory_bounds 
aud_single_cat2_index <- category_confidencedata$stim_type == "Auditory" & 
  category_confidencedata$task_type == "Single Boundary" &
  category_confidencedata$stim_orientation > taskA_auditory_bounds 

category_confidencedata$evidence[aud_single_cat1_index] <- dnorm(category_confidencedata$stim_orientation[aud_single_cat1_index], 3100, 475)/
  (dnorm(category_confidencedata$stim_orientation[aud_single_cat1_index], 2300, 475) + 
     dnorm(category_confidencedata$stim_orientation[aud_single_cat1_index], 3100, 475))

category_confidencedata$evidence[aud_single_cat2_index] <- dnorm(category_confidencedata$stim_orientation[aud_single_cat2_index], 3100, 475)/
  (dnorm(category_confidencedata$stim_orientation[aud_single_cat2_index], 2300, 475) + 
     dnorm(category_confidencedata$stim_orientation[aud_single_cat2_index], 3100, 475))


### CALCULATE EVIDENCE FOR DUAL BOUNDARY VISUAL #####
vis_dual_cat1_index <- category_confidencedata$stim_type == "Visual" & 
  category_confidencedata$task_type == "Dual Boundary" &
  abs(category_confidencedata$stim_orientation) < taskB_visual_bounds[2]
vis_dual_cat2_index <- category_confidencedata$stim_type == "Visual" & 
  category_confidencedata$task_type == "Dual Boundary" &
  abs(category_confidencedata$stim_orientation) > taskB_visual_bounds[2]

category_confidencedata$evidence[vis_dual_cat1_index] <- dnorm(category_confidencedata$stim_orientation[vis_dual_cat1_index], 0, 12)/
  (dnorm(category_confidencedata$stim_orientation[vis_dual_cat1_index], 0, 3) + 
     dnorm(category_confidencedata$stim_orientation[vis_dual_cat1_index], 0, 12))

category_confidencedata$evidence[vis_dual_cat2_index] <- dnorm(category_confidencedata$stim_orientation[vis_dual_cat2_index], 0, 12)/
  (dnorm(category_confidencedata$stim_orientation[vis_dual_cat2_index], 0, 3) + 
     dnorm(category_confidencedata$stim_orientation[vis_dual_cat2_index], 0, 12))


### CALCULATE EVIDENCE FOR DUAL BOUNDARY AUDITORY #####
aud_dual_cat1_index <- category_confidencedata$stim_type == "Auditory" & 
  category_confidencedata$task_type == "Dual Boundary" &
  category_confidencedata$stim_orientation > taskB_auditory_bounds[1] & 
  category_confidencedata$stim_orientation < taskB_auditory_bounds[2]
aud_dual_cat2_index <- category_confidencedata$stim_type == "Auditory" & 
  category_confidencedata$task_type == "Dual Boundary" &
  category_confidencedata$stim_orientation < taskB_auditory_bounds[1] | 
  category_confidencedata$stim_orientation > taskB_auditory_bounds[2]

category_confidencedata$evidence[aud_dual_cat1_index] <- dnorm(category_confidencedata$stim_orientation[aud_dual_cat1_index], 2700, 500)/
  (dnorm(category_confidencedata$stim_orientation[aud_dual_cat1_index], 2700, 125) + 
     dnorm(category_confidencedata$stim_orientation[aud_dual_cat1_index], 2700, 500))

category_confidencedata$evidence[aud_dual_cat2_index] <- dnorm(category_confidencedata$stim_orientation[aud_dual_cat2_index], 2700, 500)/
  (dnorm(category_confidencedata$stim_orientation[aud_dual_cat2_index], 2700, 125) + 
     dnorm(category_confidencedata$stim_orientation[aud_dual_cat2_index], 2700, 500))


category_confidencedata$subject_name <- as.factor(category_confidencedata$subject_name)
bin_size = 36
bin_numbers = 160

## MODEL SINGLE BOUDARY VISUAL ###
category_confidencedata_vis_single <- category_confidencedata %>%
  filter(task_type == "Single Boundary", stim_type == "Visual")

category_confidencedata_vis_single <- category_confidencedata_vis_single %>%
  mutate(m_stim_rel = scale(stim_reliability),
         m_evidence = scale(evidence), 
         int = m_stim_rel*m_evidence)

model_vis_single = glmer((resp_category - 1) ~ 1 + m_stim_rel + m_evidence + int + (1 + m_stim_rel + m_evidence + int | subject_name), 
                        data = category_confidencedata_vis_single, family =  "binomial")

category_confidencedata_vis_single <- data.frame(category_confidencedata_vis_single, 
                                        pred_group = predict(model_vis_single, re.form = NA, type = "response"),
                                        pred_individ = predict(model_vis_single, type = "response"))

category_confidencedata_vis_single_binned <- category_confidencedata_vis_single %>%
  group_by(stim_reliability_level) %>%
  arrange(evidence, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(order = rep(1:bin_numbers, each = bin_size)) %>%
  group_by(stim_reliability_level, order) %>%
  mutate(bin = mean(evidence),
         mean_category = mean(resp_category - 1),
         mean_category = mean(resp_category),
         mean_resp = mean(resp_buttonid),
         no = n(),
         std_error = sd(resp_buttonid)/sqrt(8))


## MODEL SINGLE BOUDARY AUDITORY ###
category_confidencedata_aud_single <- category_confidencedata %>%
  filter(task_type == "Single Boundary", stim_type == "Auditory")

category_confidencedata_aud_single <- category_confidencedata_aud_single %>%
  mutate(m_stim_rel = scale(stim_reliability),
         m_evidence = scale(evidence), 
         int = m_stim_rel*m_evidence)

model_aud_single = glmer((resp_category - 1) ~ 1 + m_stim_rel + m_evidence + int + (1 + m_stim_rel + m_evidence + int | subject_name), 
                        data = category_confidencedata_aud_single, family =  "binomial")

category_confidencedata_aud_single <- data.frame(category_confidencedata_aud_single, 
                                        pred_group = predict(model_aud_single, re.form = NA, type = "response"),
                                        pred_individ = predict(model_aud_single, type = "response"))

category_confidencedata_aud_single_binned <- category_confidencedata_aud_single %>%
  group_by(stim_reliability_level) %>%
  arrange(evidence, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(order = rep(1:bin_numbers, each = bin_size)) %>%
  group_by(stim_reliability_level, order) %>%
  mutate(bin = mean(evidence),
         mean_category = mean(resp_category - 1),
         mean_category = mean(resp_category),
         mean_resp = mean(resp_buttonid),
         no = n(),
         std_error = sd(resp_buttonid)/sqrt(8))


## MODEL DUAL BOUDARY VISUAL ###
category_confidencedata_vis_dual <- category_confidencedata %>%
  filter(task_type == "Dual Boundary", stim_type == "Visual")

category_confidencedata_vis_dual <- category_confidencedata_vis_dual %>%
  mutate(m_stim_rel = scale(stim_reliability),
         m_evidence = scale(evidence), 
         int = m_stim_rel*m_evidence)

model_vis_dual = glmer((resp_category - 1) ~ 1 + m_stim_rel + m_evidence + int + (1 + m_stim_rel + m_evidence + int | subject_name), 
                      data = category_confidencedata_vis_dual, family =  "binomial")

category_confidencedata_vis_dual <- data.frame(category_confidencedata_vis_dual, 
                                      pred_group = predict(model_vis_dual, re.form = NA, type = "response"),
                                      pred_individ = predict(model_vis_dual, type = "response"))

category_confidencedata_vis_dual_binned <- category_confidencedata_vis_dual %>%
  group_by(stim_reliability_level) %>%
  arrange(evidence, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(order = rep(1:bin_numbers, each = bin_size)) %>%
  group_by(stim_reliability_level, order) %>%
  mutate(bin = mean(evidence),
         mean_category = mean(resp_category - 1),
         mean_category = mean(resp_category),
         mean_resp = mean(resp_buttonid),
         no = n(),
         std_error = sd(resp_buttonid)/sqrt(8))


## MODEL DUAL BOUDARY AUDITORY ###
category_confidencedata_aud_dual <- category_confidencedata %>%
  filter(task_type == "Dual Boundary", stim_type == "Auditory")

category_confidencedata_aud_dual <- category_confidencedata_aud_dual %>%
  mutate(m_stim_rel = scale(stim_reliability),
         m_evidence = scale(evidence), 
         int = m_stim_rel*m_evidence)

model_aud_dual = glmer((resp_category - 1) ~ 1 + m_stim_rel + m_evidence + int + (1 + m_stim_rel + m_evidence + int | subject_name), 
                      data = category_confidencedata_aud_dual, family =  "binomial")

category_confidencedata_aud_dual <- data.frame(category_confidencedata_aud_dual, 
                                      pred_group = predict(model_aud_dual, re.form = NA, type = "response"),
                                      pred_individ = predict(model_aud_dual, type = "response"))

category_confidencedata_aud_dual_binned <- category_confidencedata_aud_dual %>%
  group_by(stim_reliability_level) %>%
  arrange(evidence, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(order = rep(1:bin_numbers, each = bin_size)) %>%
  group_by(stim_reliability_level, order) %>%
  mutate(bin = mean(evidence),
         mean_category = mean(resp_category - 1),
         mean_category = mean(resp_category),
         mean_resp = mean(resp_buttonid),
         no = n(),
         std_error = sd(resp_buttonid)/sqrt(8))


category_all_binned_data <- rbind(category_confidencedata_vis_single_binned,
                         category_confidencedata_aud_single_binned,
                         category_confidencedata_vis_dual_binned,
                         category_confidencedata_aud_dual_binned)

category_all_binned_data$stim_type <- factor(category_all_binned_data$stim_type, levels = c("Visual", "Auditory"))
category_all_binned_data$task_type <- factor(category_all_binned_data$task_type, levels = c("Single Boundary", "Dual Boundary"))

category_plots <- ggplot(data = category_all_binned_data) +
  geom_point(aes(x = bin, y = mean_category, col =as.factor(stim_reliability_level)), size = 1.5,
             alpha = .8) + 
  geom_line(size = 2, aes(x = evidence, y = pred_group + 1, col = as.factor(stim_reliability_level))) + 
  viridis::scale_color_viridis(discrete = TRUE, begin = 0.95, end = 0) +
  facet_grid(~task_type + stim_type, scales = 'free') +
  xlab("Evidence for Category 2") + ylab("Category") +
  labs(colour = "Stimulus Reliability") +
  theme(axis.text=element_text(size=15), axis.title=element_text(size = 18), strip.text =element_text(size= 0.5),
        legend.position = "none") + ylim(1, 2)
