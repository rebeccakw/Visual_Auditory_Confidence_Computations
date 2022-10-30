setwd("/Users/s4323621/Dropbox/Documents/multimodal_confidence/recovery_data")
files = list.files("/Users/s4323621/Dropbox/Documents/multimodal_confidence/recovery_data")
files = files[files != "correlations"]
for (f in files) {
  data = read.csv(f)
  result = data %>%
    group_by(parameter) %>%
    summarise(cor(generating, fitted))
  write.csv(result, paste0("correlations_", f))
}
