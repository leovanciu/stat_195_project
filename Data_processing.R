library(data.table)

STATLOG_PATH <- ""
stat <- fread(STATLOG_PATH, header = FALSE)
setnames(stat, c(
  'age','sex','cp','trestbps','chol','fbs','restecg',
  'thalach','exang','oldpeak','slope','ca','thal','raw_target'
))
stat[, target := as.integer(raw_target > 1)]
stat[, raw_target := NULL]
OUT_PATH <- ""
fwrite(stat, OUT_PATH)