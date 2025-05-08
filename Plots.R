library(data.table); setDTthreads(4)
library(ggplot2)
library(scales)
library(patchwork)
library(dplyr)
library(yardstick)

bne_all  <- readRDS("bne_fits_and_calibrated_probs_ones.rds")$calib
iwrf_all <- readRDS("iwrf_baseline_and_pcs_oof.rds")$baseline
hs_all   <- readRDS("horseshoe.rds")

# Read results
bne_dt <- rbindlist(lapply(names(bne_all), function(ds) {
  dt <- as.data.table(bne_all[[ds]])
  dt[, `:=`(
    dataset = ds,
    method  = "BNE",
    pred    = p_cal,
    lo      = lo95,
    hi      = hi95,
    truth   = y
  )][, .(dataset, id, pred, lo, hi, truth, method)]
}))

iwrf_dt <- rbindlist(lapply(names(iwrf_all), function(ds) {
  data.table(
    dataset = ds,
    id      = seq_along(iwrf_all[[ds]]$p),
    pred    = iwrf_all[[ds]]$p,
    truth   = iwrf_all[[ds]]$truth,
    method  = "IWRF"
  )
}))

hs_dt <- rbindlist(lapply(names(hs_all), function(ds) {
  data.table(
    dataset = ds,
    id      = seq_along(hs_all[[ds]]$p),
    pred    = hs_all[[ds]]$p,
    truth   = hs_all[[ds]]$truth,
    method  = "Horseshoe"
  )
}))

# Metrics boxplots
all_dt <- rbindlist(list(iwrf_dt, bne_dt, hs_dt))
all_dt[, truth_f := factor(truth, c(0,1), c("No","Yes"))]
all_dt[, ds_f    := factor(dataset,
                           c("Statlog","HeartFailure"),
                           c("STATLOG","Heart-Failure"))]

met <- all_dt %>%
  group_by(ds_f, method) %>%
  summarise(
    acc  = accuracy_vec(truth_f, pred>0.5),
    pre  = precision_vec(truth_f, pred>0.5),
    rec  = recall_vec(truth_f, pred>0.5),
    spec = specificity_vec(truth_f, pred>0.5),
    f1   = f_meas_vec(truth_f, pred>0.5),
    auc  = roc_auc_vec(truth_f, pred),
    .groups="drop"
  ) %>%
  mutate(across(acc:auc, ~round(.x,3)))
print(met)

p1 <- ggplot(all_dt, aes(truth_f, pred)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.4, width = 0.2) +
  facet_grid(method ~ ds_f) +
  ylab("Predicted risk") + xlab("True event?") +
  theme_minimal() +
  theme(strip.text = element_text(face="bold"))
ggsave("boxplots_3x2.pdf", p1, width=8, height=10)

# Plot probability intervals
plots <- list()
for (ds in unique(all_dt$dataset)) {
  # BNE vs underlying RF
  tmp <- bne_dt[dataset==ds][order(pred)]
  tmp[, rank := .I]
  tmp[, rf_pred := iwrf_dt[dataset==ds]$pred[id]]
  p_a <- ggplot(tmp, aes(rank)) +
    geom_point(aes(y=rf_pred), shape=3, alpha=0.6) +
    geom_linerange(aes(ymin=lo, ymax=hi, colour=factor(truth))) +
    geom_point(aes(y=pred, colour=factor(truth))) +
    scale_colour_manual(values=c("0"="blue","1"="red")) +
    scale_y_continuous(labels=percent) +
    ggtitle(paste("IWRF + BNE —", ds)) +
    theme_minimal() + theme(legend.position="top")
  
  # Horseshoe
  hh <- hs_dt[dataset==ds][order(pred)]
  hh[, rank := .I]
  p_b <- ggplot(hh, aes(rank)) +
    geom_linerange(aes(ymin=lo, ymax=hi, colour=factor(truth))) +
    geom_point(aes(y=pred, colour=factor(truth))) +
    scale_colour_manual(values=c("0"="blue","1"="red")) +
    scale_y_continuous(labels=percent) +
    ggtitle(paste("Horseshoe —", ds)) +
    theme_minimal() + theme(legend.position="none")
  
  plots[[ds]] <- p_a | p_b
}

final <- wrap_plots(plots, ncol=1, guides="collect") &
  theme(legend.position="bottom")
ggsave("calib_intervals.pdf", final, width=10, height=6)
