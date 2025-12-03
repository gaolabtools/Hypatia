library(tidyverse)
library(precrec)

load("data/simulation/run23_results.rda")
head(sim_data)

data_list <- list()
for (grp in unique(sim_data$cts.grp)) {
  
  grp_data <- sim_data %>% filter(cts.grp == grp)
  
  ## Pvalue    
  mm <- precrec::mmdata(scores = grp_data$X2.pval_log10, labels = grp_data$response)
  res <- precrec::evalmod(mm)
  prc_pval <- data.frame(Recall = res$prcs[[1]]$x,
                         Precision = res$prcs[[1]]$y,
                         auc = attr(res$prcs[[1]], "auc"),
                         Criteria = "P-value")
  head(prc_pval)
  
  ## Cramers V
  mm <- precrec::mmdata(scores = grp_data$X2.v, labels = grp_data$response)
  res <- precrec::evalmod(mm)
  prc_v <- data.frame(Recall = res$prcs[[1]]$x,
                      Precision = res$prcs[[1]]$y,
                      auc = attr(res$prcs[[1]], "auc"),
                      Criteria = "Cramer's V")
  head(prc_v)
  
  ## Pvalue (+ delta threshold)
  delta_threshold <- 0.10
  filt_data <- grp_data %>%
    mutate(delta.categorical = ifelse(delta >= delta_threshold, 1, 0)) %>% 
    filter(delta.categorical == 1) %>% 
    select(cts.grp, response, X2.pval_log10)
  mm <- precrec::mmdata(scores = filt_data$X2.pval_log10, labels = filt_data$response)
  res <- precrec::evalmod(mm)
  prc_pval_delta <- data.frame(Recall = res$prcs[[1]]$x,
                               Precision = res$prcs[[1]]$y,
                               auc = attr(res$prcs[[1]], "auc"),
                               Criteria = "P-value (Delta)")
  
  head(prc_pval_delta)
  
  ## Pvalue (+ Cramer's V threshold)
  v_threshold <- 0.20
  filt_data <- grp_data %>%
    mutate(v.categorical = ifelse(X2.v >= v_threshold, 1, 0)) %>% 
    filter(v.categorical == 1) %>% 
    select(cts.grp, response, X2.pval_log10)
  mm <- precrec::mmdata(scores = filt_data$X2.pval_log10, labels = filt_data$response)
  res <- precrec::evalmod(mm)
  prc_pval_v <- data.frame(Recall = res$prcs[[1]]$x,
                           Precision = res$prcs[[1]]$y,
                           auc = attr(res$prcs[[1]], "auc"),
                           Criteria = "P-value (Cramer's V)")
  
  head(prc_pval_v)
  
  ## combine
  data_list[[grp]] <- purrr::reduce(
    list(prc_pval, prc_v, prc_pval_delta, prc_pval_v), rbind
    ) %>%
    mutate(cts.grp = grp)
}

str(data_list)

plot_data <- purrr::reduce(data_list, rbind) %>%
  mutate(cts.grp = factor(cts.grp, levels = c("Very low", "Low", "Moderate", "High")),
         Criteria = factor(Criteria, 
                           levels = c("P-value", "Cramer's V", "P-value (Delta)", "P-value (Cramer's V)"))
  )

# Plot PRC
plot_data %>%
  ggplot() +
  geom_line(aes(x = Recall, y = Precision, color = Criteria),
            linewidth = 0.2) +
  facet_wrap(~ cts.grp, nrow = 1) +
  scale_color_manual(values = paletteer::paletteer_d("ggthemr::flat")[c(5, 2, 3, 4)], name = "Selection method") +
  guides(color = guide_legend(override.aes = list(linewidth = 0.5))) +
  theme_linedraw(base_size = 5) +
  theme(text = element_text(size = 5),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        strip.text = element_text(size = 5, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.key.size = unit(2, "mm")
        )

# AUCPR
plot_data %>%
  distinct(cts.grp, auc, Criteria)
