library(tidyverse)
library(patchwork)

results <- read.table("data/simulation/run23_results.tsv")

## diversity functions
div.shan <- function(x) {
  x <- head(sort(x, decreasing = TRUE), 2)
  -sum(x[x > 0] * log(x[x > 0]))
}

div.renyi <- function(x, renyi.alpha = 2) {
  x <- head(sort(x, decreasing = TRUE), 2)
  (1 / (1 - renyi.alpha)) * log( sum( (x[x > 0])^renyi.alpha ) )
}

div.gini.simp <- function(x) {
  1 - sum( (x[x > 0])^2 )
}

div.inv.simp <- function(x) {
  1 / sum( (x[x > 0])^2 )
}

div.tsallis <- function(x, tsallis.q = 3) {
  (1 - sum(x[x > 0]^tsallis.q)) / (tsallis.q - 1)
}

div_res <- results %>%
  select(testID, prop1, prop2) %>%
  group_by(testID) %>%
  mutate(n_iso = n(),
         n_iso = factor(n_iso, levels = c(1:10)),
         mdt1 = max(prop1),
         mdt2 = max(prop2),
         div1_sh = div.shan(prop1),
         div2_sh = div.shan(prop2),
         div1_ren = div.renyi(prop1),
         div2_ren = div.renyi(prop2),
         div1_gs = div.gini.simp(prop1),
         div2_gs = div.gini.simp(prop2),
         div1_is = div.inv.simp(prop1),
         div2_is = div.inv.simp(prop2),
         div1_ts = div.tsallis(prop1),
         div2_ts = div.tsallis(prop2)
         ) %>%
  ungroup() %>%
  select(-c(prop1, prop2)) %>%
  distinct()

## Cutoffs
div_res %>% filter(mdt1 == 0.8) %>% pull(div1_sh) %>% max()
div_res %>% filter(mdt1 == 0.8) %>% pull(div1_ren) %>% max()
div_res %>% filter(mdt1 == 0.8) %>% pull(div1_gs) %>% max()
div_res %>% filter(mdt1 == 0.8) %>% pull(div1_is) %>% max()
div_res %>% filter(mdt1 == 0.8) %>% pull(div1_ts) %>% max()

## Prepare plot data
sh_mdt <- bind_rows(div_res %>% select(x = mdt1, div = div1_sh, n_iso),
                    div_res %>% select(x = mdt2, div = div2_sh, n_iso)) %>%
  mutate("Class" = ifelse(div <= 0.5, "Monoform", "Polyform"),
         "index" = "Shannon",
         "cutoff" = 0.500)

ren_mdt <- bind_rows(div_res %>% select(x = mdt1, div = div1_ren, n_iso),
                     div_res %>% select(x = mdt2, div = div2_ren, n_iso)) %>%
  mutate("Class" = ifelse(div <= .435, "Monoform", "Polyform"),
         "index" = "Renyi",
         "cutoff" = .435)

gs_mdt <- bind_rows(div_res %>% select(x = mdt1, div = div1_gs, n_iso),
                    div_res %>% select(x = mdt2, div = div2_gs, n_iso)) %>%
  mutate("Class" = ifelse(div <= 0.348, "Monoform", "Polyform"),
         "index" = "Gini-Simpson",
         "cutoff" = 0.348)

is_mdt <- bind_rows(div_res %>% select(x = mdt1, div = div1_is, n_iso),
                    div_res %>% select(x = mdt2, div = div2_is, n_iso)) %>%
  mutate("Class" = ifelse(div <= 1.533, "Monoform", "Polyform"),
         "index" = "Inverse-Simpson",
         "cutoff" = 1.533)

ts_mdt <- bind_rows(div_res %>% select(x = mdt1, div = div1_ts, n_iso),
                    div_res %>% select(x = mdt2, div = div2_ts, n_iso)) %>%
  mutate("Class" = ifelse(div <= 0.243, "Monoform", "Polyform"),
         "index" = "Tsallis",
         "cutoff" = 0.243)

mdtdata <- rbind(sh_mdt, ren_mdt, gs_mdt, is_mdt, ts_mdt)

## Correlation
r_sh <- sh_mdt %>% 
  group_by(Class) %>%
  summarize(r = cor.test(x = 1 - x, y = div, method = "pearson")$estimate) %>%
  rbind(data.frame("Class" = "All",
                   "r" = cor.test(x = 1 - sh_mdt$x, y = sh_mdt$div, method = "pearson")$estimate)); r_sh

r_ren <- ren_mdt %>% 
  group_by(Class) %>%
  summarize(r = cor.test(x = 1 - x, y = div, method = "pearson")$estimate) %>%
  rbind(data.frame("Class" = "All",
                   "r" = cor.test(x = 1 - ren_mdt$x, y = ren_mdt$div, method = "pearson")$estimate)); r_ren

r_gs <- gs_mdt %>% 
  group_by(Class) %>%
  summarize(r = cor.test(x = 1 - x, y = div, method = "pearson")$estimate) %>%
  rbind(data.frame("Class" = "All",
                   "r" = cor.test(x = 1 - gs_mdt$x, y = gs_mdt$div, method = "pearson")$estimate)); r_gs

r_is <- is_mdt %>% 
  group_by(Class) %>%
  summarize(r = cor.test(x = 1 - x, y = div, method = "pearson")$estimate) %>%
  rbind(data.frame("Class" = "All",
                   "r" = cor.test(x = 1 - is_mdt$x, y = is_mdt$div, method = "pearson")$estimate)); r_is

r_ts <- ts_mdt %>% 
  group_by(Class) %>%
  summarize(r = cor.test(x = 1 - x, y = div, method = "pearson")$estimate) %>%
  rbind(data.frame("Class" = "All",
                   "r" = cor.test(x = 1 - ts_mdt$x, y = ts_mdt$div, method = "pearson")$estimate)); r_ts


## Plot
mdt_theme <- list(
  scale_x_reverse(limits = c(1, 0), breaks = seq(0, 1, 0.2)),
  scale_color_manual(values = paletteer::paletteer_d("RColorBrewer::Set2"), name = "Diversity\nclass."), # monoform, polyform
  xlab(expression(Prop.[MDT])),
  ylab("Diversity"),
  xlim(c(1, 0.2)),
  guides(color = guide_legend(reverse = TRUE)),
  theme_linedraw(base_size = 5, base_family = "Arial"),
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 5, hjust = 0.5),
        text = element_text(size = 5),
        axis.text = element_text(size = 5),
        legend.text = element_text(size = 5),
        strip.text = element_text(color = "black", size = 5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = NA, color = NA),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(2, "mm"))
)

mdtdata %>%
  mutate("index" = factor(index, 
                          levels = c("Tsallis", "Shannon", "Gini-Simpson", "Renyi", "Inverse-Simpson"),
                          labels = c("Tsallis", "Shannon", "Gini-Simpson", "Renyi", "Inverse-Simpson"))) %>%
  group_by(index) %>%
  sample_n(5000) %>%
  ggplot(aes(x = x, y = div)) +
  
  geom_point(aes(x = x, y = div, color = Class),
             pch = ".", stroke = 0, alpha = 1, size = 0) +
  
  facet_wrap(~ index, nrow = 1, scales = "free") +
  
  geom_hline(aes(yintercept = cutoff), color = "#C16622", linetype = "solid", 
             linewidth = 0.15, alpha = 1) +
  geom_vline(xintercept = 0.8, linewidth = 0.15, color = "#C16622", linetype = "solid") +
  
  mdt_theme
