# ==============================================================================
# 19_batch1_sensitivity_pt2.R
# run the edgeR models on batch1 alone
# ==============================================================================




# add RUV from script 18 and re-run edgeR on batch1 alone
# ------------------------------------------------------------------------------

designmat_batch1 <-
  model.matrix(~ Group + Age + Sex + Obese + ruv_correction_batch1$W,
               data = metadata_filtered[filter_samples_batch1, ])
edgeR_obj_batch1 <- estimateDisp(edgeR_obj_batch1, designmat_batch1)
edgeR_obj_batch1 <- calcNormFactors(edgeR_obj_batch1)
edgeR_fit_batch1 <- glmQLFit(edgeR_obj_batch1, designmat_batch1)



# between-group comparisons
# ------------------------------------------------------------------------------

# logFC > 0 means increase in marijuana group
sensitivity_batch_MvT <-
  extractPairwiseContrasts(3, edgeR_fit_batch1) %>%
  mutate(logFC = logFC * -1,
         logCPM = logCPM * -1)

# logFC > 0 means increase in marijuana group
sensitivity_batch_MvC <-
  extractPairwiseContrasts(2, edgeR_fit_batch1) %>%
  mutate(logFC = logFC * -1,
         logCPM = logCPM * -1)

# logFC > 0 makes increase in tobacco group
sensitivity_batch_TvC <-
  extractPairwiseContrasts(c(0, -1, 1, 0, 0, 0, 0, 0, 0),
                           edgeR_fit_batch1,
                           method = "contrast")






# logFC > 0 means increase in marijuana group
sensitivity_batch_MvT_all <-
  extractPairwiseContrasts(3, edgeR_fit_batch1, fdr_thresh = 1.1) %>%
  mutate(logFC = logFC * -1,
         logCPM = logCPM * -1) %>%
  full_join(edgeR_test_results_MvT_all, by = "ID") %>%
  mutate(NegLogP.x = -log10(PValue.x),
         NegLogP.y = -log10(PValue.y))

# logFC > 0 means increase in marijuana group
sensitivity_batch_MvC_all <-
  extractPairwiseContrasts(2, edgeR_fit_batch1, fdr_thresh = 1.1) %>%
  mutate(logFC = logFC * -1,
         logCPM = logCPM * -1) %>%
  full_join(edgeR_test_results_MvC_all, by = "ID") %>%
  mutate(NegLogP.x = -log10(PValue.x),
         NegLogP.y = -log10(PValue.y))

# logFC > 0 makes increase in tobacco group
sensitivity_batch_TvC_all <-
  extractPairwiseContrasts(c(0, -1, 1, 0, 0, 0, 0, 0, 0),
                           edgeR_fit_batch1,
                           fdr_thresh = 1.1,
                           method = "contrast") %>%
  full_join(edgeR_test_results_TvC_all, by = "ID") %>%
  mutate(NegLogP.x = -log10(PValue.x),
         NegLogP.y = -log10(PValue.y))





examine_signif <- function(joined_table) {
  tmp <-
    joined_table %>%
  transmute(Significance = 
              case_when(fdr.x > 0.05 & fdr.y <= 0.05 ~ "Signif Primary",
                        fdr.x <= 0.05 & fdr.y > 0.05 ~ "Signif Batch1-Only",
                        fdr.x <= 0.05 & fdr.y <= 0.05 ~ "Signif Both",
                        fdr.x > 0.05 & fdr.y > 0.05 ~ "Signif Neither")) %>%
              .$Significance %>%
              table
  # tmp["Signif Both"] / sum(tmp[c("Signif Both", "Signif Primary")])
}




 # get simple stats & figures for suppl file
# ------------------------------------------------------------------------------

# number significant in both
list(sensitivity_batch_MvC_all,
     sensitivity_batch_MvT_all,
     sensitivity_batch_TvC_all) %>%
  sapply(examine_signif)


# plot p-values
plot_grid(
  
  ggplot(data = sensitivity_batch_MvC_all,
         aes(x = NegLogP.x, y = NegLogP.y)) +
    geom_point(alpha = 0.5, color = "#7a7a7a") +
    geom_abline(intercept = 0, slope = 1) +
    theme_classic() +
    scale_x_continuous(expression("MvC"~-log[10]~P~", Batch 1 Only"),
                       limits = c(0, 14), expand = c(0, 0)) +
    scale_y_continuous(expression("MvC"~-log[10]~P~", Primary Analysis"),
                       limits = c(0, 14), expand = c(0, 0)) +
    annotate("text", x = 4, y = 13.5,
             label = paste0("Pearson r = ",
                            cor(sensitivity_batch_MvC_all$NegLogP.x,
                                sensitivity_batch_MvC_all$NegLogP.y, use = "complete") %>%
                              format(digits = 2))),
  
  ggplot(data = sensitivity_batch_MvT_all,
         aes(x = NegLogP.x, y = NegLogP.y)) +
    geom_point(alpha = 0.5, color = "#7a7a7a") +
    geom_abline(intercept = 0, slope = 1) +
    theme_classic() +
    scale_x_continuous(expression("MvT"~-log[10]~P~", Batch 1 Only"),
                       limits = c(0, 14), expand = c(0, 0)) +
    scale_y_continuous(expression("MvT"~-log[10]~P~", Primary Analysis"),
                       limits = c(0, 14), expand = c(0, 0)) +
    annotate("text", x = 4, y = 13.5,
             label = paste0("Pearson r = ",
                            cor(sensitivity_batch_MvT_all$NegLogP.x,
                                sensitivity_batch_MvT_all$NegLogP.y, use = "complete") %>%
                              format(digits = 2))),
  
  ggplot(data = sensitivity_batch_TvC_all,
         aes(x = NegLogP.x, y = NegLogP.y)) +
    geom_point(alpha = 0.5, color = "#7a7a7a") +
    geom_abline(intercept = 0, slope = 1) +
    theme_classic() +
    scale_x_continuous(expression("TvC"~-log[10]~P~", Batch 1 Only"),
                       limits = c(0, 14), expand = c(0, 0)) +
    scale_y_continuous(expression("TvC"~-log[10]~P~", Primary Analysis"),
                       limits = c(0, 14), expand = c(0, 0)) +
    annotate("text", x = 4, y = 13.5,
             label = paste0("Pearson r = ",
                            cor(sensitivity_batch_TvC_all$NegLogP.x,
                                sensitivity_batch_TvC_all$NegLogP.y, use = "complete") %>%
                              format(digits = 2))),
  
  nrow = 1)

ggsave("./FiguresTables/FigS3A.png", width = 9, height = 3)
system('open "./FiguresTables/FigS3A.png"')

plot_grid(
  
  ggplot(data = sensitivity_batch_MvC_all,
         aes(x = logFC.x, y = logFC.y)) +
    geom_point(alpha = 0.3, aes(color = fdr.y < 0.05, shape = fdr.y < 0.05)) +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = 0, lty = 3) + geom_hline(yintercept = 0, lty = 3) +
    theme_classic() + theme(legend.position = "none") +
    scale_color_manual(values = c("#7a7a7a", "#67001f")) +
    scale_shape_manual(values = c(1, 19)) +
    scale_x_continuous(expression("MvC"~logFC~", Batch 1 Only"),
                       limits = c(-10, 10), expand = c(0, 0)) +
    scale_y_continuous(expression("MvC"~logFC~", Primary Analysis"),
                       limits = c(-10, 10), expand = c(0, 0)) +
    annotate("text", x = -3.5, y = 9.5,
             label = paste0("Pearson r = ",
                            cor(sensitivity_batch_MvC_all$logFC.x,
                                sensitivity_batch_MvC_all$logFC.y, use = "complete") %>%
                              format(digits = 2))),
  
  
  ggplot(data = sensitivity_batch_MvT_all,
         aes(x = logFC.x, y = logFC.y)) +
    geom_point(alpha = 0.3, aes(color = fdr.y < 0.05, shape = fdr.y < 0.05)) +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = 0, lty = 3) + geom_hline(yintercept = 0, lty = 3) +
    theme_classic() + theme(legend.position = "none") +
    scale_color_manual(values = c("#7a7a7a", "#67001f")) +
    scale_shape_manual(values = c(1, 19)) +
    scale_x_continuous(expression(" MvT"~logFC~", Batch 1 Only"),
                       limits = c(-10, 10), expand = c(0, 0)) +
    scale_y_continuous(expression(" MvT"~logFC~", Primary Analysis"),
                       limits = c(-10, 10), expand = c(0, 0)) +
    annotate("text", x = -3.5, y = 9.5,
             label = paste0("Pearson r = ",
                            cor(sensitivity_batch_MvT_all$logFC.x,
                                sensitivity_batch_MvT_all$logFC.y, use = "complete") %>%
                              format(digits = 2))),
  
  
  ggplot(data = sensitivity_batch_TvC_all,
         aes(x = logFC.x, y = logFC.y)) +
    geom_point(alpha = 0.3, aes(color = fdr.y < 0.05, shape = fdr.y < 0.05)) +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept = 0, lty = 3) + geom_hline(yintercept = 0, lty = 3) +
    theme_classic() + theme(legend.position = "none") +
    scale_color_manual(values = c("#7a7a7a", "#67001f")) +
    scale_shape_manual(values = c(1, 19)) +
    scale_x_continuous(expression("TvC"~logFC~", Batch 1 Only"),
                       limits = c(-10, 10), expand = c(0, 0)) +
    scale_y_continuous(expression("TvC"~logFC~", Primary Analysis"),
                       limits = c(-10, 10), expand = c(0, 0)) +
    annotate("text", x = -3.5, y = 9.5,
             label = paste0("Pearson r = ",
                            cor(sensitivity_batch_TvC_all$logFC.x,
                                sensitivity_batch_TvC_all$logFC.y, use = "complete") %>%
                              format(digits = 2))),
  nrow = 1)

ggsave("./FiguresTables/FigS3B.png", width = 9, height = 3)
system('open "./FiguresTables/FigS3B.png"')




# examples of MvT genes significant in combined analysis only
# ------------------------------------------------------------------------------
sensitivity_batch_MvT_all %>%
  filter(fdr.x < 0.05 & fdr.y > 0.05) %>%
  arrange(-logFC.x) %>%
  .$ID %>% .[1] %>%
plot_gene_by_group(
  ensg = ., residuals = T) +
  facet_grid(~ metadata_filtered$Batch)

sensitivity_batch_MvT_all %>%
  filter(fdr.x < 0.05 & fdr.y > 0.05) %>%
  arrange(logFC.x) %>%
  .$ID %>% .[1] %>%
  plot_gene_by_group(
    ensg = ., residuals = T) +
  facet_grid(~ metadata_filtered$Batch)


# examples of MvT genes exclusively significant in Batch1-only analysis 
# ------------------------------------------------------------------------------
sensitivity_batch_MvT_all %>%
  filter(fdr.x > 0.05 & fdr.y < 0.05) %>%
  arrange(-logFC.x) %>%
  .$ID %>% .[1] %>%
  plot_gene_by_group(
    ensg = ., residuals = T) +
  facet_grid(~ metadata_filtered$Batch)

