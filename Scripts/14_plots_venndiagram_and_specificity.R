# ==============================================================================
# 14_plots_venndiagram_and_specificity.R
# "specificity" = examine overlap between pair-wise comparisons 
# ==============================================================================




# generating venn diagram (using eulerr package)
# ------------------------------------------------------------------------------

# ultimately created venn externally for better graphical cntl
# (20210217-VennDiagram.png = AJRCMB letter Figure 2)
# but numbers in the figure were generated here + checked manually
euler(
  list(`Marijuana vs Non-Smoker (MvC)` = edgeR_test_results_MvC$ID,
       `Marijuana vs Tobacco (MvT)` = edgeR_test_results_MvT$ID,
       `Tobacco vs Control (TvC)` = edgeR_test_results_TvC$ID),
  shape = "ellipse") %>%
  plot(labels = T, quantities = T)




# exploring degree of set overlap
# (i.e., is there a shared vs specific tobacco signature)
# ------------------------------------------------------------------------------

# number of significant results in each
c(MvC = nrow(edgeR_test_results_MvC),
  MvT = nrow(edgeR_test_results_MvT),
  TvC = nrow(edgeR_test_results_TvC))

# percent of TvC genes also differentially expressed in MvC
intersect(edgeR_test_results_MvC$external_gene_name,
          edgeR_test_results_TvC$external_gene_name) %>%
  length %>%
  `/`(., edgeR_test_results_TvC$external_gene_name %>% length)

MvC_TvC_specificity$dirsame %>% table



# join MvC and TvC effect estimates
# ------------------------------------------------------------------------------
MvC_TvC_specificity <-
  left_join(edgeR_test_results_MvC_all,
            edgeR_test_results_TvC_all,
            by = c("ID", "external_gene_name"), suffix = c(".MvC", ".TvC")) %>%
  mutate(Significance = case_when(
    fdr.MvC < 0.05 & fdr.TvC < 0.05 ~ "MvC & TvC",
    fdr.MvC < 0.05 & fdr.TvC > 0.05 ~ "MvC",
    fdr.MvC > 0.05 & fdr.TvC < 0.05 ~ "TvC",
    fdr.MvC > 0.05 & fdr.TvC > 0.05 ~ "neither") %>%
      fct_relevel(., c("MvC", "TvC", "MvC & TvC", "neither"))) %>%
  mutate(dirsame = logFC.TvC * logFC.MvC > 0) %>%
  mutate(text_label =
           if_else((Significance != "neither" & dirsame & (abs(logFC.MvC) > 4 | abs(logFC.TvC) > 4)) |
                     (Significance != "neither" & !dirsame & (abs(logFC.MvC) > 1.5 | abs(logFC.TvC) > 1.5)), external_gene_name, "")) %>%
  mutate(text_label = if_else(duplicated(text_label), "", text_label))



# filter above for a given significance group for plotting/corr
# (below demonstrates for either signif in MvC or TvC)
# ------------------------------------------------------------------------------

filter_for_correlation_plot <- T # all 
# filter_for_correlation_plot <- MvC_TvC_specificity$Significance != "neither" # signif in at least one
# filter_for_correlation_plot <- MvC_TvC_specificity$Significance == "MvC & TvC" # signif in both

# Pearson correlation
cor(MvC_TvC_specificity$logFC.MvC[filter_for_correlation_plot],
    MvC_TvC_specificity$logFC.TvC[filter_for_correlation_plot],
    method = "pearson", use = "complete")

# how many effects have same direction?
( count_dir_same <-
  MvC_TvC_specificity$dirsame[filter_for_correlation_plot] %>% table )
count_dir_same["FALSE"] / (count_dir_same["FALSE"] + count_dir_same["TRUE"])

# same direction, in genes significant in either comparison only
MvC_TvC_specificity$dirsame[filter_for_correlation_plot] %>% table

# plot each gene

palette_subspecific_color <-
  c(neither = "black", `MvC` = "#1b9e77", `TvC` = "#7570b3", `MvC & TvC` = "#980043")
palette_subspecific_shape <-
  c(neither = 1, `MvC` = 19, `TvC` = 19, `MvC & TvC` = 19)
palette_subspecific_alpha <-
  c(neither = 0.2, `MvC` = 0.5, `TvC` = 0.5, `MvC & TvC` = 0.5)

specificity_plot_A <-
  ggplot(data = MvC_TvC_specificity,
       aes(x = logFC.MvC, y = logFC.TvC,
           color = Significance, shape = Significance)) +
  geom_point(aes(alpha = Significance), size = 1) +
  geom_text_repel(aes(label = text_label, color = Significance),
                  size = 2.5, segment.size = 0.1,
                  box.padding = unit(0.01, "lines"), show.legend = F) +
  theme_few(base_size = 9) +
  geom_abline(slope = 1, intercept = 0, size = 0.25, color = "#555555", lty = 2) +
  geom_hline(yintercept = 0, size = 0.25, color = "#555555") +
  geom_vline(xintercept = 0, size = 0.25, color = "#555555") +
  scale_color_manual(values = palette_subspecific_color) +
  scale_shape_manual(values = palette_subspecific_shape) +
  scale_alpha_manual(values = palette_subspecific_alpha, guide = F) +
  scale_x_continuous(expression(log[2]~"FC: Marijuana vs. Non-Smoker (MvC)"),
                     limits = c(-8, 8), breaks = seq(-8, 8, 4)) +
  scale_y_continuous(expression(log[2]~"FC: Tobacco vs. Non-Smoker (TvC)"),
                     limits = c(-8, 8), breaks = seq(-8, 8, 4)) +
  annotate(x = 3.5, y = -7, geom = "text", size = 3, hjust = 0,
           label = "Pearson r = 0.56\n67% same directionality") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(1, "pt"))


ggsave("./FiguresTables/FigS4.png", width = 5, height = 5, dpi = 300)
system('open "./FiguresTables/FigS4.png"')

