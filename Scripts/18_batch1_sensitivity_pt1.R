# ==============================================================================
# 18_batch1_sensitivity_pt1.R
# re-do pre-filtering, RUVSeq k on batch 1 only
# ==============================================================================



# prep edgeR object & metadata
# ==============================================================================

filter_samples_batch1 <- metadata_filtered$Batch == 1

# do filtering
edgeR_obj_batch1 <- edgeR_obj[ , filter_samples_batch1]


# determine number of RUV components (RUVr based on first-pass)
# ==============================================================================

dm_firstpass_batch1 <-
  model.matrix(~ Group + Age + Sex + Obese,
               data = metadata_filtered[filter_samples_batch1, ])
edgeR_obj_firstpass_batch1 <- estimateDisp(edgeR_obj_batch1, dm_firstpass_batch1)
edgeR_obj_firstpass_batch1 <- calcNormFactors(edgeR_obj_firstpass_batch1)
edgeR_fit_batch1 <- glmQLFit(edgeR_obj_firstpass_batch1, dm_firstpass_batch1)
resid_firstpass_batch1 <- residuals(edgeR_fit_batch1, type = "deviance")


evaluate_btwnparticipant_by_k_batch1 <-
  function (k, variance_function = iqr,
            resid_in = resid_firstpass, dm_in = dm_firstpass) {
    
    ruv_results <-
      RUVr(
        edgeR_obj_batch1$counts,
        T,
        k = k,
        resid_in)
    
    ruv_results$normalizedCounts %>% # 17602 x 41
      apply(., 1, variance_function) %>% # within-subject variability
      (function(x) { max(x) - min(x) }) # between-subjects
    
  }

plot(1:15,
     sapply(1:15, evaluate_btwnparticipant_by_k_batch1,
            variance_function = iqr,
            resid_in = resid_firstpass_batch1,
            dm_in = dm_firstpass_batch1))

cv_filter_batch1 <-
  apply(edgeR_obj_batch1$counts, 1,
        function(x) { sd(x) / mean(x) }) %>%
  order(decreasing = T) %>%
  .[1:2500]

dist_by_sample_batch1 <-
  log(edgeR_obj$counts[cv_filter_batch1, filter_samples_batch1] + 1) %>% t %>% dist

evaluate_ruv_k_adonis_batch1 <- function(k_to_try) {
  ruv_results <-
    RUVr(
      edgeR_obj_batch1$counts,
      T,
      k = k_to_try,
      resid_firstpass_batch1
    )
  
  adonis(dist_by_sample_batch1 ~ ruv_results$W,
         data = NULL, permutations = 1000) %>%
    .$aov.tab %>% .$R2 %>% .[1]
  
}


plot(1:15, sapply(1:15, evaluate_ruv_k_adonis_batch1))









# finally, k = 3 RUVSeq correction
# save counts to new edgeR object (edgeR_obj_RUV)
# ==============================================================================
ruv_correction_batch1 <-
  RUVr(
    edgeR_obj_batch1$counts,
    T,
    k = 3,
    resid_firstpass_batch1)

