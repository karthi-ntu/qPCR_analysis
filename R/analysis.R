# qPCR analysis pure functions

compute_delta_ct <- function(df) {
  df$delta_ct <- df$Ct_target - df$Ct_reference
  df
}

compute_fold_change <- function(df, control_group) {
  if (!control_group %in% df$Group) {
    stop(sprintf("control_group '%s' not found in df$Group", control_group))
  }
  genes <- unique(df$Gene)
  result_list <- lapply(genes, function(gene) {
    sub <- df[df$Gene == gene, ]
    control_mean <- mean(sub$delta_ct[sub$Group == control_group], na.rm = TRUE)
    sub$delta_delta_ct <- sub$delta_ct - control_mean
    sub$fold_change <- 2^(-sub$delta_delta_ct)
    sub
  })
  do.call(rbind, result_list)
}
