# qPCR analysis pure functions

compute_delta_ct <- function(df) {
  df$delta_ct <- df$Ct_target - df$Ct_reference
  df
}

compute_fold_change <- function(df, control_group) {
  control_dct <- df$delta_ct[df$Group == control_group]
  control_mean <- mean(control_dct, na.rm = TRUE)
  df$delta_delta_ct <- df$delta_ct - control_mean
  df$fold_change <- 2^(-df$delta_delta_ct)
  df
}
