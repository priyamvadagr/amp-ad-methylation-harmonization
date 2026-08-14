library(dplyr)


acc_df <- read.table('/home/ec2-user/data/methyl_harmonization/ROSMAP/Cortical_clock/ROSMAP_clockID.csv', header = TRUE, sep = ',')
ind_metadata <- read.table('/home/ec2-user/data/methyl_harmonization/ROSMAP/metadata/ROSMAP_individual_metadata_processed.csv', header = TRUE, sep = ',')
comb <- merge(acc_df, ind_metadata, by = 'individualID')

comb$diagnosis <- "OTHER"
ad <- comb$Braak %in% c("Stage IV","Stage V","Stage VI") &
      comb$amyCerad %in% c("Frequent/Definite/C3","Moderate/Probable/C2") &
      comb$cogdx %in% 4
comb$diagnosis[ad] <- "AD"

ct <- comb$Braak %in% c("Stage III","Stage II","Stage I","None") &
      comb$amyCerad %in% c("Sparse/Possible/C1","None/No AD/C0") &
      comb$cogdx %in% 1
comb$diagnosis[ct] <- "CT"

table(comb$diagnosis, comb$cohort, useNA = "always")

# diag2: pathology-only (no cogdx) -- harmonizable across all four cohorts
comb$diag2 <- "OTHER2"

ad2 <- comb$Braak %in% c("Stage IV", "Stage V", "Stage VI") &
       comb$amyCerad %in% c("Frequent/Definite/C3", "Moderate/Probable/C2")
comb$diag2[ad2] <- "AD2"

ct2 <- comb$Braak %in% c("Stage III", "Stage II", "Stage I", "None") &
       comb$amyCerad %in% c("Sparse/Possible/C1", "None/No AD/C0")
comb$diag2[ct2] <- "CT2"

table(comb$diag2, useNA = "always")
table(comb$diag2, comb$sex,    useNA = "always")
table(comb$diag2, comb$cohort, useNA = "always")

library(dplyr)
library(ggplot2)

yvar   <- "age_acceleration_cat"       # continuous acceleration variable
grpvar <- "diagnosis"     # AD vs CT (use "diag2" for AD2/CT2)

df <- comb %>%
  filter(.data[[grpvar]] %in% c("AD", "CT"),
         !is.na(.data[[yvar]])) %>%
  mutate(grp = factor(.data[[grpvar]], levels = c("CT", "AD")))

# group comparisons
p_wilcox <- signif(wilcox.test(df[[yvar]] ~ df$grp)$p.value, 3)   # nonparametric
p_ttest  <- signif(t.test(df[[yvar]] ~ df$grp)$p.value, 3)        # Welch t-test

ggplot(df, aes(x = grp, y = .data[[yvar]], fill = grp)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.08, size = 1, alpha = 0.4) +
  scale_fill_manual(values = c(CT = "#2C7FB8", AD = "#C0392B")) +
  labs(x = NULL, y = yvar,
       title = "Cortical clock acceleration: AD vs CT",
       subtitle = paste0("t-test p = ", p_ttest, "   |   Wilcoxon p = ", p_wilcox)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggsave("/home/ec2-user/AMP-AD_methylation_harmonization/Results/Cortical_clock/ROSMAP/acc_AD_vs_CT_violin.pdf", width = 5, height = 5)

df_90 <- df[df$age_cat >= 90,]
p_wilcox_df90 <- signif(wilcox.test(df_90[[yvar]] ~ df_90$grp)$p.value, 3)   # nonparametric
p_ttest_df90  <- signif(t.test(df_90[[yvar]] ~ df_90$grp)$p.value, 3)        # Welch t-test
ggplot(df_90, aes(x = grp, y = .data[[yvar]], fill = grp)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.08, size = 1, alpha = 0.4) +
  scale_fill_manual(values = c(CT = "#2C7FB8", AD = "#C0392B")) +
  labs(x = NULL, y = yvar,
       title = "Cortical clock acceleration: AD vs CT for age >= 90",
       subtitle = paste0("t-test p = ", p_ttest_df90, "   |   Wilcoxon p = ", p_wilcox_df90)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
ggsave("/home/ec2-user/AMP-AD_methylation_harmonization/Results/Cortical_clock/ROSMAP/acc_AD_vs_CT_violin_over_90.pdf", width = 10, height = 5)


df_89 <- df[df$age_cat < 90,]
p_wilcox_df89 <- signif(wilcox.test(df_89[[yvar]] ~ df_89$grp)$p.value, 3)   # nonparametric
p_ttest_df89  <- signif(t.test(df_89[[yvar]] ~ df_89$grp)$p.value, 3)        # Welch t-test
ggplot(df_89, aes(x = grp, y = .data[[yvar]], fill = grp)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.08, size = 1, alpha = 0.4) +
  scale_fill_manual(values = c(CT = "#2C7FB8", AD = "#C0392B")) +
  labs(x = NULL, y = yvar,
       title = "Cortical clock acceleration: AD vs CT for age < 90",
       subtitle = paste0("t-test p = ", p_ttest_df89, "   |   Wilcoxon p = ", p_wilcox_df89)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
ggsave("/home/ec2-user/AMP-AD_methylation_harmonization/Results/Cortical_clock/ROSMAP/acc_AD_vs_CT_violin_under_90.pdf", width = 10, height = 5)
summary(acc_df$age_cat)

