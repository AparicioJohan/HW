library(tidyverse)
library(agriutilities)
library(ASRgenomics)
library(ggpubr)
library(openxlsx)
library(broom)
pBreaks <- c(0, 0.001, 0.01, 0.05, Inf)
pLabels <- c("***", "**", "*", "ns")

# Reading Files
dt_1 <- read.delim(
  file = "genotypicMeans/genotypicMeans_exp1.txt",
  sep = " "
) %>%
  rownames_to_column("Genotype") %>%
  mutate(Exp = "Exp_1", .before = Genotype)
dt_2 <- read.delim(
  file = "genotypicMeans/genotypicMeans_exp2.txt",
  sep = " "
) %>%
  rownames_to_column("Genotype") %>%
  mutate(Exp = "Exp_2", .before = Genotype)
dt_3 <- read.delim(
  file = "genotypicMeans/genotypicMeans_exp3.txt",
  sep = " "
) %>%
  rownames_to_column("Genotype") %>%
  mutate(Exp = "Exp_3", .before = Genotype)

# Stacking Data
dt_all <- rbind.data.frame(dt_1, dt_2, dt_3) %>%
  mutate(plot = 1:n(), .before = Exp)
head(dt_all)

dt_all %>%
  gather(key = "locus", "Count", -c(1:5)) %>%
  gather(key = "trait", "response", -c(1:3, 6, 7)) %>%
  ggplot(
    aes(x = Count, y = response)
  ) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_regline_equation() +
  facet_grid(trait ~ Exp + locus) +
  theme_bw()

# -------------------------------------------------------------------------
# Point 1 -----------------------------------------------------------------
# -------------------------------------------------------------------------

dt_1 %>%
  gather(key = "locus", "Count", -c(1:4)) %>%
  gather(key = "trait", "response", -c(1:2, 5, 6)) %>%
  ggplot(
    aes(x = Count, y = response)
  ) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_regline_equation(label.y = 10, label.x = 0) +
  facet_grid(trait ~ locus) +
  theme_bw()

# Minor Allele Frequency
M <- dt_1 %>%
  select(Genotype, locus1:locus4) %>%
  column_to_rownames("Genotype") %>%
  as.matrix()
M

maf <- colMeans(M, na.rm = TRUE) / 2
maf <- apply(cbind(maf, 1 - maf), 1, min)
maf

M_clean <- qc.filtering(M = M, maf = 0.05)
M_clean

# Analysis Plant Height ---------------------------------------------------

res_1_PHT <- dt_1 %>%
  gather("Locus", "Count", -c(1:4)) %>%
  group_by(Locus, Count) %>%
  summarise(n = n(), Value = mean(PlantHeight)) %>%
  # Gene and Genotypic Frequencies
  group_by(Locus) %>%
  mutate(
    Genotype_F = n / sum(n),
    p = Genotype_F[3] + 1 / 2 * Genotype_F[2],
    q = Genotype_F[1] + 1 / 2 * Genotype_F[2]
  ) %>%
  mutate(Genotype = c("aa", "Aa", "AA"), .before = n) %>%
  # Expected Genotypic Frequencies
  group_by(Locus, Count) %>%
  mutate(
    Exp_Genotype_F = ifelse(
      test = Genotype == "aa",
      yes = q^2,
      no = ifelse(
        test = Genotype == "AA",
        yes = p^2,
        no = 2 * p * q
      )
    ),
    .before = p
  ) %>%
  mutate(Exp_n = Exp_Genotype_F * 36, .before = Genotype_F) %>%
  # Chi-Square Test (Yate's Correction)
  group_by(Locus) %>%
  mutate(
    chi.sqr = sum((abs(n - Exp_n) - 0.5)^2 / Exp_n),
    pvalue = pchisq(chi.sqr, df = 1, lower.tail = FALSE)
  ) %>%
  # Additive value
  mutate(
    a = ifelse(
      test = Genotype == "aa",
      yes = -(Value[3] - Value[1]) / 2,
      no = ifelse(
        test = Genotype == "AA",
        yes = (Value[3] - Value[1]) / 2,
        no = Value[2] - (Value[3] + Value[1]) / 2
      )
    ),
    .before = p
  ) %>%
  # Gene-Substitution
  mutate(alpha = (a[3] + a[2] * (q - p))) %>%
  # Average Effect
  mutate(alpha_1 = q * alpha, alpha_2 = -p * alpha) %>%
  # Breeding Value
  mutate(
    BV = ifelse(
      test = Genotype == "aa",
      yes = -2 * p * alpha,
      no = ifelse(
        test = Genotype == "AA",
        yes = 2 * q * alpha,
        no = (q - p) * alpha
      )
    ),
    .after = alpha_2
  ) %>%
  # Dominance Deviation
  mutate(Dom = Value - BV)
res_1_PHT

res_1_PHT <- dt_1 %>%
  gather("Locus", "Count", -c(1:4)) %>%
  group_by(Locus, Count) %>%
  summarise(n = n(), Value = mean(PlantHeight)) %>%
  # Gene and Genotypic Frequencies
  group_by(Locus) %>%
  mutate(
    Genotype_F = n / sum(n),
    p = Genotype_F[3] + 1 / 2 * Genotype_F[2],
    q = Genotype_F[1] + 1 / 2 * Genotype_F[2]
  ) %>%
  mutate(Genotype = c("aa", "Aa", "AA"), .before = n) %>%
  # Expected Genotypic Frequencies
  group_by(Locus, Count) %>%
  mutate(
    Exp_Genotype_F = ifelse(
      test = Genotype == "aa",
      yes = q^2,
      no = ifelse(
        test = Genotype == "AA",
        yes = p^2,
        no = 2 * p * q
      )
    ),
    .before = p
  ) %>%
  mutate(Exp_n = Exp_Genotype_F * 36, .before = Genotype_F) %>%
  # Chi-Square Test (Yate's Correction)
  group_by(Locus) %>%
  mutate(
    chi.sqr = sum((abs(n - Exp_n) - 0.5)^2 / Exp_n),
    pvalue = pchisq(chi.sqr, df = 1, lower.tail = FALSE)
  ) %>%
  # Additive value
  mutate(
    a = ifelse(
      test = Genotype == "aa",
      yes = -(Value[3] - Value[1]) / 2,
      no = ifelse(
        test = Genotype == "AA",
        yes = (Value[3] - Value[1]) / 2,
        no = Value[2] - (Value[3] + Value[1]) / 2
      )
    ),
    .before = p
  ) %>%
  # Gene-Substitution
  mutate(alpha = (a[3] + a[2] * (q - p))) %>%
  # Average Effect
  mutate(alpha_1 = q * alpha, alpha_2 = -p * alpha) %>%
  # Breeding Value
  mutate(
    BV = ifelse(
      test = Genotype == "aa",
      yes = -2 * p * alpha,
      no = ifelse(
        test = Genotype == "AA",
        yes = 2 * q * alpha,
        no = (q - p) * alpha
      )
    ),
    .after = alpha_2
  ) %>%
  # Dominance Deviation
  mutate(Dom = Value - BV)
res_1_PHT

# Analysis Ear Height -----------------------------------------------------

res_1_EHT <- dt_1 %>%
  gather("Locus", "Count", -c(1:4)) %>%
  group_by(Locus, Count) %>%
  summarise(n = n(), Value = mean(EarHeight)) %>%
  # Gene and Genotypic Frequencies
  group_by(Locus) %>%
  mutate(
    Genotype_F = n / sum(n),
    p = Genotype_F[3] + 1 / 2 * Genotype_F[2],
    q = Genotype_F[1] + 1 / 2 * Genotype_F[2]
  ) %>%
  mutate(Genotype = c("aa", "Aa", "AA"), .before = n) %>%
  # Expected Genotypic Frequencies
  group_by(Locus, Count) %>%
  mutate(
    Exp_Genotype_F = ifelse(
      test = Genotype == "aa",
      yes = q^2,
      no = ifelse(
        test = Genotype == "AA",
        yes = p^2,
        no = 2 * p * q
      )
    ),
    .before = p
  ) %>%
  mutate(Exp_n = Exp_Genotype_F * 36, .before = Genotype_F) %>%
  # Chi-Square Test (Yate's Correction)
  group_by(Locus) %>%
  mutate(
    chi.sqr = sum((abs(n - Exp_n) - 0.5)^2 / Exp_n),
    pvalue = pchisq(chi.sqr, df = 1, lower.tail = FALSE)
  ) %>%
  # Additive value
  mutate(
    a = ifelse(
      test = Genotype == "aa",
      yes = -(Value[3] - Value[1]) / 2,
      no = ifelse(
        test = Genotype == "AA",
        yes = (Value[3] - Value[1]) / 2,
        no = Value[2] - (Value[3] + Value[1]) / 2
      )
    ),
    .before = p
  ) %>%
  # Gene-Substitution
  mutate(alpha = (a[3] + a[2] * (q - p))) %>%
  # Average Effect
  mutate(alpha_1 = q * alpha, alpha_2 = -p * alpha) %>%
  # Breeding Value
  mutate(
    BV = ifelse(
      test = Genotype == "aa",
      yes = -2 * p * alpha,
      no = ifelse(
        test = Genotype == "AA",
        yes = 2 * q * alpha,
        no = (q - p) * alpha
      )
    ),
    .after = alpha_2
  ) %>%
  # Dominance Deviation
  mutate(Dom = Value - BV)
res_1_EHT

# -------------------------------------------------------------------------

# Graphical Representation
p0 <- dt_1 %>%
  gather(key = "Locus", "Count", -c(1:4)) %>%
  gather(key = "trait", "response", -c(1:2, 5, 6)) %>%
  filter(trait %in% "PlantHeight" & !Locus %in% "locus3") %>%
  ggplot(
    aes(x = Count, y = response)
  ) +
  facet_grid(trait ~ Locus) +
  geom_point(color = "grey") +
  geom_point(
    data = res_1_PHT %>% filter(!Locus %in% "locus3"),
    aes(x = Count, y = Value),
    color = "black",
    fill = "blue",
    shape = 21,
    size = 3,
    alpha = 0.7
  ) +
  geom_text(
    data = res_1_PHT %>% filter(!Locus %in% "locus3"),
    aes(
      x = 1,
      y = -16,
      label = paste0("y = ", round(-2 * p * alpha, 2), " + ", round(alpha, 2), "x")
    ),
    size = 5,
    alpha = 0.4
  ) +
  geom_text(
    data = res_1_PHT %>% filter(!Locus %in% "locus3"),
    aes(
      x = 1,
      y = -18,
      label = paste0("(p = ", round(p, 2), ", q = ", round(q, 2), ")")
    ),
    alpha = 0.3
  ) +
  theme_bw(base_size = 15) +
  geom_abline(
    data = res_1_PHT %>% filter(!Locus %in% "locus3"),
    aes(slope = alpha, intercept = -2 * p * alpha)
  ) +
  geom_point(
    data = res_1_PHT %>% filter(!Locus %in% "locus3"),
    aes(x = Count, y = BV),
    color = "black",
    shape = 21,
    fill = "red",
    size = 3,
    alpha = 0.7
  ) +
  labs(
    y = NULL,
    x = "Number of A genes in the genotype",
    caption =
      "Red points represent Breeding Values while Blue points Genotypic Values,
    of the genotypes for a locus with two alleles, A and a, at frequencies p and q.
    The red dashed line represents the Dominance Deviation."
  ) +
  geom_segment(
    data = res_1_PHT %>% filter(!Locus %in% "locus3"),
    aes(
      xend = Count,
      yend = BV,
      y = Value
    ),
    col = "red",
    lty = "dashed"
  )
p0

ggsave(
  filename = "images/plot_point_1.png",
  plot = p0,
  units = "in",
  dpi = 300,
  width = 10,
  height = 5
)

# Linear Model with Raw Data (Wald's Test)
wald_exp_1 <- dt_1 %>%
  gather(key = "locus", "Count", -c(1:4)) %>%
  gather(key = "trait", "response", -c(1:2, 5, 6)) %>%
  split(f = ~ .$locus + .$trait, sep = "_") %>%
  map(
    .f = \(x) lm(response ~ Count, data = x) %>%
    summary() %>%
    .$coef %>%
    data.frame(coef = row.names(.), ., check.names = TRUE)
  ) %>%
  data.table::rbindlist(idcol = "locus_trait") %>%
  mutate(
    signi = cut(
      x = Pr...t..,
      right = FALSE,
      breaks = pBreaks,
      labels = pLabels
    )
  ) %>%
  mutate_if(is.numeric, round, 2)
wald_exp_1

# Variance Explained
ve_1 <- dt_1 %>%
  gather(key = "locus", "Count", -c(1:4)) %>%
  gather(key = "trait", "response", -c(1:2, 5, 6)) %>%
  split(f = ~ .$locus + .$trait, sep = "_") %>%
  map(.f = \(x) lm(response ~ Count, data = x) %>%
    glance() %>%
    data.frame(.)) %>%
  data.table::rbindlist(idcol = "trait") %>%
  mutate_if(is.numeric, round, 2)
ve_1

ve_all_1 <- dt_1 %>%
  gather(key = "trait", "response", -c(1:2, 5:8)) %>%
  split(f = ~ .$trait) %>%
  map(
    .f = \(x) lm(response ~ locus1 + locus2 + locus3 + locus4, data = x) %>%
      glance() %>%
      data.frame(.)
  ) %>%
  data.table::rbindlist(idcol = "trait") %>%
  mutate_if(is.numeric, round, 2)
ve_all_1

ve_all_1 <- ve_1 %>% rbind.data.frame(ve_all_1)
ve_all_1

# -------------------------------------------------------------------------
# Experiment 3 ------------------------------------------------------------
# -------------------------------------------------------------------------

dt_3 %>%
  gather(key = "locus", "Count", -c(1:4)) %>%
  gather(key = "trait", "response", -c(1:2, 5, 6)) %>%
  ggplot(
    aes(x = Count, y = response)
  ) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_regline_equation(label.y = 10, label.x = 0) +
  facet_grid(trait ~ locus) +
  theme_bw()

# Minor Allele Frequency
M <- dt_3 %>%
  select(Genotype, locus1:locus4) %>%
  column_to_rownames("Genotype") %>%
  as.matrix()
M

maf <- colMeans(M, na.rm = TRUE) / 2
maf <- apply(cbind(maf, 1 - maf), 1, min)
maf

M_clean <- qc.filtering(M = M, maf = 0.05)
M_clean

# Analysis
res_3_PHT <- dt_3 %>%
  gather("Locus", "Count", -c(1:4)) %>%
  group_by(Locus, Count) %>%
  summarise(n = n(), Value = mean(PlantHeight)) %>%
  # filter(!Locus %in% "locus2") %>%
  ungroup() %>%
  add_row(
    Locus = c("locus1", "locus2", "locus3", "locus4"),
    n = rep(0, 4),
    Count = rep(1, 4),
    Value = rep(0, 4)
  ) %>%
  arrange(Locus, Count) %>%
  # Gene and Genotypic Frequencies
  group_by(Locus) %>%
  mutate(
    Genotype_F = n / sum(n),
    p = Genotype_F[3] + 1 / 2 * Genotype_F[2],
    q = Genotype_F[1] + 1 / 2 * Genotype_F[2]
  ) %>%
  mutate(Genotype = c("aa", "Aa", "AA"), .before = n) %>%
  # Expected Genotypic Frequencies
  group_by(Locus, Count) %>%
  mutate(
    Exp_Genotype_F = ifelse(
      test = Genotype == "aa",
      yes = q^2,
      no = ifelse(
        test = Genotype == "AA",
        yes = p^2,
        no = 2 * p * q
      )
    ),
    .before = p
  ) %>%
  mutate(Exp_n = Exp_Genotype_F * 36, .before = Genotype_F) %>%
  # Chi-Square Test (Yate's Correction)
  group_by(Locus) %>%
  mutate(
    chi.sqr = sum((abs(n - Exp_n) - 0.5)^2 / Exp_n),
    pvalue = pchisq(chi.sqr, df = 1, lower.tail = FALSE)
  ) %>%
  # Additive value
  mutate(
    a = ifelse(
      test = Genotype == "aa",
      yes = -(Value[3] - Value[1]) / 2,
      no = ifelse(
        test = Genotype == "AA",
        yes = (Value[3] - Value[1]) / 2,
        no = Value[2] - (Value[3] + Value[1]) / 2
      )
    ),
    .before = p
  ) %>%
  # Gene-Substitution
  mutate(alpha = (a[3] + a[2] * (q - p))) %>%
  # Average Effect
  mutate(alpha_1 = q * alpha, alpha_2 = -p * alpha) %>%
  # Breeding Value
  mutate(
    BV = ifelse(
      test = Genotype == "aa",
      yes = -2 * p * alpha,
      no = ifelse(
        test = Genotype == "AA",
        yes = 2 * q * alpha,
        no = (q - p) * alpha
      )
    ),
    .after = alpha_2
  ) %>%
  # Dominance Deviation
  mutate(Dom = Value - BV)
res_3_PHT

# Graphical Representation
p0 <- dt_3 %>%
  gather(key = "Locus", "Count", -c(1:4)) %>%
  gather(key = "trait", "response", -c(1:2, 5, 6)) %>%
  filter(trait %in% "PlantHeight" & !Locus %in% "locus2") %>%
  ggplot(
    aes(x = Count, y = response)
  ) +
  facet_grid(trait ~ Locus) +
  geom_point(color = "grey") +
  geom_point(
    data = res_3_PHT %>% filter(!Locus %in% "locus2"),
    aes(x = Count, y = Value),
    color = "black",
    fill = "blue",
    shape = 21,
    size = 3,
    alpha = 0.7
  ) +
  geom_text(
    data = res_3_PHT %>% filter(!Locus %in% "locus2"),
    aes(
      x = 1,
      y = -8,
      label = paste0("y = ", round(-2 * p * alpha, 2), " + ", round(alpha, 2), "x")
    ),
    size = 5,
    alpha = 0.4
  ) +
  geom_text(
    data = res_3_PHT %>% filter(!Locus %in% "locus2"),
    aes(
      x = 1,
      y = -10,
      label = paste0("(p = ", round(p, 2), ", q = ", round(q, 2), ")")
    ),
    alpha = 0.3
  ) +
  theme_bw(base_size = 15) +
  geom_abline(
    data = res_3_PHT %>% filter(!Locus %in% "locus2"),
    aes(slope = alpha, intercept = -2 * p * alpha)
  ) +
  geom_point(
    data = res_3_PHT %>% filter(!Locus %in% "locus2"),
    aes(x = Count, y = BV),
    color = "black",
    shape = 21,
    fill = "red",
    size = 3,
    alpha = 0.7
  ) +
  labs(
    y = NULL,
    x = "Number of A genes in the genotype",
    caption =
      "Red points represent Breeding Values while Blue points Genotypic Values,
    of the genotypes for a locus with two alleles, A and a, at frequencies p and q.
    The red dashed line represents the Dominance Deviation."
  ) +
  geom_segment(
    data = res_3_PHT %>% filter(!Locus %in% "locus2"),
    aes(
      xend = Count,
      yend = BV,
      y = Value
    ),
    col = "red",
    lty = "dashed"
  )
p0

ggsave(
  filename = "images/plot_point_4.png",
  plot = p0,
  units = "in",
  dpi = 300,
  width = 10,
  height = 5
)

# Linear Model with Raw Data (Wald's Test)
wald_exp_3 <- dt_3 %>%
  gather(key = "locus", "Count", -c(1:4)) %>%
  gather(key = "trait", "response", -c(1:2, 5, 6)) %>%
  split(f = ~ .$locus + .$trait, sep = "_") %>%
  map(.f = \(x) lm(response ~ Count, data = x) %>%
    summary() %>%
    .$coef %>%
    data.frame(coef = row.names(.), ., check.names = TRUE))%>%
  data.table::rbindlist(idcol = "locus_trait") %>%
  mutate(
    signi = cut(
      x = Pr...t..,
      right = FALSE,
      breaks = pBreaks,
      labels = pLabels
    )
  ) %>%
  mutate_if(is.numeric, round, 2)
wald_exp_3

# -------------------------------------------------------------------------
# Point 6 -----------------------------------------------------------------
# -------------------------------------------------------------------------

# Analysis
res_3_EHT <- dt_3 %>%
  gather("Locus", "Count", -c(1:4)) %>%
  group_by(Locus, Count) %>%
  summarise(n = n(), Value = mean(EarHeight)) %>%
  ungroup() %>%
  add_row(
    Locus = c("locus1", "locus2","locus3", "locus4"),
    n = rep(0, 4),
    Count = rep(1, 4),
    Value = rep(0, 4)
  ) %>%
  arrange(Locus, Count) %>%
  # Gene and Genotypic Frequencies
  group_by(Locus) %>%
  mutate(
    Genotype_F = n / sum(n),
    p = Genotype_F[3] + 1 / 2 * Genotype_F[2],
    q = Genotype_F[1] + 1 / 2 * Genotype_F[2]
  ) %>%
  mutate(Genotype = c("aa", "Aa", "AA"), .before = n) %>%
  # Expected Genotypic Frequencies
  group_by(Locus, Count) %>%
  mutate(
    Exp_Genotype_F = ifelse(
      test = Genotype == "aa",
      yes = q^2,
      no = ifelse(
        test = Genotype == "AA",
        yes = p^2,
        no = 2 * p * q
      )
    ),
    .before = p
  ) %>%
  mutate(Exp_n = Exp_Genotype_F * 36, .before = Genotype_F) %>%
  # Chi-Square Test (Yate's Correction)
  group_by(Locus) %>%
  mutate(
    chi.sqr = sum((abs(n - Exp_n) - 0.5)^2 / Exp_n),
    pvalue = pchisq(chi.sqr, df = 1, lower.tail = FALSE)
  ) %>%
  # Additive value
  mutate(
    a = ifelse(
      test = Genotype == "aa",
      yes = -(Value[3] - Value[1]) / 2,
      no = ifelse(
        test = Genotype == "AA",
        yes = (Value[3] - Value[1]) / 2,
        no = Value[2] - (Value[3] + Value[1]) / 2
      )
    ),
    .before = p
  ) %>%
  # Gene-Substitution
  mutate(alpha = (a[3] + a[2] * (q - p))) %>%
  # Average Effect
  mutate(alpha_1 = q * alpha, alpha_2 = -p * alpha) %>%
  # Breeding Value
  mutate(
    BV = ifelse(
      test = Genotype == "aa",
      yes = -2 * p * alpha,
      no = ifelse(
        test = Genotype == "AA",
        yes = 2 * q * alpha,
        no = (q - p) * alpha
      )
    ),
    .after = alpha_2
  ) %>%
  # Dominance Deviation
  mutate(Dom = Value - BV)
res_3_EHT

# Graphical Representation
p0 <- dt_3 %>%
  gather(key = "Locus", "Count", -c(1:4)) %>%
  gather(key = "trait", "response", -c(1:2, 5, 6)) %>%
  filter(trait %in% "EarHeight" & !Locus %in% "locus2") %>%
  ggplot(
    aes(x = Count, y = response)
  ) +
  facet_grid(trait ~ Locus) +
  geom_point(color = "grey") +
  geom_point(
    data = res_3_EHT %>% filter(!Locus %in% "locus2"),
    aes(x = Count, y = Value),
    color = "black",
    fill = "blue",
    shape = 21,
    size = 3,
    alpha = 0.7
  ) +
  geom_text(
    data = res_3_EHT %>% filter(!Locus %in% "locus2"),
    aes(
      x = 1,
      y = -8,
      label = paste0("y = ", round(-2 * p * alpha, 2), " + ", round(alpha, 2), "x")
    ),
    size = 5,
    alpha = 0.4
  ) +
  geom_text(
    data = res_3_EHT %>% filter(!Locus %in% "locus2"),
    aes(
      x = 1,
      y = -10,
      label = paste0("(p = ", round(p, 2), ", q = ", round(q, 2), ")")
    ),
    alpha = 0.3
  ) +
  theme_bw(base_size = 15) +
  geom_abline(
    data = res_3_EHT %>% filter(!Locus %in% "locus2"),
    aes(slope = alpha, intercept = -2 * p * alpha)
  ) +
  geom_point(
    data = res_3_EHT %>% filter(!Locus %in% "locus2"),
    aes(x = Count, y = BV),
    color = "black",
    shape = 21,
    fill = "red",
    size = 3,
    alpha = 0.7
  ) +
  labs(
    y = NULL,
    x = "Number of A genes in the genotype",
    caption =
      "Red points represent Breeding Values while Blue points Genotypic Values,
    of the genotypes for a locus with two alleles, A and a, at frequencies p and q.
    The red dashed line represents the Dominance Deviation."
  ) +
  geom_segment(
    data = res_3_EHT %>% filter(!Locus %in% "locus2"),
    aes(
      xend = Count,
      yend = BV,
      y = Value
    ),
    col = "red",
    lty = "dashed"
  )
p0

ggsave(
  filename = "images/plot_point_6.png",
  plot = p0,
  units = "in",
  dpi = 300,
  width = 10,
  height = 5
)

# Variance Explained
ve_3 <- dt_3 %>%
  gather(key = "locus", "Count", -c(1:4)) %>%
  gather(key = "trait", "response", -c(1:2, 5, 6)) %>%
  split(f = ~ .$locus + .$trait, sep = "_") %>%
  map(
    .f = \(x) lm(response ~ Count, data = x) %>%
      glance() %>%
      data.frame(.)
  ) %>%
  data.table::rbindlist(idcol = "trait") %>%
  mutate_if(is.numeric, round, 2)
ve_3

ve_all_3 <- dt_3 %>%
  gather(key = "trait", "response", -c(1:2, 5:8)) %>%
  split(f = ~ .$trait) %>%
  map(
    .f = \(x) lm(response ~ locus1 + locus2 + locus4, data = x) %>%
      glance() %>%
      data.frame(.)
  ) %>%
  data.table::rbindlist(idcol = "trait") %>%
  mutate_if(is.numeric, round, 2)
ve_all_3

ve_all_3 <- ve_3 %>% rbind.data.frame(ve_all_3)
ve_all_3


# -------------------------------------------------------------------------
# Stacking Tables ---------------------------------------------------------
# -------------------------------------------------------------------------

res_1 <- rbind.data.frame(
  res_1_PHT %>% mutate(trait = "PlantHeight", .before = Locus),
  res_1_EHT %>% mutate(trait = "EarHeight", .before = Locus)
)

res_3 <- rbind.data.frame(
  res_3_PHT %>% mutate(trait = "PlantHeight", .before = Locus),
  res_3_EHT %>% mutate(trait = "EarHeight", .before = Locus)
)

var_exp <- rbind.data.frame(
  ve_all_1 %>% mutate(exp = "Exp_1", .before = trait),
  ve_all_3 %>% mutate(exp = "Exp_3", .before = trait)
)

# -------------------------------------------------------------------------
# Excel Tables ------------------------------------------------------------
# -------------------------------------------------------------------------

OUT <- createWorkbook()
addWorksheet(OUT, "summary_res_1")
addWorksheet(OUT, "wald_exp_1")
addWorksheet(OUT, "summary_res_3")
addWorksheet(OUT, "wald_exp_3")
addWorksheet(OUT, "variance_exp")
writeData(OUT, sheet = "summary_res_1", x = res_1)
writeData(OUT, sheet = "wald_exp_1", x = wald_exp_1)
writeData(OUT, sheet = "summary_res_3", x = res_3)
writeData(OUT, sheet = "wald_exp_3", x = wald_exp_3)
writeData(OUT, sheet = "variance_exp", x = var_exp)
saveWorkbook(OUT, "tables/results_assigment_2.xlsx", overwrite = TRUE)
