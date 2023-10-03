library(tidyverse)
library(agriutilities)
library(ggsci)
library(ggpubr)
library(stringr)
library(lmerTest)
library(openxlsx)
library(car)
library(performance)
library(asreml)
library(desplot)
library(broom)
library(lmtest)
library(emmeans)
emm_options(rg.limit = 100000)

# Reading Files
dt_1 <- read.delim(file = "data/exp1.txt") %>%
  mutate(Exp = "Exp_1", .before = Block) %>%
  mutate(Loc = paste0("Loc_", Loc))
dt_2 <- read.delim(file = "data/exp2.txt") %>%
  mutate(Exp = "Exp_2", .before = Block) %>%
  mutate(Loc = paste0("Loc_", Loc))
dt_3 <- read.delim(file = "data/exp3.txt") %>%
  mutate(Exp = "Exp_3", .before = Block) %>%
  mutate(Loc = paste0("Loc_", Loc))

# Stacking Data
dt_all <- rbind.data.frame(dt_1, dt_2, dt_3) %>%
  mutate(Trial = paste(Exp, Loc, sep = " - "), .before = row) %>%
  mutate(Block = as.factor(Block)) %>%
  mutate(plot = 1:n(), .before = Exp)
head(dt_all)

dt_all %>%
  mutate(diff = PlantHeight - EarHeight) %>%
  arrange(diff)

# -------------------------------------------------------------------------
# Exploring Trial Layout - Exp Design -------------------------------------
# -------------------------------------------------------------------------

p0 <- dt_all %>%
  mutate(Genotype = str_remove(pattern = "Geno_", Genotype)) %>%
  ggplot(aes(x = col, y = row)) +
  geom_tile(aes(fill = Block), color = "white", alpha = 0.8) +
  geom_tileborder(
    mapping = aes(group = 1, grp = Block),
    lineend = "round",
    color = "black",
    lwd = 0.6,
    linetype = 2
  ) +
  facet_grid(Exp ~ Loc, switch = "y") +
  geom_text(aes(label = Genotype), size = 3) +
  theme_bw(base_size = 15) +
  theme(legend.position = "top") +
  fill_palette(palette = "jco") +
  labs(
    title = "Spatial Layout",
    subtitle = "3 Experiments in 3 Locations",
    caption = "(RCBD - 36 Genotypes - 4 Complete Blocks)"
  )
p0

ggsave(
  filename = "images/experimental_design.png",
  plot = p0,
  units = "in",
  dpi = 300,
  width = 10,
  height = 8
)

# Summary by Exp by Trait
summary_table <- dt_all %>%
  dplyr::select(Exp, Block, Loc, Trial, PlantHeight, EarHeight) %>%
  gather(key = "Trait", value = "Value", -c(1:4)) %>%
  group_by(Exp, Trait, Loc) %>%
  summarise(
    min = min(Value, na.rm = TRUE),
    mean = mean(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    cv = sd / mean * 100 %>% round(., 2)
  )
summary_table

# Boxplot by Trial
p1 <- dt_all %>%
  dplyr::select(Exp, Block, Loc, Trial, PlantHeight, EarHeight) %>%
  gather(key = "Trait", value = "Value", -c(1:4)) %>%
  na.omit() %>%
  ggplot(
    aes(x = Exp, y = Value, fill = Loc)
  ) +
  geom_boxplot() +
  facet_wrap(~Trait, scales = "free") +
  theme_bw(base_size = 15) +
  theme(legend.position = "top") +
  fill_palette(palette = "jco")
p1

# Removing the most atypical Values
boxplot.stats(dt_all$PlantHeight)$out

dt_filtered <- dt_all %>%
  filter(!plot %in% c("972", "154", "651")) %>%
  droplevels()

# Boxplot Without Outliers
p1 <- dt_filtered %>%
  dplyr::select(Exp, Block, Loc, Trial, PlantHeight, EarHeight) %>%
  gather(key = "Trait", value = "Value", -c(1:4)) %>%
  ggplot(
    aes(x = Exp, y = Value, fill = Loc)
  ) +
  geom_boxplot() +
  facet_wrap(~Trait, scales = "free") +
  theme_bw(base_size = 15) +
  theme(legend.position = "top") +
  fill_palette(palette = "jco")
p1

ggsave(
  filename = "images/boxplot_raw_data.png",
  plot = p1,
  units = "in",
  dpi = 300,
  width = 7,
  height = 4
)

# Dot Plot Genotypes
tt_means <- dt_filtered %>%
  dplyr::select(Exp, Genotype, Block, Loc, Trial, PlantHeight, EarHeight) %>%
  gather(key = "Trait", value = "Value", -c(1:5)) %>%
  group_by(Exp, Genotype, Trait, Loc) %>%
  summarise(mean = mean(Value, na.rm = TRUE))

tmp_objt <- list()
for (i in c("Exp_1", "Exp_2", "Exp_3")) {
  gen_plot <- dt_filtered %>%
    dplyr::select(Exp, Genotype, Block, Loc, Trial, PlantHeight, EarHeight) %>%
    gather(key = "Trait", value = "Value", -c(1:5)) %>%
    filter(Exp %in% i) %>%
    droplevels() %>%
    ggplot(
      aes(x = Genotype, y = Value, color = Loc)
    ) +
    geom_point(alpha = 0.6) +
    geom_point(
      data = tt_means %>%
        filter(Exp %in% i) %>%
        droplevels(), aes(x = Genotype, y = mean), color = "black", shape = 3
    ) +
    geom_line(
      data = tt_means %>% filter(Exp %in% i) %>%
        droplevels(), aes(x = Genotype, y = mean, group = paste(Exp, Trait, Loc))
    ) +
    facet_grid(Trait ~ Exp, scales = "free", switch = "y") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    ) +
    color_palette(palette = "jco") +
    labs(y = NULL, x = NULL)
  gen_plot

  tmp_objt[[i]] <- gen_plot

  ggsave(
    filename = paste0("images/gen_by_loc_rawdata_", i, ".png"),
    plot = gen_plot,
    units = "in",
    dpi = 300,
    width = 8,
    height = 6
  )
}

p_tmp <- ggarrange(plotlist = tmp_objt, common.legend = TRUE, ncol = 1)
p_tmp

ggsave(
  filename = paste0("images/gen_by_loc_rawdata_all", ".png"),
  plot = p_tmp,
  units = "in",
  dpi = 300,
  width = 9,
  height = 12
)

# Correlation Between Traits
p2 <- dt_filtered %>%
  ggplot(
    aes(x = PlantHeight, y = EarHeight)
  ) +
  geom_point(aes(color = Loc), alpha = 0.5, size = 2) +
  facet_wrap(~Trial, scales = "free") +
  theme_bw(base_size = 15) +
  stat_cor(size = 3, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "grey30", linewidth = 0.8) +
  color_palette(palette = "jco") +
  theme(legend.position = "top")
p2

ggsave(
  filename = "images/pheno_corr_raw_data.png",
  plot = p2,
  units = "in",
  dpi = 300,
  width = 9,
  height = 8
)

# Correlation Between Locations by Experiment
p3 <- dt_filtered %>%
  filter(Exp %in% "Exp_1") %>%
  gather(key = "Trait", value = "Value", -c(1:8)) %>%
  group_by(Trial, Trait, Genotype) %>%
  summarise(mean = mean(Value, na.rm = TRUE)) %>%
  mutate(Trial = str_remove(Trial, "Exp_1 - ")) %>%
  filter(Trait %in% "EarHeight") %>%
  pivot_wider(names_from = c(Trial), values_from = mean) %>%
  select_if(is.numeric) %>%
  gg_cor(label_size = 4) +
  ggtitle(
    label = "Ear Height",
    subtitle = "Experiment 1"
  )
p3

p3.1 <- dt_filtered %>%
  filter(Exp %in% "Exp_1") %>%
  gather(key = "Trait", value = "Value", -c(1:8)) %>%
  group_by(Trial, Trait, Genotype) %>%
  summarise(mean = mean(Value, na.rm = TRUE)) %>%
  mutate(Trial = str_remove(Trial, "Exp_1 - ")) %>%
  filter(Trait %in% "PlantHeight") %>%
  pivot_wider(names_from = c(Trial), values_from = mean) %>%
  select_if(is.numeric) %>%
  gg_cor(label_size = 4) +
  ggtitle(
    label = "Plant Height",
    subtitle = "Experiment 1"
  )
p3.1

p4 <- dt_filtered %>%
  filter(Exp %in% "Exp_2") %>%
  gather(key = "Trait", value = "Value", -c(1:8)) %>%
  group_by(Trial, Trait, Genotype) %>%
  summarise(mean = mean(Value, na.rm = TRUE)) %>%
  mutate(Trial = str_remove(Trial, "Exp_2 - ")) %>%
  filter(Trait %in% "EarHeight") %>%
  pivot_wider(names_from = c(Trial), values_from = mean) %>%
  select_if(is.numeric) %>%
  gg_cor(label_size = 4) +
  ggtitle(
    label = "Ear Height",
    subtitle = "Experiment 2"
  )
p4
p4.1 <- dt_filtered %>%
  filter(Exp %in% "Exp_2") %>%
  gather(key = "Trait", value = "Value", -c(1:8)) %>%
  group_by(Trial, Trait, Genotype) %>%
  summarise(mean = mean(Value, na.rm = TRUE)) %>%
  mutate(Trial = str_remove(Trial, "Exp_2 - ")) %>%
  filter(Trait %in% "PlantHeight") %>%
  pivot_wider(names_from = c(Trial), values_from = mean) %>%
  select_if(is.numeric) %>%
  gg_cor(label_size = 4) +
  ggtitle(
    label = "Plant Height",
    subtitle = "Experiment 2"
  )
p4.1

p5 <- dt_filtered %>%
  filter(Exp %in% "Exp_3") %>%
  gather(key = "Trait", value = "Value", -c(1:8)) %>%
  group_by(Trial, Trait, Genotype) %>%
  summarise(mean = mean(Value, na.rm = TRUE)) %>%
  mutate(Trial = str_remove(Trial, "Exp_3 - ")) %>%
  filter(Trait %in% "EarHeight") %>%
  pivot_wider(names_from = c(Trial), values_from = mean) %>%
  select_if(is.numeric) %>%
  gg_cor(label_size = 4) +
  ggtitle(
    label = "Ear Height",
    subtitle = "Experiment 3"
  )
p5
p5.1 <- dt_filtered %>%
  filter(Exp %in% "Exp_3") %>%
  gather(key = "Trait", value = "Value", -c(1:8)) %>%
  group_by(Trial, Trait, Genotype) %>%
  summarise(mean = mean(Value, na.rm = TRUE)) %>%
  mutate(Trial = str_remove(Trial, "Exp_3 - ")) %>%
  filter(Trait %in% "PlantHeight") %>%
  pivot_wider(names_from = c(Trial), values_from = mean) %>%
  select_if(is.numeric) %>%
  gg_cor(label_size = 4) +
  ggtitle(
    label = "Plant Height",
    subtitle = "Experiment 3"
  )
p5.1

a <- ggarrange(p3, p4, p5, p3.1, p4.1, p5.1, ncol = 3, nrow = 2)
a

ggsave(
  filename = "images/pheno_corr_loc_raw_data.png",
  plot = a,
  units = "in",
  dpi = 300,
  width = 7,
  height = 5
)

# -------------------------------------------------------------------------
# Spatial Trend -----------------------------------------------------------
# -------------------------------------------------------------------------

p6 <- dt_filtered %>%
  mutate(Genotype = str_remove(pattern = "Geno_", Genotype)) %>%
  ggplot(aes(x = col, y = row)) +
  geom_tile(aes(fill = EarHeight), color = "white", alpha = 0.8) +
  geom_tileborder(
    mapping = aes(group = 1, grp = Block),
    lineend = "round",
    color = "black",
    lwd = 0.5,
    linetype = 2
  ) +
  scale_fill_gradientn(colors = topo.colors(100), na.value = "white") +
  facet_grid(Exp ~ Loc, switch = "y") +
  geom_text(aes(label = Genotype, ), size = 3) +
  theme_bw(base_size = 15) +
  theme(
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  labs(title = "Spatial Trend (Ear Height)")
p6

ggsave(
  filename = "images/spatial_trend_earheight.png",
  plot = p6,
  units = "in",
  dpi = 300,
  width = 12,
  height = 8
)

p7 <- dt_filtered %>%
  mutate(Genotype = str_remove(pattern = "Geno_", Genotype)) %>%
  ggplot(aes(x = col, y = row)) +
  geom_tile(aes(fill = PlantHeight), color = "white", alpha = 0.8) +
  geom_tileborder(
    mapping = aes(group = 1, grp = Block),
    lineend = "round",
    color = "black",
    lwd = 0.5,
    linetype = 2
  ) +
  scale_fill_gradientn(colors = topo.colors(100), na.value = "white") +
  facet_grid(Exp ~ Loc, switch = "y") +
  geom_text(aes(label = Genotype, ), size = 3) +
  theme_bw(base_size = 15) +
  theme(
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  labs(title = "Spatial Trend (Plant Height)")
p7

ggsave(
  filename = "images/spatial_trend_plantheight.png",
  plot = p7,
  units = "in",
  dpi = 300,
  width = 12,
  height = 8
)

# -------------------------------------------------------------------------
# Problem 2 ---------------------------------------------------------------
# -------------------------------------------------------------------------

pBreaks <- c(0, 0.001, 0.01, 0.05, Inf)
pLabels <- c("***", "**", "*", "ns")

dt_filtered <- dt_filtered %>%
  mutate(col_f = as.factor(col), row_f = as.factor(row))

traits <- c("PlantHeight", "EarHeight")
exps <- unique(dt_filtered$Exp)

varcomp <- aov <- res_table <- BLUEs_loc <- list()

i <- 1
for (trait in traits) {
  for (exp in exps) {
    # Subset Experiment
    data_tmp <- dt_filtered %>%
      filter(Exp %in% exp) %>%
      droplevels()
    # Equation
    equation <- reformulate(
      termlabels = c(
        "Genotype",
        "Loc",
        "Genotype:Loc",
        "Loc:Block",
        "Loc:Block:row_f",
        "Loc:Block:col_f"
      ),
      response = trait
    )
    # Modelling
    model <- lm(formula = equation, data = data_tmp)
    # Variance Components
    varcomp[[i]] <- tidy(anova(model)) %>%
      as.data.frame() %>%
      dplyr::select(term, sumsq) %>%
      mutate(percent = sumsq / sum(sumsq) * 100) %>%
      mutate_if(is.numeric, round, 2) %>%
      mutate(exp = exp, trait = trait, .before = term)
    # ANOVA Fixed Effects
    aov[[i]] <- broom::tidy(anova(model)) %>%
      mutate(exp = exp, trait = trait, .before = term) %>%
      mutate(
        signi = cut(
          x = p.value,
          right = FALSE,
          breaks = pBreaks,
          labels = pLabels
        )
      )
    # Residuals
    res_table[[i]] <- data.frame(
      exp = exp,
      trait = trait,
      residuals = rstandard(model)
    )
    # Assumptions
    png(
      filename = paste0("images/assumptions/", exp, "_", trait, ".png"),
      width = 7,
      height = 8,
      units = "in",
      res = 300
    )
    print(
      check_model(
        x = model,
        check = c("qq", "normality", "linearity", "outliers", "homogeneity"),
        detrend = FALSE
      )
    )
    dev.off()
    # BLUEs
    mean_comparisons <- model %>%
      emmeans(pairwise ~ "Loc", adjust = "tukey") %>%
      pluck("emmeans") %>%
      multcomp::cld(details = TRUE, Letters = letters)
    BLUEs_loc[[i]] <- mean_comparisons$emmeans %>%
      data.frame() %>%
      mutate(exp = exp, trait = trait, .before = Loc)
    i <- i + 1
  }
}

# -------------------------------------------------------------------------
# Saving Results ----------------------------------------------------------
# -------------------------------------------------------------------------

varcomp <- data.table::rbindlist(varcomp)
aov <- data.table::rbindlist(aov)
res_table <- data.table::rbindlist(res_table)
BLUEs_loc <- data.table::rbindlist(BLUEs_loc)

p_anova <- aov %>%
  filter(!term %in% "Residuals") %>%
  ggplot(aes(x = term, y = -log10(p.value), fill = term)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
  facet_grid(trait ~ exp, scales = "free") +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "bottom"
  ) +
  fill_palette(palette = "jco") +
  geom_text(
    aes(label = signi),
    position = position_dodge2(width = 0.9, preserve = "single"),
    vjust = -0.1,
    size = 5
  ) +
  ylim(c(0, 170)) +
  labs(
    x = NULL,
    fill = "Source",
    title = "ANOVA - Fixed Effects",
    subtitle = "Significance"
  )
p_anova

ggsave(
  filename = "images/anova_fixed_significance.png",
  plot = p_anova,
  units = "in",
  dpi = 300,
  width = 8,
  height = 7
)

# Nomality
k <- res_table %>%
  select(exp, trait) %>%
  unique.data.frame() %>%
  mutate(id = paste(exp, trait, sep = " - ")) %>%
  pull(id)
k

norm_table <- list()
for (i in k) {
  tmp_dt <- res_table %>%
    mutate(id = paste(exp, trait, sep = " - ")) %>%
    filter(id %in% i) %>%
    droplevels()

  ShapiroWilk <- shapiro.test(tmp_dt$residuals)$p.value
  AndersonDarling <- nortest::ad.test(tmp_dt$residuals)$p.value
  Lilliefors <- nortest::lillie.test(tmp_dt$residuals)$p.value
  KS <- ks.test(tmp_dt$residuals, "pnorm", alternative = "two.sided")$p.value

  norm_table[[i]] <- data.frame(
    "Exp_Trait" = i,
    ShapiroWilk = ShapiroWilk,
    AndersonDarling = AndersonDarling,
    Lilliefors = Lilliefors,
    KolmogorovS = KS
  )
}
norm_table <- data.table::rbindlist(norm_table) %>%
  mutate_if(is.numeric, round, 3)

# -------------------------------------------------------------------------
# Problem 3 ---------------------------------------------------------------
# -------------------------------------------------------------------------

dt_filter_1 <- dt_filtered %>%
  filter(Exp %in% "Exp_1") %>%
  droplevels()

# Model without blocking
model_1 <- lm(formula = PlantHeight ~ Genotype + Loc, data = dt_filter_1)
res_mod_1 <- data.frame(dt_filter_1, res = rstandard(model_1))

# Model blocking
model_2 <- lm(
  formula = PlantHeight ~ Genotype + Loc + Loc:Block,
  data = dt_filter_1
)
res_mod_2 <- data.frame(dt_filter_1, res = rstandard(model_2))

# Comparison
anova(model_1, model_2)
mod_comp <- rbind.data.frame(
  data.frame(model = "Model_null", glance(model_1)),
  data.frame(model = "Model_full", glance(model_2))
) %>%
  dplyr::select(-statistic, -p.value)
mod_comp

p6 <- res_mod_1 %>%
  mutate(Genotype = str_remove(pattern = "Geno_", Genotype)) %>%
  ggplot(aes(x = col, y = row)) +
  geom_tile(aes(fill = res), color = "white", alpha = 0.8) +
  geom_tileborder(
    mapping = aes(group = 1, grp = Block),
    lineend = "round",
    color = "black",
    lwd = 0.5,
    linetype = 2
  ) +
  scale_fill_gradientn(colors = topo.colors(100), na.value = "white") +
  facet_grid(Exp ~ Loc, switch = "y") +
  geom_text(aes(label = Genotype, ), size = 3) +
  theme_bw(base_size = 15) +
  theme(
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  labs(title = "Residuals Spatial Trend (Null Model)")
p6

ggsave(
  filename = "images/spatial_trend_exp1_m0.png",
  plot = p6,
  units = "in",
  dpi = 300,
  width = 12,
  height = 5
)

p7 <- res_mod_2 %>%
  mutate(Genotype = str_remove(pattern = "Geno_", Genotype)) %>%
  ggplot(aes(x = col, y = row)) +
  geom_tile(aes(fill = res), color = "white", alpha = 0.8) +
  geom_tileborder(
    mapping = aes(group = 1, grp = Block),
    lineend = "round",
    color = "black",
    lwd = 0.5,
    linetype = 2
  ) +
  scale_fill_gradientn(colors = topo.colors(100), na.value = "white") +
  facet_grid(Exp ~ Loc, switch = "y") +
  geom_text(aes(label = Genotype, ), size = 3) +
  theme_bw(base_size = 15) +
  theme(
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  labs(title = "Residuals Spatial Trend (Full Model)")
p7

ggsave(
  filename = "images/spatial_trend_exp1_m1.png",
  plot = p7,
  units = "in",
  dpi = 300,
  width = 12,
  height = 5
)

# -------------------------------------------------------------------------
# Problem 4 ---------------------------------------------------------------
# -------------------------------------------------------------------------

# MET Analysis
traits <- c("PlantHeight", "EarHeight", "PlantHeight + EarHeight")

var_table <- aov_mixed <- h2_table <- resid_mixed <- BLUPs <- raov <- list()

i <- 1
for (trait in traits) {
  for (exp in exps) {
    # Subset Experiment
    data_tmp <- dt_filtered %>%
      filter(Exp %in% exp) %>%
      droplevels()
    # Equation
    equation <- reformulate(
      termlabels = c(
        "(1|Loc)",
        "(1|Loc:Block)",
        "(1|Genotype)",
        "(1|Genotype:Loc)",
        "(1|Loc:Block:row)",
        "(1|Loc:Block:col)"
      ),
      response = trait
    )
    # Modelling
    model <- lmer(formula = equation, data = data_tmp)
    # Variance Components
    var_table[[i]] <- VarCorr(model) %>%
      as.data.frame() %>%
      dplyr::select(grp, vcov) %>%
      mutate(percent = vcov / sum(vcov) * 100) %>%
      mutate(exp = exp, trait = trait, .before = grp)
    # ANOVA Fixed Effects
    aov_mixed[[i]] <- broom.mixed::tidy(anova(model)) %>%
      mutate(exp = exp, trait = trait, .before = term)
    # Heritability
    h2 <- MrBean:::h.cullis(model, gen = "Genotype")
    h2_table[[i]] <- data.frame(exp, trait, h2)
    # Residuals
    resid_mixed[[i]] <- data.frame(
      exp = exp,
      trait = trait,
      residuals = residuals(model, scale = TRUE)
    )
    # BLUPs
    BLUPs[[i]] <- broom.mixed::augment(ranef(model)) %>%
      filter(grp %in% "Genotype") %>%
      dplyr::select(level, estimate) %>%
      mutate(exp = exp, trait = trait, .before = level) %>%
      merge(
        y = broom.mixed::augment(ranef(model)) %>%
          filter(grp %in% "Genotype:Loc") %>%
          dplyr::select(level, estimate) %>%
          mutate(trait = trait, .before = level) %>%
          separate(level, into = c("level", "loc"), sep = ":") %>%
          group_by(trait, level) %>%
          summarise(gxe = mean(estimate)),
        by = c("trait", "level")
      ) %>%
      transmute(exp, trait, level, estimate = estimate + gxe)
    
    # RANOVA
    raov[[i]] <- ranova(model) %>% 
      data.frame(comp = rownames(.), .) %>% 
      mutate(exp = exp, trait = trait, .before = comp) %>% 
      .[-1,] %>% 
      mutate(
        signi = cut(
          x = Pr..Chisq.,
          right = FALSE,
          breaks = pBreaks,
          labels = pLabels
        )
      )
    i <- i + 1
  }
}

# -------------------------------------------------------------------------
# Saving Results ----------------------------------------------------------
# -------------------------------------------------------------------------

var_table <- data.table::rbindlist(var_table)
aov_mixed <- data.table::rbindlist(aov_mixed)
raov <- data.table::rbindlist(raov) %>% 
  select(exp, trait, comp, signi) %>% 
  spread(comp, value = signi) %>% 
  .[, c(1, 2, 4, 3, 8, 7, 5, 6)] 

h2_table <- data.table::rbindlist(h2_table)
resid_mixed <- data.table::rbindlist(resid_mixed)
BLUPs <- data.table::rbindlist(BLUPs) %>%
  droplevels() %>%
  pivot_wider(names_from = trait, values_from = estimate)

var_table_reformate <- var_table %>%
  mutate_if(is.numeric, round, 1) %>%
  # mutate(value = paste0(vcov, " (", round(percent, 1), "%)")) %>%
  dplyr::select(-percent) %>%
  pivot_wider(names_from = grp, values_from = vcov) %>%
  .[, c(1, 2, 6, 3, 8, 7, 4, 5, 9)] %>%
  merge(h2_table, by = c("exp", "trait"))

# -------------------------------------------------------------------------
# Phenotypic Correlation --------------------------------------------------
# -------------------------------------------------------------------------

# Correlation Between BLUPs
cor_plot <- BLUPs %>%
  filter(!exp %in% "Exp_2") %>%
  mutate(exp = ifelse(exp == "Exp_1", "Experiment 1", "Experiment 3")) %>%
  ggplot(
    aes(x = PlantHeight, y = EarHeight)
  ) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "grey30") +
  stat_cor() +
  theme_bw(base_size = 15) +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  labs(
    title = NULL,
    x = "Plant Height (BLUPs)",
    y = "Ear Height (BLUPs)"
  ) +
  facet_wrap(~exp) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    strip.background = element_rect(
      color = "black", fill = "#00000000"
    )
  ) +
  coord_fixed()
cor_plot

ggsave(
  filename = "images/corr_blups_earheight_plantheight.png",
  plot = cor_plot,
  units = "in",
  dpi = 300,
  width = 9,
  height = 3
)

cor(BLUPs %>% select_if(is.numeric))

pheno_geno_corr <- function(varcomps) {
  v1 <- varcomps

  # Phenotypic Correlation
  Px <- v1["PlantHeight - Genotype"] +
    v1["PlantHeight - Genotype:Loc"] / 3 +
    v1["PlantHeight - Residual"] / (3 * 4)
  names(Px) <- "Px"

  Py <- v1["EarHeight - Genotype"] +
    v1["EarHeight - Genotype:Loc"] / 3 +
    v1["EarHeight - Residual"] / (3 * 4)
  names(Py) <- "Py"

  Pz <- v1["PlantHeight + EarHeight - Genotype"] +
    v1["PlantHeight + EarHeight - Genotype:Loc"] / 3 +
    v1["PlantHeight + EarHeight - Residual"] / (3 * 4)
  names(Pz) <- "Pz"

  cov_xy <- (Pz - Px - Py) / 2
  names(cov_xy) <- "Cov_XY"
  P_cor_xy <- cov_xy / sqrt(Px * Py)
  names(P_cor_xy) <- "P_corr_XY"

  # Genotypic
  Gx <- v1["PlantHeight - Genotype"]
  names(Gx) <- "Gx"
  Gy <- v1["EarHeight - Genotype"]
  names(Gy) <- "Gy"
  Gz <- v1["PlantHeight + EarHeight - Genotype"]
  cov_g_xy <- (Gz - Gx - Gy) / 2
  names(cov_g_xy) <- "Cov_g_XY"
  corr_g_xy <- cov_g_xy / sqrt(Gx * Gy)
  names(corr_g_xy) <- "Corr_g_XY"

  # Outputs
  out <- data.frame(
    id = c("V(X)", "V(Y)", "V(Z)", "Cov(X, Y)", "Corr(X, Y)"),
    Phenotypic = c(Px, Py, Pz, cov_xy, P_cor_xy),
    Genotypic = c(Gx, Gy, Gz, cov_g_xy, corr_g_xy),
    row.names = NULL
  )
  return(out)
}

# Experiment 1
pg_vcov_table <- pheno_geno_corr(
  varcomps = var_table[var_table$exp %in% "Exp_1", ] %>%
    transmute(id = paste0(trait, " - ", grp), vcov) %>%
    pull(vcov, name = id)
)
pg_vcov_table

p10 <- var_table %>%
  mutate_if(is.numeric, round, 1) %>%
  mutate(perc = paste0(percent, "%")) %>%
  filter(!trait %in% "PlantHeight + EarHeight") %>%
  ggplot(
    aes(x = grp, y = percent, label = perc)
  ) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_text(nudge_y = 5, size = 5, angle = 0) +
  facet_grid(exp ~ trait, switch = "y") +
  theme_bw(base_size = 15) +
  fill_palette(palette = "jco") +
  labs(x = "Source of Variation", y = "Percentage of Variance Explained (%)") +
  ylim(c(0, 110)) +
  geom_label(
    data = h2_table %>% filter(!trait %in% "PlantHeight + EarHeight"),
    mapping = aes(y = 95, x = 4, label = paste0("H2 = ", h2)),
    color = "red",
    size = 5
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(title = "Variance Components (Experiment 1)")
p10

ggsave(
  filename = "images/var_components_exp1.png",
  plot = p10,
  units = "in",
  dpi = 300,
  width = 8,
  height = 5.5
)

pheno_VCOV <- matrix(
  data = c(
    pg_vcov_table[1, "Phenotypic"],
    pg_vcov_table[4, "Phenotypic"],
    pg_vcov_table[4, "Phenotypic"],
    pg_vcov_table[2, "Phenotypic"]
  ),
  nrow = 2,
  ncol = 2,
  byrow = TRUE,
  dimnames = list(c("PlantHeight", "EarHeight"), c("PlantHeight", "EarHeight"))
)
geno_VCOV <- matrix(
  data = c(
    pg_vcov_table[1, "Genotypic"],
    pg_vcov_table[4, "Genotypic"],
    pg_vcov_table[4, "Genotypic"],
    pg_vcov_table[2, "Genotypic"]
  ),
  nrow = 2,
  ncol = 2,
  byrow = TRUE,
  dimnames = list(c("PlantHeight", "EarHeight"), c("PlantHeight", "EarHeight"))
)
geno_VCOV


a1 <- ggarrange(
  agriutilities::covcor_heat(
    matrix = pheno_VCOV,
    corr = FALSE,
    digits = 2,
    size = 5
  ) +
    theme(legend.position = "none") +
    ggtitle(label = "VCOV"),
  agriutilities::covcor_heat(
    matrix = cov2cor(pheno_VCOV),
    corr = TRUE,
    digits = 2,
    size = 5
  ) +
    theme(legend.position = "none") +
    ggtitle(label = "CORR")
)
a1

ggsave(
  filename = "images/pheno_vcov.png",
  plot = a1,
  units = "in",
  dpi = 300,
  width = 6,
  height = 3
)


a2 <- ggarrange(
  agriutilities::covcor_heat(
    matrix = geno_VCOV,
    corr = FALSE,
    digits = 2,
    size = 5
  ) +
    theme(legend.position = "none") +
    ggtitle(label = "VCOV"),
  agriutilities::covcor_heat(
    matrix = cov2cor(geno_VCOV),
    corr = TRUE,
    digits = 1,
    size = 5
  ) +
    theme(legend.position = "none") +
    ggtitle(label = "CORR")
)
a2

ggsave(
  filename = "images/geno_vcov.png",
  plot = a2,
  units = "in",
  dpi = 300,
  width = 6,
  height = 3
)

# Experiment 3
pg_vcov_table_3 <- pheno_geno_corr(
  varcomps = var_table[var_table$exp %in% "Exp_3", ] %>%
    transmute(id = paste0(trait, " - ", grp), vcov) %>%
    pull(vcov, name = id)
)
pg_vcov_table_3

pg_vcov_table <- rbind.data.frame(
  pg_vcov_table %>% mutate(Exp = "Exp_1", .before = id),
  pg_vcov_table_3 %>% mutate(Exp = "Exp_3", .before = id)
)

# -------------------------------------------------------------------------
# Excel Tables ------------------------------------------------------------
# -------------------------------------------------------------------------

OUT <- createWorkbook()
addWorksheet(OUT, "summary_table")
addWorksheet(OUT, "aov_fixed")
addWorksheet(OUT, "BLUEs_loc")
addWorksheet(OUT, "norm_table")
addWorksheet(OUT, "mod_comp")
addWorksheet(OUT, "varcomp_mixed")
addWorksheet(OUT, "varcomp_mixed_2")
addWorksheet(OUT, "covariance")
addWorksheet(OUT, "BLUPs")
addWorksheet(OUT, "ranova")
writeData(OUT, sheet = "summary_table", x = summary_table)
writeData(OUT, sheet = "aov_fixed", x = aov)
writeData(OUT, sheet = "BLUEs_loc", x = BLUEs_loc)
writeData(OUT, sheet = "norm_table", x = norm_table)
writeData(OUT, sheet = "mod_comp", x = mod_comp)
writeData(OUT, sheet = "varcomp_mixed", x = var_table)
writeData(OUT, sheet = "varcomp_mixed_2", x = var_table_reformate)
writeData(OUT, sheet = "BLUPs", x = BLUPs)
writeData(OUT, sheet = "covariance", x = pg_vcov_table)
writeData(OUT, sheet = "ranova", x = raov)
saveWorkbook(OUT, "tables/summary_results_assigment.xlsx", overwrite = TRUE)

# -------------------------------------------------------------------------
# Multi Trait Model -------------------------------------------------------
# -------------------------------------------------------------------------

dt_filter_1 <- dt_filtered %>%
  filter(Exp %in% "Exp_1") %>%
  mutate(Genotype = as.factor(Genotype), Loc = as.factor(Loc)) %>%
  droplevels()

mod_mtraits <- asreml(
  fixed = cbind(PlantHeight, EarHeight) ~ trait + trait:Loc + trait:Loc:Block,
  random = ~
    corgh(trait):Genotype +
      diag(trait):Genotype:Loc +
      diag(trait):Block:col_f +
      diag(trait):Block:row_f,
  residual = ~ id(units):diag(trait),
  data = dt_filter_1,
  maxit = 100
)
lucid::vc(mod_mtraits)
#
# dt_filter_long <- dt_filtered %>%
#   filter(Exp %in% "Exp_1") %>%
#   mutate(Genotype = as.factor(Genotype), Loc = as.factor(Loc)) %>%
#   droplevels() %>%
#   gather(key = "Trait", value = "Response", -c(1:8, 11, 12)) %>%
#   mutate(Trait = as.factor(Trait))
#
# mod_mtraits <- asreml(
#   fixed = Response ~ Trait + Trait:Loc + Trait:Loc:Block,
#   random = ~
#     us(Trait):Genotype +
#       diag(Trait):Loc:Genotype +
#       diag(Trait):Block:col_f +
#       diag(Trait):Block:row_f,
#   residual = ~ dsum(~ units | Trait),
#   data = dt_filter_long,
#   maxit = 100,
#   trace = 0
# )
# lucid::vc(mod_mtraits)
