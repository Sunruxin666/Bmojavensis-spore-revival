# ============================================================
# 最终规范版代码：Figure 3 - 辐射与模拟微重力对种子萌发的影响
# 布局：第一行1张图（A），第二行3张图（B、C、D）
# 图例放在正上方
# ANOVA结果保留在正文/统计审计中，不放入主图
# 修改：统一术语为"Gravity condition"和"Radiation dose"
# 修改：星号统一标识在图内最上方水平线
# 修改：输出文件名为Fig3_Revised_V2.pdf
# 修改：Panel D标签已修正（fit_d对应t50，fit_e对应rel3h）
# 修改：Panel A纵坐标：Δ Germination rate [h⁻¹]
# 修改：Panel C纵坐标：Relative germination rate
# 新增：导出所有绘图数据为Fig3_Data.xlsx或CSV
# ============================================================

# ============================================================
# 1. 载入所需包
# ============================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(car)
  library(stringr)
  library(showtext)
})

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
input_file <- file.path(project_root, "data", "phenotype", "cleaned_germination_data.csv")
outdir <- file.path(project_root, "results")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. 字体设置（Arial）
# ============================================================
font_regular <- if (Sys.info()[["sysname"]] == "Darwin") {
  "/System/Library/Fonts/Supplemental/Arial.ttf"
} else {
  "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"
}
font_bold <- if (Sys.info()[["sysname"]] == "Darwin") {
  "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
} else {
  "/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"
}
font_add("Arial", regular = font_regular, bold = font_bold)
showtext_auto()
showtext_opts(dpi = 600)

# ============================================================
# 3. 设置全局对比编码
# ============================================================
options(contrasts = c('contr.sum', 'contr.poly'))

# ============================================================
# 4. 统一配色方案
# ============================================================
col_gravity <- c('1g' = '#0072B5', 'SMG' = '#BC3C29')
col_effects <- c('Radiation dose' = '#20854E', 'Gravity condition' = '#0072B5', 
                 'Radiation dose×Gravity condition' = '#BC3C29')
col_stage <- c('0-1 h' = '#E18727', '1-2 h' = '#7876B1', '2-3 h' = '#7E6148')

# ============================================================
# 5. 辅助函数
# ============================================================

fmt_p <- function(p){
  if(!is.finite(p)) return('p = NA')
  if(p < 0.001) return('p < 0.001')
  sprintf('p = %.3f', p)
}

# ANOVA结果格式化（用于图注）- 修改术语
anova_caption <- function(fit, metric_name){
  tab   <- car::Anova(fit, type = 3)
  p_dose <- tab['Dose_f', 'Pr(>F)']
  p_cond <- tab['Condition', 'Pr(>F)']
  p_int  <- tab['Dose_f:Condition', 'Pr(>F)']
  
  paste0(metric_name, ': Radiation dose (', fmt_p(p_dose), '), Gravity condition (',
         fmt_p(p_cond), '), Radiation dose×Gravity condition (', fmt_p(p_int), ')')
}

calc_t50 <- function(r1, r2, r3){
  t <- c(0, 1, 2, 3)
  y <- c(0, r1, r2, r3)
  if(max(y) < 0.5) return(NA_real_)
  i <- which(y >= 0.5)[1]
  if(i == 1) return(0)
  x0 <- t[i-1]; x1 <- t[i]
  y0 <- y[i-1]; y1 <- y[i]
  if(isTRUE(all.equal(y0, y1))) return(x1)
  x0 + (0.5 - y0) / (y1 - y0) * (x1 - x0)
}

decomp_terms <- function(fit, metric){
  tb      <- as.data.frame(car::Anova(fit, type = 3))
  tb$term <- rownames(tb)
  keep    <- tb %>% filter(term %in% c('Dose_f', 'Condition', 'Dose_f:Condition'))
  ss_sum  <- sum(keep$`Sum Sq`)
  keep %>% transmute(
    metric = metric,
    term   = dplyr::recode(term,
                           Dose_f             = 'Radiation dose',
                           Condition          = 'Gravity condition',
                           `Dose_f:Condition` = 'Radiation dose×Gravity condition'),
    prop   = 100 * `Sum Sq` / ss_sum
  )
}

# ============================================================
# 6. 读取数据
# ============================================================
g <- read_csv(input_file, show_col_types = FALSE)

g <- g %>%
  mutate(
    Dose      = as.numeric(Dose),
    Dose_f    = factor(Dose, levels = c(0, 6.5, 13, 19.5)),
    Condition = factor(Condition, levels = c('1g', 'SMG')),
    t50       = mapply(calc_t50, Rate_1h, Rate_2h, Rate_3h),
    d01       = Rate_1h,
    d12       = Rate_2h - Rate_1h,
    d23       = Rate_3h - Rate_2h
  )

base_tbl <- g %>%
  filter(Dose == 0) %>%
  group_by(Condition) %>%
  summarise(base = mean(Rate_3h), .groups = 'drop')

g <- g %>%
  left_join(base_tbl, by = 'Condition') %>%
  mutate(rel3h = Rate_3h / base)

# ============================================================
# 7. 阶段速率数据框
# ============================================================
stage <- g %>%
  select(Condition, Dose, Dose_f, d01, d12, d23) %>%
  pivot_longer(cols = c(d01, d12, d23),
               names_to = 'Stage',
               values_to = 'StageRate') %>%
  mutate(
    Stage = dplyr::recode(Stage,
                          d01 = '0-1 h',
                          d12 = '1-2 h',
                          d23 = '2-3 h'),
    Stage = factor(Stage, levels = c('0-1 h', '1-2 h', '2-3 h'))
  )

# ============================================================
# 8. Panel A 效应量（无横杠误差棒）
# ============================================================
stage_summary <- stage %>%
  group_by(Stage, Dose, Condition) %>%
  summarise(
    mean = mean(StageRate, na.rm = TRUE),
    se = sd(StageRate, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

stage_wide <- stage_summary %>%
  pivot_wider(id_cols = c(Stage, Dose),
              names_from = Condition,
              values_from = c(mean, se))

stage_eff <- stage_wide %>%
  mutate(
    effect = mean_SMG - mean_1g,
    se_effect = sqrt(se_SMG^2 + se_1g^2),
    lwr = effect - se_effect,
    upr = effect + se_effect
  ) %>%
  select(Stage, Dose, effect, lwr, upr, se_effect)

stage_raw_diff <- stage %>%
  group_by(Stage, Dose, Condition) %>%
  mutate(obs_id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(Stage, Dose, obs_id),
              names_from = Condition,
              values_from = StageRate) %>%
  mutate(diff = SMG - `1g`)

# ============================================================
# 9. Panel B和C汇总统计（无横杠误差棒）
# ============================================================
t50_sum <- g %>%
  group_by(Condition, Dose, Dose_f) %>%
  summarise(
    mean = mean(t50, na.rm = TRUE),
    se = sd(t50, na.rm = TRUE) / sqrt(n()),
    lwr = mean - se,
    upr = mean + se,
    .groups = 'drop'
  )

rel_sum <- g %>%
  group_by(Condition, Dose, Dose_f) %>%
  summarise(
    mean = mean(rel3h, na.rm = TRUE),
    se = sd(rel3h, na.rm = TRUE) / sqrt(n()),
    lwr = mean - se,
    upr = mean + se,
    .groups = 'drop'
  )

# ============================================================
# 10. 显著性星号（统一标识在图内最上方水平线）
# ============================================================
# Panel B - 计算统一的最高位置
t50_max_all <- max(t50_sum$mean + t50_sum$se, na.rm = TRUE)
t50_range <- diff(range(c(t50_sum$mean - t50_sum$se, t50_sum$mean + t50_sum$se), na.rm = TRUE))
t50_star_y <- t50_max_all + t50_range * 0.25

stars_t50 <- lapply(sort(unique(g$Dose)), function(d){
  a <- g %>% filter(Condition == '1g', Dose == d) %>% pull(t50)
  b <- g %>% filter(Condition == 'SMG', Dose == d) %>% pull(t50)
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if(length(a) < 2 || length(b) < 2) return(data.frame(Dose = d, label = '', y = NA))
  p <- t.test(a, b)$p.value
  s <- if(p < 0.01) '**' else if(p < 0.05) '*' else ''
  data.frame(Dose = d, label = s, y = t50_star_y)
}) %>% bind_rows() %>% filter(label != '' & !is.na(y))

# Panel C - 计算统一的最高位置
rel_max_all <- max(rel_sum$mean + rel_sum$se, na.rm = TRUE)
rel_range <- diff(range(c(rel_sum$mean - rel_sum$se, rel_sum$mean + rel_sum$se), na.rm = TRUE))
rel_star_y <- rel_max_all + rel_range * 0.25

stars_rel <- lapply(sort(unique(g$Dose)), function(d){
  a <- g %>% filter(Condition == '1g', Dose == d) %>% pull(rel3h)
  b <- g %>% filter(Condition == 'SMG', Dose == d) %>% pull(rel3h)
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if(length(a) < 2 || length(b) < 2) return(data.frame(Dose = d, label = '', y = NA))
  p <- t.test(a, b)$p.value
  s <- if(p < 0.01) '**' else if(p < 0.05) '*' else ''
  data.frame(Dose = d, label = s, y = rel_star_y)
}) %>% bind_rows() %>% filter(label != '' & !is.na(y))

# ============================================================
# 11. 线性模型
# ============================================================
fit_d <- lm(t50 ~ Dose_f * Condition, data = g)           # fit_d 是 t50
fit_e <- lm(rel3h ~ Dose_f * Condition, data = g)         # fit_e 是 rel3h

# ANOVA结果用于图注
anova_caption_t50 <- anova_caption(fit_d, 't50')
anova_caption_rel <- anova_caption(fit_e, 'Relative germination fraction at 3 h')

# ============================================================
# 12. 方差分解（修正：标签与模型对应）
# ============================================================
var_df <- bind_rows(
  decomp_terms(fit_d, 't50'),                    # fit_d 对应 t50
  decomp_terms(fit_e, 'Final germination (3 h)') # fit_e 对应 rel3h/相对萌发率
)

# ============================================================
# 12.5 导出所有绘图数据为CSV（新增功能）
# ============================================================

# 创建导出数据列表
export_data <- list()

# 1. 原始处理数据（包含t50, rel3h, 阶段速率等）
export_data$raw_data <- g

# 2. Panel A: 阶段特异性耦合效应（SMG-1g）及误差
export_data$panel_A_stage_effect <- stage_eff

# 3. Panel A原始差值（个体观测）
export_data$panel_A_raw_diff <- stage_raw_diff

# 4. Panel B: t50汇总统计（均值±SE）
export_data$panel_B_t50_summary <- t50_sum

# 5. Panel B: t50显著性星号
export_data$panel_B_t50_stars <- stars_t50

# 6. Panel C: 相对萌发率汇总统计（均值±SE）
export_data$panel_C_rel_summary <- rel_sum

# 7. Panel C: 相对萌发率显著性星号
export_data$panel_C_rel_stars <- stars_rel

# 8. Panel D: 方差分解结果（解释比例%）
export_data$panel_D_variance_decomposition <- var_df

# 9. ANOVA完整结果（t50模型）
export_data$anova_t50_full <- as.data.frame(car::Anova(fit_d, type = 3)) %>%
  mutate(term = rownames(.))

# 10. ANOVA完整结果（相对萌发率模型）
export_data$anova_rel_full <- as.data.frame(car::Anova(fit_e, type = 3)) %>%
  mutate(term = rownames(.))

# 11. 逐剂量t检验结果（t50）
export_data$ttest_t50_per_dose <- lapply(sort(unique(g$Dose)), function(d){
  a <- g %>% filter(Condition == '1g', Dose == d) %>% pull(t50)
  b <- g %>% filter(Condition == 'SMG', Dose == d) %>% pull(t50)
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if(length(a) < 2 || length(b) < 2) {
    return(data.frame(Dose = d, p_value = NA, t_statistic = NA, 
                      mean_1g = mean(a), mean_SMG = mean(b), n_1g = length(a), n_SMG = length(b)))
  }
  test <- t.test(a, b)
  data.frame(Dose = d, p_value = test$p.value, t_statistic = test$statistic,
             mean_1g = mean(a), mean_SMG = mean(b), n_1g = length(a), n_SMG = length(b))
}) %>% bind_rows()

# 12. 逐剂量t检验结果（相对萌发率）
export_data$ttest_rel_per_dose <- lapply(sort(unique(g$Dose)), function(d){
  a <- g %>% filter(Condition == '1g', Dose == d) %>% pull(rel3h)
  b <- g %>% filter(Condition == 'SMG', Dose == d) %>% pull(rel3h)
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if(length(a) < 2 || length(b) < 2) {
    return(data.frame(Dose = d, p_value = NA, t_statistic = NA,
                      mean_1g = mean(a), mean_SMG = mean(b), n_1g = length(a), n_SMG = length(b)))
  }
  test <- t.test(a, b)
  data.frame(Dose = d, p_value = test$p.value, t_statistic = test$statistic,
             mean_1g = mean(a), mean_SMG = mean(b), n_1g = length(a), n_SMG = length(b))
}) %>% bind_rows()

# 将所有数据写入一个Excel文件（使用openxlsx，如果没有安装则保存为多个CSV）
# 方法1：保存为单个Excel文件（推荐）
if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  wb <- createWorkbook()
  for (name in names(export_data)) {
    addWorksheet(wb, name)
    writeData(wb, name, export_data[[name]])
  }
  saveWorkbook(wb, file.path(outdir, "Figure3_source_data_v16.xlsx"), overwrite = TRUE)
  cat('数据已保存为: Fig3_Data.xlsx\n')
} else {
  # 方法2：如果没有openxlsx，保存为多个CSV文件
  fig3_data_dir <- file.path(outdir, "Figure3_source_data_v16")
  dir.create(fig3_data_dir, showWarnings = FALSE)
  for (name in names(export_data)) {
    write_csv(export_data[[name]], file.path(fig3_data_dir, paste0(name, ".csv")))
  }
  cat('数据已保存至文件夹: ./Fig3_Data/\n')
}

# ============================================================
# 13. 统一主题设置（图例在正上方，字体使用pt单位）
# ============================================================
base_theme <- theme_bw(base_size = 8) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_blank(),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7, color = 'black'),
    axis.line = element_line(linewidth = 0.35, color = 'black'),
    axis.ticks = element_line(linewidth = 0.35, color = 'black'),
    legend.position = 'top',
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(6, 'pt'),
    legend.spacing.x = unit(3, 'pt'),
    legend.margin = margin(b = 2, t = 2),
    strip.background = element_rect(fill = 'grey95', color = NA, linewidth = 0),
    strip.text = element_text(face = 'bold', size = 8, margin = margin(3, 0, 3, 0)),
    panel.border = element_rect(linewidth = 0.35, color = 'black', fill = NA),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 5, 5),
    panel.spacing = unit(6, 'pt')
  )

# ============================================================
# 14. Panel A（修改纵坐标标签：Δ Germination rate [h⁻¹]）
# ============================================================
pA <- ggplot(stage_eff, aes(x = Dose, y = effect)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                width = 0, linewidth = 0.5, color = 'grey40') +
  geom_hline(yintercept = 0, linetype = 'dashed', 
             linewidth = 0.35, color = 'grey50') +
  geom_point(data = stage_raw_diff,
             aes(x = Dose, y = diff),
             alpha = 0.4, size = 1.2, color = 'grey50',
             position = position_jitter(width = 0.15, height = 0, seed = 3001)) +
  geom_line(linewidth = 0.5, color = '#3C5488FF') +
  geom_point(size = 1.8, color = '#3C5488FF', shape = 16) +
  facet_wrap(~ Stage, nrow = 1) +
  scale_x_continuous(breaks = c(0, 6.5, 13, 19.5), 
                     labels = c('0', '6.5', '13', '19.5')) +
  labs(
    x = 'Radiation dose [Gy]',
    y = expression(Delta~Germination~rate~'[h'^-1*']')
  ) +
  base_theme

# ============================================================
# 15. Panel B（星号在最上方）
# ============================================================
y_max_t50 <- max(t50_sum$mean + t50_sum$se, na.rm = TRUE)
y_min_t50 <- min(t50_sum$mean - t50_sum$se, na.rm = TRUE)
y_range_t50 <- y_max_t50 - y_min_t50

pB <- ggplot() +
  geom_point(data = g, aes(Dose, t50, color = Condition),
             alpha = 0.4, size = 1.1,
             position = position_jitter(width = 0.12, height = 0, seed = 3002)) +
  geom_errorbar(data = t50_sum, 
                aes(Dose, ymin = lwr, ymax = upr, color = Condition),
                position = position_dodge(width = 0.5), 
                width = 0, linewidth = 0.5, alpha = 0.8) +
  geom_line(data = t50_sum, 
            aes(Dose, mean, color = Condition, group = Condition),
            linewidth = 0.5) +
  geom_point(data = t50_sum, aes(Dose, mean, color = Condition),
             size = 2, shape = 16) +
  geom_text(data = stars_t50, aes(Dose, y, label = label),
            inherit.aes = FALSE, size = 8/.pt, fontface = 'bold', vjust = 0) +
  scale_color_manual(values = col_gravity, labels = c('1 g', 'SMG')) +
  scale_x_continuous(breaks = c(0, 6.5, 13, 19.5), 
                     labels = c('0', '6.5', '13', '19.5')) +
  coord_cartesian(ylim = c(y_min_t50 - y_range_t50 * 0.05, 
                           t50_star_y + y_range_t50 * 0.05)) +
  labs(x = 'Radiation dose [Gy]', 
       y = expression(italic(t)[50] ~ '[h]')) +
  base_theme

# ============================================================
# 16. Panel C（纵坐标：Relative germination rate）
# ============================================================
y_max_rel <- max(rel_sum$mean + rel_sum$se, na.rm = TRUE)
y_min_rel <- min(rel_sum$mean - rel_sum$se, na.rm = TRUE)
y_range_rel <- y_max_rel - y_min_rel

pC <- ggplot() +
  geom_hline(yintercept = 1.0, linetype = 'dashed', 
             linewidth = 0.35, color = 'grey60') +
  geom_point(data = g, aes(Dose, rel3h, color = Condition),
             alpha = 0.4, size = 1.1,
             position = position_jitter(width = 0.12, height = 0, seed = 3003)) +
  geom_errorbar(data = rel_sum, 
                aes(Dose, ymin = lwr, ymax = upr, color = Condition),
                position = position_dodge(width = 0.5), 
                width = 0, linewidth = 0.5, alpha = 0.8) +
  geom_line(data = rel_sum, 
            aes(Dose, mean, color = Condition, group = Condition),
            linewidth = 0.5) +
  geom_point(data = rel_sum, aes(Dose, mean, color = Condition),
             size = 2, shape = 16) +
  geom_text(data = stars_rel, aes(Dose, y, label = label),
            inherit.aes = FALSE, size = 8/.pt, fontface = 'bold', vjust = 0) +
  scale_color_manual(values = col_gravity, labels = c('1 g', 'SMG')) +
  scale_x_continuous(breaks = c(0, 6.5, 13, 19.5), 
                     labels = c('0', '6.5', '13', '19.5')) +
  coord_cartesian(ylim = c(y_min_rel - y_range_rel * 0.05,
                           rel_star_y + y_range_rel * 0.05)) +
  labs(x = 'Radiation dose [Gy]', 
       y = 'Relative germination fraction at 3 h') +
  base_theme

# ============================================================
# 17. Panel D（图例两行）
# ============================================================
pD <- ggplot(var_df, aes(x = metric, y = prop, fill = term)) +
  geom_col(width = 0.65, color = 'white', linewidth = 0.25) +
  geom_text(aes(label = sprintf('%.1f%%', prop)), 
            position = position_stack(vjust = 0.5), 
            size = 7/.pt, color = 'white', fontface = 'bold', family = 'Arial') +
  scale_fill_manual(values = col_effects,
                    breaks = c('Radiation dose', 'Gravity condition', 'Radiation dose×Gravity condition')) +
  labs(x = NULL, 
       y = 'Explained proportion among fixed terms [%]') +
  base_theme +
  theme(legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7)) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

# ============================================================
# 18. 拼合：第一行1张图（A），第二行3张图（B、C、D）
#     统计检验结果由正文报告，不重复占用图面
# ============================================================
fig3_final <- (pA) / (pB | pC | pD) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 9, face = 'bold', family = 'Arial'))

print(fig3_final)

# ============================================================
# 19. 导出双栏图（宽度 19 cm），版本 V2
# ============================================================
ggsave(file.path(outdir, 'Figure3_submission_v16.pdf'), fig3_final,
       width = 18.3, height = 15.4, units = 'cm', device = cairo_pdf)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(outdir, 'Figure3_submission_v16.svg'), fig3_final,
         width = 18.3, height = 15.4, units = 'cm', device = svglite::svglite)
}
ggsave(file.path(outdir, 'Figure3_submission_v16.tiff'), fig3_final,
       width = 18.3, height = 15.4, units = 'cm', dpi = 600,
       device = ragg::agg_tiff, compression = 'lzw')
ggsave(file.path(outdir, 'Figure3_submission_v16.png'), fig3_final,
       width = 18.3, height = 15.4, units = 'cm', dpi = 600,
       device = ragg::agg_png)

cat('\n========== 完成 ==========\n')
cat('图片已保存为: Fig3_Revised_V2.pdf (双栏, 19 cm width)\n')
cat('数据已保存为: Fig3_Data.xlsx (若已安装openxlsx) 或 ./Fig3_Data/ 文件夹 (多个CSV)\n')
cat('布局：第一行 Panel A，第二行 Panel B | Panel C | Panel D\n')
cat('图例位置：正上方\n')
cat('术语统一：Gravity condition 和 Radiation dose\n')
cat('星号位置：统一标识在图内最上方水平线\n')
cat('ANOVA结果：由正文报告，未重复放入主图\n')
cat('Panel D图例：两行显示\n')
cat('误差棒：无横杠（width=0）\n')
cat('原始数据点：已添加（半透明，带jitter）\n')
cat('坐标轴标签：Panel A = Δ Germination rate [h⁻¹], Panel C = Relative germination rate\n')
