# ============================================================
# 修改版代码：Figure 4 - 芽孢杆菌群体生长分析
# 版本：V6
# 布局：2x2布局（A、B、C、D）
# 图例放在正上方
# 修改要点：
# 1. Panel A：改用误差棒（geom_errorbar无横杠）+ 原始数据点 + position_dodge，纵坐标为ΔOD600
# 2. Panel B：Condition决定颜色（2条线），Dose决定形状
# 3. 字体统一使用.pt单位
# 4. 配色统一
# 5. ANOVA结果保留在正文/统计审计中，不放入主图
# 6. 图CD的星号都在图上方的一条水平线上
# 新增：导出所有绘图数据为Excel或CSV文件
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
  library(minpack.lm)
  library(emmeans)
})

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
input_file <- file.path(project_root, "data", "phenotype", "cleaned_OD600_data.csv")
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
col_dose <- c('0'    = '#0072B5',
              '6.5'  = '#20854E',
              '13'   = '#E18727',
              '19.5' = '#7876B1')
col_cond <- c('1g'  = '#0072B5',
              'SMG' = '#BC3C29')

# ============================================================
# 5. 辅助函数
# ============================================================

fmt_p <- function(p){
  if(!is.finite(p)) return('p = NA')
  if(p < 0.001) return('p < 0.001')
  sprintf('p = %.3f', p)
}

# ANOVA结果格式化（用于图注）
anova_caption <- function(fit, metric_name){
  tab   <- car::Anova(fit, type = 3)
  p_dose <- tab['Dose_f', 'Pr(>F)']
  p_cond <- tab['Condition', 'Pr(>F)']
  p_int  <- tab['Dose_f:Condition', 'Pr(>F)']
  
  paste0(metric_name, ': Radiation dose (', fmt_p(p_dose), '), Gravity condition (',
         fmt_p(p_cond), '), Radiation dose×Gravity condition (', fmt_p(p_int), ')')
}

# 梯形积分
trapz_local <- function(x, y){
  o <- order(x); x <- x[o]; y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# 4参数Logistic函数
logistic4 <- function(t, A, K, r, t0){
  A + (K - A) / (1.0 + exp(-r * (t - t0)))
}

# 计算均值 ± 标准误（SEM）
ci_sem <- function(x){
  x <- x[is.finite(x)]
  n <- length(x)
  m <- mean(x)
  if(n <= 1) return(c(mean=m, se=NA, lwr=NA, upr=NA))
  se <- sd(x) / sqrt(n)
  c(mean=m, se=se, lwr=m-se, upr=m+se)
}

# ============================================================
# 6. 读取数据
# ============================================================
od <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    Dose      = as.numeric(Dose),
    Dose_f    = factor(Dose, levels = c(0, 6.5, 13, 19.5)),
    Condition = factor(Condition, levels = c('1g', 'SMG'))
  )

# ============================================================
# 7. Logistic拟合参数设置
# ============================================================
fit_t    <- c(0, 1, 2, 3, 4.5, 6, 7.5, 9, 10)
fit_cols <- c("0h", "1h", "2h", "3h", "4.5h", "6h", "7.5h", "9h", "10h")

lb <- c(A = 0.0,  K = 0.05, r = 0.01, t0 = -5.0)
ub <- c(A = 0.2,  K = 0.35, r = 3.0,  t0 = 20.0)

# ============================================================
# 8. 单行拟合函数
# ============================================================
fit_one_row <- function(row_vals){
  y   <- as.numeric(row_vals[fit_cols])
  A0  <- min(y, na.rm = TRUE)
  K0  <- max(y, na.rm = TRUE)
  r0  <- 0.8
  t00 <- 5.0
  
  fit_ok    <- TRUE
  hit_bound <- FALSE
  A <- K <- r <- t0 <- NA_real_
  
  tryCatch({
    fit <- nlsLM(
      y ~ logistic4(t, A, K, r, t0),
      data    = data.frame(t = fit_t, y = y),
      start   = c(A = A0, K = K0, r = r0, t0 = t00),
      lower   = lb,
      upper   = ub,
      control = nls.lm.control(maxiter = 500)
    )
    p  <- coef(fit)
    A  <- p["A"]; K <- p["K"]; r <- p["r"]; t0 <- p["t0"]
    
    eps <- 1e-5
    hit_bound <- any(
      abs(A  - lb["A"])  < eps, abs(A  - ub["A"])  < eps,
      abs(K  - lb["K"])  < eps, abs(K  - ub["K"])  < eps,
      abs(r  - lb["r"])  < eps, abs(r  - ub["r"])  < eps,
      abs(t0 - lb["t0"]) < eps, abs(t0 - ub["t0"]) < eps
    )
  }, error = function(e){
    fit_ok <<- FALSE
  })
  
  mu_max        <- if(fit_ok) r else NA_real_
  lag_time      <- if(fit_ok && !is.na(r) && r > 0) t0 - 2.0 / r else NA_real_
  max_slope_fit <- if(fit_ok) (K - A) * r / 4.0 else NA_real_
  
  t_4_10  <- c(4.5, 6.0, 7.5, 9.0, 10.0)
  y_4_10  <- as.numeric(row_vals[c("4.5h", "6h", "7.5h", "9h", "10h")])
  auc_4_10   <- trapz_local(t_4_10, y_4_10)
  odmax_4_10 <- max(y_4_10, na.rm = TRUE)
  
  tibble(
    SampleID      = row_vals[["SampleID"]],
    Condition     = row_vals[["Condition"]],
    Dose          = as.numeric(row_vals[["Dose"]]),
    A_fit         = A,
    K_fit         = K,
    r_fit         = r,
    t0_fit        = t0,
    mu_max        = mu_max,
    lag_time      = lag_time,
    max_slope_fit = max_slope_fit,
    AUC_4_10      = auc_4_10,
    ODmax_4_10    = odmax_4_10,
    fit_ok        = fit_ok,
    hit_bound     = hit_bound
  )
}

# ============================================================
# 9. 逐行拟合
# ============================================================
met <- od %>%
  rowwise() %>%
  group_map(~ fit_one_row(.x), .keep = TRUE) %>%
  bind_rows() %>%
  mutate(
    Dose_f    = factor(Dose, levels = c(0, 6.5, 13, 19.5)),
    Condition = factor(Condition, levels = c('1g', 'SMG')),
    inv_lag   = ifelse(lag_time > 0, 1 / lag_time, NA_real_)
  )

# ============================================================
# 10. 基线中心化 OD 动力学（ΔOD600）
# ============================================================
long <- od %>%
  pivot_longer(
    cols      = c(`0h`,`1h`,`2h`,`3h`,`4.5h`,`6h`,`7.5h`,`9h`,`10h`),
    names_to  = 'TimeLabel',
    values_to = 'OD600'
  ) %>%
  mutate(Time_h = as.numeric(gsub('h', '', TimeLabel)))

long <- long %>%
  group_by(SampleID) %>%
  mutate(
    OD0      = OD600[Time_h == 0][1],
    OD_delta = OD600 - OD0   # ΔOD600
  ) %>%
  ungroup()

end_delta <- long %>%
  filter(Time_h >= 4.5, Time_h <= 10) %>%
  group_by(SampleID, Condition, Dose, Dose_f) %>%
  summarise(
    AUCd_4_10 = trapz_local(Time_h, OD_delta),
    dOD_10    = OD_delta[Time_h == 10][1],
    .groups   = 'drop'
  )

met2 <- met %>%
  left_join(end_delta, by = c('SampleID', 'Condition', 'Dose', 'Dose_f'))

# ============================================================
# 11. 线性模型（用于图注）
# ============================================================
fit_mu    <- lm(mu_max   ~ Dose_f * Condition, data = met2)
fit_lag   <- lm(lag_time ~ Dose_f * Condition, data = met2)

# ANOVA结果用于图注
anova_caption_mu  <- anova_caption(fit_mu, 'μmax')
anova_caption_lag <- anova_caption(fit_lag, 'Lag time')

# emmeans成对比较（用于星号）
pair_mu <- emmeans(fit_mu, ~ Condition | Dose_f) %>%
  contrast(method = 'pairwise') %>% as.data.frame()
pair_lag <- emmeans(fit_lag, ~ Condition | Dose_f) %>%
  contrast(method = 'pairwise') %>% as.data.frame()

star <- function(p) ifelse(p < 0.01, '**', ifelse(p < 0.05, '*', ''))
pair_mu$star <- star(pair_mu$p.value)
pair_lag$star <- star(pair_lag$p.value)

# ============================================================
# 12. Panel A 汇总统计（纵坐标为ΔOD600，使用无横杠误差棒）
# ============================================================
traj_summary <- long %>%
  group_by(Condition, Dose_f, Time_h) %>%
  summarise(
    mean = ci_sem(OD_delta)['mean'],
    se = ci_sem(OD_delta)['se'],
    lwr = ci_sem(OD_delta)['lwr'],
    upr = ci_sem(OD_delta)['upr'],
    .groups = 'drop'
  )

# 原始数据点（每个时间点3个重复）
raw_points <- long %>%
  group_by(Condition, Dose_f, Time_h) %>%
  mutate(rep_id = row_number()) %>%
  ungroup()

# ============================================================
# 13. Panel B 数据准备：Condition决定颜色，Dose决定形状
# ============================================================
# 移除NA值
trade_data <- met2 %>% filter(!is.na(mu_max) & !is.na(inv_lag))

# 按Condition分别拟合回归线
fit_1g <- lm(inv_lag ~ mu_max, data = trade_data %>% filter(Condition == '1g'))
fit_smg <- lm(inv_lag ~ mu_max, data = trade_data %>% filter(Condition == 'SMG'))

# 生成预测线数据
pred_lines <- bind_rows(
  data.frame(
    mu_max = seq(min(trade_data$mu_max[trade_data$Condition == '1g'], na.rm = TRUE),
                 max(trade_data$mu_max[trade_data$Condition == '1g'], na.rm = TRUE),
                 length = 100),
    Condition = '1g'
  ) %>% mutate(inv_lag = predict(fit_1g, newdata = data.frame(mu_max = mu_max))),
  data.frame(
    mu_max = seq(min(trade_data$mu_max[trade_data$Condition == 'SMG'], na.rm = TRUE),
                 max(trade_data$mu_max[trade_data$Condition == 'SMG'], na.rm = TRUE),
                 length = 100),
    Condition = 'SMG'
  ) %>% mutate(inv_lag = predict(fit_smg, newdata = data.frame(mu_max = mu_max)))
)

# ============================================================
# 14. 统一主题设置（图例在正上方）
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
# 15. Panel A：基线中心化OD动力学曲线（ΔOD600，误差棒无横杠 + 原始数据点）
# ============================================================
pA <- ggplot(traj_summary, aes(Time_h, mean, color = Dose_f, fill = Dose_f, group = Dose_f)) +
  geom_hline(yintercept = 0, linetype = 'dashed', 
             linewidth = 0.35, color = 'grey50') +
  # 无横杠误差棒
  geom_errorbar(aes(ymin = lwr, ymax = upr),
                width = 0, linewidth = 0.4, alpha = 0.6) +
  # 原始数据点（半透明，带轻微jitter）
  geom_point(data = raw_points, aes(Time_h, OD_delta, color = Dose_f),
             alpha = 0.3, size = 0.8, 
             position = position_jitter(width = 0.08, height = 0, seed = 4001)) +
  geom_line(linewidth = 0.5) +
  geom_point(aes(y = mean), size = 1.5) +
  geom_vline(xintercept = 3, linetype = 'dashed', 
             linewidth = 0.35, color = 'grey60') +
  facet_wrap(~ Condition, nrow = 1,
             labeller = labeller(Condition = c('1g' = '1 g', 'SMG' = 'SMG'))) +
  scale_color_manual(values = col_dose, name = 'Radiation dose [Gy]') +
  scale_fill_manual(values = col_dose, guide = 'none') +
  labs(
    x = 'Time [h]',
    y = expression(Delta * OD[600])  # ΔOD600
  ) +
  base_theme

# ============================================================
# 16. Panel B：延滞期-增长率权衡散点图
# Condition决定颜色（1g=蓝色，SMG=红色）
# Dose决定形状（不同剂量不同形状）
# ============================================================
pB <- ggplot(trade_data, aes(mu_max, inv_lag)) +
  # 散点：Condition决定颜色，Dose决定形状
  geom_point(aes(color = Condition, shape = Dose_f), size = 1.8, alpha = 0.7) +
  # 回归线：Condition决定颜色（1g=蓝色，SMG=红色）
  geom_line(data = pred_lines, 
            aes(mu_max, inv_lag, color = Condition),
            linewidth = 0.6) +
  # 配色：Condition用cond配色
  scale_color_manual(values = col_cond,
                     labels = c('1 g', 'SMG')) +
  # 形状：Dose用不同形状
  scale_shape_manual(values = c('0' = 16, '6.5' = 17, '13' = 15, '19.5' = 18),
                     labels = c('0 Gy', '6.5 Gy', '13 Gy', '19.5 Gy')) +
  labs(
    x = expression(mu[max] ~ '[h'^-1~']'),
    y = expression('1/lag [h'^-1~']')
  ) +
  guides(
    color = guide_legend(title = 'Condition', order = 1),
    shape = guide_legend(title = 'Radiation dose', order = 2)
  ) +
  base_theme +
  theme(legend.position = 'top')

# ============================================================
# 17. Panel C：mu_max端点图（星号在图上方的一条水平线上）
# ============================================================
sum_mu <- met2 %>%
  group_by(Condition, Dose, Dose_f) %>%
  summarise(
    mean = ci_sem(mu_max)['mean'],
    se = ci_sem(mu_max)['se'],
    lwr = ci_sem(mu_max)['lwr'],
    upr = ci_sem(mu_max)['upr'],
    .groups = 'drop'
  )

# 计算星号位置：图上方的一条水平线
y_max_all_mu <- max(met2$mu_max, na.rm = TRUE)
y_range_mu <- y_max_all_mu - min(met2$mu_max, na.rm = TRUE)
star_y_mu <- y_max_all_mu + y_range_mu * 0.15

mu_star <- pair_mu %>%
  transmute(
    Dose = as.numeric(as.character(Dose_f)),
    star = star(p.value)
  ) %>%
  filter(star != '')

pC <- ggplot() +
  geom_point(data = met2, aes(Dose, mu_max, color = Condition),
             alpha = 0.4, size = 1.1,
             position = position_jitter(width = 0.1, height = 0, seed = 4002)) +
  geom_errorbar(data = sum_mu, 
                aes(Dose, ymin = lwr, ymax = upr, color = Condition),
                position = position_dodge(width = 0.5), 
                width = 0, linewidth = 0.5, alpha = 0.8) +
  geom_line(data = sum_mu, 
            aes(Dose, mean, color = Condition, group = Condition),
            linewidth = 0.5) +
  geom_point(data = sum_mu, aes(Dose, mean, color = Condition),
             size = 2, shape = 16) +
  geom_text(data = mu_star, aes(Dose, star_y_mu, label = star),
            inherit.aes = FALSE, size = 8/.pt, fontface = 'bold', vjust = 0.5) +
  scale_color_manual(values = col_cond, labels = c('1 g', 'SMG')) +
  scale_x_continuous(breaks = c(0, 6.5, 13, 19.5), 
                     labels = c('0', '6.5', '13', '19.5')) +
  coord_cartesian(ylim = c(min(met2$mu_max, na.rm = TRUE) - y_range_mu * 0.05,
                           star_y_mu + y_range_mu * 0.1)) +
  labs(x = 'Radiation dose [Gy]', 
       y = expression(mu[max] ~ '[h'^-1~']')) +
  base_theme

# ============================================================
# 18. Panel D：lag_time端点图（星号在图上方的一条水平线上）
# ============================================================
sum_lag <- met2 %>%
  group_by(Condition, Dose, Dose_f) %>%
  summarise(
    mean = ci_sem(lag_time)['mean'],
    se = ci_sem(lag_time)['se'],
    lwr = ci_sem(lag_time)['lwr'],
    upr = ci_sem(lag_time)['upr'],
    .groups = 'drop'
  )

y_max_all_lag <- max(met2$lag_time, na.rm = TRUE)
y_range_lag <- y_max_all_lag - min(met2$lag_time, na.rm = TRUE)
star_y_lag <- y_max_all_lag + y_range_lag * 0.15

lag_star <- pair_lag %>%
  transmute(
    Dose = as.numeric(as.character(Dose_f)),
    star = star(p.value)
  ) %>%
  filter(star != '')

pD <- ggplot() +
  geom_point(data = met2, aes(Dose, lag_time, color = Condition),
             alpha = 0.4, size = 1.1,
             position = position_jitter(width = 0.1, height = 0, seed = 4003)) +
  geom_errorbar(data = sum_lag, 
                aes(Dose, ymin = lwr, ymax = upr, color = Condition),
                position = position_dodge(width = 0.5), 
                width = 0, linewidth = 0.5, alpha = 0.8) +
  geom_line(data = sum_lag, 
            aes(Dose, mean, color = Condition, group = Condition),
            linewidth = 0.5) +
  geom_point(data = sum_lag, aes(Dose, mean, color = Condition),
             size = 2, shape = 16) +
  geom_text(data = lag_star, aes(Dose, star_y_lag, label = star),
            inherit.aes = FALSE, size = 8/.pt, fontface = 'bold', vjust = 0.5) +
  scale_color_manual(values = col_cond, labels = c('1 g', 'SMG')) +
  scale_x_continuous(breaks = c(0, 6.5, 13, 19.5), 
                     labels = c('0', '6.5', '13', '19.5')) +
  coord_cartesian(ylim = c(min(met2$lag_time, na.rm = TRUE) - y_range_lag * 0.05,
                           star_y_lag + y_range_lag * 0.1)) +
  labs(x = 'Radiation dose [Gy]', 
       y = 'Lag time [h]') +
  base_theme

# ============================================================
# 19. 拼合：2x2布局（A、B、C、D）；统计结果由正文报告
# ============================================================
fig4_revised <- (pA | pB) / (pC | pD) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 9, face = 'bold', family = 'Arial'))

print(fig4_revised)

# ============================================================
# 20. 导出双栏图（宽度 19 cm）
# ============================================================
ggsave(file.path(outdir, 'Figure4_submission_v16.pdf'), fig4_revised,
       width = 18.3, height = 15.4, units = 'cm', device = cairo_pdf)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(outdir, 'Figure4_submission_v16.svg'), fig4_revised,
         width = 18.3, height = 15.4, units = 'cm', device = svglite::svglite)
}
ggsave(file.path(outdir, 'Figure4_submission_v16.tiff'), fig4_revised,
       width = 18.3, height = 15.4, units = 'cm', dpi = 600,
       device = ragg::agg_tiff, compression = 'lzw')
ggsave(file.path(outdir, 'Figure4_submission_v16.png'), fig4_revised,
       width = 18.3, height = 15.4, units = 'cm', dpi = 600,
       device = ragg::agg_png)

# ============================================================
# 21. 导出所有绘图数据
# ============================================================

export_data_fig4 <- list()

# 1. 原始OD数据（基线中心化后的长格式，ΔOD600）
export_data_fig4$raw_OD_delta <- long

# 2. Panel A汇总统计（轨迹均值±SEM）
export_data_fig4$panel_A_trajectory_summary <- traj_summary

# 3. Panel B数据（mu_max和inv_lag）
export_data_fig4$panel_B_trade_data <- trade_data

# 4. Panel B回归线预测数据
export_data_fig4$panel_B_regression_lines <- pred_lines

# 5. Panel C mu_max汇总统计
export_data_fig4$panel_C_mu_summary <- sum_mu

# 6. Panel C mu_max显著性星号
export_data_fig4$panel_C_mu_stars <- mu_star

# 7. Panel D lag_time汇总统计
export_data_fig4$panel_D_lag_summary <- sum_lag

# 8. Panel D lag_time显著性星号
export_data_fig4$panel_D_lag_stars <- lag_star

# 9. 拟合指标完整数据
export_data_fig4$fitted_metrics <- met2

# 10. ANOVA完整结果
export_data_fig4$anova_mu_full <- as.data.frame(car::Anova(fit_mu, type = 3)) %>%
  mutate(term = rownames(.))
export_data_fig4$anova_lag_full <- as.data.frame(car::Anova(fit_lag, type = 3)) %>%
  mutate(term = rownames(.))

# 11. emmeans成对比较结果
export_data_fig4$emmeans_mu <- pair_mu
export_data_fig4$emmeans_lag <- pair_lag

# 保存为Excel文件
if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  wb_fig4 <- createWorkbook()
  for (name in names(export_data_fig4)) {
    addWorksheet(wb_fig4, name)
    writeData(wb_fig4, name, export_data_fig4[[name]])
  }
  saveWorkbook(wb_fig4, file.path(outdir, "Fig4_Data_V6.xlsx"), overwrite = TRUE)
  cat('数据已保存为: Fig4_Data_V6.xlsx\n')
} else {
  data_dir <- file.path(outdir, "Fig4_Data_V6")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(export_data_fig4)) {
    write_csv(export_data_fig4[[name]], file.path(data_dir, paste0(name, ".csv")))
  }
  cat('数据已保存至文件夹:', data_dir, '\n')
}

# ============================================================
# 22. 完成信息
# ============================================================
cat('\n========== 完成 ==========\n')
cat('版本: V6\n')
cat('图片已保存为: Fig4_Revised_V6.pdf (双栏, 19 cm width)\n')
cat('纵坐标: Panel A 使用 ΔOD600\n')
cat('布局：2x2布局 (A | B) / (C | D)\n')
cat('图例位置：正上方\n')
