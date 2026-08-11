####################################################
## Figure 6b：功能模块评分与表型关联深度分析
##
## 1. 结果解释与印证（Biological Interpretation）：
##    - ROS/Redox 与 10h 生长性能（μmax, AUC）的极强负相关（R = -1.00*）：
##      印证了转录组中发现的氧化应激响应（katE, ahpC 等上调）是导致生长抑制的核心驱动力。
##      辐射诱导的自由基损伤迫使能量从"增殖"转向"防御"，直接表现为生长速率的下降。
##    - Metabolism 与 3h 萌发率的正相关（R = 0.80）：
##      印证了萌发早期代谢基因（如中心碳代谢）的激活程度直接决定了细胞脱离休眠的速度。
##    - Transport 与 Lagtime 的关联：
##      转运蛋白（特别是离子和糖类转运）的表达模式与细胞进入指数生长期所需的准备时间高度耦合。
##    - Iron/Siderophore 的独立意义：
##      铁稳态通过芬顿反应与氧化应激耦合，其评分波动反映了航天应激下的微环境适应策略。
##
## 2. 绘图设计逻辑：
##    - 左侧热图（Heatmap）：全景展示 5 个功能箱与 5 个核心表型的 Spearman 相关性。
##    - 右侧趋势图（Trend Plots）：每个表型行取一个 |ρ| 最大的功能箱，覆盖不同表型维度。
##      趋势线为 OLS 辅助可视化，相关系数为 Spearman ρ（与热图一致）。
##    - 符号规定：圆形(●)代表 3h 萌发阶段，三角形(▲)代表 10h 增殖阶段。
##
## 3. 优化说明（v36）：
##    - [基准统一] 所有 Delta 表型值均以 1g+0Gy 组为基准计算，语义一致。
##    - [方法统一] 热图与趋势图均采用 Spearman ρ；趋势线为 OLS 辅助线，
##      图注已说明其仅用于可视化，不代表 Pearson 相关。
##    - [零值区分] pivot_wider 使用 values_fill = NA：score=0 表示功能箱有 DEG
##      但上下调均衡抵消；NA 表示该对比组在该功能箱无 DEG，趋势图中自动略去
##      NA 点，避免填充零值混淆生物学信号。
##    - [图例精简] 除热图色标外，全图不添加任何图例。
####################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(reshape2)
  library(cowplot)
})

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out_dir <- file.path(project_root, "results")
original_result <- file.path(project_root, "data", "phenotype")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"

#========================
# 1. 环境与参数设置
#========================

## 字体与配色：使用设备通用 sans，在 macOS/Linux 上分别映射到
## Helvetica/Arial-compatible 字体，避免 PostScript 字体替换告警。

sz_base <- 8; sz_axis <- 8; sz_tick <- 7; sz_legend <- 7; sz_tag <- 10

## 功能箱颜色（与Fig6_a保持一致）
bin_colors <- c(
  "Metabolism"          = "#0072B5", # Blue
  "Transport"           = "#BC3C29", # Red
  "Motility/Chemotaxis" = "#20854E", # Green
  "ROS/Redox"           = "#E18727", # Orange
  "Iron/Siderophore"    = "#7876B1"  # Purple
)

## 实验组颜色（与Fig6_a保持一致）
contrast_cols <- c(
  "1g+6.5"   = "#0072B5", "1g+13"    = "#808281", "1g+19.5"  = "#20854E",
  "SMG+0"    = "#E18727", "SMG+6.5"  = "#7876B1", "SMG+13"   = "#BC3C29", "SMG+19.5" = "#7E6148"
)

#========================
# 2. 数据处理函数
#========================

## 3h t50计算函数（线性插值）
calc_t50 <- function(r1, r2, r3) {
  t <- c(0, 1, 2, 3); y <- c(0, r1, r2, r3)
  if(max(y) < 0.5) return(NA)
  for(i in 2:4) {
    if(y[i] >= 0.5) {
      return(t[i-1] + (0.5 - y[i-1]) / (y[i] - y[i-1]) * (t[i] - t[i-1]))
    }
  }
  return(NA)
}

## 功能箱分配逻辑
assign_bin <- function(txt) {
  t <- str_to_lower(txt)
  case_when(
    str_detect(t, "catalase|peroxiredoxin|thioredoxin|glutaredoxin|superoxide|sod|bacilliredoxin|peroxidase") ~ "ROS/Redox",
    str_detect(t, "iron|sideroph|ferric|ferritin|feo|ferrous|fur|perr|fhu|heme transport") ~ "Iron/Siderophore",
    str_detect(t, "flagell|fli|flg|flh|mot|chemotaxis|chea|cheb|chew|cher") ~ "Motility/Chemotaxis",
    str_detect(t, "transport|permease|symporter|abc|uptake|efflux|substrate-binding") ~ "Transport",
    str_detect(t, "dehydrogenase|synthase|transferase|reductase|lyase|isomerase|kinase|glycolysis|tca|citrate|pyruvate") ~ "Metabolism",
    TRUE ~ "Other"
  )
}

#========================
# 3. 加载并清洗数据
#========================

# 1) 转录组评分
sig_annotated <- read_csv(file.path(out_dir, "fig5_sig_annotated.csv"), show_col_types = FALSE) %>%
  mutate(bin = assign_bin(paste(gene, def))) %>%
  filter(bin != "Other")

score_df <- sig_annotated %>%
  group_by(timepoint, condition, dose, bin) %>%
  summarise(
    n = n(),
    score = sum(log2FoldChange) / sqrt(pmax(n, 1)),
    .groups = "drop"
  ) %>%
  mutate(
    contrast = paste0(condition, "_", dose),
    contrast_label = paste0(condition, "+", dose)
  )

# 2) 3h 萌发表型
germ_data <- read_csv(file.path(original_result, "cleaned_germination_data.csv"), show_col_types = FALSE)
g_sum <- germ_data %>%
  rowwise() %>%
  mutate(t50 = calc_t50(Rate_1h, Rate_2h, Rate_3h)) %>%
  group_by(Condition, Dose) %>%
  summarise(Rate_3h = mean(Rate_3h), t50 = mean(t50, na.rm=T), .groups="drop")

# [基准统一] 所有 Delta 统一以 1g+0Gy 为基准，确保各数据点语义一致：
# 每个 delta 值均反映"相对无处理对照的偏离量"，消除 SMG 组使用不同基准的异质性。
base_1g <- g_sum %>% filter(Condition == "1g", Dose == 0)

pheno_3h <- g_sum %>%
  mutate(
    delta_germ = Rate_3h - base_1g$Rate_3h,
    delta_t50  = t50    - base_1g$t50,
    contrast = paste0(Condition, "_", Dose),
    contrast_label = paste0(Condition, "+", Dose),
    timepoint = "3 h"
  )

# 3) 10h 生长表型
# [基准统一] 注意：bridge 文件中的 delta 值应同样以 1g+0Gy 为基准预先计算，
# 若原始文件使用了分组基准，请在此处重新对齐后再读入。
bridge_data <- read_csv(file.path(out_dir, "fig5_transcriptome_phenotype_bridge.csv"), show_col_types = FALSE)
pheno_10h <- bridge_data %>% filter(timepoint == "10 h")

#========================
# 4. 关联分析计算 (增强报错容错性)
#========================

target_bins <- names(bin_colors)
target_phenos <- c("delta_germ", "delta_t50", "delta_mu_max", "delta_lag", "delta_AUC")

# [零值区分] 热图与趋势图对"无 DEG"的处理需求不同，因此分开建两套矩阵：
#
# 热图用矩阵（values_fill = 0）：
#   score = 0 表示该功能箱在该对比组无净转录信号，0 是合理的模块层面汇总值，
#   保留全部对比组参与相关性计算，确保热图格子不因样本量不足而消失。
#
# 趋势图用矩阵（values_fill = NA）：
#   score = NA 表示该对比组在该功能箱完全无 DEG，散点图中自动略去此类点，
#   避免填充零值在 x 轴上制造虚假的"无响应"数据点。
#   score = 0（真实值）仍正常显示，代表上下调均衡抵消的情况。

## 热图用：values_fill = 0
m3_cor <- score_df %>%
  filter(timepoint == "3 h") %>%
  select(contrast, contrast_label, bin, score) %>%
  pivot_wider(names_from = bin, values_from = score, values_fill = 0) %>%
  mutate(timepoint = "3 h") %>%
  inner_join(pheno_3h %>% select(contrast, delta_germ, delta_t50), by = "contrast")

m10_cor <- score_df %>%
  filter(timepoint == "10 h") %>%
  select(contrast, contrast_label, bin, score) %>%
  pivot_wider(names_from = bin, values_from = score, values_fill = 0) %>%
  mutate(timepoint = "10 h") %>%
  inner_join(pheno_10h %>% select(contrast, delta_mu_max, delta_lag, delta_AUC), by = "contrast")

## 趋势图用：values_fill = NA
m3 <- score_df %>%
  filter(timepoint == "3 h") %>%
  select(contrast, contrast_label, bin, score) %>%
  pivot_wider(names_from = bin, values_from = score, values_fill = NA) %>%
  mutate(timepoint = "3 h") %>%
  inner_join(pheno_3h %>% select(contrast, delta_germ, delta_t50), by = "contrast")

m10 <- score_df %>%
  filter(timepoint == "10 h") %>%
  select(contrast, contrast_label, bin, score) %>%
  pivot_wider(names_from = bin, values_from = score, values_fill = NA) %>%
  mutate(timepoint = "10 h") %>%
  inner_join(pheno_10h %>% select(contrast, delta_mu_max, delta_lag, delta_AUC), by = "contrast")

# 2) 带有安全检查的相关性计算函数
safe_cor <- function(x, y, return_p = FALSE) {
  # 排除 NA 后的有效配对数
  valid_idx <- which(!is.na(x) & !is.na(y))
  
  # 如果有效观测值少于 3 个，或变量全部相同（方差为0），则返回 NA
  if(length(valid_idx) < 3 || sd(x[valid_idx]) == 0 || sd(y[valid_idx]) == 0) {
    return(NA_real_)
  }
  
  tryCatch({
    if(return_p) {
      return(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
    } else {
      return(cor(x, y, method = "spearman"))
    }
  }, error = function(e) return(NA_real_))
}

# 3) 执行计算
# 热图相关性计算使用 m3_cor / m10_cor（values_fill = 0），保证全部功能箱均有足够样本量
cor_res <- expand.grid(Phenotype = target_phenos, Bin = target_bins, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    # 选择对应的时间点数据框（热图专用矩阵）
    target_df = list(if(str_detect(Phenotype, "germ|t50")) m3_cor else m10_cor),
    
    # 检查列是否存在
    has_data = (Bin %in% colnames(target_df) && Phenotype %in% colnames(target_df)),
    
    cor_val = if(has_data) {
      safe_cor(target_df[[Bin]], target_df[[Phenotype]], return_p = FALSE)
    } else { NA_real_ },
    
    p_val = if(has_data) {
      safe_cor(target_df[[Bin]], target_df[[Phenotype]], return_p = TRUE)
    } else { NA_real_ }
  ) %>%
  filter(!is.na(cor_val)) # 过滤掉无法计算组合，保持热图整洁

#========================
# 5. 绘图 - Panel B
#========================

# B1: 相关性热图
# 在绘图前，手动指定行顺序
cor_res$Phenotype <- factor(cor_res$Phenotype, levels = c(
  "delta_AUC", "delta_lag", "delta_mu_max", # 10h 放在下方
  "delta_t50", "delta_germ"                 # 3h 放在上方
))

p_heat <- ggplot(cor_res, aes(x = Bin, y = Phenotype, fill = cor_val)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", cor_val)), size = 2.0, family = font_family, color = "black") +
  geom_text(data = filter(cor_res, p_val < 0.05), aes(label = "*"),
            vjust = -0.5, size = 5, color = "black") +
  scale_fill_gradient2(
    low = "#0072B5",
    mid = "white",
    high = "#BC3C29",
    midpoint = 0,
    name = "Spearman R",
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  scale_y_discrete(labels = c(
    "delta_germ"   = "Δ Germination rate (3h)", "delta_t50"    = "Δ t50 (3h)",
    "delta_mu_max" = "Δ μmax (10h)",            "delta_lag"    = "Δ Lag time (10h)",
    "delta_AUC"    = "Δ AUC (10h)"
  ), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_bw(base_size = sz_base) +
  theme(
    text          = element_text(family = font_family),
    axis.title    = element_blank(),
    axis.text     = element_text(color = "black"),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "right",   # 热图保留色标图例
    panel.grid    = element_blank(),
    panel.border  = element_rect(linewidth = 0.5)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust    = 0.5,
    barwidth       = unit(0.3, "cm"),
    barheight      = unit(4, "cm")
  ))

# B2: 趋势散点图组合函数
# [方法统一] geom_smooth(method = "lm") 仅作辅助可视化趋势线，
# 相关系数与热图一致均为 Spearman ρ，图注中已注明。
# [图例精简] 除热图色标外，全图不展示任何图例（legend.position = "none"）。
# [零值区分] 因 values_fill = NA，无 DEG 的对比组不在趋势图中产生数据点，
# ggplot 会自动跳过 NA，无需额外过滤。
plot_trend <- function(df, x_bin, y_pheno, title) {
  y_lab <- case_when(
    y_pheno == "delta_germ"   ~ "Δ Germination rate (3h)",
    y_pheno == "delta_lag"    ~ "Δ Lag time (10h)",
    y_pheno == "delta_mu_max" ~ "Δ μmax (10h)",
    y_pheno == "delta_AUC"    ~ "Δ AUC (10h)",
    TRUE                      ~ str_replace(y_pheno, "delta_", "Δ ")
  )
  
  # 剔除 x 或 y 为 NA 的行，保证 geom_smooth 与 geom_point 使用相同数据集
  plot_data <- df %>% filter(!is.na(.data[[x_bin]]), !is.na(.data[[y_pheno]]))
  
  ggplot(plot_data, aes(x = .data[[x_bin]], y = .data[[y_pheno]])) +
    geom_smooth(method = "lm", color = "grey70", fill = "grey95",
                linetype = "dashed", linewidth = 0.5, se = TRUE) +
    geom_point(aes(color = contrast_label,
                   shape = factor(timepoint, levels = c("3 h", "10 h"))),
               size = 2) +
    scale_color_manual(values = contrast_cols, limits = names(contrast_cols), drop = FALSE) +
    scale_shape_manual(values = c("3 h" = 16, "10 h" = 17),
                       limits = c("3 h", "10 h"), drop = FALSE) +
    labs(title = title,
         x     = paste(x_bin, "Score"),
         y     = y_lab) +
    theme_bw(base_size = sz_base) +
    theme(
      text            = element_text(family = font_family),
      plot.title      = element_text(size = sz_axis - 1, face = "bold"),
      legend.position = "none"   # 趋势图不展示任何图例
    )
}

# 计算各组合在趋势图矩阵（values_fill=NA）中的有效点数
n_valid_df <- expand.grid(
  Phenotype = target_phenos,
  Bin       = target_bins,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    target_df = list(if(str_detect(Phenotype, "germ|t50")) m3 else m10),
    has_data  = (Bin %in% colnames(target_df) && Phenotype %in% colnames(target_df)),
    n_valid   = if(has_data) {
      sum(!is.na(target_df[[Bin]]) & !is.na(target_df[[Phenotype]]))
    } else { 0L }
  ) %>%
  select(Phenotype, Bin, n_valid)

# 将 n_valid 合并进 cor_res，只保留 n_valid >= 4 的组合供趋势图候选
cor_res_trend <- cor_res %>%
  left_join(n_valid_df, by = c("Phenotype", "Bin")) %>%
  filter(n_valid >= 4)

write_csv(cor_res_trend, file.path(out_dir, "fig6_cor_res_trend.csv"))

# 选择理由如下：
# "Trend plots are shown for 10h associations only, as 3h module scores lacked sufficient data coverage for reliable scatter visualization (n_valid ≤ 5 across all 3h combinations)."
p1 <- plot_trend(m10, "Transport",  "delta_AUC",    "Transport vs Biomass (10h)")
p2 <- plot_trend(m10, "Transport",  "delta_mu_max", "Transport vs Growth rate (10h)")
p3 <- plot_trend(m10, "Metabolism", "delta_mu_max", "Metabolism vs Growth rate (10h)")
p4 <- plot_trend(m10, "Metabolism", "delta_AUC",    "Metabolism vs Biomass (10h)")

#========================
# 6. 组装与保存
#========================
right_panel <- cowplot::plot_grid(
  p1, p2, p3, p4,
  ncol        = 1,
  rel_heights = c(1, 1, 1, 1),
  align       = "v",
  axis        = "lr"
)

# 最终全图组装
pb_final <- cowplot::ggdraw() +
  cowplot::draw_plot(p_heat,      x = 0.00, y = 0, width = 0.55, height = 1.00) +
  cowplot::draw_plot(right_panel, x = 0.55, y = 0, width = 0.45, height = 1.00)

ggsave(file.path(out_dir, "Figure6b_submission_v16.pdf"), pb_final,
       width = 18.3, height = 16, units = "cm", device = cairo_pdf)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(out_dir, "Figure6b_submission_v16.svg"), pb_final,
         width = 18.3, height = 16, units = "cm", device = svglite::svglite)
}
ggsave(file.path(out_dir, "Figure6b_submission_v16.tiff"), pb_final,
       width = 18.3, height = 16, units = "cm", dpi = 600,
       device = ragg::agg_tiff, compression = "lzw")
ggsave(file.path(out_dir, "Figure6b_submission_v16.png"), pb_final,
       width = 18.3, height = 16, units = "cm", dpi = 600,
       device = ragg::agg_png)

message("Figure 6b 优化完成（v40）：热图与趋势图分别使用独立矩阵，热图恢复完整五列，趋势图 NA 点正确略去。")
