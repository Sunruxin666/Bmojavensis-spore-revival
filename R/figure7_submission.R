####################################################
## Figure 7：第二代芽孢消杀表型 × 转录组耦合可视化
##
## 目的：
## 1）展示二代孢子无消杀剂基础生长（应激历史记忆），并突出 Ctrl 与 SMG+13 组；
## 2）量化 H₂O₂ 和庆大霉素对各组二代孢子的消杀动力学；
## 3）通过 AUC 条形图横向比较消杀效力；
## 4）将消杀 AUC 与10h转录组4个功能模块耦合，揭示分子机制。
##
##
##
## 前置文件（均放工作目录）：
##   control_1.xlsx / control_2.xlsx     ← 二代孢子纯生长基线
##   h2o2.xlsx / gentamicin.xlsx         ← 正文消杀实验（主图）
##   DDAC.xlsx / soda_1.xlsx             ← 补充图消杀实验
##   arial.ttf                           ← Arial 字体文件
##   Result/fig5_sig_annotated.csv       ← Fig5.R 输出的10h DEG 注释表
##
## 输出：
##   Result/Figure7.pdf       双栏，19 cm，正文主图（Panel A–D）
##   Result/Figure7_S1.pdf    双栏，19 cm，补充图（DDAC + NaHCO₃）
##   Result/fig7_auc_table.csv           净消杀 AUC 数值汇总
####################################################

#========================================================================
# Panel D 双热图说明（左：功能箱 右：精选通路）
#========================================================================
# 
# [背景与目的]
# Panel D 旨在验证 10 h 转录模块/通路活性是否能预示二代孢子对不同
# 杀菌剂的耐受性。通过分别对“辐射”(1g+剂量) 和“微重力+辐射”
# (SMG+剂量) 两种预处理类型内部进行跨组相关性分析，定量评估功能模块
# 或代谢通路活性 (score_norm) 与净杀伤 (net_kill) 的关联强度。
#
# [左图：功能箱热图]
# 1. 数据来源：Fig5 输出的显著 DEG 注释表 (fig5_sig_annotated.csv)，
#    使用与 Fig6 完全相同的 assign_bin 函数，将每个显著基因划入 5 个
#    功能箱 (Metabolism, Transport, Motility/Chemotaxis, ROS/Redox,
#    Iron/Siderophore)，并计算 10 h 各箱的 score_norm。
# 2. 筛选：仅保留 Metabolism, Transport, Iron/Siderophore 三个模块
#    (删去 Motility 和 ROS)，以聚焦能量/底物转运/铁稳态。
# 3. 分类处理：去掉 “SMG alone” (仅1组，无法计算相关性)，仅保留
#    “Radiation” (1g+6.5/13/19.5) 和 “SMG + Radiation”
#    (SMG+6.5/13/19.5)，每组 N=3。
# 4. 关联计算：对每个模块×杀菌剂×处理类型组合，计算 score_norm 与
#    net_kill 的 Pearson 相关系数 r 及 p 值。若某组合内变量无变化则
#    记为 NA (显示为灰色 “–”)。
# 5. 热图绘制：x 轴为 “杀菌剂\n处理类型”，y 轴为功能箱，填充色为 r
#    (蓝: 负相关，红: 正相关，白: 0)，格内标注 r 值、显著性 (*p<0.05,
#    **p<0.01) 和样本量 (N=3)。abs(r) > 0.6 时文字换白色以保证可读。
#
# [右图：通路热图]
# 1. 数据来源：Fig5 输出的 KEGG 通路评分表 (fig56_kegg_pathway_scores.csv)，
#    提取 10 h 所有对比组的 score_norm，并按实验组统一为 group 名。
#    与 net_all (消杀净杀伤表) 合并。
# 2. 筛选通路：仅保留 Three-component system (双组分系统)、
#    Carbon metabolism (碳代谢)、PTS (磷酸转移酶系统) 三条与感知/
#    能量/糖摄取密切相关的通路。
# 3. 分类处理、关联计算和绘图参数与左图完全一致，确保可比性。
#
# [辅助说明]
# - Safe correlation: 自定义 safe_cor 函数先用 unique 检验变量方差，
#   再用 tryCatch 捕获任何错误，避免某组合因常数导致脚本中断。
# - Ctrl 组没有 score_norm (logFC 为 0)，因此不参与相关性计算，但
#   其在两种杀菌剂下的 net_kill 值可通过 AUC 表获知，作为基线参考。
# - 两图共用同一个色标图例，放在右侧，避免重复。

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)   # 用于提取图例
  library(scales)
})

set.seed(42)

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
source_dir <- file.path(project_root, "data", "secondary_spores")
out_dir <- file.path(project_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"

#========================
# 1. 全局参数（与原脚本一致）
#========================

fig_w_double <- 18.3
fig_res      <- 600
sz_base      <- 8
sz_axis      <- 8
sz_tick      <- 7
sz_legend    <- 7
sz_tag       <- 9
sz_annot     <- 6.5

# 主图 A–C 展示的 5 组
main_groups <- c("Ctrl", "1g+13", "1g+19.5", "SMG+6.5", "SMG+13")
main_colors <- c(
  "Ctrl"    = "#7BA3C7",
  "1g+13"   = "#4A7BA7",
  "1g+19.5" = "#2E5A8B",
  "SMG+6.5" = "#E8956F",
  "SMG+13"  = "#B91C1C"
)
main_linetypes <- c(
  "Ctrl"    = "solid",
  "1g+13"   = "dashed",
  "1g+19.5" = "dotted",
  "SMG+6.5" = "dashed",
  "SMG+13"  = "solid"
)
main_linewidths <- c(
  "Ctrl"    = 0.65,
  "1g+13"   = 0.50,
  "1g+19.5" = 0.50,
  "SMG+6.5" = 0.50,
  "SMG+13"  = 0.75
)
main_labels <- c(
  "Ctrl"    = "Ctrl (1g)",
  "1g+13"   = "1g+13Gy",
  "1g+19.5" = "1g+19.5Gy",
  "SMG+6.5" = "SMG+6.5Gy",
  "SMG+13"  = "SMG+13Gy"
)

# 全组配色（用于 AUC bar 和热图）
all_group_colors <- c(
  "Ctrl"     = "#7BA3C7",
  "1g+6.5"   = "#5A8DB8",
  "1g+13"    = "#4A7BA7",
  "1g+19.5"  = "#2E5A8B",
  "SMG"      = "#F0A070",
  "SMG+6.5"  = "#E8956F",
  "SMG+13"   = "#B91C1C",
  "SMG+19.5" = "#A02800"
)
all_group_levels <- names(all_group_colors)

# 统一主题（与 Fig5 / Fig6 保持一致）
theme_journal <- function() {
  theme_bw(base_size = sz_base, base_family = font_family) +
    theme(
      axis.title         = element_text(size = sz_axis, family = font_family),
      axis.text          = element_text(size = sz_tick, family = font_family),
      axis.line          = element_line(linewidth = 0.4),
      plot.margin = margin(2, 2, 2, 2, "pt"),
      panel.border       = element_rect(linewidth = 0.4, fill = NA),
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.text        = element_text(size = sz_legend, family = font_family),
      legend.title       = element_text(size = sz_legend, family = font_family),
      legend.key.size    = unit(0.35, "cm"),
      strip.text         = element_text(size = sz_tick, family = font_family),
      strip.background   = element_rect(fill = "grey95", color = "grey70", linewidth = 0.3)
    )
}

#========================
# 2. 路径验证
#========================

req <- c("control_1.xlsx", "control_2.xlsx", "h2o2.xlsx", "gentamicin.xlsx",
         "DDAC.xlsx", "soda_1.xlsx")
for (f in req) {
  if (!file.exists(file.path(source_dir, f))) stop("缺失文件：", f)
}
for (f in c("fig5_sig_annotated.csv", "fig56_kegg_pathway_scores.csv")) {
  if (!file.exists(file.path(out_dir, f))) stop("缺失文件：", f)
}
message("路径验证通过！")

#========================
# 3. 样本列名元数据
#========================

group_cols <- list(
  "Ctrl"     = c("Ctrl-1", "Ctrl-2", "Ctrl-3"),
  "SMG"      = c("SMG-1", "SMG-2", "SMG-3"),
  "1g+6.5"   = c("SR（6.5 Gy）-1", "SR（6.5 Gy）-2", "SR（6.5 Gy）-3"),
  "1g+13"    = c("SR（13 Gy）-1", "SR（13 Gy）-2", "SR（13 Gy）-3"),
  "1g+19.5"  = c("SR（19.5 Gy）-1", "SR（19.5 Gy）-2", "SR（19.5 Gy）-3"),
  "SMG+6.5"  = c("SMG+SR（6.5 Gy）-1", "SMG+SR（6.5 Gy）-2", "SMG+SR（6.5 Gy）-3"),
  "SMG+13"   = c("SMG+SR（13 Gy）-1", "SMG+SR（13 Gy）-2", "SMG+SR（13 Gy）-3"),
  "SMG+19.5" = c("SMG+SR（19.5 Gy）-1", "SMG+SR（19.5 Gy）-2", "SMG+SR（19.5 Gy）-3")
)

#========================
# 4. 读取与处理 OD 数据（Panel A–C 所用）
#========================

read_od <- function(path) {
  read_excel(path) %>% mutate(across(-time, as.numeric))
}

wide_to_rel_long <- function(df_wide, grp_cols, treatment) {
  time_vec <- df_wide$time
  bind_rows(lapply(names(grp_cols), function(grp) {
    cols <- intersect(grp_cols[[grp]], colnames(df_wide))
    if (length(cols) == 0) return(NULL)
    sub <- df_wide[, cols, drop = FALSE]
    od0 <- as.numeric(sub[1, ])
    rel <- sweep(sub, 2, od0, `-`) |> sweep(2, od0, `/`)
    rel |>
      mutate(time = time_vec) |>
      pivot_longer(-time, names_to = "rep", values_to = "rel_od") |>
      mutate(group = grp, treatment = treatment)
  }))
}

summarise_rel <- function(long_df) {
  long_df |>
    group_by(treatment, group, time) |>
    summarise(mean_rel = mean(rel_od, na.rm = TRUE),
              sd_rel   = sd(rel_od, na.rm = TRUE),
              .groups  = "drop") |>
    mutate(time_h = time / 60,
           group  = factor(group, levels = all_group_levels))
}

calc_auc_rep <- function(long_df) {
  long_df |>
    group_by(treatment, group, rep) |>
    arrange(time) |>
    summarise(
      auc = sum(diff(time) * (head(rel_od, -1) + tail(rel_od, -1)) / 2) /
        (max(time) - min(time)),
      .groups = "drop"
    )
}

calc_auc_group_mean <- function(sum_df) {
  sum_df |>
    group_by(treatment, group) |>
    arrange(time) |>
    summarise(
      auc = sum(diff(time) * (head(mean_rel, -1) + tail(mean_rel, -1)) / 2) /
        (max(time) - min(time)),
      .groups = "drop"
    )
}

message("\n===== 读取 OD600 数据 =====")
rl_c1   <- wide_to_rel_long(read_od(file.path(source_dir, "control_1.xlsx")), group_cols, "ctrl1")
rl_c2   <- wide_to_rel_long(read_od(file.path(source_dir, "control_2.xlsx")), group_cols, "ctrl2")
rl_h2o2 <- wide_to_rel_long(read_od(file.path(source_dir, "h2o2.xlsx")), group_cols, "H2O2")
rl_gent <- wide_to_rel_long(read_od(file.path(source_dir, "gentamicin.xlsx")), group_cols, "Gent")
rl_ddac <- wide_to_rel_long(read_od(file.path(source_dir, "DDAC.xlsx")), group_cols, "DDAC")
rl_soda <- wide_to_rel_long(read_od(file.path(source_dir, "soda_1.xlsx")), group_cols, "NaHCO3")

# Author decision: treat control_1 and control_2 as one combined baseline
# dataset. The six traces are measurement traces, not six independent
# biological replicates. The baseline AUC is therefore calculated from the
# combined mean trajectory, matching the original analysis convention.
rl_base <- bind_rows(rl_c1, rl_c2) %>% mutate(treatment = "Baseline")

sum_base <- summarise_rel(rl_base)
sum_h2o2 <- summarise_rel(rl_h2o2)
sum_gent <- summarise_rel(rl_gent)
sum_ddac <- summarise_rel(rl_ddac)
sum_soda <- summarise_rel(rl_soda)

auc_base <- calc_auc_group_mean(sum_base) |> rename(auc_base = auc)
auc_h2o2 <- calc_auc_rep(rl_h2o2) |> rename(auc_biocide = auc)
auc_gent <- calc_auc_rep(rl_gent) |> rename(auc_biocide = auc)
auc_ddac <- calc_auc_rep(rl_ddac) |> rename(auc_biocide = auc)
auc_soda <- calc_auc_rep(rl_soda) |> rename(auc_biocide = auc)

calc_net <- function(auc_treat, biocide_name) {
  auc_base |>
    select(group, auc_base) |>
    left_join(auc_treat |> select(group, rep, auc_biocide), by = "group") |>
    mutate(net_kill = auc_base - auc_biocide,
           biocide   = biocide_name)
}

net_h2o2 <- calc_net(auc_h2o2, "H2O2")
net_gent <- calc_net(auc_gent, "Gentamicin")
net_ddac <- calc_net(auc_ddac, "DDAC")
net_soda <- calc_net(auc_soda, "NaHCO3")

net_all_rep <- bind_rows(net_h2o2, net_gent, net_ddac, net_soda)

net_all_rep |>
  select(biocide, group, rep, auc_base, auc_biocide, net_kill) |>
  write_csv(file.path(out_dir, "Figure7_replicate_net_kill_v18.csv"))

net_all_rep |>
  group_by(biocide, group) |>
  summarise(
    n = sum(is.finite(net_kill)),
    auc_base = mean(auc_base, na.rm = TRUE),
    auc_biocide = mean(auc_biocide, na.rm = TRUE),
    net_kill_sd = sd(net_kill, na.rm = TRUE),
    net_kill = mean(net_kill, na.rm = TRUE),
    .groups = "drop"
  ) |>
  write_csv(file.path(out_dir, "fig7_auc_table.csv"))

bind_rows(
  auc_base %>% transmute(treatment = "Combined baseline", group, rep = "combined mean", auc = auc_base),
  auc_h2o2 %>% transmute(treatment = "H2O2", group, rep, auc = auc_biocide),
  auc_gent %>% transmute(treatment = "Gentamicin", group, rep, auc = auc_biocide),
  auc_ddac %>% transmute(treatment = "DDAC", group, rep, auc = auc_biocide),
  auc_soda %>% transmute(treatment = "NaHCO3", group, rep, auc = auc_biocide)
) %>%
  write_csv(file.path(out_dir, "Figure7_replicate_AUC_v18.csv"))
bind_rows(sum_base, sum_h2o2, sum_gent) %>%
  mutate(group = as.character(group)) %>%
  write_csv(file.path(out_dir, "Figure7_kinetics_summary_v18.csv"))
message("AUC 汇总已保存：Result/fig7_auc_table.csv")

#========================
# 5. 读取转录组数据并计算功能箱评分（用于新 Panel D）
#========================

message("\n===== 准备 Panel D 转录组数据 =====")

# 5.1 读取 Fig5 输出的显著 DEG 注释表
sig <- read_csv(file.path(out_dir, "fig5_sig_annotated.csv"), show_col_types = FALSE)

# 5.2 功能箱分配函数（完全复用 Fig6 的 assign_bin 逻辑）
assign_bin <- function(txt) {
  t <- tolower(coalesce(txt, ""))
  case_when(
    # ROS/Redox：精确抗氧化酶和氧化还原蛋白
    str_detect(t, paste0("catalase|peroxiredoxin|thioredoxin|",
                         "glutaredoxin|superoxide|sod|bacilliredoxin|",
                         "lipoyl-dependent peroxidase"))
    ~ "ROS/Redox",
    # Iron/Siderophore：铁获取与铁相关氧化应激调控
    str_detect(t, paste0("iron|sideroph|ferric|ferritin|feo|",
                         "ferrous|ferrichrome|hydroxamate|",
                         "\\bfur\\b|\\bperr\\b|heme transport|",
                         "fhu[abcd]|efeob"))
    ~ "Iron/Siderophore",
    # Motility/Chemotaxis：鞭毛组装与趋化
    str_detect(t, paste0("flagell|\\bfli[a-z]\\b|\\bflg[a-z]\\b|",
                         "\\bflh[a-z]\\b|\\bmot[ab]\\b|",
                         "chemotaxis|methyl-accepting|\\bchea\\b|",
                         "\\bcheb\\b|\\bchew\\b|\\bcher\\b"))
    ~ "Motility/Chemotaxis",
    # Transport：营养物质和离子转运（排除铁转运蛋白）
    str_detect(t, paste0("transport|permease|symporter|",
                         "\\babc\\b.*protein|uptake|efflux|",
                         "substrate-binding protein|atp-binding protein"))
    ~ "Transport",
    # Metabolism：中心碳代谢与氨基酸代谢
    str_detect(t, paste0("dehydrogenase|synthase|transferase|",
                         "reductase|lyase|isomerase|kinase|",
                         "glycolysis|gluconate|amino acid|",
                         "\\btca\\b|citrate|pyruvate|acetyl"))
    ~ "Metabolism",
    TRUE ~ "Other"
  )
}

# 5.3 分配功能箱
sig_bin <- sig %>%
  mutate(def = coalesce(def, ""),
         bin = assign_bin(paste(gene, def)))

# 5.4 定义保留的功能箱（与 Fig6 一致）
focus_bins <- c("Metabolism", "Transport",
                "Motility/Chemotaxis", "ROS/Redox",
                "Iron/Siderophore")

# 5.5 计算 10 h 功能箱评分
mod10 <- sig_bin %>%
  filter(bin %in% focus_bins, timepoint == "10 h") %>%
  group_by(contrast, bin) %>%
  summarise(
    n          = n(),
    signed_sum = sum(log2FoldChange, na.rm = TRUE),
    score_norm = signed_sum / sqrt(pmax(n, 1)),
    mean_lfc   = mean(log2FoldChange, na.rm = TRUE),
    n_sig      = n_distinct(gene),
    .groups    = "drop"
  )

# 5.6 将内部 contrast 映射为展示标签，并与 full group 名对齐
contrast_map <- c(
  "1g_6.5"   = "1g+6.5",
  "1g_13"    = "1g+13",
  "1g_19.5"  = "1g+19.5",
  "SMG_0"    = "SMG",
  "SMG_6.5"  = "SMG+6.5",
  "SMG_13"   = "SMG+13",
  "SMG_19.5" = "SMG+19.5"
)

mod10 <- mod10 %>%
  mutate(
    group = contrast_map[contrast],
    group = factor(group, levels = all_group_levels),
    bin   = factor(bin, levels = focus_bins)
  )

# 5.7 读取消杀 net_kill 表
net_all <- read_csv(file.path(out_dir, "fig7_auc_table.csv"), show_col_types = FALSE) %>%
  mutate(group = factor(group, levels = all_group_levels))

# 5.8 合并：每个 (group, bin, biocide) 一行，携带 score_norm 和 net_kill
panelD_data <- net_all %>%
  select(biocide, group, net_kill) %>%
  left_join(mod10, by = "group", relationship = "many-to-many") %>%
  filter(!is.na(bin))

message("Panel D 数据合并完成，共 ", nrow(panelD_data), " 行")

#========================
# 6. 绘图函数（Panel A–C 与辅助）
#========================

plot_kinetics <- function(sum_df, title_str, show_ylab = FALSE, annotate = NULL) {
  df <- sum_df %>%
    filter(group %in% main_groups) %>%
    mutate(group = factor(group, levels = main_groups))
  
  p <- ggplot(df, aes(x = time_h, y = mean_rel, color = group, fill = group, linetype = group)) +
    geom_ribbon(aes(ymin = mean_rel - sd_rel, ymax = mean_rel + sd_rel), alpha = 0.15, color = NA) +
    geom_line(aes(linewidth = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
    scale_color_manual(name = "Experimental Group", values = main_colors, labels = main_labels) +
    scale_fill_manual(name = "Experimental Group", values = main_colors, labels = main_labels) +
    scale_linetype_manual(name = "Experimental Group", values = main_linetypes, labels = main_labels) +
    scale_linewidth_manual(name = "Experimental Group", values = main_linewidths, labels = main_labels) +
    labs(
      title = title_str,
      x     = "Time [h]",
      y     = if (show_ylab) expression(Relative~OD[600]~change) else NULL
    ) +
    theme_journal() +
    theme(
      plot.title      = element_text(size = 9, face = "plain", hjust = 0.5),
      legend.position = "none",
      axis.title.y    = if (show_ylab) element_text(size = sz_axis, family = font_family) else element_blank(),
      axis.text.y     = element_text(size = sz_tick, family = font_family)
    )
  
  if (!is.null(annotate)) {
    for (a in annotate) {
      p <- p + annotate(
        "text",
        x = a$xt,
        y = a$yt,
        label = a$txt,
        family = font_family,
        size = 2.6,
        fontface = ifelse(is.null(a$face), "plain", a$face),
        hjust = ifelse(is.null(a$hjust), 0.5, a$hjust),
        vjust = ifelse(is.null(a$vjust), 0.5, a$vjust)
      )
    }
  }
  
  p
}

create_shared_legend <- function() {
  legend_df <- expand.grid(group = factor(main_groups, levels = main_groups), x = c(1, 2)) %>%
    mutate(y = 1)
  p_dummy <- ggplot(legend_df,
                    aes(x = x, y = y, color = group, linetype = group, fill = group)) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(aes(ymin = 0.9, ymax = 1.1), alpha = 0.15, color = NA) +
    scale_color_manual(name = "Experimental Group", values = main_colors, labels = main_labels) +
    scale_fill_manual(name = "Experimental Group", values = main_colors, labels = main_labels) +
    scale_linetype_manual(name = "Experimental Group", values = main_linetypes, labels = main_labels) +
    theme_journal() +
    theme(
      legend.position = "right",
      legend.key = element_blank(),
      legend.key.width = unit(0.5, "cm"),
      legend.background = element_blank(),
      legend.margin = margin(l = 10)
    ) +
    guides(
      color = guide_legend(
        override.aes = list(alpha = 1, linewidth = 0.8, fill = alpha(main_colors, 0.15), linetype = main_linetypes)
      ),
      fill = "none",
      linetype = "none"
    )
  cowplot::get_legend(p_dummy)
}

plot_auc_replicates <- function(net_df, x_col, title_str, xlab_str,
                                show_ylab = FALSE) {
  df_rep <- net_df |>
    rename(x_val = all_of(x_col)) |>
    mutate(
      group = factor(group, levels = all_group_levels)
    )

  df_sum <- df_rep %>%
    group_by(group) %>%
    summarise(
      n = sum(is.finite(x_val)),
      mean = mean(x_val, na.rm = TRUE),
      sd = ifelse(n > 1, sd(x_val, na.rm = TRUE), 0),
      .groups = "drop"
    )

  p <- ggplot(df_sum, aes(x = mean, y = group, color = group)) +
    geom_errorbar(aes(xmin = mean - sd, xmax = mean + sd),
                  orientation = "y", width = 0.18,
                  linewidth = 0.45, show.legend = FALSE) +
    geom_point(data = df_rep, aes(x = x_val, y = group, fill = group),
               inherit.aes = FALSE, shape = 21, size = 1.8, stroke = 0.35,
               position = position_jitter(height = 0.08, width = 0),
               show.legend = FALSE) +
    geom_point(shape = 18, size = 2.5, show.legend = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    geom_hline(yintercept = 4.5, linetype = "dashed", linewidth = 0.4, color = "grey55") +
    scale_color_manual(values = all_group_colors) +
    scale_fill_manual(values = all_group_colors) +
    scale_x_continuous(name = xlab_str, expand = expansion(mult = c(0.08, 0.08))) +
    labs(
      title = title_str,
      y = if (show_ylab) "Experimental Group" else NULL
    ) +
    theme_journal() +
    theme(
      plot.title      = element_text(size = 9, face = "plain", family = font_family, hjust = 0.5),
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.2, color = "grey92"),
      axis.title.y    = if (show_ylab) element_text(size = sz_axis, family = font_family) else element_blank(),
      axis.text.y     = if (show_ylab) element_text(size = sz_tick, family = font_family) else element_blank()
    )
  p
}

#========================
# 7. 绘制主图 A–C
#========================

message("\n===== 绘制 Panel A-C =====")

pA_kin <- plot_kinetics(sum_base, "No biocide (growth baseline)", show_ylab = TRUE)
pB_kin <- plot_kinetics(sum_h2o2, expression(paste("3% H"[2], "O"[2], " (Oxidative agent)")), show_ylab = FALSE)
pC_kin <- plot_kinetics(sum_gent, "4 μg/mL Gentamicin\n(Protein synthesis inhibitor)", show_ylab = FALSE)

shared_legend <- create_shared_legend()

pA_bar <- plot_auc_replicates(
  auc_base |> rename(auc = auc_base), "auc", "Baseline growth [AUC]",
  expression("Relative OD"[600]~"AUC"), show_ylab = TRUE
)
pB_bar <- plot_auc_replicates(
  net_h2o2, "net_kill",
  expression(paste("H"[2], "O"[2], " net killing [AUC]")),
  expression("Net killing AUC")
)
pC_bar <- plot_auc_replicates(
  net_gent, "net_kill", "Gentamicin net killing [AUC]", "Net killing AUC"
)

pA_bar_tagged <- pA_bar + labs(tag = "a")
pB_bar_tagged <- pB_bar
pC_bar_tagged <- pC_bar

pA_kin_tagged <- pA_kin + labs(tag = "b")
pB_kin_tagged <- pB_kin
pC_kin_tagged <- pC_kin

#========================
# 8. 新 Panel D：双热图（功能箱 + 精选通路）并排，共用图例
#========================

message("\n===== 绘制新 Panel D (双热图并排) =====")

# 8.1 通用参数 --------------------------------------------------------------
selected_biocides <- c("H2O2", "Gentamicin")

# ---- 功能箱数据准备 ----
pd_fun <- panelD_data %>%
  filter(biocide %in% selected_biocides) %>%
  mutate(
    biocide = factor(biocide, levels = selected_biocides),
    bin     = factor(bin, levels = focus_bins),
    treat_type = case_when(
      str_detect(group, "^1g\\+")          ~ "Radiation",
      str_detect(group, "^SMG\\+")         ~ "SMG + Radiation",
      group == "SMG"                       ~ "SMG alone",
      TRUE                                 ~ NA_character_
    )
  ) %>%
  filter(!is.na(treat_type), treat_type != "SMG alone")   # 删除 SMG alone

# 仅保留三个功能箱（删除 ROS/Redox 和 Motility/Chemotaxis）
keep_bins <- c("Metabolism", "Transport", "Iron/Siderophore")
pd_fun <- pd_fun %>% filter(bin %in% keep_bins) %>%
  mutate(bin = factor(bin, levels = keep_bins))

# 安全性相关系数函数
safe_cor <- function(x, y) {
  if (length(unique(x)) < 2 || length(unique(y)) < 2) return(list(r = NA_real_, p = NA_real_))
  res <- tryCatch(cor.test(x, y, method = "pearson"), error = function(e) NULL)
  if (is.null(res)) return(list(r = NA_real_, p = NA_real_))
  list(r = as.numeric(res$estimate), p = res$p.value)
}

# 计算功能箱相关性
cor_fun <- pd_fun %>%
  filter(treat_type %in% c("Radiation", "SMG + Radiation")) %>%
  group_by(biocide, bin, treat_type) %>%
  summarise(
    n = n(),
    r = safe_cor(score_norm, net_kill)$r,
    p = safe_cor(score_norm, net_kill)$p,
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(!is.na(r),
                   paste0(sprintf("%.2f", r), "\n(n=", n, ")"),
                   "–"),
    sublabel = ifelse(!is.na(r), paste0("(N=", n, ")"), ""),
    text_col = ifelse(abs(r) > 0.6 & !is.na(r), "white", "#111111")
  )

# 构建完整因子网格（三种处理类型，但 SMG alone 已移除，只保留两种）
fun_grid <- expand.grid(
  biocide    = factor(selected_biocides, levels = selected_biocides),
  bin        = factor(keep_bins, levels = keep_bins),
  treat_type = c("Radiation", "SMG + Radiation"),
  stringsAsFactors = FALSE
) %>% as_tibble()

cor_fun_mat <- fun_grid %>%
  left_join(cor_fun, by = c("biocide", "bin", "treat_type")) %>%
  mutate(
    label      = coalesce(label, "–"),
    sublabel   = coalesce(sublabel, ""),
    text_col   = coalesce(text_col, "#111111"),
    x_label    = paste(biocide, treat_type, sep = "\n")
  )

# ---- 通路数据准备 ----
path_sc <- read_csv(file.path(out_dir, "fig56_kegg_pathway_scores.csv"), show_col_types = FALSE)

contrast_map <- c(
  "1g_6.5"   = "1g+6.5",
  "1g_13"    = "1g+13",
  "1g_19.5"  = "1g+19.5",
  "SMG_0"    = "SMG",
  "SMG_6.5"  = "SMG+6.5",
  "SMG_13"   = "SMG+13",
  "SMG_19.5" = "SMG+19.5"
)

path10 <- path_sc %>%
  filter(timepoint == "10 h") %>%
  mutate(group = factor(contrast_map[contrast], levels = all_group_levels)) %>%
  filter(!is.na(group))

path_pd <- path10 %>%
  left_join(
    net_all %>% filter(biocide %in% selected_biocides) %>% select(biocide, group, net_kill),
    by = "group",
    relationship = "many-to-many"
  ) %>%
  filter(!is.na(net_kill), !is.na(biocide)) %>%
  mutate(
    biocide = factor(biocide, levels = selected_biocides),
    treat_type = case_when(
      str_detect(group, "^1g\\+")          ~ "Radiation",
      str_detect(group, "^SMG\\+")         ~ "SMG + Radiation",
      group == "SMG"                       ~ "SMG alone",
      TRUE                                 ~ NA_character_
    )
  ) %>%
  filter(!is.na(treat_type), treat_type != "SMG alone")

# 仅保留三个通路
keep_paths <- c("Two-component system", "Carbon metabolism", "Phosphotransferase system (PTS)")
path_pd <- path_pd %>%
  filter(pathway %in% keep_paths) %>%
  mutate(pathway = factor(pathway, levels = keep_paths))

# 计算通路相关性
cor_path <- path_pd %>%
  filter(treat_type %in% c("Radiation", "SMG + Radiation")) %>%
  group_by(biocide, pathway, treat_type) %>%
  summarise(
    n = n(),
    r = safe_cor(score_norm, net_kill)$r,
    p = safe_cor(score_norm, net_kill)$p,
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(!is.na(r),
                   paste0(sprintf("%.2f", r), "\n(n=", n, ")"),
                   "–"),
    sublabel = ifelse(!is.na(r), paste0("(N=", n, ")"), ""),
    text_col = ifelse(abs(r) > 0.6 & !is.na(r), "white", "#111111")
  )

bind_rows(
  cor_fun %>% transmute(feature_type = "functional module", feature = as.character(bin),
                        biocide = as.character(biocide), treat_type, n, r, p),
  cor_path %>% transmute(feature_type = "KEGG pathway", feature = as.character(pathway),
                         biocide = as.character(biocide), treat_type, n, r, p)
) %>%
  write_csv(file.path(out_dir, "Figure7_exploratory_correlations_v18.csv"))

path_grid <- expand.grid(
  biocide    = factor(selected_biocides, levels = selected_biocides),
  pathway    = factor(keep_paths, levels = keep_paths),
  treat_type = c("Radiation", "SMG + Radiation"),
  stringsAsFactors = FALSE
) %>% as_tibble()

cor_path_mat <- path_grid %>%
  left_join(cor_path, by = c("biocide", "pathway", "treat_type")) %>%
  mutate(
    label      = coalesce(label, "–"),
    sublabel   = coalesce(sublabel, ""),
    text_col   = coalesce(text_col, "#111111"),
    x_label    = paste(biocide, treat_type, sep = "\n")
  )

# ---- 绘制两个热图 ----
make_heatmap <- function(data, x, y, fill, label_col, subtitle) {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(
      aes(label = .data[[label_col]], color = text_col),
      size = 2.2,
      family = font_family,
      show.legend = FALSE,
      lineheight = 0.85
    ) +
    scale_fill_gradient2(
      low      = "#0072B5",
      mid      = "white",
      high     = "#BC3C29",
      midpoint = 0,
      limits   = c(-1, 1),
      na.value = "grey85",
      name     = expression("Pearson " * italic(r))
    ) +
    scale_color_identity() +
    labs(x = NULL, y = NULL, subtitle = subtitle) +
    theme_journal() +
    theme(
      plot.subtitle = element_text(size = 8, face = "plain", hjust = 0.5),
      axis.text.x  = element_text(angle = 35, hjust = 1, size = sz_tick, margin = margin(t = 2)),
      panel.grid   = element_blank(),
      legend.position = "none"   # 稍后提取共享图例
    )
}

p_left <- make_heatmap(
  cor_fun_mat,
  x = "x_label", y = "bin", fill = "r",
  label_col = "label", subtitle = "Exploratory functional-module correlations"
) + labs(tag = "c")

p_right <- make_heatmap(
  cor_path_mat,
  x = "x_label", y = "pathway", fill = "r",
  label_col = "label", subtitle = "Exploratory KEGG-pathway correlations"
)

# 合并左右两图，提取共享图例（从任一热图中抓取）
shared_legend_pD <- cowplot::get_legend(
  p_left + theme(legend.position = "right")
)

# 最终 Row3：两图并排 + 右侧共享图例
row3 <- wrap_plots(
  p_left, p_right, wrap_elements(full = shared_legend_pD),
  nrow = 1, widths = c(1, 1, 0.2)
)

#========================
# 9. 拼图并保存（总高微调）
#========================

message("\n===== 拼合最终 Figure 7 =====")

row1 <- wrap_plots(pA_bar_tagged, pB_bar_tagged, pC_bar_tagged, plot_spacer(),
                   nrow = 1, widths = c(1, 1, 1, 0.2))
row2 <- wrap_plots(pA_kin_tagged, pB_kin_tagged, pC_kin_tagged,
                   wrap_elements(full = shared_legend),  # 动力学曲线共享图例
                   nrow = 1, widths = c(1, 1, 1, 0.5))

fig7 <- (row1 / row2 / row3) +
  plot_layout(heights = c(1, 0.95, 1.25)) &   # 双热图行高适当缩小
  theme(plot.tag = element_text(size = sz_tag, face = "bold", family = font_family))

ggsave(
  file.path(out_dir, "Figure7_submission_v18.pdf"),
  fig7,
  width  = fig_w_double,
  height = 24,   # 总高微调以容纳双热图
  units  = "cm",
  device = cairo_pdf
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(out_dir, "Figure7_submission_v18.svg"), fig7,
         width = fig_w_double, height = 24, units = "cm", device = svglite::svglite)
}
ggsave(file.path(out_dir, "Figure7_submission_v18.tiff"), fig7,
       width = fig_w_double, height = 24, units = "cm", dpi = fig_res,
       device = ragg::agg_tiff, compression = "lzw")
ggsave(file.path(out_dir, "Figure7_submission_v18.png"), fig7,
       width = fig_w_double, height = 24, units = "cm", dpi = fig_res,
       device = ragg::agg_png)

message("Figure 7 replicate-level AUC version generated: Figure7_submission_v18.pdf")
message("\n===== Figure 7 全部完成 =====")

