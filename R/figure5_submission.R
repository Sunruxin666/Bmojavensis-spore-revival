####################################################
## Figure 5：转录组DEG负担分析 + KEGG通路重塑可视化
##
## 目的：
## 1）量化各对比组显著差异基因数量（DEG负担）；
## 2）构建转录组-表型桥接数据（与生长动力学关联）；
## 3）基于KEGG通路评分展示耦合应激下的通路重编程；
## 4）读取第一脚本（RDA）输出，绘制PERMANOVA因子重要性图。
##
## 核心思路：
## 1）扫描14组DEG结果文件，解析实验元数据；
## 2）读取6个all_genes.kegg3 + KEGG官方表，构建基因注释；
## 3）筛选显著DEG（padj<0.05，|log2FC|≥1），关联功能注释；
## 4）在KEGG通路层面计算方向性评分（score_norm）；
## 5）读取生长动力学数据，计算各对比组相对基线的表型偏移；
## 6）绘制主图（PERMANOVA柱状图 + KEGG通路热图）；
## 7）绘制补充图（SMG重写指数 + DEG负担图）。
##
## 前置要求：
## 需先运行第一脚本（RDA + PERMANOVA），生成：
##   Result/Fig5A_permanova_3h.csv
##   Result/Fig5A_permanova_10h.csv
## Result/目录下需有 growth_kinetics_group_summary.csv（Fig4输出）；
## 工作目录下需有已下载的KEGG对照表：
##   ko_definition.tsv / ko_to_pathway.tsv / pathway_names.tsv
##
## Repository data structure:
##   data/
##   ├── annotations/
##   │   ├── ko_definition.tsv
##   │   ├── ko_to_pathway.tsv
##   │   └── pathway_names.tsv
##   ├── phenotype/
##   ├── transcriptomics/3h/
##   │   ├── Ctrl_all/
##   │   │   ├── all_genes.kegg3
##   │   │   └── results/DEG_full_*_vs_Ctrl.csv
##   │   ├── Ctrl_vs_SMG/
##   │   │   ├── all_genes.kegg3
##   │   │   └── results/
##   │   └── Ctrl_vs_SMG+SR/
##   │       ├── all_genes.kegg3
##   │       └── results/
##   └── transcriptomics/10h/       (same structure)
##
## 输出内容：
## Figure5.pdf（主图：PERMANOVA柱状图 + KEGG通路热图）
## Figure5_S1_rewriting.pdf（补充图1：SMG重写指数热图）
## Figure5_S2_DEGburden.pdf（补充图2：DEG负担三联图）
## fig5_sig_annotated.csv（注释后的显著DEG，供Fig6.R使用）
## fig5_transcriptome_phenotype_bridge.csv（表型桥接数据，供Fig6.R使用）
## fig56_kegg_pathway_scores.csv（KEGG通路评分中间数据）
## fig5_deg_burden.csv（各组DEG统计汇总）
####################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(showtext)
})

set.seed(123)

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
annotation_dir <- file.path(project_root, "data", "annotations")
original_result <- file.path(project_root, "data", "phenotype")

#========================
# 1. 参数设置
#========================

## 数据根目录：工作目录往上一级即 dataprocess/
data_root <- file.path(project_root, "data", "transcriptomics")

## 输出目录
out_dir <- file.path(project_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ── 字体设置（期刊要求：Arial）──────────────────────────────
## arial.ttf 需放置在工作目录下
## Windows系统通常可在 C:/Windows/Fonts/arial.ttf 找到并复制过来
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
showtext_auto()   # 自动对所有图形设备启用showtext渲染
showtext_opts(dpi = 600)  # 必须与高分辨率输出一致，避免字体按96 dpi渲染后过小

## ── 期刊尺寸参数 ────────────────────────────────────────────
## 双栏图宽度固定19cm；高度根据内容自定义（见各ggsave）
fig_width_cm  <- 18.3    # 双栏图宽度（cm）
fig_res       <- 600     # 分辨率（dpi）

## ── 颜色方案（期刊配色）─────────────────────────────────────
col_blue <- "#0072B5"    # 蓝色：下调/负值/10h
col_red  <- "#BC3C29"    # 红色：上调/正值/显著/3h
col_grey <- "#BDBDBD"    # 灰色：不显著

## ── 期刊字号参数 ─────────────────────────────────────────────
sz_base   <- 9   # 基础字号（pt）；保证缩放至双栏宽度后仍清晰
sz_axis   <- 8.5 # 坐标轴标题
sz_tick   <- 8   # 坐标轴刻度与通路名称
sz_legend <- 8   # 图例文字与p值标签
sz_tag    <- 9.5 # 分图标号（a、b）

## ── 统一主题（所有子图共用）────────────────────────────────
## 基于期刊要求：theme_bw去灰底，Arial字体，统一字号，细线框
theme_journal <- function() {
  theme_bw(base_size = sz_base, base_family = "Arial") +
    theme(
      ## 坐标轴
      axis.title       = element_text(size = sz_axis,   family = "Arial"),
      axis.text        = element_text(size = sz_tick,   family = "Arial"),
      axis.line        = element_line(linewidth = 0.4),
      ## 图框边线
      panel.border     = element_rect(linewidth = 0.4, fill = NA),
      ## 网格线（保留水平主网格，去除次网格和垂直网格）
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      ## 图例
      legend.text      = element_text(size = sz_legend, family = "Arial"),
      legend.title     = element_text(size = sz_legend, family = "Arial"),
      legend.key.size  = unit(0.35, "cm"),
      ## 分面标签
      strip.text       = element_text(size = sz_tick,   family = "Arial"),
      strip.background = element_rect(fill = "grey95",  color = "grey70",
                                      linewidth = 0.3),
      ## 副标题
      plot.subtitle    = element_text(size = 7.5, family = "Arial",
                                      color = "grey40")
    )
}

## 路径验证：启动前确认所有关键文件可访问
stopifnot(dir.exists(file.path(data_root, "3h")))
stopifnot(dir.exists(file.path(data_root, "10h")))
stopifnot(file.exists(font_regular), file.exists(font_bold))
stopifnot(file.exists(file.path(annotation_dir, "ko_definition.tsv")))
stopifnot(file.exists(file.path(annotation_dir, "ko_to_pathway.tsv")))
stopifnot(file.exists(file.path(annotation_dir, "pathway_names.tsv")))
stopifnot(file.exists(file.path(original_result, "Fig5A_permanova_3h.csv")))
stopifnot(file.exists(file.path(original_result, "Fig5A_permanova_10h.csv")))
stopifnot(file.exists(file.path(original_result, "growth_kinetics_group_summary.csv")))
message("路径验证通过！")

## 差异基因筛选阈值
fdr_cutoff <- 0.05   # FDR阈值
fc_cutoff  <- 1      # |log2FC|阈值

## ── KEGG主图展示通路（按生物学主题分组选定，共12条）──────────
## 筛选原则：
##   1. 排除 Metabolic pathways / Microbial metabolism 等超级聚合通路；
##   2. score_norm 绝对值 ≥ 3.2（高信号强度）；
##   3. 覆盖与 B. mojavensis 芽孢复苏相关的5个核心主题：
##      信号感知 / 能量代谢 / 碳源摄取 / 运动能力 / 铁获取。
## 排列顺序：由上到下反映生物学逻辑（感知→能量→运动→铁）
path_keep <- c(
  ## 信号感知与转导
  "Two-component system",
  "Inositol phosphate metabolism",
  ## 碳代谢与能量供给
  "Carbon metabolism",
  "Glycolysis / Gluconeogenesis",
  "Citrate cycle (TCA cycle)",
  "Glyoxylate and dicarboxylate metabolism",
  ## 糖摄取（芽孢萌发关键）
  "Phosphotransferase system (PTS)",
  ## 氧化还原（辐射应激响应）
  "Pentose phosphate pathway",
  ## 氨基酸代谢
  "Alanine, aspartate and glutamate metabolism",
  ## 运动能力（SMG敏感）
  "Bacterial chemotaxis",
  "Flagellar assembly",
  ## 铁获取（SMG环境下铁代谢重排）
  "Biosynthesis of siderophore group nonribosomal peptides"
)

## 通路显示名称映射：仅对超长名称做简化，其余保持原名
## 简化原则：保留生物学核心词，去除冗余修饰语
pathway_display_map <- c(
  "Biosynthesis of siderophore group nonribosomal peptides" =
    "Siderophore biosynthesis"
)

## 生成用于因子水平的显示顺序（path_keep中的名称替换后）
path_display <- ifelse(
  path_keep %in% names(pathway_display_map),
  pathway_display_map[path_keep],
  path_keep
)

## 补充图自动选取的最大通路数（排除超级通路后按score排名）
top_n_pathway <- 8

## 全局对比组标签映射：内部编码 → 图中展示标签
contrast_map <- c(
  "1g_6.5"   = "1 g + 6.5",
  "1g_13"    = "1 g + 13",
  "1g_19.5"  = "1 g + 19.5",
  "SMG_0"    = "SMG + 0",
  "SMG_6.5"  = "SMG + 6.5",
  "SMG_13"   = "SMG + 13",
  "SMG_19.5" = "SMG + 19.5"
)

contrast_levels <- unname(contrast_map)   # 图例从左到右的因子水平

## PERMANOVA效应项标签映射
term_map <- c(
  "condition"        = "Gravity condition",
  "dose_f"           = "Radiation dose",
  "condition:dose_f" = "Gravity condition × Radiation dose"
)

#========================
# 2. 扫描DEG文件并解析元数据
#========================

message("\n===== 扫描DEG文件 =====")

## 递归扫描两个时间点下所有results/子目录中的DEG文件
## 路径规则：dataprocess/A/B/results/DEG_full_*_vs_Ctrl.csv
files <- unique(c(
  list.files(file.path(data_root, "3h"),
             pattern = "DEG_full_.*_vs_Ctrl\\.csv$",
             recursive = TRUE, full.names = TRUE),
  list.files(file.path(data_root, "10h"),
             pattern = "DEG_full_.*_vs_Ctrl\\.csv$",
             recursive = TRUE, full.names = TRUE)
))

if (length(files) == 0) stop("未找到任何DEG文件，请检查目录结构")
message("共扫描到 ", length(files), " 个DEG文件")

## ── 辅助函数：从文件路径解析实验元数据 ────────────────────
## 输入：单个DEG文件路径
## 输出：tibble（timepoint / condition / dose / contrast）
parse_meta <- function(path) {
  fn        <- basename(path)
  timepoint <- ifelse(str_detect(path, paste0("[/\\\\]3h[/\\\\]")), "3 h", "10 h")
  
  if (str_detect(fn, "DEG_full_SR[0-9.]+Gy_vs_Ctrl")) {
    ## 纯辐射组（正常重力 + 辐射）：提取剂量数值
    dose      <- as.numeric(str_match(fn, "SR([0-9.]+)Gy_vs_Ctrl")[, 2])
    condition <- "1g"
    contrast  <- paste0("1g_", dose)
    
  } else if (str_detect(fn, "DEG_full_SMG_vs_Ctrl")) {
    ## 纯SMG组（无辐射）：剂量设为0
    dose      <- 0
    condition <- "SMG"
    contrast  <- "SMG_0"
    
  } else if (str_detect(fn, "DEG_full_SMG\\+SR[0-9.]+Gy_vs_Ctrl")) {
    ## SMG + 辐射联合组：提取剂量数值
    dose      <- as.numeric(str_match(fn, "SMG\\+SR([0-9.]+)Gy_vs_Ctrl")[, 2])
    condition <- "SMG"
    contrast  <- paste0("SMG_", dose)
    
  } else {
    return(NULL)   # 无法识别的文件名，跳过
  }
  
  tibble(path = path, timepoint = timepoint,
         condition = condition, dose = dose, contrast = contrast)
}

## 合并元数据，同一时间点 × 对比组仅保留一条
meta <- bind_rows(lapply(files, parse_meta)) %>%
  distinct(timepoint, contrast, .keep_all = TRUE)

message("解析到有效对比组：", nrow(meta), " 个")

#========================
# 3. 读取KEGG官方对照表
#========================

message("\n===== 读取KEGG官方对照表 =====")

## 1️⃣ KO功能描述表（KO → def）
## 已手动从 https://rest.kegg.jp/list/ko 下载到工作目录
k2def <- read_tsv(
  file.path(annotation_dir, "ko_definition.tsv"),
  col_names = c("KO", "def"),
  show_col_types = FALSE
) %>%
  mutate(KO = str_replace(KO, "^ko:", "")) %>%   # 去除"ko:"前缀
  distinct(KO, .keep_all = TRUE)

message("  KO功能描述条目数：", nrow(k2def))

## 2️⃣ KO → 通路映射表
## 已从 https://rest.kegg.jp/link/pathway/ko 下载
k2p <- read_tsv(
  file.path(annotation_dir, "ko_to_pathway.tsv"),
  col_names = c("ko_raw", "path_raw"),
  show_col_types = FALSE
) %>%
  transmute(
    KO   = str_replace(ko_raw,   "^ko:",   ""),
    path = str_replace(path_raw, "^path:", "")
  ) %>%
  filter(str_detect(path, "^ko\\d+")) %>%   # 仅保留KEGG参考通路
  distinct(KO, path)

message("  KO→通路映射条目数：", nrow(k2p))

## 3️⃣ 通路名称表
## 已从 https://rest.kegg.jp/list/pathway 下载
pnames <- read_tsv(
  file.path(annotation_dir, "pathway_names.tsv"),
  col_names = c("mapid", "pathway"),
  show_col_types = FALSE
) %>%
  transmute(
    path    = str_replace(mapid, "^map", "ko"),   # map00010 → ko00010
    pathway = str_replace(pathway, " - .*$", "")  # 去除物种后缀（如" - Escherichia coli"）
  )

message("  通路名称条目数：", nrow(pnames))

#========================
# 4. 构建gene→KO→def映射
#========================

message("\n===== 构建基因注释映射表 =====")

## 读取6个all_genes.kegg3（2时间点 × 3子目录），合并构建gene→KO映射
## 路径规则：dataprocess/A/B/all_genes.kegg3
## CountMatrix覆盖全部检测基因，故6文件合并后天然覆盖14组所有DEG
kegg3_files <- c(
  file.path(data_root, "3h/Ctrl_all/all_genes.kegg3"),
  file.path(data_root, "3h/Ctrl_vs_SMG/all_genes.kegg3"),
  file.path(data_root, "3h/Ctrl_vs_SMG+SR/all_genes.kegg3"),
  file.path(data_root, "10h/Ctrl_all/all_genes.kegg3"),
  file.path(data_root, "10h/Ctrl_vs_SMG/all_genes.kegg3"),
  file.path(data_root, "10h/Ctrl_vs_SMG+SR/all_genes.kegg3")
)

## 缺失文件给出警告（不终止，避免因单文件缺失中断整个流程）
for (f in kegg3_files) {
  if (!file.exists(f)) message("  警告：文件不存在 → ", f)
}

## 合并6个kegg3文件，过滤非规范KO，每个gene保留唯一KO条目
g2k <- bind_rows(lapply(kegg3_files, function(f) {
  if (!file.exists(f)) return(tibble(gene = character(), KO = character()))
  read.table(f, col.names = c("gene", "KO", "score"),
             stringsAsFactors = FALSE) %>%
    select(gene, KO)
})) %>%
  mutate(
    gene = str_replace_all(as.character(gene), '"', ''),
    KO   = str_replace_all(as.character(KO),   '"', '')
  ) %>%
  filter(str_detect(KO, "^K\\d+")) %>%   # 仅保留K+数字格式的KO编号
  distinct(gene, KO)

message("gene→KO映射条目数：", nrow(g2k))

## 构建ann_ko（gene → KO + def）
## g2k关联KEGG官方k2def，每个基因只保留第一个有效KO和功能描述
ann_ko <- g2k %>%
  left_join(k2def, by = "KO") %>%
  group_by(gene) %>%
  summarise(
    KO  = first(KO[!is.na(KO)]),
    def = first(def[!is.na(def)]),
    .groups = "drop"
  )

message("ann_ko构建完成，覆盖基因数：", nrow(ann_ko),
        "；有def注释：", sum(!is.na(ann_ko$def)))

#========================
# 5. 统计各对比组DEG负担
#========================

message("\n===== 统计DEG负担 =====")

## 逐组读取DEG结果，统计显著基因总数及上下调方向
## 筛选标准：padj < fdr_cutoff 且 |log2FC| ≥ fc_cutoff
deg_burden <- bind_rows(lapply(seq_len(nrow(meta)), function(i) {
  m   <- meta[i, ]
  dat <- read_csv(m$path, show_col_types = FALSE) %>%
    mutate(padj           = as.numeric(padj),
           log2FoldChange = as.numeric(log2FoldChange))
  
  sig_i <- dat %>%
    filter(!is.na(padj), padj < fdr_cutoff, abs(log2FoldChange) >= fc_cutoff)
  
  tibble(
    timepoint        = m$timepoint,
    condition        = m$condition,
    dose             = m$dose,
    contrast         = m$contrast,
    n_tested         = nrow(dat),
    n_sig            = nrow(sig_i),
    n_up             = sum(sig_i$log2FoldChange > 0, na.rm = TRUE),
    n_down           = sum(sig_i$log2FoldChange < 0, na.rm = TRUE),
    mean_abs_lfc_sig = ifelse(nrow(sig_i) > 0,
                              mean(abs(sig_i$log2FoldChange), na.rm = TRUE),
                              NA_real_)
  )
})) %>%
  mutate(
    contrast       = factor(contrast, levels = names(contrast_map)),
    contrast_label = factor(contrast_map[as.character(contrast)],
                            levels = contrast_levels)
  ) %>%
  arrange(timepoint, contrast)

message("DEG负担统计完成，共 ", nrow(deg_burden), " 个对比组")

#========================
# 6. 读取生长动力学数据，构建表型桥接
#========================

message("\n===== 构建转录组-表型桥接数据 =====")

## 读取Fig4脚本输出的生长动力学分组汇总数据
growth <- read_csv(file.path(original_result, "growth_kinetics_group_summary.csv"),
                   show_col_types = FALSE) %>%
  mutate(dose = as.numeric(Dose), condition = Condition)

## 提取基线（1g 0 Gy）用于计算各组相对偏移量
baseline <- growth %>% filter(condition == "1g", dose == 0) %>% slice(1)

if (nrow(baseline) != 1) {
  stop("基线（1g, 0 Gy）在growth_kinetics_group_summary.csv中不唯一，请核查数据")
}

## 合并DEG统计与生长动力学，计算各指标相对基线的偏移量
bridge <- deg_burden %>%
  mutate(condition = as.character(condition), dose = as.numeric(dose)) %>%
  left_join(
    growth %>% transmute(condition, dose,
                         mu_max = mu_max_mean, lag = lag_mean, AUC = AUC_mean),
    by = c("condition", "dose")
  ) %>%
  mutate(
    delta_mu_max = mu_max - baseline$mu_max_mean,   # 最大增长率偏移
    delta_lag    = lag    - baseline$lag_mean,       # 延滞期偏移
    delta_AUC    = AUC    - baseline$AUC_mean        # 净增长AUC偏移
  ) %>%
  select(timepoint, condition, dose, contrast, contrast_label,
         n_sig, delta_mu_max, delta_lag, delta_AUC)

message("表型桥接数据构建完成，共 ", nrow(bridge), " 条记录")

#========================
# 7. 合并显著DEG并关联注释
#========================

message("\n===== 合并显著差异基因 =====")

## 逐组筛选显著DEG并关联KO编号和功能描述
## 该数据集保存为CSV供Fig6.R直接读取，无需重复运行注释流程
sig <- bind_rows(lapply(seq_len(nrow(meta)), function(i) {
  m <- meta[i, ]
  message("  读取：", m$timepoint, " | ", m$contrast)
  
  # 1️⃣ 读取DEG结果，筛选显著基因
  read_csv(m$path, show_col_types = FALSE) %>%
    transmute(
      gene           = str_replace_all(as.character(gene), '"', ''),
      log2FoldChange = as.numeric(log2FoldChange),
      padj           = as.numeric(padj)
    ) %>%
    filter(!is.na(padj), padj < fdr_cutoff, abs(log2FoldChange) >= fc_cutoff) %>%
    
    # 2️⃣ 关联KO编号和功能描述（来自g2k和ann_ko）
    left_join(g2k,                          by = "gene") %>%
    left_join(ann_ko %>% select(gene, def), by = "gene") %>%
    
    # 3️⃣ 添加分组信息；def缺失时填为空字符串
    mutate(
      def       = coalesce(def, ""),
      timepoint = m$timepoint,
      condition = m$condition,
      dose      = m$dose,
      contrast  = m$contrast
    )
}))

message("显著DEG记录总数：", nrow(sig))

#========================
# 8. 计算KEGG通路层面评分
#========================

message("\n===== 计算KEGG通路评分 =====")

## 逐通路计算方向性评分（score_norm = signed_sum / sqrt(k)）
## 分母sqrt(k)归一化基因数量，分子signed_sum反映整体方向与幅度
path_sc <- sig %>%
  filter(!is.na(KO), KO != "") %>%
  left_join(k2p,    by = "KO",   relationship = "many-to-many") %>%
  filter(!is.na(path)) %>%
  left_join(pnames, by = "path") %>%
  filter(!is.na(pathway)) %>%
  group_by(timepoint, condition, dose, contrast, path, pathway) %>%
  summarise(
    k          = n_distinct(gene),
    n_sig      = k,
    up_n       = sum(log2FoldChange > 0),
    down_n     = sum(log2FoldChange < 0),
    signed_sum = sum(log2FoldChange),
    mean_lfc   = mean(log2FoldChange),
    score_norm = signed_sum / sqrt(pmax(k, 1)),
    .groups    = "drop"
  )

message("通路评分计算完成，共 ", nrow(path_sc), " 条通路-对比组记录")

## 排除超级聚合通路后，按max|score_norm|自动选取补充图代表性通路
super_pathways <- c("Metabolic pathways",
                    "Microbial metabolism in diverse environments",
                    "Biosynthesis of secondary metabolites")

path_top_auto <- path_sc %>%
  filter(!pathway %in% super_pathways) %>%
  group_by(pathway) %>%
  summarise(max_abs = max(abs(score_norm), na.rm = TRUE),
            tot_n   = sum(n_sig), .groups = "drop") %>%
  arrange(desc(max_abs), desc(tot_n)) %>%
  slice_head(n = top_n_pathway) %>%
  pull(pathway)

#========================
# 9. 绘制Figure 5主图
#========================

message("\n===== 绘制Figure 5主图 =====")

## ── Panel a：PERMANOVA因子重要性柱状图 ──────────────────────

# 1️⃣ 读取第一脚本输出的PERMANOVA结果（Result/目录下）
perm <- bind_rows(
  read_csv(file.path(original_result, "Fig5A_permanova_3h.csv"),
           show_col_types = FALSE) %>% mutate(timepoint = "3 h"),
  read_csv(file.path(original_result, "Fig5A_permanova_10h.csv"),
           show_col_types = FALSE) %>% mutate(timepoint = "10 h")
) %>%
  filter(term %in% c("condition", "dose_f", "condition:dose_f")) %>%
  transmute(
    timepoint,
    term      = term_map[term],   # 命名向量直接索引，转为可读标签
    R2,
    p         = `Pr(>F)`
  ) %>%
  mutate(
    term      = factor(term, levels = unname(term_map)),
    timepoint = factor(timepoint, levels = c("3 h", "10 h"))
  )

# 2️⃣ 绘制横向柱状图
## 显著（p<0.05）= col_red；不显著 = col_grey；p值标注在柱右侧
pa <- ggplot(perm, aes(x = R2, y = term, fill = p < 0.05)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))),
    hjust = -0.1, size = sz_legend / .pt,   # .pt将pt转换为ggplot内部单位
    family = "Arial"
  ) +
  facet_wrap(~ timepoint, nrow = 1) +
  scale_fill_manual(
    values = c(`TRUE` = col_red, `FALSE` = col_grey),
    labels = c(`TRUE` = "Significant", `FALSE` = "Not significant")
  ) +
  expand_limits(x = max(perm$R2, na.rm = TRUE) * 1.3) +
  labs(x = expression("PERMANOVA" ~ R^2),
       y = NULL,
       fill = "p < 0.05") +
  theme_journal() +
  theme(
    legend.position  = "top",
    ## 去掉水平方向网格（横向柱状图不需要）
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2, color = "grey90")
  )

## ── Panel b：KEGG通路活性热图 ─────────────────────────────

# 1️⃣ 筛选path_keep通路，补全缺失的pathway-contrast组合（填0）
## 关键改动：
##   a) contrast_map已去"Gy"，标签更简洁
##   b) 应用pathway_display_map缩短超长通路名
##   c) complete()补全空缺组合，缺失格显示为白色（score=0）
pdat <- path_sc %>%
  filter(pathway %in% path_keep) %>%
  mutate(
    contrast_label = contrast_map[as.character(contrast)],
    ## 应用显示名称映射（仅对名称过长的通路生效）
    pathway = ifelse(pathway %in% names(pathway_display_map),
                     pathway_display_map[pathway],
                     pathway)
  ) %>%
  select(timepoint, contrast_label, pathway, score_norm, n_sig) %>%
  complete(
    timepoint      = c("3 h", "10 h"),
    contrast_label = contrast_levels,
    ## complete()中也使用显示名称
    pathway        = path_display,
    fill           = list(score_norm = 0, n_sig = 0)
  ) %>%
  mutate(
    ## rev(path_display)：热图纵轴从上到下与path_keep排列一致
    pathway        = factor(pathway, levels = rev(path_display)),
    contrast_label = factor(contrast_label, levels = contrast_levels),
    timepoint      = factor(timepoint, levels = c("3 h", "10 h"))
  )

# 2️⃣ 绘制热图
## 颜色填充 = score_norm（col_blue=下调，白色=0，col_red=上调）
## 主要改动：
##   a) 评分公式和阈值仅在图例中定义，避免图内重复方法文字
##   b) x轴：45°角 + 标题标注"dose / Gy"
##   c) 色阶：对称limits + squish处理极端值，两时间点颜色一致
##   d) 图例：竖向堆叠，节省宽度
pb <- ggplot(pdat,
             aes(x = contrast_label, y = pathway, fill = score_norm)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_point(aes(size = n_sig),
             shape = 21, fill = NA, color = "black", stroke = 0.3) +
  facet_wrap(~ timepoint, nrow = 1) +
  scale_fill_gradient2(
    low      = col_blue,
    mid      = "#F7F7F7",
    high     = col_red,
    midpoint = 0,
    ## 对称色阶范围：确保白色严格对应0，两个时间点颜色可横向对比
    ## squish：超出范围的值截断到边界颜色，不显示为NA
    limits   = c(-5.5, 5.5),
    oob      = scales::squish,
    name     = "Pathway\nscore"
  ) +
  scale_size_continuous(
    range  = c(0.3, 3.5),
    breaks = c(0, 2, 5, 10),
    ## n_sig用斜体下标表示，与公式中的n对应
    name   = expression(italic(n)[sig])
  ) +
  labs(
    ## 评分公式、DEG阈值和k的定义由Figure 5 legend统一说明
    ## x轴标题标注Gy单位，不在每个标签中重复
    x        = "Gravity condition + Radiation dose [Gy]",
    y        = NULL
  ) +
  theme_journal() +
  theme(
    ## 45°角：在19cm宽度下比35°更紧凑，不重叠
    axis.text.x     = element_text(angle = 45, hjust = 1, size = sz_tick),
    axis.text.y     = element_text(size  = sz_tick),
    axis.title.x    = element_text(size  = sz_axis),
    panel.grid      = element_blank(),
    ## 图例竖向堆叠（fill在上，size在下），节省横向空间
    legend.position = "right",
    legend.box      = "vertical",
    legend.spacing  = unit(0.1, "cm")
  )

## ── 拼合主图并保存 ─────────────────────────────────────────
## heights比例：Panel a（3行）约占Panel b（12行）的1/3
fig5 <- (pa / pb) +
  plot_layout(heights = c(1, 2.2)) +
  plot_annotation(
    tag_levels = 'a',
    theme = theme(
      plot.tag = element_text(size   = sz_tag,
                              face   = "bold",
                              family = "Arial")
    )
  )

ggsave(
  file.path(out_dir, "Figure5_submission_v16.pdf"),
  fig5,
  width  = fig_width_cm,
  height = 14,
  units  = "cm",
  device = cairo_pdf,
  dpi    = fig_res
)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(out_dir, "Figure5_submission_v16.svg"), fig5,
         width = fig_width_cm, height = 14, units = "cm", device = svglite::svglite)
}
ggsave(file.path(out_dir, "Figure5_submission_v16.tiff"), fig5,
       width = fig_width_cm, height = 14, units = "cm", dpi = fig_res,
       device = ragg::agg_tiff, compression = "lzw")
ggsave(file.path(out_dir, "Figure5_submission_v16.png"), fig5,
       width = fig_width_cm, height = 14, units = "cm", dpi = fig_res,
       device = ragg::agg_png)
message("Figure 5主图已保存：Result/Figure5-v4.pdf")

#========================
# 10. 绘制Figure 5补充图1：SMG重写指数
#========================

message("\n===== 绘制补充图1：SMG重写指数 =====")

## 计算重写指数：同剂量下 SMG组score_norm − 1g组score_norm
## 正值 = SMG相对1g更上调；负值 = SMG相对1g更下调
rw <- path_sc %>%
  filter(pathway %in% path_top_auto, dose > 0) %>%
  select(timepoint, pathway, dose, condition, score_norm) %>%
  group_by(timepoint, pathway, dose, condition) %>%
  summarise(score_norm = mean(score_norm), .groups = "drop") %>%
  pivot_wider(names_from  = condition,
              values_from = score_norm,
              values_fill = 0) %>%
  mutate(rewrite = SMG - `1g`) %>%
  group_by(pathway) %>%
  mutate(max_abs = max(abs(rewrite), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(max_abs)) %>%
  slice_head(n = top_n_pathway) %>%
  mutate(
    pathway   = fct_reorder(pathway, max_abs, .desc = TRUE),
    dose_f    = factor(dose, levels = c(6.5, 13, 19.5),
                       labels = c("6.5 Gy", "13 Gy", "19.5 Gy")),
    timepoint = factor(timepoint, levels = c("3 h", "10 h"))
  )

ps1 <- ggplot(rw, aes(x = dose_f, y = pathway, fill = rewrite)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~ timepoint, nrow = 1) +
  scale_fill_gradient2(
    low      = col_blue,
    mid      = "#F7F7F7",
    high     = col_red,
    midpoint = 0,
    name     = "SMG − 1g"
  ) +
  labs(x = "Dose (Gy)", y = NULL) +
  theme_journal() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid  = element_blank()
  )

ggsave(
  file.path(out_dir, "Figure5_S1_rewriting_v16.pdf"),
  ps1,
  width  = fig_width_cm,
  height = 7,
  units  = "cm",
  device = cairo_pdf,
  dpi    = fig_res
)
message("补充图1已保存：Result/Figure5_S1_rewriting-v2.pdf")

#========================
# 11. 绘制Figure 5补充图2：DEG负担三联图
#========================

message("\n===== 绘制补充图2：DEG负担分析 =====")

## ── Panel a：DEG负担热图 ─────────────────────────────────────
## 颜色深浅反映各组显著DEG数量（浅蓝→深蓝），格内标注具体数值
ps2a <- ggplot(deg_burden,
               aes(x = contrast_label, y = timepoint, fill = n_sig)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = n_sig),
            size = sz_tick / .pt, family = "Arial") +
  scale_fill_gradient(
    low  = "#C6DBF0",   # 浅蓝
    high = col_blue,    # 深蓝（#0072B5）
    name = "Sig. DEGs"
  ) +
  labs(x = NULL, y = NULL) +
  theme_journal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid  = element_blank())

## ── Panel b：上调/下调基因数量柱状图 ─────────────────────────
## 上调（col_red，正值向上）/ 下调（col_blue，负值向下）
## 零线分隔方向，按时间点分面
bar_dat <- deg_burden %>%
  transmute(timepoint, contrast_label,
            Up   =  n_up,
            Down = -n_down) %>%   # 下调取负使柱向下延伸
  pivot_longer(cols = c("Up", "Down"),
               names_to = "direction", values_to = "n")

ps2b <- ggplot(bar_dat,
               aes(x = contrast_label, y = n, fill = direction)) +
  geom_col(width = 0.72, linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
  facet_wrap(~ timepoint, nrow = 1) +
  scale_fill_manual(values = c(Up = col_red, Down = col_blue)) +
  labs(x = NULL, y = "DEG count (up / down)", fill = NULL) +
  theme_journal() +
  theme(axis.text.x     = element_text(angle = 35, hjust = 1),
        legend.position = "top")

## ── Panel c：DEG负担 - 表型偏移散点图 ───────────────────────
## x = 显著DEG数量；y = 各表型指标相对基线偏移量
## 时间点用颜色区分（3h = col_blue，10h = col_red）
## 趋势线不含CI，仅展示方向
metric_labs <- c(
  "delta_mu_max" = expression(Delta * mu[max]),
  "delta_lag"    = expression(Delta * lag),
  "delta_AUC"    = expression(Delta * AUC)
)

bridge_long <- bridge %>%
  pivot_longer(cols      = c(delta_mu_max, delta_lag, delta_AUC),
               names_to  = "metric",
               values_to = "delta_value")

ps2c <- ggplot(bridge_long,
               aes(x = n_sig, y = delta_value,
                   color = timepoint, shape = condition)) +
  geom_hline(yintercept = 0, linetype = 2,
             linewidth = 0.3, color = "grey45") +
  geom_point(size = 1.5, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1,
             labeller = as_labeller(metric_labs, default = label_parsed)) +
  scale_color_manual(
    values = c(`3 h` = col_blue, `10 h` = col_red),
    name   = "RNA time"
  ) +
  scale_shape_manual(
    values = c(`1g` = 16, `SMG` = 17),
    name   = "Condition"
  ) +
  labs(x = "Significant DEG count",
       y = "Phenotypic shift vs 1g 0 Gy") +
  theme_journal() +
  theme(legend.position  = "top",
        panel.grid.minor = element_blank())

## ── 拼合补充图2并保存 ────────────────────────────────────────
fig5_s2 <- (ps2a / ps2b / ps2c) +
  plot_annotation(
    tag_levels = 'a',
    theme = theme(
      plot.tag = element_text(size   = sz_tag,
                              face   = "bold",
                              family = "Arial")
    )
  )

ggsave(
  file.path(out_dir, "Figure5_S2_DEGburden_v16.pdf"),
  fig5_s2,
  width  = fig_width_cm,
  height = 18,
  units  = "cm",
  device = cairo_pdf,
  dpi    = fig_res
)
message("补充图2已保存：Result/Figure5_S2_DEGburden-v2.pdf")

#========================
# 12. 保存所有中间数据CSV
#========================

message("\n===== 保存中间数据CSV =====")

## sig：注释后的显著DEG全集（Fig6.R直接读取，无需重复注释）
write_csv(sig,        file.path(out_dir, "fig5_sig_annotated.csv"))

## bridge：转录组-表型桥接数据（Fig6.R Panel b需要）
write_csv(bridge,     file.path(out_dir, "fig5_transcriptome_phenotype_bridge.csv"))

## path_sc：KEGG通路评分矩阵（供独立分析或Fig6.R使用）
write_csv(path_sc,    file.path(out_dir, "fig56_kegg_pathway_scores.csv"))

## deg_burden：各对比组DEG统计汇总
write_csv(deg_burden, file.path(out_dir, "fig5_deg_burden.csv"))

## rw：SMG重写指数数据（补充图1的数据来源）
write_csv(rw,         file.path(out_dir, "fig56_kegg_rewriting.csv"))

message("所有中间数据已保存至 Result/ 目录")
message("\n===== Figure 5 分析完成 =====")
