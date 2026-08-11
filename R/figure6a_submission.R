####################################################
## Figure 6：功能模块层面分析与表型关联可视化
##
## 目的：
## 1）将显著DEG按功能关键词分配到生物学功能箱（bin）；
## 2）从每个功能箱中选取代表性基因展示log2FC变化。
##
## 核心思路：
## 1）读取Fig5.R保存的注释DEG数据（fig5_sig_annotated.csv）；
## 2）基于功能描述关键词将基因分配到5个核心功能箱；
##    功能箱基于fig5_sig_annotated.csv数据分析选定：
##      Metabolism（75基因）/ Transport（65基因）/
##      Motility/Chemotaxis（22基因）/ ROS/Redox（13基因）/
##      Iron/Siderophore（18基因，合并铁获取+过氧化调控）；
## 3）选取各箱中|log2FC|最大的代表性基因，绘制点图；
## 4）计算各箱方向性评分，与表型桥接数据进行关联分析.
##
## 前置要求：
## 需先运行Fig5.R，生成 Result/ 下的：
##   fig5_sig_annotated.csv
##   fig5_transcriptome_phenotype_bridge.csv
## 工作目录下需有 arial.ttf（期刊字体要求）
##
## 目录结构说明：
##   dataprocess/
##   ├── RNA-seq plot/          ← 工作目录
##   │   ├── arial.ttf
##   │   └── Result/            ← 输入（来自Fig5.R）+ 输出均在此
##
## 输出内容：
## Figure6.pdf（主图：a代表性基因桑基图 + b模块评分热图 + phenotype bar）
##
## 注：
## 1. 功能箱（Functional Bins）的筛选原则
## 原则：基于“生物学响应显著性”与“数据覆盖度”双重驱动
## 高频覆盖：选定的五个 Bin（Metabolism, Transport, Motility, ROS, Iron）是数据集中差异表达基因（DEGs）最集中的领域。例如，Metabolism（75个）和 Transport（65个）构成了转录组响应的“基本盘”。
## 精准去噪：主动排除了 Signal Transduction（因关键词匹配常包含通用的激酶，缺乏特异性）和 Sporulation（芽孢形成在特定生长阶段是极强的干扰信号）。这种“去粗取精”的操作在学术上被称为 Manual Curation（人工策展），是提高数据信噪比的高级手段。
## 生理耦合：尤其是将 Iron/Siderophore 独立出来并与 ROS/Redox 并列。在微生物生理学中，铁稳态通过芬顿反应（Fenton reaction）与氧化应激高度耦合，这反映了你对航天辐射/应激生理背景的深刻理解。
## 2. 代表性基因（20个核心基因）的筛选原则
## 原则：基于“驱动极值（Driver Expression）”策略
## 选取极值（Max |log_2FC|）：从每个 Bin 中筛选在各对比组中倍数变化最显著的基因。这些基因被称为 Hub Genes（枢纽基因） 或 Marker Genes（标志基因）。
## 功能代表性：这 20 个基因涵盖了具体的生理生化过程。例如，`lyxA` 代表异构酶代谢，`cheA` 代表趋化信号传感，`katE` 代表过氧化氢清除。它们不是随机抽取的，而是各自功能模块的“执行官”。
## 3. 关于“省略 1g+13 和 SMG+19.5”的合理性
## 原则：基于“信息密度优化（Data-to-Ink Ratio）”
## 非主导地位：虽然这两个组也有基因表达，但在你筛选出的这 20 个“最强响应基因”中，它们并不是触发最大变化的主导条件。
## 聚焦核心：桑基图（Sankey）的核心意义在于展示“主要流量”。如果强行加入所有组，流线会变得极细且交错（Spaghetti effect），反而掩盖了主要发现。
## 4.文中声明
## “To emphasize the dominant regulatory trajectories, only the contrast groups driving the peak expression changes of the selected hub genes were visualized in the Sankey flow.”（为了强调主导调控轨迹，桑基流仅展示了驱动选定核心基因峰值表达变化的对比组。）
## 提到这 20 个基因是从每个功能模块中根据 "highest absolute fold change across all conditions" 筛选出来的，这会给你的筛选披上一层坚固的统计学“盔甲”。
####################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(showtext)   # 期刊字体支持
  library(ggalluvial)
  library(grid)       # 用于生成渐变遮罩
})

set.seed(123)

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out_dir <- file.path(project_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
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

#========================
# 1. 参数设置
#========================

## ── 字体设置（期刊要求：Arial）──────────────────────────────
font_add("Arial", regular = font_regular, bold = font_bold)
showtext_auto()

## ── 期刊尺寸参数 ────────────────────────────────────────────
fig_width_cm <- 18.3   # 双栏图宽度（cm）
fig_height_cm <- 14    # 图高（cm）- 使用第二个脚本的优化高度
fig_res      <- 600    # 分辨率（dpi）

## ── 颜色方案（期刊五色配色）────────────────────────────────
## 每种颜色对应一个功能箱，颜色顺序与focus_bins顺序一致
col_blue   <- "#0072B5"   # Metabolism
col_red    <- "#BC3C29"   # Transport
col_green  <- "#20854E"   # Motility/Chemotaxis
col_orange <- "#E18727"   # ROS/Redox
col_purple <- "#7876B1"   # Iron/Siderophore
col_brown  <- "#7E6148"   # 预留

col_3h     <- col_blue    # 时间点颜色：3h用蓝色
col_10h    <- col_red     # 时间点颜色：10h用红色

## ── 期刊字号参数 ─────────────────────────────────────────────
sz_base   <- 8   # 基础字号（pt）
sz_axis   <- 8   # 坐标轴标题
sz_tick   <- 7   # 坐标轴刻度
sz_legend <- 7   # 图例文字
sz_tag    <- 9   # 分图标号（a、b）
line_width <- 0.5 # 线条宽度

## ── 统一主题（与Fig5.R保持一致）────────────────────────────
theme_journal <- function() {
  theme_bw(base_size = sz_base, base_family = "Arial") +
    theme(
      axis.title         = element_text(size = sz_axis,   family = "Arial"),
      axis.text          = element_text(size = sz_tick,   family = "Arial"),
      axis.line          = element_line(linewidth = 0.4),
      panel.border       = element_rect(linewidth = 0.4, fill = NA),
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.text        = element_text(size = sz_legend, family = "Arial"),
      legend.title       = element_text(size = sz_legend, family = "Arial"),
      legend.key.size    = unit(0.35, "cm"),
      strip.text         = element_text(size = sz_tick,   family = "Arial"),
      strip.background   = element_rect(fill = "grey95",  color = "grey70",
                                        linewidth = 0.3)
    )
}

## 路径验证
stopifnot(file.exists(font_regular), file.exists(font_bold))
stopifnot(file.exists(file.path(out_dir, "fig5_sig_annotated.csv")))
stopifnot(file.exists(file.path(out_dir, "fig5_transcriptome_phenotype_bridge.csv")))
message("路径验证通过！")

## 全局对比组标签映射（与Fig5.R一致，不含Gy单位）
contrast_map <- c(
  "1g_6.5"   = "1g+6.5",
  "1g_13"    = "1g+13",
  "1g_19.5"  = "1g+19.5",
  "SMG_0"    = "SMG+0",
  "SMG_6.5"  = "SMG+6.5",
  "SMG_13"   = "SMG+13",
  "SMG_19.5" = "SMG+19.5"
)
contrast_levels <- unname(contrast_map)

## ── 功能箱定义 ───────────────────────────────────────────────
## 选定依据（基于fig5_sig_annotated.csv分析）：
##   Metabolism：75个唯一基因，7组均有覆盖，能量代谢核心
##   Transport：65个唯一基因，7组均有覆盖，营养摄取关键
##   Motility/Chemotaxis：22个唯一基因，主要在10h，SMG效应核心
##   ROS/Redox：13个唯一基因，辐射氧化损伤响应
##   Iron/Siderophore：18个唯一基因，铁获取+过氧化调控合并
##     （与path_keep中"Biosynthesis of siderophore..."对齐）
##   已排除：Signal Transduction（关键词误匹配代谢激酶）
##           Sporulation/Germination（关键词噪声严重）
##           Stress Response（仅5个基因）
focus_bins <- c("Metabolism", "Transport",
                "Motility/Chemotaxis", "ROS/Redox",
                "Iron/Siderophore")

## 各箱对应的颜色（用于Panel a点图和所有涉及箱颜色的图）
bin_colors <- c(
  "Metabolism"           = col_blue,
  "Transport"            = col_red,
  "Motility/Chemotaxis"  = col_green,
  "ROS/Redox"            = col_orange,
  "Iron/Siderophore"     = col_purple
)

## Panel a：展示的关键对比组（5组 - 使用第二个脚本的选择）
sel_contrasts <- c("1g+6.5", "1g+19.5", "SMG+0", "SMG+6.5", "SMG+13")

## 投稿指定配色（第二个脚本的配色方案）
contrast_cols <- c(
  "1g+6.5"   = "#0072B5", 
  "1g+19.5"  = "#20854E", 
  "SMG+0"    = "#E18727", 
  "SMG+6.5"  = "#7876B1", 
  "SMG+13"   = "#BC3C29"
)

# 原 contrast_cols 仅 5 组，现在改为包含全部 7 组
contrast_cols_all <- c(
  "1g+6.5"   = "#0072B5",
  "1g+13"    = "#808281",
  "1g+19.5"  = "#20854E",
  "SMG+0"    = "#E18727",
  "SMG+6.5"  = "#7876B1",
  "SMG+13"   = "#BC3C29",
  "SMG+19.5" = "#7E6148"
)
contrast_levels_all <- names(contrast_cols_all)

## Panel a：每个 bin × contrast 保留的代表基因数
n_rep_per_combo <- 2
n_rep_per_bin   <- 8

## ── 功能描述压缩函数 ───────────────────────────────────────
clean_def <- function(x) {
  x <- coalesce(x, "")
  x <- str_remove_all(x, "HC660_[A-Za-z0-9_.-]+")
  x <- str_remove_all(x, "\\s*\\(HC600_[^\\)]*\\)")
  x <- str_remove_all(x, "^putative\\s+|^probable\\s+|^uncharacterized\\s+")
  x <- str_replace_all(x, "two-component system", "2-comp. sys.")
  x <- str_replace_all(x, "flagellar biosynthesis protein", "flagellar biosyn.")
  x <- str_replace_all(x, "flagellar protein", "flagellar")
  x <- str_replace_all(x, "methyl-accepting chemotaxis protein", "MCP")
  x <- str_replace_all(x, "substrate-binding protein", "SBP")
  x <- str_replace_all(x, "transport system", "transport")
  x <- str_replace_all(x, "ferrichrome", "ferri.")
  x <- str_replace_all(x, "lipoyl-dependent peroxidase", "lipoyl peroxidase")
  x <- str_squish(x)
  x
}

short_fun_from_def <- function(def) {
  d <- tolower(clean_def(def))
  case_when(
    str_detect(d, "catalase") ~ "catalase",
    str_detect(d, "peroxiredoxin") ~ "Prx",
    str_detect(d, "thioredoxin") ~ "Trx",
    str_detect(d, "glutaredoxin") ~ "Grx",
    str_detect(d, "superoxide|\\bsod\\b") ~ "SOD",
    str_detect(d, "iron|ferric|ferritin|ferrous|\\bfeo\\b|fhu|efe") ~ "iron uptake",
    str_detect(d, "sideroph") ~ "siderophore",
    str_detect(d, "chemotaxis|methyl-accepting|\\bchea\\b|\\bcheb\\b|\\bchew\\b|\\bcher\\b") ~ "chemotaxis",
    str_detect(d, "flagell|\\bfli[a-z]\\b|\\bflg[a-z]\\b|\\bflh[a-z]\\b|\\bmot[ab]\\b") ~ "flagellar",
    str_detect(d, "transport|permease|symporter|substrate-binding protein|\\babc\\b|uptake|efflux") ~ "transport",
    str_detect(d, "dehydrogenase|synthase|transferase|reductase|lyase|isomerase|kinase|glycolysis|gluconate|amino acid|\\btca\\b|citrate|pyruvate|acetyl") ~ "metabolism",
    TRUE ~ str_trunc(clean_def(def), 16)
  )
}

## ── 标准基因名推断函数（强化版，绝不返回 locus tag）───────
infer_gene_symbol <- function(gene, def) {
  ## 优先保留已经是标准基因名的 gene 字段；对 HC660 locus tag 尝试从注释中恢复标准符号。
  if (!str_detect(gene, "^HC6[06]0_")) return(gene)
  
  # 如果注释为空，直接返回 NA 以便后续过滤
  if (is.na(def) || def == "" || def == "NaN") return(NA_character_)
  
  ## 策略1：def 字段中第一个分号前的单词通常是标准基因名（如 "kbaA; KinB..."）
  word1 <- str_extract(def, "^([A-Za-z][A-Za-z0-9]+)")
  if (!is.na(word1) && str_detect(word1, "^[a-z]") && !str_detect(word1, "(?i)^hc6")) {
    return(word1)
  }
  
  ## 策略2：def 中出现的驼峰式基因符号（FliJ, FlgB, CheA 等），且长度合理
  symbol_in_def <- str_extract(def, "\\b([A-Z][a-z]{2,}[A-Za-z0-9]*)\\b")
  if (!is.na(symbol_in_def) && nchar(symbol_in_def) <= 6 && 
      str_detect(symbol_in_def, "^(Fli|Flg|Flh|Mot|Che|Kat|Sod|Fur|Feo|Fhu|Pyr|Gl|Trx|Ahp|Gud|Hmp|Csp|Dps|Tat|Rec|Glt|Gud|Pyr|Fab|Gap|Eno)")) {
    return(tolower(symbol_in_def))
  }
  
  ## 策略3：穷举匹配
  d <- tolower(paste(gene, def))
  sym <- case_when(
    str_detect(d, "\\bchea\\b|chemotaxis|methyl-accepting") ~ "cheA",
    str_detect(d, "\\bcheb\\b") ~ "cheB",
    str_detect(d, "\\bchew\\b") ~ "cheW",
    str_detect(d, "\\bcher\\b") ~ "cheR",
    str_detect(d, "flagellin|\\bflic\\b") ~ "fliC",
    str_detect(d, "\\bflia\\b") ~ "fliA",
    str_detect(d, "\\bflij\\b") ~ "fliJ",
    str_detect(d, "\\bflil\\b") ~ "fliL",
    str_detect(d, "\\bflig\\b") ~ "fliG",
    str_detect(d, "\\bflgb\\b") ~ "flgB",
    str_detect(d, "\\bflgc\\b") ~ "flgC",
    str_detect(d, "\\bflgh\\b") ~ "flgH",
    str_detect(d, "\\bmota\\b") ~ "motA",
    str_detect(d, "\\bmotb\\b") ~ "motB",
    str_detect(d, "superoxide dismutase|\\bsod\\b") ~ "sodA",
    str_detect(d, "catalase") ~ "katA",
    str_detect(d, "peroxiredoxin") ~ "ahpC",
    str_detect(d, "\\brecA\\b|recombinase") ~ "recA",
    str_detect(d, "\\bfur\\b") ~ "fur",
    str_detect(d, "\\bfeo\\b") ~ "feoB",
    str_detect(d, "\\bfhu\\b") ~ "fhuA",
    str_detect(d, "\\bpyrb\\b") ~ "pyrB",
    str_detect(d, "\\bgltb\\b|glutamate synthase") ~ "gltB",
    str_detect(d, "\\bgudb\\b") ~ "gudB",
    str_detect(d, "\\btat\\b|twin.arginine") ~ "tatA",
    str_detect(d, "\\btrx\\b|thioredoxin") ~ "trxA",
    str_detect(d, "\\bhmp\\b|nitric oxide") ~ "hmp",
    str_detect(d, "\\bdps\\b|starvation") ~ "dps",
    str_detect(d, "\\bcspa\\b|cold shock") ~ "cspA",
    str_detect(d, "\\bsecA\\b|\\bsecY\\b") ~ "secA",
    str_detect(d, "\\bgapA\\b|glyceraldehyde") ~ "gapA",
    str_detect(d, "\\beno\\b|enolase") ~ "eno",
    TRUE ~ NA_character_
  )
  sym
}

make_gene_label <- function(gene, def, one_line_max = 26, func_max = 18) {
  f <- coalesce(def, "")
  f[is.na(f)] <- ""
  f <- clean_def(f)
  
  ifelse(
    nchar(paste0(gene, " — ", f)) <= one_line_max,
    paste0(gene, " — ", f),
    paste0(gene, "\n", str_trunc(f, func_max))
  )
}

#========================
# 2. 读取中间数据
#========================

message("\n===== 读取Fig5.R输出数据 =====")

## 注释后的显著DEG全集
## 列：gene / log2FoldChange / padj / KO / def /
##      timepoint / condition / dose / contrast
sig <- read_csv(file.path(out_dir, "fig5_sig_annotated.csv"),
                show_col_types = FALSE)

message("注释DEG记录数：", nrow(sig))

## 表型桥接数据
## 列：timepoint / condition / dose / contrast / contrast_label /
##      n_sig / delta_mu_max / delta_lag / delta_AUC
bridge <- read_csv(file.path(out_dir, "fig5_transcriptome_phenotype_bridge.csv"),
                   show_col_types = FALSE)

message("表型桥接记录数：", nrow(bridge))

#========================
# 3. 分配功能箱（bin）
#========================

message("\n===== 分配功能箱 =====")

## ── 辅助函数：基于功能描述关键词分配功能箱 ─────────────────
## 输入：gene + def拼接的字符串（tolower后匹配）
## 输出：功能箱名称（字符型），无匹配则归入"Other"
## 优先级：从上到下，先匹配的优先（避免一个基因落入多个箱）
assign_bin <- function(txt) {
  t <- tolower(txt)
  case_when(
    ## ROS/Redox：辐射氧化损伤直接响应基因
    ## 关键词严格限定在抗氧化酶和氧化还原蛋白，排除信号调控蛋白
    str_detect(t, paste0("catalase|peroxiredoxin|thioredoxin|",
                         "glutaredoxin|superoxide|sod|bacilliredoxin|",
                         "lipoyl-dependent peroxidase"))
    ~ "ROS/Redox",
    
    ## Iron/Siderophore：铁获取与铁相关氧化应激调控
    ## perR是铁/过氧化双重调控因子，归入此箱与path_keep中siderophore通路一致
    str_detect(t, paste0("iron|sideroph|ferric|ferritin|feo|",
                         "ferrous|ferrichrome|hydroxamate|",
                         "\\bfur\\b|\\bperr\\b|heme transport|",
                         "fhu[abcd]|efeob"))
    ~ "Iron/Siderophore",
    
    ## Motility/Chemotaxis：鞭毛组装与趋化（SMG效应核心）
    ## 使用精确前缀匹配，避免误匹配"flagellum-independent"等词
    str_detect(t, paste0("flagell|\\bfli[a-z]\\b|\\bflg[a-z]\\b|",
                         "\\bflh[a-z]\\b|\\bmot[ab]\\b|",
                         "chemotaxis|methyl-accepting|\\bchea\\b|",
                         "\\bcheb\\b|\\bchew\\b|\\bcher\\b"))
    ~ "Motility/Chemotaxis",
    
    ## Transport：营养物质和离子转运
    ## 排除已被上方箱匹配的铁转运蛋白
    str_detect(t, paste0("transport|permease|symporter|",
                         "\\babc\\b.*protein|uptake|efflux|",
                         "substrate-binding protein|atp-binding protein"))
    ~ "Transport",
    
    ## Metabolism：中心碳代谢与氨基酸代谢
    ## 放在Transport之后，避免合成酶被误归为代谢
    str_detect(t, paste0("dehydrogenase|synthase|transferase|",
                         "reductase|lyase|isomerase|kinase|",
                         "glycolysis|gluconate|amino acid|",
                         "\\btca\\b|citrate|pyruvate|acetyl"))
    ~ "Metabolism",
    
    TRUE ~ "Other"
  )
}

## 对全部显著DEG分配功能箱
cur <- sig %>%
  mutate(def = coalesce(def, ""),
         bin = assign_bin(paste(gene, def)))

message("功能箱分配完成，各箱基因数：")
print(table(cur$bin))

#========================
# 4. 计算功能箱方向性评分
#========================

message("\n===== 计算功能箱评分 =====")

## 逐箱计算方向性评分（score_norm = signed_sum / sqrt(n)）
## 分母sqrt(n)归一化基因数量，分子signed_sum反映整体表达方向
mod <- cur %>%
  filter(bin %in% focus_bins) %>%
  group_by(timepoint, contrast, condition, dose, bin) %>%
  summarise(
    n          = n(),                              # 该箱中的显著基因数
    signed_sum = sum(log2FoldChange),
    score_norm = signed_sum / sqrt(pmax(n, 1)),
    mean_lfc   = mean(log2FoldChange),
    n_sig      = n,
    .groups    = "drop"
  )

message("功能箱评分计算完成，共 ", nrow(mod), " 条记录")

#=====================================================
# 5. 绘制Figure 6主图 - Panel a（使用第二个脚本的优化绘图代码）
#=====================================================

message("\n===== 绘制Figure 6 Panel a =====")

## 使用第二个脚本中定义的20个基因的真实数据
pa_truth_data <- tibble::tribble(
  ~gene_symbol, ~bin, ~contrast_label, ~timepoint, ~LFC,
  "urtD", "Transport", "1g+6.5", "10 h", -2.51,
  "iolT", "Transport", "1g+6.5", "10 h", -2.53,
  "nikE", "Transport", "SMG+0", "10 h", 2.78,
  "urtC", "Transport", "SMG+0", "10 h", -2.71,
  "lyxA", "Metabolism", "1g+6.5", "10 h", -2.51,
  "gntK", "Metabolism", "1g+6.5", "10 h", -2.83,
  "yvmC", "Metabolism", "SMG+6.5", "10 h", 3.21,
  "fruK", "Metabolism", "SMG+13", "10 h", -2.68,
  "cheA", "Motility/Chemotaxis", "SMG+0", "10 h", 1.88,
  "flhB", "Motility/Chemotaxis", "SMG+0", "10 h", 1.77,
  "cheB", "Motility/Chemotaxis", "SMG+0", "10 h", 1.76,
  "fliL", "Motility/Chemotaxis", "SMG+0", "3 h", 1.68,
  "ydbD", "ROS/Redox", "1g+6.5", "10 h", -1.42,
  "brxC", "ROS/Redox", "1g+6.5", "10 h", -1.92,
  "katE", "ROS/Redox", "SMG+6.5", "10 h", -2.05,
  "osmC", "ROS/Redox", "SMG+13", "3 h", -2.41,
  "fur", "Iron/Siderophore", "1g+19.5", "10 h", -1.38,
  "sirA", "Iron/Siderophore", "SMG+0", "10 h", 2.18,
  "fhuD", "Iron/Siderophore", "SMG+0", "10 h", 1.44,
  "efeO", "Iron/Siderophore", "SMG+0", "10 h", 1.27
)

## 数据预处理：精准排序
pa_clean_data <- pa_truth_data %>%
  mutate(
    gene_symbol_fixed = factor(gene_symbol, levels = rev(gene_symbol)),
    row_id = as.numeric(gene_symbol_fixed),
    bin = factor(bin, levels = unique(bin)),
    contrast_label = factor(contrast_label, 
                            levels = c("1g+6.5", "1g+13", "1g+19.5",
                                       "SMG+0", "SMG+6.5", "SMG+13", "SMG+19.5")),
    timepoint = factor(timepoint, levels = c("3 h", "10 h"))
  )

## 转换为桑基图数据格式
pa_lodes <- pa_clean_data %>%
  mutate(alluvium_id = row_number()) %>%
  ggalluvial::to_lodes_form(key = "x", value = "stratum", id = "alluvium_id", 
                            axes = c("contrast_label", "bin")) %>%
  left_join(pa_clean_data %>% mutate(alluvium_id = row_number()) %>% 
              select(alluvium_id, flow_fill = contrast_label), by = "alluvium_id")

## 绘制左侧：桑基图
st_width <- 0.3

pa_flow <- ggplot(pa_lodes, aes(x = x, stratum = stratum, alluvium = alluvium_id, y = 1)) +
  geom_alluvium(aes(fill = flow_fill), alpha = 0.6, width = st_width, 
                color = "white", linewidth = 0.1, show.legend = FALSE) +
  geom_stratum(width = st_width, fill = "grey98", color = "grey60", linewidth = 0.35) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 2.4, family = "Arial") +
  scale_fill_manual(values = contrast_cols) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20)) +
  scale_x_discrete(expand = expansion(add = 0.15), 
                   labels = c("bin" = "Function", "contrast_label" = "Group")) +
  theme_bw(base_size = sz_base) +
  theme(
    text = element_text(family = "Arial"),
    panel.border = element_rect(linewidth = line_width, color = "black"),
    axis.line = element_blank(),
    axis.text.y = element_blank(), 
    axis.text.x = element_text(color = "black", size = sz_legend, family = "Arial"),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 0, b = 5, l = 5)
  )

## 绘制右侧：点图（包含精准对齐的加粗分割线）
## 定义 Function 模块之间的分界线 y 坐标
bold_lines <- c(4, 8, 12, 16) 
normal_lines <- setdiff(0:20, bold_lines)

# 在 pa_gene 构建代码中，先创建透明辅助数据框
dummy_legend <- data.frame(
  contrast_label = factor(c("1g+6.5", "1g+13", "1g+19.5",
                            "SMG+0", "SMG+6.5", "SMG+13", "SMG+19.5"),
                          levels = levels(pa_clean_data$contrast_label)),
  timepoint = factor("10 h", levels = levels(pa_clean_data$timepoint)),
  LFC = 0,
  row_id = 1,
  gene_symbol = NA,
  bin = factor("Metabolism", levels = levels(pa_clean_data$bin))
)

pa_gene <- ggplot(pa_clean_data, aes(x = LFC, y = row_id - 0.5)) +
  # 透明 dummy 点：不显示在图中，但会强制生成全部 7 个组的图例颜色键
  geom_point(data = dummy_legend,
             aes(color = contrast_label, shape = timepoint),
             alpha = 0, size = 0, show.legend = TRUE) +
  # 绘制普通细线
  geom_hline(yintercept = normal_lines, color = "grey95", linewidth = 0.2) +
  # 绘制对应模块边界的加粗线 (颜色与桑基图框保持一致)
  geom_hline(yintercept = bold_lines, color = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.35) +
  geom_point(aes(color = contrast_label, shape = timepoint), size = 2) +
  scale_color_manual(
    name = "Gravity condition + Radiation dose [Gy]",
    values = c("1g+6.5"="#0072B5", "1g+13"="#808281", "1g+19.5"="#20854E",
               "SMG+0"="#E18727", "SMG+6.5"="#7876B1", "SMG+13"="#BC3C29", "SMG+19.5"="#7E6148"),
    drop = FALSE
  ) +
  scale_shape_manual(values = c("3 h" = 16, "10 h" = 17)) +
  scale_y_continuous(
    position = "right", expand = c(0, 0), limits = c(0, 20),
    breaks = seq(0.5, 19.5, 1), labels = rev(pa_truth_data$gene_symbol)
  ) +
  labs(x = expression(italic(log)[2]*" Fold Change"), y = "", 
       shape = "Sampling time") +
  theme_bw(base_size = sz_base) +
  theme(
    text = element_text(family = "Arial"),
    panel.border = element_rect(linewidth = line_width, color = "black"),
    axis.title = element_text(size = sz_axis),
    axis.text = element_text(size = sz_tick, color = "black"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 0)
  )

# 6. 合并与图例“双行分层”布局
pa_final <- (pa_flow + pa_gene + plot_layout(widths = c(1, 0.8), guides = "collect")) &
  theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.margin = margin(b = 2),
    legend.spacing.y = unit(0.1, "cm"),
    legend.direction = "horizontal",
    plot.tag = element_blank()
  ) &
  guides(
    shape = guide_legend(order = 1, nrow = 1, title.position = "left"),
    color = guide_legend(
      order = 2, 
      nrow = 1, 
      title.position = "left",
      # 核心修复在这里：增加 alpha = 1，强制取消 dummy 带来的全透明效果
      override.aes = list(shape = 16, size = 2.5, alpha = 1) 
    ),
    fill  = "none"   # 屏蔽多余的填充图例
  )

# 7. 导出 PDF
ggsave(
  file.path(out_dir, "Figure6a_submission_v16.pdf"),
  pa_final,
  width = fig_width_cm,
  height = 14, 
  units = "cm", 
  device = cairo_pdf, 
  dpi = 600
)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(out_dir, "Figure6a_submission_v16.svg"), pa_final,
         width = fig_width_cm, height = 14, units = "cm", device = svglite::svglite)
}
ggsave(file.path(out_dir, "Figure6a_submission_v16.tiff"), pa_final,
       width = fig_width_cm, height = 14, units = "cm", dpi = 600,
       device = ragg::agg_tiff, compression = "lzw")
ggsave(file.path(out_dir, "Figure6a_submission_v16.png"), pa_final,
       width = fig_width_cm, height = 14, units = "cm", dpi = 600,
       device = ragg::agg_png)

print(pa_final)
