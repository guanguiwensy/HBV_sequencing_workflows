````markdown
# RNA-seq 上下游分析完整流程 (含HBV整合分析)

## 简介
欢迎使用本RNA-seq数据分析流程！该流程旨在指导您完成从原始测序数据（FASTQ文件）到下游差异表达分析和可变剪接分析的全过程。特别地，本流程考虑了样本中可能存在的HBV（乙型肝炎病毒）序列，并整合了相应的分析步骤。

## 流程概览:

**上游分析 (命令行环境):**
- 原始数据质量控制 (QC)
- 数据清洗 (Trimming)
- HBV 基因型鉴定
- 参考基因组与注释文件准备
- 参考基因组索引构建
- 序列比对
- 基因/转录本表达定量

**下游分析 (主要使用R语言):**
- 差异表达基因分析 (DGE)
- 可变剪接分析
- HBV剪接事件坐标还原
- HBV RNA 分析
  - HBV pgRNA 表达分析
  - HBV转录本覆盖度可视化及坐标还原

## 1. 上游分析

### 1.1 原始数据质量控制 (QC)
使用 `FastQC` 对原始 FASTQ 文件进行质量评估。

```bash
# 安装 FastQC (如果尚未安装)
# sudo apt-get install fastqc # Debian/Ubuntu
# conda install -c bioconda fastqc # Conda

# 为每个 FASTQ 文件运行 FastQC (示例)
fastqc sample1_R1.fastq.gz sample1_R2.fastq.gz
fastqc sample2_R1.fastq.gz sample2_R2.fastq.gz
# ... 为您的所有样本运行

# 推荐使用 MultiQC 汇总所有 FastQC 报告
# pip install multiqc
multiqc . # 在包含 FastQC 输出的目录中运行
````

`FastQC`的报告将帮助您了解测序数据的质量，例如碱基质量分布、接头含量等。

### 1.2 数据清洗 (Trimming)

使用 `Trimmomatic` 或 `Cutadapt` 去除测序接头序列和低质量碱基。

```bash
# 安装 Trimmomatic (如果尚未安装)
# sudo apt-get install trimmomatic # Debian/Ubuntu
# conda install -c bioconda trimmomatic # Conda

# 示例：使用 Trimmomatic PE (Paired-End) 模式
java -jar /path/to/trimmomatic-0.39.jar PE \
    -threads 4 \
    sample1_R1.fastq.gz sample1_R2.fastq.gz \
    sample1_R1_paired_trimmed.fastq.gz sample1_R1_unpaired_trimmed.fastq.gz \
    sample1_R2_paired_trimmed.fastq.gz sample1_R2_unpaired_trimmed.fastq.gz \
    ILLUMINACLIP:/path/to/adapters.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# 注意:
# - /path/to/trimmomatic-0.39.jar: 替换为您的Trimmomatic jar文件路径。
# - /path/to/adapters.fa: 替换为包含您测序所用接头序列的FASTA文件路径。
# - 对所有样本重复此操作。
```

数据清洗是确保下游分析准确性的关键步骤。

### 1.3 HBV 基因型鉴定

由于样本中可能存在HBV，且其基因型对于后续分析可能很重要，建议先进行基因型鉴定。

**策略：**

1.  **提取疑似HBV reads：** 将一部分清洗后的reads比对到一个包含多种已知HBV参考基因组的数据库。
2.  **BLAST比对：** 将上一步提取出的HBV reads，或通过对这些reads进行de novo组装得到的HBV contigs，与NCBI的HBV序列数据库（如nt库）进行BLAST比对，以确定最相似的基因型。

**示例命令 (概念性):**

```bash
# 1. 准备一个包含多种HBV基因型参考序列的FASTA文件 (例如: hbv_genotypes_db.fasta)
#    您可以从 NCBI Viral Genomes Resource ([https://www.ncbi.nlm.nih.gov/genome/viruses/](https://www.ncbi.nlm.nih.gov/genome/viruses/)) 下载。

# 2. 使用 bowtie2 (或类似比对工具) 将一部分trimmed reads比对到此HBV数据库
# conda install -c bioconda bowtie2 samtools
bowtie2-build hbv_genotypes_db.fasta hbv_genotypes_db_idx
bowtie2 -x hbv_genotypes_db_idx \
        -1 sample1_R1_paired_trimmed.fastq.gz \
        -2 sample1_R2_paired_trimmed.fastq.gz \
        -S temp_hbv_align.sam --al-conc-gz mapped_to_hbv_%.fastq.gz --threads 4

# mapped_to_hbv_1.fastq.gz 和 mapped_to_hbv_2.fastq.gz 将包含映射到HBV数据库的reads

# 3. (可选) De novo 组装这些reads得到HBV contigs
# 例如使用 SPAdes:
# conda install -c bioconda spades
# spades.py --rna -1 mapped_to_hbv_1.fastq.gz -2 mapped_to_hbv_2.fastq.gz -o hbv_assembly

# 4. 使用 blastn 将提取的reads或contigs与NCBI nt数据库比对
# conda install -c bioconda blast
# blastn -query mapped_to_hbv_1.fastq.gz -db nt -remote -out hbv_blast_results.txt -max_target_seqs 10 -html
# 或者使用组装的contigs:
# blastn -query hbv_assembly/contigs.fasta -db nt -remote -out hbv_contigs_blast_results.txt -max_target_seqs 10 -html

# 分析 blastn 的输出 (hbv_blast_results.txt 或 hbv_contigs_blast_results.txt)
# 以确定样本中最主要的HBV基因型及其最相似的参考序列。
```

**注意:** 即便后续分析统一使用HBV基因型B的参考序列，此步骤仍有助于了解样本中实际存在的HBV类型，或确认是否存在预期的基因型B。

### 1.4 参考基因组与注释文件准备 (整合人与HBV)

为了同时分析宿主（人）和病毒（HBV）的转录本，需要准备一个包含两者序列的联合参考基因组和相应的基因注释（GTF）文件。

#### 1.4.1 人类参考基因组和GTF

建议从 GENCODE 或 Ensembl 下载最新版的人类参考基因组 (FASTA) 和基因注释文件 (GTF)。

```bash
# 示例下载链接 (请访问官网获取最新版本链接)
# GRCh38 参考基因组
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_XX/GRCh38.primary_assembly.genome.fa.gz
# GENCODE 全面基因注释
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_XX/gencode.vXX.annotation.gtf.gz

gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.vXX.annotation.gtf.gz

HUMAN_GENOME_FA="GRCh38.primary_assembly.genome.fa" # 替换为您的文件名
HUMAN_GTF="gencode.vXX.annotation.gtf"         # 替换为您的文件名
```

#### 1.4.2 HBV Genotype B 参考基因组和GTF (用户提供)

本流程假设您将提供一个针对HBV基因型B的FASTA参考基因组文件 (`HBV_B.fasta`) 和相应的GTF注释文件 (`HBV_B.gtf`)。

**重要说明:**
您提供的 `HBV_B.fasta` 文件应具有特定的线性化结构，以准确分析可能跨越HBV环状基因组复制起点的转录本（如pgRNA）。该结构描述如下：

  * 原始HBV基因型B全长为3215bp。
  * `HBV_B.fasta` 中的序列由两部分拼接而成：
    1.  第一部分：取原始HBV基因型B序列的 1600nt 至 3215nt。
    2.  第二部分：取原始HBV基因型B序列的 1nt 至 1950nt。
  * 这两部分序列直接拼接，形成一个总长度为 (3215 - 1600 + 1) + 1950 = 1616 + 1950 = 3566bp 的线性参考序列。
  * 您提供的 `HBV_B.gtf` 文件应包含HBV基因（特别是pgRNA）在该线性化 `HBV_B.fasta` 参考序列上的坐标注释。

请将您准备好的文件命名为：

```bash
HBV_REF_FASTA="HBV_B.fasta"
HBV_REF_GTF="HBV_B.gtf"
```

并确保这两个文件位于您的工作目录中或提供正确路径。

#### 1.4.3 合并参考基因组和GTF

将人类和HBV的参考基因组及GTF文件分别合并。

```bash
# 合并FASTA文件
COMBINED_GENOME_FA="human_hbv_combined.fasta"
cat ${HUMAN_GENOME_FA} ${HBV_REF_FASTA} > ${COMBINED_GENOME_FA}

# 合并GTF文件
COMBINED_GTF="human_hbv_combined.gtf"
cat ${HUMAN_GTF} ${HBV_REF_GTF} > ${COMBINED_GTF}
```

### 1.5 参考基因组索引构建

使用 `STAR` (推荐) 或 `HISAT2` 等比对工具为合并后的参考基因组构建索引。

```bash
# 安装 STAR (如果尚未安装)
# conda install -c bioconda star

# STAR 索引构建
GENOME_DIR="star_index_human_hbv"
mkdir -p ${GENOME_DIR}

STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir ${GENOME_DIR} \
     --genomeFastaFiles ${COMBINED_GENOME_FA} \
     --sjdbGTFfile ${COMBINED_GTF} \
     --sjdbOverhang 99 # 通常设置为 ReadLength - 1，例如读长100bp则设为99
```

`--sjdbOverhang` 参数应根据您的测序读长（Read Length）进行调整。

### 1.6 序列比对

使用 `STAR` 将清洗后的RNA-seq reads比对到构建好的联合参考基因组索引。

```bash
# STAR 比对 (示例为单个样本)
OUTPUT_DIR_SAMPLE1="star_output_sample1" # 为每个样本创建独立的输出目录
mkdir -p ${OUTPUT_DIR_SAMPLE1}

STAR --runThreadN 8 \
     --genomeDir ${GENOME_DIR} \
     --readFilesIn sample1_R1_paired_trimmed.fastq.gz sample1_R2_paired_trimmed.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix ${OUTPUT_DIR_SAMPLE1}/sample1_ \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --quantMode GeneCounts TranscriptomeSAM # 可选，STAR可直接输出基因计数

# 对所有样本重复此比对操作。
# 主要输出文件: sample1_Aligned.sortedByCoord.out.bam (用于后续定量和可视化)
# 如果使用 --quantMode GeneCounts, 还会生成 sample1_ReadsPerGene.out.tab
```

`--quantMode GeneCounts` 选项会让STAR直接输出一个基因计数文件，可以作为后续`featureCounts`步骤的替代或补充。

### 1.7 基因/转录本表达定量

如果未使用 STAR 的 `--quantMode GeneCounts`，或者需要更灵活的定量选项（例如，在exon水平定量），可以使用 [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) (Subread包的一部分)。

```bash
# 安装 Subread 包 (包含 featureCounts)
# conda install -c bioconda subread

# featureCounts 示例 (通常一次处理所有BAM文件以生成统一的计数矩阵)
# 假设您的BAM文件位于各自的STAR输出目录中
# 列出所有BAM文件
BAM_FILES=$(ls star_output_*/sample*_Aligned.sortedByCoord.out.bam)

featureCounts -T 8 -p \
    -a ${COMBINED_GTF} \
    -o all_samples_featureCounts.txt \
    ${BAM_FILES}

# 参数说明:
# -T: 线程数
# -p: 输入数据是paired-end。如果单端测序，则不使用此参数。
# -a: 指定GTF注释文件。
# -o: 输出的计数文件。

# all_samples_featureCounts.txt 将包含所有样本的基因水平read计数，
# 是下游差异表达分析的主要输入。
```

## 2\. 下游分析 (R 语言)

下游分析主要在R环境中进行，利用Bioconductor等包进行差异表达、可变剪接等分析。

```r
# R 环境准备: 安装必要的包
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install(c("DESeq2", "pheatmap", "ggplot2", "EnhancedVolcano", "dplyr"))
#
# # 对于通路富集分析 (可选)
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db") # 人类注释数据库

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)
```

### 2.1 差异表达基因 (DGE) 分析

使用 `DESeq2` 包进行差异表达基因分析。

```r
# 1. 读取计数矩阵
# all_samples_featureCounts.txt 是 featureCounts 的输出
countData_raw <- read.table("all_samples_featureCounts.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE, comment.char="#")

# 清理 featureCounts 输出:
# featureCounts的输出通常第1-5列是注释(Chr,Start,End,Strand,Length)，之后是样本的counts
# 保留基因ID (行名) 和实际的count列
countData <- countData_raw[, (colnames(countData_raw) %in% c("Chr","Start","End","Strand","Length")) == FALSE]

# 清理样本名 (如果需要)
# colnames(countData) <- gsub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(countData))
# colnames(countData) <- gsub("^star_output_sample\\d+/", "", colnames(countData)) # 示例清理

# 2. 准备样本信息表 (colData)
# colData 是一个DataFrame，行名与countData的列名(样本名)匹配，并包含实验条件等信息。
# 请根据您的实验设计创建此文件。
# 示例：假设有两组条件 "control" 和 "treatment"，各有3个重复样本
# 样本名应与countData中的列名完全一致且顺序对应
sample_names <- colnames(countData) # 确保与countData列名一致
# 示例条件，请务必根据您的实际情况修改
conditions <- factor(c(rep("control", 3), rep("treatment", 3)))
colData <- data.frame(row.names=sample_names, condition=conditions)

# 确保 countData 和 colData 的样本名和顺序完全一致
stopifnot(all(rownames(colData) == colnames(countData)))

# 3. 创建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition) # design formula 根据您的实验设计

# 预过滤 (去除在所有样本中表达量都极低的基因，可选但推荐)
keep <- rowSums(counts(dds)) >= 10 # 例如，至少在所有样本中总共有10个reads
dds <- dds[keep,]

# 4. 运行 DESeq2 分析
dds <- DESeq(dds)

# 5. 获取结果
# contrast 参数指定比较的组别，格式为 c("factorName", "level1_treatment", "level2_control")
res <- results(dds, contrast=c("condition", "treatment", "control"))
summary(res)

# 按调整后的p-value (padj) 排序
resOrdered <- res[order(res$padj),]

# 提取显著差异表达的基因 (例如: padj < 0.05 和 |log2FoldChange| > 1)
significant_res <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(significant_res), file="DGE_significant_results.csv")
print(paste("Found", nrow(significant_res), "significant DEGs."))

# 6. 可视化
# Volcano plot
volcano_plot <- EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj', # 使用调整后的p值
    title = 'Treatment vs Control',
    pCutoff = 0.05,
    FCcutoff = 1.0, # log2FC > 1 or < -1
    pointSize = 2.0,
    labSize = 3.0,
    colAlpha = 0.5,
    legendPosition = 'right',
    drawConnectors = TRUE,
    widthConnectors = 0.2,
    colConnectors = 'grey30')
print(volcano_plot)
ggsave("volcano_plot_DEGs.png", plot=volcano_plot, width=10, height=8)

# Heatmap of top N significant genes (e.g., top 50)
# 需要先对数据进行转换，例如 rlog (regularized log) 或 vst (variance stabilizing transformation)
vsd <- vst(dds, blind=FALSE) # blind=FALSE因为我们已经进行了差异分析
topN_genes <- head(rownames(significant_res), 50) # 取padj最小的前50个基因
mat <- assay(vsd)[topN_genes,]
mat <- mat - rowMeans(mat) # Z-score scaling per gene (centering)

heatmap_plot <- pheatmap(mat, annotation_col=colData, scale="row",
                         main="Top 50 DEGs Heatmap",
                         fontsize_row=8, fontsize_col=10)
# ggsave("heatmap_top50_DEGs.png", plot=heatmap_plot, width=8, height=10) # pheatmap直接保存或在RStudio中导出

# MA plot
plotMA(res, ylim=c(-3,3), main="MA Plot")
```

### 2.2 可变剪接分析

使用 `rMATS` (replicate Multivariate Analysis of Transcript Splicing) 进行可变剪接分析。`rMATS` 直接在比对产生的BAM文件上运行，并需要GTF注释文件。

#### 2.2.1 运行 rMATS (命令行)

`rMATS` (特别是 `rMATS.turbo`) 通常在Linux命令行环境运行。

```bash
# 准备输入文件列表 (BAM文件路径)
# b1.txt: 包含第一组条件 (e.g., control) 的BAM文件路径列表，逗号分隔，无空格
# echo "star_output_control1/control1_Aligned.sortedByCoord.out.bam,star_output_control2/control2_Aligned.sortedByCoord.out.bam" > b1.txt
# b2.txt: 包含第二组条件 (e.g., treatment) 的BAM文件路径列表，逗号分隔，无空格
# echo "star_output_treatment1/treatment1_Aligned.sortedByCoord.out.bam,star_output_treatment2/treatment2_Aligned.sortedByCoord.out.bam" > b2.txt

# 运行 rMATS (示例使用 rmats.py, 即rMATS v4.x.x, Python 3兼容)
# conda install -c bioconda rmats # 如果通过conda安装
RMATS_PY_PATH="/path/to/your/rmats.py" # 或者直接用 rmats.py 如果在PATH中
RMATS_OUTPUT_DIR="rmats_output_human_hbv"
mkdir -p ${RMATS_OUTPUT_DIR}

python ${RMATS_PY_PATH} \
    --b1 b1.txt \
    --b2 b2.txt \
    --gtf ${COMBINED_GTF} \
    -t paired \
    --readLength 99 \
    --nthread 8 \
    --od ${RMATS_OUTPUT_DIR} \
    --tmp ${RMATS_OUTPUT_DIR}/temp \
    --statoff # 如果每组样本数少于3个，可能需要此选项，否则通常不需要

# 参数说明:
# --gtf: 使用我们合并的 ${COMBINED_GTF}
# --readLength: 您的测序读长 (通常是 --sjdbOverhang + 1)
# --statoff: 当每组重复数小于3时，rMATS的统计检验可能无法进行，此参数关闭统计。
# 详细参数请查阅rMATS官方文档。
```

#### 2.2.2 rMATS 输出解读 (R环境)

`rMATS` 会为每种主要的可变剪接事件类型（SE: Skipped Exon, RI: Retained Intron, A5SS: Alternative 5' Splice Site, A3SS: Alternative 3' Splice Site, MXE: Mutually Exclusive Exons）生成一个 `MATS.JCEC.txt` (Junction Counts and Exon Counts) 文件。

```r
# 读取rMATS输出，例如SE (Skipped Exon) 事件
# 路径应指向您的rMATS输出目录
rmats_output_folder <- "rmats_output_human_hbv"
se_results_file <- file.path(rmats_output_folder, "SE.MATS.JCEC.txt")

if (file.exists(se_results_file)) {
  se_results <- read.table(se_results_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

  # 筛选显著差异剪接事件
  # IncLevelDifference: PSI (Percent Spliced In) 或 PS (Percent Spliced) 值的差异
  # 阈值可根据研究需要调整
  significant_splicing_se <- subset(se_results, FDR < 0.05 & abs(IncLevelDifference) >= 0.1)
  
  if (nrow(significant_splicing_se) > 0) {
    print(paste("Found", nrow(significant_splicing_se), "significant SE events."))
    write.csv(significant_splicing_se, "significant_SE_events.csv", row.names=FALSE)
  } else {
    print("No significant SE events found with the current criteria.")
  }
} else {
  print(paste("rMATS SE results file not found:", se_results_file))
}

# 对其他剪接类型 (RI, A5SS, A3SS, MXE) 重复类似操作。
```

#### 2.2.3 HBV剪接事件坐标还原

您提供的 `HBV_B.fasta` 是一个线性化的参考序列，其坐标系统与原始的环状HBV基因组不同。如果`rMATS`或其他剪接分析工具在 `HBV_B.fasta` (染色体名可能为 `HBV_B` 或您在FASTA中定义的名字) 上检测到剪接事件，其报告的坐标（如外显子起始终止位点）是基于这个线性参考的。为了更好地理解这些事件在原始HBV基因组上的位置，需要进行坐标还原。

**HBV线性参考结构回顾:**

  * 原始HBV基因型B全长 ($L\_{orig}$) = 3215 bp。
  * 线性参考 `HBV_B.fasta` 构成:
      * 第一段 (Part A): 原始HBV的 1600bp - 3215bp。在线性参考上为 1bp - 1616bp ($L\_A$ = 1616bp)。
      * 第二段 (Part B): 原始HBV的 1bp - 1950bp。在线性参考上为 1617bp - 3566bp ($L\_B$ = 1950bp)。

**坐标转换逻辑 (1-based):**
设 `coord_linear` 是在 `HBV_B.fasta` 上的坐标。$L\_A = 1616$ (第一段的长度)

  * 如果 $1 \\le \\text{coord\_linear} \\le L\_A$ (即 $\\text{coord\_linear} \\le 1616$):
    $\\text{coord\_original\_circular} = \\text{coord\_linear} + 1600 - 1$
  * 如果 $L\_A + 1 \\le \\text{coord\_linear} \\le L\_A + L\_B$ (即 $1617 \\le \\text{coord\_linear} \\le 3566$):
    $\\text{coord\_original\_circular} = \\text{coord\_linear} - L\_A$ (即 $\\text{coord\_linear} - 1616$)

**应用于rMATS结果:**
对于rMATS输出的 `*.MATS.JCEC.txt` 文件中涉及HBV染色体的坐标列 (如 `exonStart_0base`, `exonEnd`, `upstreamES`, `upstreamEE`, `downstreamES`, `downstreamEE` 等)，如果这些事件发生在HBV上，您需要应用上述转换逻辑。注意rMATS的 `*_0base` 坐标是0-based，转换为1-based后再应用上述逻辑，或调整逻辑以适应0-based。

例如，在R中处理一个包含HBV剪接事件的数据框 `hbv_splicing_events_df`，其中有一列 `linear_exon_start` (假设已转为1-based):

```r
# 假设 hbv_splicing_events_df 包含来自rMATS的、发生在HBV上的事件
# 并且 linear_exon_start 是1-based的线性坐标
L_A <- 1616 # Length of Part A in the linear HBV reference

hbv_splicing_events_df <- hbv_splicing_events_df %>%
  mutate(
    original_exon_start = ifelse(linear_exon_start <= L_A,
                                 linear_exon_start + 1600 - 1,
                                 linear_exon_start - L_A)
  )
# 对其他相关坐标列 (exonEnd等) 应用相同逻辑。
```

这种坐标还原对于解释HBV内部的剪接模式至关重要。

### 2.3 HBV RNA 分析

#### 2.3.1 HBV pgRNA 表达分析

HBV pgRNA (pregenomic RNA) 的表达量可以直接从 `DESeq2` 分析产生的归一化计数或差异表达结果 (`res` 对象) 中提取。您需要知道在您提供的 `HBV_B.gtf` 中，pgRNA的基因ID是什么。

```r
# 假设 dds 对象已经创建并运行了 DESeq
# 获取归一化计数
normalized_counts <- counts(dds, normalized=TRUE)

# 替换为您的HBV pgRNA在GTF中的实际gene_id
hbv_pgRNA_gene_id <- "HBV_pgRNA_gene" # 示例ID，请确认

if (hbv_pgRNA_gene_id %in% rownames(normalized_counts)) {
  hbv_pgRNA_counts_matrix <- normalized_counts[hbv_pgRNA_gene_id, , drop=FALSE] # 保留矩阵形式
  
  # 转换为长格式数据框用于ggplot
  hbv_pgRNA_df_long <- data.frame(
    Sample = colnames(hbv_pgRNA_counts_matrix),
    NormalizedCounts = as.numeric(hbv_pgRNA_counts_matrix[1,]),
    Condition = colData(dds)$condition # 从colData获取条件信息
  )

  print("Normalized counts for HBV pgRNA:")
  print(hbv_pgRNA_df_long)

  # 绘制箱线图比较不同条件下的表达
  pgRNA_plot <- ggplot(hbv_pgRNA_df_long, aes(x=Condition, y=NormalizedCounts, fill=Condition)) +
    geom_boxplot(outlier.shape = NA) + # 不显示离群点，让jitter显示
    geom_jitter(shape=16, position=position_jitter(width=0.2, height=0)) +
    labs(title="HBV pgRNA Expression",
         x="Condition",
         y="Normalized Counts") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  print(pgRNA_plot)
  ggsave("HBV_pgRNA_expression_boxplot.png", plot=pgRNA_plot, width=6, height=5)

} else {
  print(paste("HBV pgRNA gene ID '", hbv_pgRNA_gene_id, "' not found in count matrix.", sep=""))
}

# 检查HBV pgRNA是否在差异表达结果中
if (exists("res") && hbv_pgRNA_gene_id %in% rownames(res)) {
  hbv_pgRNA_de_results <- res[hbv_pgRNA_gene_id,]
  print("Differential expression results for HBV pgRNA:")
  print(hbv_pgRNA_de_results)
}
```

#### 2.3.2 HBV转录本覆盖度可视化及坐标还原

为了解HBV基因组上各个位置的reads覆盖情况，可以生成覆盖度折线图。这通常涉及从BAM文件中提取每个碱基的覆盖深度，然后将这些深度值绘制出来。关键在于将绘图时的x轴坐标还原为原始的环状HBV基因组坐标。

**1. 生成每个碱基的覆盖深度 (命令行)**

首先，从比对产生的BAM文件中提取HBV染色体上的reads，然后计算每个碱基的覆盖深度。

```bash
# 假设您的HBV染色体在 .fasta 和 .gtf 中被称为 "HBV_B" (请替换为实际名称)
HBV_CHROM_NAME="HBV_B" # 替换为您的HBV染色体名

# 过滤得到只包含HBV染色体上比对的BAM文件 (可选，但推荐)
# samtools idxstats sample1_Aligned.sortedByCoord.out.bam # 查看染色体名称
samtools view -b star_output_sample1/sample1_Aligned.sortedByCoord.out.bam ${HBV_CHROM_NAME} > sample1_hbv_reads.bam
samtools index sample1_hbv_reads.bam

# 计算每个碱基的覆盖深度
# -a: 输出所有位置 (包括深度为0的)
# -d 0: 无最大深度限制
samtools depth -a -d 0 sample1_hbv_reads.bam > sample1_hbv_coverage.txt
# 对每个样本重复此操作
```

`sample1_hbv_coverage.txt` 文件将包含三列：染色体名、位置 (1-based)、覆盖深度。

**2. 在R中读取、转换坐标并绘图**

```r
# 假设为 sample1 处理
coverage_file <- "sample1_hbv_coverage.txt"
if (file.exists(coverage_file)) {
  hbv_cov_data <- read.table(coverage_file, header=FALSE, col.names=c("chromosome", "linear_pos", "depth"))

  # 坐标转换参数
  L_A <- 1616 # 线性参考中第一段 (源自原始1600-3215) 的长度
  L_orig <- 3215 # 原始HBV基因组长度

  # 应用坐标转换逻辑
  hbv_cov_data <- hbv_cov_data %>%
    mutate(
      original_pos = ifelse(linear_pos <= L_A,
                            linear_pos + 1600 - 1,
                            linear_pos - L_A)
    ) %>%
    filter(original_pos <= L_orig) # 确保坐标在原始基因组范围内

  # 绘制覆盖度折线图
  coverage_line_plot <- ggplot(hbv_cov_data, aes(x=original_pos, y=depth)) +
    geom_line() +
    labs(title=paste("HBV Coverage Plot for Sample1 (Original Coordinates)"),
         x="Position on Original HBV Genome (bp)",
         y="Read Depth") +
    theme_bw() +
    # 可以考虑添加垂直线标出重要区域，如pgRNA的起始终止
    # geom_vline(xintercept = 1701, linetype="dashed", color="blue") + # pgRNA start (original)
    # geom_vline(xintercept = 1930, linetype="dashed", color="red")   # pgRNA end (original)
    scale_x_continuous(breaks = seq(0, L_orig, by = 500)) # 调整x轴刻度
  
  print(coverage_line_plot)
  ggsave("sample1_hbv_coverage_original_coords.png", plot=coverage_line_plot, width=12, height=6)

} else {
  print(paste("Coverage file not found:", coverage_file))
}

# 您可能希望将多个样本的覆盖度绘制在同一张图上或分面显示，
# 这需要合并来自多个样本的 hbv_cov_data 数据框，并添加样本标识列。
```

## 3\. 结论与进一步分析

本流程提供了一个从原始RNA-seq数据到差异表达、可变剪接分析以及HBV整合分析的框架。

可能的进一步分析方向：

  * **通路富集分析 (Pathway Enrichment Analysis):** 对鉴定出的差异表达基因进行GO (Gene Ontology) 和 KEGG 通路富集分析，以了解其生物学功能。R包 `clusterProfiler` 是一个常用工具。
  * **HBV-宿主融合转录本检测:** 如果怀疑存在HBV与宿主基因组的融合转录本，可以使用专门的工具 (如 `STAR-Fusion`, `EricScript`) 进行检测。这通常需要特定的比对参数和后处理流程。
  * **高级可视化:** 使用 `IGV` (Integrative Genomics Viewer) 查看特定基因的reads覆盖情况、剪接事件以及HBV区域的比对细节。IGV可以直接加载BAM文件和GTF文件。

请记住，所有提供的代码和参数都需要根据您的具体数据、文件路径和计算环境进行调整。这是一个指导性的流程，实际操作中可能需要根据遇到的问题进行调试和优化。祝您分析顺利！

```
```
