#RNAseq_ldg integrated
##RNAseq_human 
#获得fastq原始数据
#仅仅扫描FASTQ文件一次，就完成比FASTQC+ cutadapt + Trimmomatic 这三个软件加起来的功能还多很多的功能，
#而且速度上比仅仅使用Trimmomatic一个软件还要快3倍左右
#进行质控：进行fastp
mkdir fastp_result
fastp -i DXY_treat_S161_R1_001.fastq.gz      -o ./fastp_result/DXY_treat_S161_R1_001_clean.fastq.gz\
      -I DXY_treat_S161_R2_001.fastq.gz     -O ./fastp_result/DXY_treat_S161_R2_001_clean.fastq.gz


#使用提前构建好的index salmon_index
#使用的是以下版本：wget -c https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz 

#Salmon 定量
#创建sbatch文档
vim CZQ-CD235-H.sh

#!/bin/bash
#SBATCH --job-name=salmon_quant   # 作业名
#SBATCH --partition=cpu   # cpu 队列
#SBATCH -n 80 # 总核数
#SBATCH --ntasks-per-node=40 # 每节点核数
#SBATCH --output=%j_salmon_quant.out #标准输出文件
#SBATCH --error=%%j_salmon_quant.err
cd /lustre/home/acct-medty/medty-a/LDG/thalidomide_thalassemia/CD235/QD0_HD85_pair/thalidomide/rawdata/nonresponse/fastp_result
source activate salmon
salmon quant -i ~/LDG/human_transciptomeIndex_salmonBased/salmon_index -l A \
         -1 CXQ-CD235-H_S47_R1_001_clean.fastq.gz \
         -2 CXQ-CD235-H_S47_R2_001_clean.fastq.gz\
         -p 8 --validateMappings --gcBias --numGibbsSample 20 \
         -o quants/CXQ-CD235-H
#Q为D0，H为D85

#提交任务
sbatch CZQ-CD235-H.sh

#整理Salmon定量文件用于DESeq2差异基因鉴定：
#教程见 https://blog.csdn.net/qazplm12_3/article/details/111056012
#找到Salmon的输出文件并压缩起来,用于下载到本地进行差异分析。
# 列出salmon的输出文件
find . -name quant.sf
# 这个压缩包下载到到本地
zip quant.sf.zip `find . -name quant.sf`

#已经构建好转录本和基因对应的关系文件：该版本Homo_sapiens.GRCh38.cdna.all.fa.gz 
GRCh38.tx2gene 文件如下：.后为版本号
Tv     Gv
ENST00000624431.2 ENSG00000279928.2
ENST00000424215.1 ENSG00000228037.1
ENST00000511072.5 ENSG00000142611.17
ENST00000607632.1 ENSG00000142611.17

#至此就完成了基于Salmon的所有样本基因和转录本的定量。下载quant.sf.zip文件到本地进行下游分析。


#下游：R语言操作:DEseq2配对T检验
rm(list=ls())
#BiocManager::install("tximport")
library("tximport")
library(DESeq2)
setwd("~/LDG/RNAseq_ldg/RNAseq_ldg/CD235_thalidomide/")
txTogene <- read.table("~/LDG/RNAseq_ldg/GRCh38_tx2gene2/GRCh38.tx2gene.csv",sep = ",",header = T)
sampleList <- c("CZQ-CD235-H","CZQ-CD235-Q","DXY-CD235-H", "DXY-CD235-Q","LFF-CD235-H","LFF-CD235-Q","LHJ-CD235-H",
                "LHJ-CD235-Q","LHY-CD235-H","LHY-CD235-Q","LMX-CD235-H","LMX-CD235-Q","LZH-CD235-H","LZH-CD235-Q",
                "WYT-CD235-H","WYT-CD235-Q","ZWW-CD235-H","ZWW-CD235-Q")
fileList <- file.path("quants", sampleList, "quant.sf")
names(fileList) <- sampleList
fileList
file.exists(fileList)
txi <- tximport(fileList, type = "salmon",tx2gene = txTogene)

head(txi$counts) # 查看矩阵
head(txi$abundance)
CD235_thalidomide_abundance_ENSMEBL_tximport<- txi$abundance
##################
#ID转换,输出symbol的矩阵
library(dplyr)
CD235 <- as.data.frame(CD235_thalidomide_abundance_ENSMEBL_tximport)
CD235 <- tibble::rownames_to_column(CD235,"ENSEMBL")
#colnames(Leucine)[1] <- "ENSEMBL"
library(stringi)##加载包
CD235_2 <- CD235
CD235_2$ENSEMBL=stri_sub(CD235_2$ENSEMBL,1,15)##保留前15位
#加载相关包
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("igraph")
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
# 查看org.Hs.eg.db 包提供的转换类型
keytypes(org.Hs.eg.db)
# 需要转换的Ensembl_ID
Ensembl_ID <-CD235_2$ENSEMBL
# 采用bitr()函数进行转换
gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
# 查看转换的结果
head(gene_symbol)
merge_matrix <- merge(gene_symbol,CD235_2,by="ENSEMBL")
merge_matrix <- merge_matrix[,-1]
write.table(merge_matrix,"CD235_thalidomide_abundance_symbol_matrix.csv", sep = ",", row.names = TRUE)

########################
#配对信息,构建DEseq2对象
colData <- read.csv("colData.csv")
#配对信息格式如下：
               condition      pair
Patient1_Pre         Pre   Patient1
Patient1_Post       Post   Patient1
Patient2_Pre         Pre   Patient2
Patient2_Post       Post   Patient2
Patient3_Pre         Pre   Patient3
Patient3_Post       Post   Patient3
Patient4_Pre         Pre   Patient4
Patient4_Post       Post   Patient4

# 将 colData 中的字符型变量转换为因子
colData$condition <- as.factor(colData$condition)
colData$pair <- as.factor(colData$pair)
str(colData)
dds <- DESeqDataSetFromTximport(txi, colData = colData, design = ~ pair + condition)


#进行差异分析
dds<-DESeq(dds)  #将数据进行标准化，必要步骤
res<-results(dds)        ##从DESeq 分析中提取结果表
resordered<-res[order(res$padj),]        ##以padj对res排序
summary(res)          #看一下结果的概要信息，p值默认小于0.1。
write.csv(resordered,'res.csv',row.names=TRUE) # 保存全部结果

##################
#ID转换
library(dplyr)
thalidomide <- read.csv("result/res.csv")
#tha <- tibble::rownames_to_column(D8_HUDEP2_thalidomide_rawcount2,"ENSEMBL")
colnames(thalidomide)[1] <- "ENSEMBL"
library(stringi)##加载包
thalidomide2 <- thalidomide
thalidomide2$ENSEMBL=stri_sub(thalidomide2$ENSEMBL,1,15)##保留前15位
#加载相关包
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)
# 查看org.Hs.eg.db 包提供的转换类型
keytypes(org.Hs.eg.db)
# 需要转换的Ensembl_ID
Ensembl_ID <-thalidomide2$ENSEMBL
# 采用bitr()函数进行转换
gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
# 查看转换的结果
head(gene_symbol)
merge <- merge(gene_symbol,thalidomide2,by="ENSEMBL")
merge <- merge[,-1]
resordered_merge<-merge[order(merge$log2FoldChange,decreasing = TRUE),]##以Basemean对res排序
write.table(resordered_merge,"CD235_thalidomdie_DEseq2_output_SYMBOL_all.csv", sep = ",")

#########################循环跑RNAseq上游###########################3
#获得fastq原始数据
#仅仅扫描FASTQ文件一次，就完成比FASTQC+ cutadapt + Trimmomatic 这三个软件加起来的功能还多很多的功能，
#而且速度上比仅仅使用Trimmomatic一个软件还要快3倍左右
#进行质控：进行fastp
#更改样本名，方便循环：
C1_LiuRan-D12_S1_L001_R1_001.fastq.gz  LiuRan-D12_S1_L001_1_R1.fastq.gz
C1_LiuRan-D12_S1_L001_R2_001.fastq.gz  LiuRan-D12_S1_L001_1_R2.fastq.gz
C2_LiuRan-D12_S1_L001_R1_001.fastq.gz  LiuRan-D12_S1_L001_2_R1.fastq.gz
C2_LiuRan-D12_S1_L001_R2_001.fastq.gz  LiuRan-D12_S1_L001_2_R2.fastq.gz
C3_LiuRan-D12_S1_L001_R1_001.fastq.gz  LiuRan-D12_S1_L001_3_R1.fastq.gz
C3_LiuRan-D12_S1_L001_R2_001.fastq.gz  LiuRan-D12_S1_L001_3_R2.fastq.gz
T1_LiuRan-D12_S1_L001_R1_001.fastq.gz  LiuRan-D12_S1_L001_4_R1.fastq.gz
T1_LiuRan-D12_S1_L001_R2_001.fastq.gz  LiuRan-D12_S1_L001_4_R2.fastq.gz
T2_LiuRan-D12_S1_L001_R1_001.fastq.gz  LiuRan-D12_S1_L001_5_R1.fastq.gz
T2_LiuRan-D12_S1_L001_R2_001.fastq.gz  LiuRan-D12_S1_L001_5_R2.fastq.gz
T3_LiuRan-D12_S1_L001_R1_001.fastq.gz  LiuRan-D12_S1_L001_6_R1.fastq.gz
T3_LiuRan-D12_S1_L001_R2_001.fastq.gz  LiuRan-D12_S1_L001_6_R2.fastq.gz

#创建系统命令
vim  fastp.sh 
#!/bin/bash
#SBATCH --job-name=fastp  # 作业名
#SBATCH --partition=cpu   # cpu 队列
#SBATCH -n 80 # 总核数
#SBATCH --ntasks-per-node=40 # 每节点核数
#SBATCH --output=%j_fastp.out #标准输出文件
#SBATCH --error=%j_fastp.err

# 进入工作目录
cd /lustre/home/acct-medty/medty-a/LDG/Thalassemia/thalassemia_RNAseq/Differentiation_in_vitro

# 激活 salmon 环境
source activate salmon

# 使用 for 循环处理从 1 到 6 的文件
for fn in {1..6}; do #fn 一般只设置变化部分
    samp="LiuRan-D12_S1_L001_${fn}"  # 构造样本名称，把可变部分和不变部分组合到一起
    echo "Processing sample ${samp}"
    # 运行 fastp 进行质控
    fastp -i ${samp}_R1.fastq.gz \
          -o ./fastp/${samp}_R1_clean.fastq.gz \
          -I ${samp}_R2.fastq.gz \
          -O ./fastp/${samp}_R2_clean.fastq.gz
done


#提交系统命令
sbatch fastp.sh 


#使用提前构建好的index salmon_index
#使用的是以下版本：wget -c https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz 

#Salmon 定量
#创建sbatch文档
vim salmon_quant.sh

#!/bin/bash
#SBATCH --job-name=salmon_quant   # 作业名
#SBATCH --partition=cpu   # cpu 队列
#SBATCH -n 80 # 总核数
#SBATCH --ntasks-per-node=40 # 每节点核数
#SBATCH --output=%j_salmon_quant.out #标准输出文件
#SBATCH --error=%j_salmon_quant.err
cd /lustre/home/acct-medty/medty-a/LDG/Thalassemia/thalassemia_RNAseq/Differentiation_in_vitro/fastp
source activate salmon
# 创建输出目录
mkdir -p quants
# 使用 for 循环处理从 1 到 6 的文件
for fn in {1..6}; do #fn 一般只设置变化部分
    samp="LiuRan-D12_S1_L001_${fn}"  # 构造样本名称,把可变部分和不变部分组合到一起
    echo "Processing sample ${samp}"
salmon quant -i ~/LDG/human_transciptomeIndex_salmonBased/salmon_index -l A \
         -1  ${samp}_R1_clean.fastq.gz\
         -2  ${samp}_R2_clean.fastq.gz\
         -p 8 --validateMappings --gcBias --numGibbsSample 20 \
         -o quants/${samp}
done

#提交任务
sbatch salmon_quant.sh

#整理Salmon定量文件用于DESeq2差异基因鉴定：
#教程见 https://blog.csdn.net/qazplm12_3/article/details/111056012
#找到Salmon的输出文件并压缩起来,用于下载到本地进行差异分析。
# 列出salmon的输出文件
find . -name quant.sf
# 这个压缩包下载到到本地
zip quant.sf.zip `find . -name quant.sf`

#已经构建好转录本和基因对应的关系文件：该版本Homo_sapiens.GRCh38.cdna.all.fa.gz 
GRCh38.tx2gene 文件如下：.后为版本号
Tv     Gv
ENST00000624431.2 ENSG00000279928.2
ENST00000424215.1 ENSG00000228037.1
ENST00000511072.5 ENSG00000142611.17
ENST00000607632.1 ENSG00000142611.17

#该对应关系已经保存在路径/lustre/home/acct-medty/medty-a/LDG/RNAseq_ldg/GRCh38_tx2gene2/GRCh38.tx2gene.csv，可直接使用：
#至此就完成了基于Salmon的所有样本基因和转录本的定量。下载quant.sf.zip文件到本地进行下游分析。


#下游：R语言操作，成组DEseq2
rm(list=ls())
#BiocManager::install("tximport")
library("tximport")
library(DESeq2)
setwd("~/LDG/RNAseq_ldg/RNAseq_ldg/CD235_thalidomide/data")
txTogene <- read.table("~/LDG/RNAseq_ldg/GRCh38_tx2gene2/GRCh38.tx2gene.csv",sep = ",",header = T)
sampleList <- c("LiuRan-D12_S1_L001_1","LiuRan-D12_S1_L001_2","LiuRan-D12_S1_L001_3", "LiuRan-D12_S1_L001_4",
                "LiuRan-D12_S1_L001_5","LiuRan-D12_S1_L001_6")
fileList <- file.path("quants", sampleList, "quant.sf")
names(fileList) <- sampleList
fileList
file.exists(fileList)
txi <- tximport(fileList, type = "salmon",tx2gene = txTogene)

head(txi$counts) # 查看矩阵
head(txi$abundance)# 查看标准化的TPM，方便导出直接比较不同样本之间的基因表达水平
CD235_thalidomide_abundance_ENSMEBL_tximport<- txi$abundance
write.table(CD235_thalidomide_abundance_ENSMEBL_tximport,"CD235_thalidomide_abundance_ENSMEBL_tximport.csv", sep = ",", row.names = TRUE)

#########################################################
#官方教程：https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts
#构建DEseq2分析的对象
library(DESeq2)
sampleNames<-c("LiuRan-D12_S1_L001_1","LiuRan-D12_S1_L001_2","LiuRan-D12_S1_L001_3", "LiuRan-D12_S1_L001_4","LiuRan-D12_S1_L001_5","LiuRan-D12_S1_L001_6")
sampleGroup<-factor(rep(c("DMSO","Thal"),each=3),level=c("DMSO","Thal"))     ##差异比较矩阵。R语言中，factor 的level默认按照字母排序排在最前面的对照组，手动设置level.
sampleTable<-data.frame(sampleName=sampleNames,type=sampleGroup)     ##分组矩阵
rownames(sampleTable)<-colnames(txi$counts)
dds<-DESeqDataSetFromTximport(txi,sampleTable,design=~type)       #构建dds矩阵

#利用层次聚类的方法评估一下数据质量（vst标准化---hclust)
#BiocManager::install('factoextra') 
library('factoextra')          #下载并加载factoextra包
vsd <- vst(dds,blind = TRUE)    #参数blind=TRUE是为了不把样本分组信息考虑在内——即以无偏的方式进行转换
sampleDists <- dist(t(assay(vsd))) 
hclust <- hcut(sampleDists, k = 2, stand = FALSE,hc_method ="average" ) 
g <- fviz_dend(hclust,
          rect_fill = T,
          # 字体大小
          cex = 1,
          # 字体颜色
          color_labels_by_k=T,
          # 平行放置
          horiz=T)
g
ggsave("vst标准化_hclust.pdf", g, width = 8, height = 5)

#利用PCA主成分分析的方法评估一下数据质量。
rld <- vst(dds, blind=FALSE)     #vst()函数效果和rlog（）一样，且速度更快。
g <- plotPCA(rld, intgroup="type",ntop=500)
g
ggsave("PCA.pdf", g, width = 8, height = 5)



#进行差异分析
dds<-DESeq(dds)  #将数据进行标准化，必要步骤
res<-results(dds)        ##从DESeq 分析中提取结果表
resordered<-res[order(res$padj),]        ##以padj对res排序
summary(res)          #看一下结果的概要信息，p值默认小于0.1。
write.csv(resordered,'LiuRan_CD34_Differentiation_in_vitro_thalidomide_ENSEMBL_DEseq2_output.csv',row.names=TRUE) # 保存全部结果s

###以下可做？？？？去重复基因

#按照BaseMean排序，保留表达量最大的第一个基因，其他重复基因删除
resordered_merge<-merge[order(merge$baseMean,na.last = TRUE,decreasing = TRUE),]        ##以Basemean对res排序
#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
keep=!duplicated(resordered_merge$SYMBOL)
#得到最后处理之后的表达谱矩阵
expr_max=resordered_merge[keep,]
expr_max
expr_max2  <- na.omit(expr_max[abs(expr_max$log2FoldChange) >1 & expr_max$padj <0.05,])
dim(expr_max2)

write.table(expr_max,"Nalm6_LD-leucine_DEseq2_output_SYMBOL_abspj0.05_log2FC1.csv", sep = ",", row.names = TRUE)
