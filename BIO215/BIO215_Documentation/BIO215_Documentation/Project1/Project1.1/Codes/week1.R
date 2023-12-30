##############################################################################
########             RNA-Seq Analysis for Breast Cancer                 ######
##############################################################################
#查看文件结构
system("tree ~/redo/BIO215")

#创建一个总结表，并且在data目录下导出为csv
coldata <- data.frame(srr_id = paste0("SRR2228734", 1:8),
                      library = rep(c("MeRIP-Seq IP","RNA-Seq"), times = 4),
                      tissue = rep(c("Breast cancer tissue","Breast tissue"), each = 4)
                     )

write.csv(coldata, "~/redo/BIO215/project1/data/colData.csv")

#查看8个FASTQ文件，并且软链接到本地data目录下
system("ls -lh /data/BIO215/*")
#只要RNA-Seq
system("ln -s /data/BIO215/SRR22287342.fastq /home/jiaming.huang21/redo/BIO215/project1/data/SRR22287342.fastq;
       ln -s /data/BIO215/SRR22287344.fastq /home/jiaming.huang21/redo/BIO215/project1/data/SRR22287344.fastq;
       ln -s /data/BIO215/SRR22287346.fastq /home/jiaming.huang21/redo/BIO215/project1/data/SRR22287346.fastq;
       ln -s /data/BIO215/SRR22287348.fastq /home/jiaming.huang21/redo/BIO215/project1/data/SRR22287348.fastq
       ")

#再次查看文件结构
system("tree ~/redo/BIO215")

#Step6: 使用genome aligner Hisat2，将FASTQ文件与hg38 map
system("hisat2 -x /data/genome_indx/genome_hg38 -U /home/jiaming.huang21/redo/BIO215/project1/data/SRR22287342.fastq -S /home/jiaming.huang21/redo/BIO215/project1/data/SRR22287342.sam;
       hisat2 -x /data/genome_indx/genome_hg38 -U /home/jiaming.huang21/redo/BIO215/project1/data/SRR22287344.fastq -S /home/jiaming.huang21/redo/BIO215/project1/data/SRR22287344.sam;
       hisat2 -x /data/genome_indx/genome_hg38 -U /home/jiaming.huang21/redo/BIO215/project1/data/SRR22287346.fastq -S /home/jiaming.huang21/redo/BIO215/project1/data/SRR22287346.sam;
       hisat2 -x /data/genome_indx/genome_hg38 -U /home/jiaming.huang21/redo/BIO215/project1/data/SRR22287348.fastq -S /home/jiaming.huang21/redo/BIO215/project1/data/SRR22287348.sam
       ")

#Step7: SAM转BAM
system("echo $HOME")
system("cd ~/redo/BIO215/project1/data")
system("pwd") #发现异常：工作目录不变，但在terminal直接运行会改变
system("samtools view -Sb SRR22287342.sam > SRR22287342.bam;
       samtools view -Sb SRR22287344.sam > SRR22287344.bam;
       samtools view -Sb SRR22287346.sam > SRR22287346.bam;
       samtools view -Sb SRR22287348.sam > SRR22287348.bam
       ")

#Step8: 利用BAM与注释信息中的基因信息进行SummarizeOverlaps，得到基因表达量

#获得基因注释
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
features <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

#重命名features的名字
library(org.Hs.eg.db)
names(features) <- AnnotationDbi::select(org.Hs.eg.db, names(features), "SYMBOL", "ENTREZID")$SYMBOL

#计算基因表达量
data_dir <- "~/redo/BIO215/project1/data"
bams <- BamFileList(file.path(data_dir,c("SRR22287342.bam",
                                         "SRR22287344.bam",
                                         "SRR22287346.bam",
                                         "SRR22287348.bam"
                                         )
                              )
                   ) #制作BamFileList对象：bam文件路径组合
count_se <- summarizeOverlaps(features, bams) #用summarizeOverlaps数bams在每个基因上的数量
head(assay(count_se)) #查看表达量
count_se <- count_se[rowSums(assay(count_se) >= 10) > 1,] #筛选显著基因
head(assay(count_se)) #查看过滤后的基因的表达量

gc() #清空内存

#########################差异基因表达分析和自由探索#############################

#方法一：
library(DESeq2)
count_se$tissue <- factor(c("Cancer", "Cancer", "Normal", "Normal"),
                          levels = c("Normal", "Cancer") #"Normal"对照， 
                          ) #colData给因子
dds <- DESeqDataSet(count_se, ~tissue) #准备DESeqDataSet对象：数据和标签
dds <- DESeq(dds) #DESeq()相当于做线性回归

res <- results(dds) #提取result
res[which(res$padj < 0.05),] #padj是adjusted p-value

DESeq2::plotDispEsts(dds)

DESeq2::plotMA(dds)

colnames(dds) <- c("Cancer-Rep1", "Cancer-Rep2", "Normal-Rep1", "Normal-Rep2")
DESeq2::plotPCA(rlog(dds), "tissue")

DESeq2::plotCounts(dds, "TP53", "tissue")
res["TP53",]

DESeq2::plotCounts(dds, "BRCA1", "tissue")
res["BRCA1",]

heatmap(assay(rlog(dds)))

pheatmap::pheatmap(assay(rlog(dds)), scale = "row", show_rownames = FALSE)

#方法二：
library(DESeq2)

#count tables from RNA-Seq data
cnts <- assay(count_se)
cond <- factor(c("Control", "Control", "Cancer", "Cancer"), levels = c("Control", "Cancer"))

#object construction
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)

#standard analysis
dds <- DESeq(dds)
res <- results(dds)
sig_genes <- rownames(res)[res$padj < 0.05]

#plot
plotMA(dds)
heatmap(assay(rlog(dds)))
pheatmap::pheatmap(assay(rlog(dds)))















