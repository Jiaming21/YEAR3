##############################################################################
########      Decoding m6A Sites Through MeRIP-Seq Analysis             ######
##############################################################################
#查看8个FASTQ文件，并且软链接到本地data目录下
system("ls -lh /data/BIO215/*")
#只要RNA-Seq
system("ln -s /data/BIO215/SRR22287341.fastq /home/jiaming.huang21/redo/BIO215/project2/data/SRR22287341.fastq;
       ln -s /data/BIO215/SRR22287343.fastq /home/jiaming.huang21/redo/BIO215/project2/data/SRR22287343.fastq;
       ln -s /data/BIO215/SRR22287345.fastq /home/jiaming.huang21/redo/BIO215/project2/data/SRR22287345.fastq;
       ln -s /data/BIO215/SRR22287347.fastq /home/jiaming.huang21/redo/BIO215/project2/data/SRR22287347.fastq
       ")
#查看文件结构
system("tree ~/redo/BIO215/project2")

#Step1: 使用genome aligner Hisat2，将FASTQ文件与hg38 map
system("hisat2 -x /data/genome_indx/genome_hg38 -U /home/jiaming.huang21/redo/BIO215/project2/data/SRR22287341.fastq -S /home/jiaming.huang21/redo/BIO215/project2/data/SRR22287341.sam;
       hisat2 -x /data/genome_indx/genome_hg38 -U /home/jiaming.huang21/redo/BIO215/project2/data/SRR22287343.fastq -S /home/jiaming.huang21/redo/BIO215/project2/data/SRR22287343.sam;
       hisat2 -x /data/genome_indx/genome_hg38 -U /home/jiaming.huang21/redo/BIO215/project2/data/SRR22287345.fastq -S /home/jiaming.huang21/redo/BIO215/project2/data/SRR22287345.sam;
       hisat2 -x /data/genome_indx/genome_hg38 -U /home/jiaming.huang21/redo/BIO215/project2/data/SRR22287347.fastq -S /home/jiaming.huang21/redo/BIO215/project2/data/SRR22287347.sam
       ")

#Step2: bam转sam
system("cd ~/redo/BIO215/project1.2/data")
system("samtools view -Sb SRR22287341.sam > SRR22287341.bam;
       samtools view -Sb SRR22287343.sam > SRR22287343.bam;
       samtools view -Sb SRR22287345.sam > SRR22287345.bam;
       samtools view -Sb SRR22287347.sam > SRR22287347.bam
       ")

#Step3: 利用peak calling比较IP与input/RNA-Seq的覆盖率

#加载包
library(exomePeak2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

original_seqlevels <- seqlevels(txdb)
print(original_seqlevels)

#TxDb取子集只包括17号染色体
target_seqlevels <- "chr17"
seqlevels(txdb) <- target_seqlevels

txdb_subset <- keepSeqlevels(txdb, target_seqlevels)


#Peak calling for cancer condition
data_dir_1 <- "~/redo/BIO215/project1.2/data"
bam_ip_1 <- file.path(data_dir_1,c(
                                  "SRR22287345.bam",
                                  "SRR22287347.bam"
                                  )
                     )
                    
data_dir <- "~/redo/BIO215/project1.1/data"
bam_input_1 <- file.path(data_dir,c(
                                  "SRR22287346.bam",
                                  "SRR22287348.bam"
                                  )
                   )
                      

peaks_cancer <- exomePeak2(#导入MeRIP-Seq (cancer) 的BAM文件
                           bam_ip = bam_ip,
                           #导入RNA-Seq (control) (cancer) 的BAM文件
                           bam_input = bam_input,
                           #转录组注释
                           txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           #参考基因组名字
                           genome = "hg38",
                           #命名输出文件
                           experiment_name = "Peaks_cancer")


data_dir_2 <- "~/redo/BIO215/project1.2/data"
bam_ip_2 <- file.path(data_dir_2,c(
                                  "SRR22287341.bam",
                                  "SRR22287343.bam"
                                  )
                     )

data_dir_2 <- "~/redo/BIO215/project1.1/data"
bam_input_2 <- file.path(data_dir_2,c(
                                     "SRR22287342.bam",
                                     "SRR22287344.bam"
                                     )
                         )

peaks_normal <- exomePeak2(#导入MeRIP-Seq (normal) 的BAM文件
                           bam_ip = bam_ip_2,
                           #导入RNA-Seq (control) (normal) 的BAM文件
                           bam_input = bam_input_2,
                           #转录组注释
                           txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           #参考基因组名字
                           genome = "hg38",
                           #命名输出文件
                           experiment_name = "Peaks_normal")
#bam要可视化需要sort再index

system("samtools sort /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287341.bam -o /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287341.bam")
system("samtools index /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287341.bam > /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287341.bam.bai")

system("samtools sort /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287343.bam -o /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287343.bam")
system("samtools index /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287343.bam > /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287343.bam.bai")

system("samtools sort /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287345.bam -o /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287345.bam")
system("samtools index /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287345.bam > /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287345.bam.bai")

system("samtools sort /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287347.bam -o /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287347.bam")
system("samtools index /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287347.bam > /home/jiaming.huang21/redo/BIO215/project1.2/data/SRR22287347.bam.bai")

#Step4: m6A peaks 功能性分析

#使用clusterProfiler去辨认出在normal和cancer中的peaks的富集的GO terms

#enrichGO: enrich gene oncology

library(clusterProfiler)
#自动识别gene_id类型，再进行转换
result <- enrichGO(unlist(strsplit(mcols(peaks_cancer)$geneID, ",")), org.Hs.eg.db, "ENTREZID")

#为什么要unlist()?
barplot(result)

#创建pie charts 和 meta-gene plots 

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene #注释文件：储存genomic的信息
utr5_ranges <- fiveUTRsByTranscript(txdb) #5’UTR
utr3_ranges <- threeUTRsByTranscript(txdb) #3’UTR
cds_ranges <- cdsBy(txdb) #CDS

#计算peaks在这些区域上的数量
peaks_cancer_sum_utr5 <- sum(peaks_cancer %over% utr5_ranges)
peaks_normal_sum_utr5 <- sum(peaks_normal %over% utr5_ranges)

peaks_cancer_sum_utr3 <- sum(peaks_cancer %over% utr3_ranges)
peaks_normal_sum_utr3 <- sum(peaks_normal %over% utr3_ranges)

peaks_cancer_sum_cds <- sum(peaks_cancer %over% cds_ranges)
peaks_normal_sum_cds <- sum(peaks_normal %over% cds_ranges)

#画饼图(显示在各个区域的占比)
region_names <- c("Cancer UTR5", "Cancer UTR3", "Cancer CDS")
region_sizes <- c(peaks_cancer_sum_utr5, 
                  peaks_cancer_sum_utr3, 
                  peaks_cancer_sum_cds
                  ) #创建一个包含占比的向量
colors <- rainbow(length(region_names)) #创建颜色向量，用于标识各个区域
pie(region_sizes, labels = region_names, col = colors, main = "Region Proportions") #绘制饼图
legend("topright", legend = region_names, fill = colors, cex = 0.8) #添加图例

region_names <- c("Normal UTR5", "Normal UTR3", "Normal CDS")
region_sizes <- c(peaks_normal_sum_utr5, 
                  peaks_normal_sum_utr3, 
                  peaks_normal_sum_cds
                  )
colors <- rainbow(length(region_names))
pie(region_sizes, labels = region_names, col = colors, main = "Region Proportions")
legend("topright", legend = region_names, fill = colors, cex = 0.8)

################################################################################

#画meta-genes plot
library(Guitar) 
#创建一个 Meta-genes plot
setwd("/home/jiaming.huang21/redo/BIO215/project1.2/fig")
meta_genes_plot_origin <- GuitarPlot(txTxdb = txdb,
                                     stGRangeLists = GRangesList("cancer" = unlist(peaks_cancer), "normal" = unlist(peaks_normal)),
                                     stGroupName = c("cancer", "normal"),
                                     pltTxType = c("mrna"),
                                     CI_interval = c(0.5, 0.5)
                                     )

meta_genes_plot_adjusted <- meta_genes_plot + theme_classic() + xlab("Coordinates")
#使用illustrator删除横坐标

sessionInfo()



