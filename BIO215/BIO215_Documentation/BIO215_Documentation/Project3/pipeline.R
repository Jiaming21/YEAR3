# 创建文件架构
system("cd /home/jiaming.huang21/redo/BIO215/project3.1/")
system("mkdir metagenomics") # Create the main working folder
system("cd metagenomics")
system("mkdir 00.data") # To store raw sequencing data & cleaned data
system("mkdir 01.assemble") # To store assembled data
system("mkdir 02.annotation") # To store annotated data

system("ll") # Check whether the four folders are created succssfully

##################### Step 1: Data retrieval ###################################
# Using Soft Links for Data Organization

system("cd /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/00.data/") 

system("ln -s /data/BIO215/metagenomics/NC001566.1.R1.fq.gz /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/00.data/NC001566.1.R1.fq.gz")
system("ln -s /data/BIO215/metagenomics/NC001566.1.R2.fq.gz /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/00.data/NC001566.1.R2.fq.gz")

# Let's check whether the two files are soft linked successfully.
system("ll") # system中ll中不可用

##################### Step 2: Quality control by fastp #########################

# 将clean data 存在 00.data下
system("cd /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/00.data/") # Enter `01.data` folder
  
system("fastp \
-i NC001566.1.R1.fq.gz \
-o NC001566.1.clean.R1.fq.gz \
-I NC001566.1.R2.fq.gz \
-O NC001566.1.clean.R2.fq.gz \
-q 20 -u 10 -w 16")

# Quality Control Criteria (Parameter `-q`, `-u`, `-w`):
# -q: Specifies the quality threshold for filtering low-quality bases (default: 15). 
# -u: Specifies the minimum length for reads to be retained (default: 35). 
# -w: Specifies the window size for quality control sliding windows (default: 4).

########################## Step 3: Assembly by MEGIHIT #########################
##########################        create contigs      ########################## 
system("cd /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/01.assemble/") # Enter "01.assmebly" folder
  
system("megahit \
-1 /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/00.data/NC001566.1.clean.R1.fq.gz \
-2 /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/00.data/NC001566.1.clean.R2.fq.gz \
--tmp-dir ./ \
-m 0.8 \
-t 1 \
--k-list 31,51,71 \
--no-mercy \
-o NC001566.1_assembly")

# Additional parameters:
# `-m`: Sets the memory usage limit to 80% (0.8) of available RAM (adjust as needed).
# `-t`: Specifies the number of threads or CPU cores to use during assembly (adjust as needed).
# `--k-list`: Specifies the list of k-mer sizes used for assembly (adjust as needed).
# `--no-mercy`: Disables "no-mercy" mode, which reduces memory consumption.

#################### Step 4: Mitogenome annotation by MitoZ ####################
############################### 标上线粒体基因 #################################
#内圈高度表示的是sequence depth: 有多少reads supporting this gene
system("conda init bash")
system("conda activate mitozEnv") # Remember to enter the MitoZ environment first.

system("cd /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/02.annotation/")

# 在assembly中寻找线粒体的scaffold
system("mitoz findmitoscaf \
--fastafile /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/01.assemble/NC001566.1_assembly/final.contigs.fa \
--fq1 /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/00.data/NC001566.1.clean.R1.fq.gz \
--fq2 /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/00.data/NC001566.1.clean.R2.fq.gz \
--outprefix NC001566.1_findMito \
--thread_number 1 \
--requiring_taxa Arthropoda") # 只寻找节肢动物的线粒体基因，缩小范围

# 开始注释上线粒体基因
system("conda activate mitozEnv") # Remember to enter the MitoZ environment first.
system("cd /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/02.annotation/")
system("mitoz annotate \
--fastafile /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/02.annotation/NC001566.1_findMito.result/NC001566.1_findMito.mitogenome.fa \
--fq1 /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/00.data/NC001566.1.clean.R1.fq.gz \
--fq2 /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/00.data/NC001566.1.clean.R2.fq.gz \
--outprefix NC001566.1_annoMito \
--thread_number 1 \
--clade Arthropoda") # 此时还提供fq1与fq2是为了增加准确性



##################### Step 5: Sequence alignment by MUSCLE #####################

system("mkdir /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/03.alignment")
system("mkdir /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/04.tree")

system("cd /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/03.alignment")
system("ln -s /data/BIO215/metagenomics/COX1-ND5-CYTB.fa /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/03.alignment/COX1-ND5-CYTB.fa")

system("cd /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/03.alignment")
system("muscle -in COX1-ND5-CYTB.fa -out COX1-ND5-CYTB.aln.fa")

system("cd /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/03.alignment")
system("less COX1-ND5-CYTB.aln.fa")

################# Step 6: Construct phylogenetic tree by FastTree ##############

system("cd /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/04.tree")
system("FastTree /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/03.alignment/COX1-ND5-CYTB.aln.fa > COX1-ND5-CYTB.tree")

####################### Step 7: Reads mapping by SOAPAligner ###################
system("mkdir /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/05.mapping")
system("cd  /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/05.mapping")

# Soft link the reference assembly to "05.mapping" directory
system("ln -s /data/BIO215/metagenomics/COX1-ND5-CYTB.fa /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/05.mapping/COX1-ND5-CYTB.fa")
system("2bwt-builder COX1-ND5-CYTB.fa")

# Mapping reads
system("cd  /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/05.mapping")

# Soft the paired-end reads files to "05.mapping" directory
system("ln -s /data/BIO215/metagenomics/SRR2001703.clean.R1.fq.gz SRR2001703.clean.R1.fq.gz")
system("ln -s /data/BIO215/metagenomics/SRR2001703.clean.R2.fq.gz SRR2001703.clean.R2.fq.gz")

system("soap \
-a SRR2001703.clean.R1.fq.gz \
-b SRR2001703.clean.R2.fq.gz \
-D COX1-ND5-CYTB.fa.index \
-o SRR2001703.aln \
-2 SRR2001703.se.aln \
-M 4 -l 30 -r 1 -v 7 -m 200 -p 1")

# Alignment summary
system("cd  /home/jiaming.huang21/redo/BIO215/project3.1/metagenomics/05.mapping")

system("soap.coverage \
-cvg \
-i SRR2001703.aln SRR2001703.se.aln \
-refsingle COX1-ND5-CYTB.fa \
-o SRR2001703.result.txt \
-p 1")







