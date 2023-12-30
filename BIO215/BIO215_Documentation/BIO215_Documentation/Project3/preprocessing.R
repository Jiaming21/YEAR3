#################### Appendix: Pre-processing SRR from NCBI #################### 
# Step 1: Download SRR from NCBI
system("wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-11/SRR002/001/SRR2001703/SRR2001703.sralite.1](https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-11/SRR002/001/SRR2001703/SRR2001703.sralite.1")
                                                                                                                 
# Change the file name
system("mv SRR2001703.sralite.1 SRR2001703")  

# Step 2: Decompression data by fastq-dump

# 1. Use fastq-dump to decompress SRR file.
system("fastq-dump --gzip --skip-technical --split-files --clip --readids -O ./ SRR2001703")

# fastq-dump runs only on one thread, so it's really time consuming (nearly one hour). You can put this process running at background. 
system("(nohup fastq-dump --gzip --skip-technical --split-files --clip --readids -O ./ SRR2001703)&")

# 2. After decompress, you can also randomly select a small amount of data to run the test pipeline. 
# In this lab session, we selected 50% of the total data to reduce computational resources usage for every user. 
# To achieve this, use seqtk command.

system("seqtk sample -s11 SRR2001703_1.fastq.gz 0.5 | gzip > SRR2001703_1.fq.gz")
system("seqtk sample -s11 SRR2001703_2.fastq.gz 0.5 | gzip > SRR2001703_2.fq.gz")

  






                                                                                                                 