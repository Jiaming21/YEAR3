# use zcat to look at the fq.gz file
# chiplot can combine diagram together
# cat SRR2001703.result.txt 可以看各个物种的mitogene的覆盖程度

# 用split只获得物种名信息
# excel：数据
# how to get the species name: 
# code node_id 
# data  test to columns
# other _
# copy 

# 拼接后两列 打入公式=C2&" "&D2
# 需要利用COX1-ND5-CYTB.tree拖入normal tree去获得 copy leaves ID 然后到excel里split
# 利用SRR2001703.result.txt获得coverage的信息，然后spilt

reads -> contigs -> scafold (k-mer)