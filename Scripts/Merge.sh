#!/bin/bash -l
#merge and crssant
#SBATCH --job-name=fastq
#SBATCH --partition=cpuPartition
#SBATCH --nodes=16                 ##作业申请的节点数，按需修改，<=16
#SBATCH --ntasks-per-node=1        ##作业申请的每个节点的任务数，按需修改，建议为1
#SBATCH --cpus-per-task=16          ##作业申请的每个任务使用的核心数，按需修改，<=128
#SBATCH --error=%j.err
#SBATCH --output=%j.out
#SBATCH --account=minjielab

# 请确保此处设置的nodes X ntasks-per-node X cpus-per-task = 自己想用的线程数，或者脚本命令中使用-t等参数指定的线程数

##############################
merge(){
	cd $WorkPath
	python $scriptPath/merger.py $merge_file'gap1_filtered.sam' $merge_file'trans.sam' $Outprefix'_merged'
}
## Crssant script
##############################
crssant() {
        cd $WorkPath
        samtools view -bS -o $file'.bam' $file'.sam'
        samtools sort -o $file'_sort.bam' $file'.bam'
        bedtools genomecov -bg -split -strand + -ibam $file'_sort.bam' > $file'_plus.bedgraph'
        bedtools genomecov -bg -split -strand - -ibam $file'_sort.bam' > $file'_minus.bedgraph'
        #python $scriptPath/get_genefeature_bed4crssant_v2.py $genebed  $file'_plus.bedgraph'  $file'_minus.bedgraph' 20  $file'_crssant_gene' #Optional
        #python $scriptPath/crssant_v3.py -cluster cliques -t_o 0.2 -covlimit 50 -out ./ $file'.sam' $file'_crssant_gene.bed' $file'_plus.bedgraph',$file'_minus.bedgraph' #Optional
        python $scriptPath/crssant_v3.py -cluster cliques -t_o 0.2 -covlimit 50 -out ./ $file'.sam' $genebed $file'_plus.bedgraph',$file'_minus.bedgraph'
        #rm $file'.bam'
} 
##############################
scriptPath=~/software/script
genebed=~/database/hg38genmaskaddAllSnord1156/hg38mask14addSNORD113456_feature_gene_KARR_7S_reverse_240709.bed
#genebed=/mnt/hpc/home/zhangminjie/database/mm10/resort_left_joined_with_mm10_Gmtransed_linked_up_10_all.bed
#Gtf=/mnt/hpc/home/zhangminjie/database/hg38pri/gencode.v33.primary_assembly.annotation.gtf
##############################
#Outprefix=PARIS_rep1_pri
#WorkPath=/mnt/data/home/liziyang/work/STAR/hg38_index
#file=PARIS_rep1
#merge
#crssant



#file=chrM_all_merged
#crssant

