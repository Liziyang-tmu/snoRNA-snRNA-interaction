#!/bin/bash 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-12:00:00
#SBATCH --mem=200G
#SBATCH --partition=gpu


module load gcc/8.3.0
module load python/3.7.6
module load samtools


################   sub-function (Don't change any parameters)   ##########################################

	## get the primary alignments
	samtools view -h $Outprefix'Aligned.sortedByCoord.out.bam' | awk '$1~/^@/ || NF<21' > $Outprefix'_nonchimeric_temp.sam'
	samtools view -h $Outprefix'Aligned.sortedByCoord.out.bam' | awk '$1!~/^@/ && NF==21 && $13=="HI:i:1"'> $Outprefix'_chimeric_temp.sam'
	samtools view -bS -F 0x900 -o $Outprefix'_nonchimeric_pri.bam' $Outprefix'_nonchimeric_temp.sam'
	samtools view -h $Outprefix'_nonchimeric_pri.bam' > $Outprefix'_nonchimeric_pri.sam'
	cat $Outprefix'_nonchimeric_pri.sam' $Outprefix'_chimeric_temp.sam' > $Outprefix'_pri.sam'
	rm -f $Outprefix'_nonchimeric_temp.sam' $Outprefix'_chimeric_temp.sam' $Outprefix'_nonchimeric_pri.bam' $Outprefix'_nonchimeric_pri.sam' 

	## divide alignments to different types
	python $ScriptPath/gaptypes.py $Outprefix'_pri.sam' $Outprefix'_pri' -1 15 1

	## remove splice junction alignments
	python $ScriptPath/gapfilter.py $Gtf $Outprefix'_prigap1.sam' $Outprefix'_prigap1_filtered.sam' 11 yes 
}
##########################################################################################################################


################ Set up the python3 environment  ############################
#export PATH=/auto/rcf-proj3/zl2/minjiez/software/conda/anaconda3/bin:$PATH  #
#unset PYTHONPATH                                                            #
#PyPath=/auto/rcf-proj3/zl2/minjiez/software/conda/anaconda3/bin             #
#ScriptPath=/auto/rcf-proj2/zl2/minjiez/software/CRSSANT-master/scriptsi     #
ScriptPath=/project/zhipengl_72/minjiez/software/CRSSANT-master/scripts      #
#############################################################################
