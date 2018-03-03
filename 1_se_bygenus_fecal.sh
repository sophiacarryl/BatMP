#PBS -N spiec-easi_bygenus_fecal
#PBS -S /bin/bash

#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=4gb

#PBS -o /group/gilbert-lab/Lutz/Batmicrobiome/batrun_finaldata_for_phyloseq_Oct2017/spiec_easi/se_bygenus_fecal.out
#PBS -e /group/gilbert-lab/Lutz/Batmicrobiome/batrun_finaldata_for_phyloseq_Oct2017/spiec_easi/se_bygenus_fecal.err


module load gcc/6.2.0
module load python/2.7.13
module load R/3.4.1

Rscript /group/gilbert-lab/Lutz/Batmicrobiome/batrun_finaldata_for_phyloseq_Oct2017/spiec_easi/1_se_bygenus_fecal.txt
