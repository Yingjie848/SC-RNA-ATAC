#BSUB -W 6:00
#BSUB -n 3
#BSUB -R "rusage[mem=100]"
#BSUB -e log/run_ArchR.MBC5_exclude_C8.err
#BSUB -o log/run_ArchR.MBC5_exclude_C8.out

source /home/zhuy1/my_apps/miniconda3/bin/activate /home/zhuy1/my_apps/env-R-4

Rscript seurat_workflow/code/ArchR/run_ArchR.MBC5_exclude_C8.R


