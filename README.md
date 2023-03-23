### scRNA-seq-Skin-microvascular-project
##scRNA-seq Skin microvascular project- single nuclei RNA-seq
##running cellranger for each sample
##cellranger job

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --job-name=cellranger_GFP_Y3
#SBATCH --output=jlog_cellranger_GFP_Y3

date

~/projects/rrg-glettre/programs/CellRanger/cellranger-7.0.1/cellranger count --id sample_Y3_GFP --sample Y3_02KK48_11258 --fastqs ~/projects/rrg-glettre/sequencing_datastore/raw_data/scRNAseq_skin-senescense_20220928 --transcriptome ~/projects/rrg-glettre/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/align_EGFP/mouse_GFP --localcores 4 --localmem 40

date
