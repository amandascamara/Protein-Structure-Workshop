#!/bin/bash

#SBATCH --auks=yes
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100000
#SBATCH --output=4pc8.out

dbdir=/filer-dg/agruppen/DG/camara/AlphaFold_databases

. /etc/profile.d/modules.sh
module load alphafold
module load cuda-toolkit

export OPENMM_CUDA_COMPILER=/opt/Bio/cuda-toolkit/11.6/bin/nvcc <8>
export SINGULARITY_BIND="/filer,/filer-dg/agruppen/DG/camara/AlphaFold_databases/"

singularity run -B /home/camara/AlphaFold/alphafold/modified_run_alphafold.py:/app/alphafold/run_alphafold.py /opt/Bio/alphafold/2.30.0/bin/alphafold \
 --max_template_date=2023-03-16 \
 --uniref90_database_path=$dbdir/uniref90/uniref90.fasta \
 --mgnify_database_path=$dbdir/mgnify/mgy_clusters_2022_05.fa \
 --template_mmcif_dir=$dbdir/pdb_mmcif/mmcif_files \
 --obsolete_pdbs_path=$dbdir/pdb_mmcif/obsolete.dat \
 --use_gpu_relax=true \
 --bfd_database_path=$dbdir/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
 --uniref30_database_path=$dbdir/uniref30/UniRef30_2021_03 \
 --pdb70_database_path=$dbdir/pdb70/pdb70 \
 --fasta_paths=/home/camara/AlphaFold/workshop/4pc8.fasta \
 --model_preset=monomer_ptm \
 --data_dir=$dbdir \
 --output_dir=/filer-dg/agruppen/DG/camara/AlphaFold_models/workshop/
