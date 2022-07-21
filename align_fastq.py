#!/usr/bin/env python3
# coding: utf-8
# conda activate samtools

"""
1. Align metaT to reference genome with bwa

   Semidán Robaina Estévez (srobaina@ull.edu.es)
   Python >= 3.6
"""
import os
from metaT2genome.src.align import has_bwa_index, makeBWAindex, bwaAlign
from metaT2genome.src.utils import getFastqPairedFiles
from metaT2genome.src.filtersam import extractSegmentsWithMDtag

work_dir = os.getcwd()
fasta_file = 'MIT9301/Prochlorococcus_marinus_str_MIT_9301.fasta'
data_dir = "/prochlorococcus/98_filtered_fastq/"
sam_output_dir = "/prochlorococcus/sam_files"
n_threads = 20


if not has_bwa_index(work_dir):
    print('Building bwa index\n')
    makeBWAindex(fasta_file)

# Iterate over conditions
failed_conditions = []
paired_fastqs = getFastqPairedFiles(data_dir) 
n_conds = len(paired_fastqs)

for n, (condition, (fastq_1_file, fastq_2_file)) in enumerate(paired_fastqs.items()):
    
    print(f'Processing condition ({n + 1}/{n_conds}): {condition}')

    try:
        print('\tAligning metaT to reference genome')
        bwaAlign(fasta_file, os.path.join(data_dir, fastq_1_file),
                 os.path.join(data_dir, fastq_2_file), 
                 n_threads=n_threads, output_dir=f'{sam_output_dir}/{condition}.sam')
        
        extractSegmentsWithMDtag(f'{sam_output_dir}/{condition}.sam',
                                 output_dir=f'{sam_output_dir}/{condition}_MD_tags.sam')
        
        os.remove(f'{sam_output_dir}/{condition}.sam')

    except Exception:
        failed_conditions.append(condition)
        print(f'Failed condition: {condition} with exception {Exception}')
