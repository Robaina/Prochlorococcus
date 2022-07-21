#!/usr/bin/env python3
# coding: utf-8

"""
1. Filter SAM file by percent identity cutoff and percent matched sequence
2. Sort paired-end entries by name
3. Count reads with htseq-count
4. TPM-normalize counts

   Semidán Robaina Estévez (srobaina@ull.edu.es)
   Python >= 3.6
"""
import os
from metaT2genome.src.count import (htseqCount, tpmNormalizeHtseqOutput,
                                    aggregateTPMresults, aggregateCountsresults)
from metaT2genome.src.filtersam import filterSAMbyIdentity, filterSAMbyPercentMatched, filterSAMbyReadLength
from metaT2genome.src.utils import deleteTemporaryFiles, sortSAMbyName

work_dir = os.getcwd()
gtf_file = 'MIT9301/MIT9301.gtf'
gbk_file = 'MIT9301/MIT9301.gb'
samfiles_dir = "/prochlorococcus/sam_files"
PM_samfiles_dir = "/prochlorococcus/PM_sam_files"
os.makedirs(PM_samfiles_dir, exist_ok=True)
os.chmod(PM_samfiles_dir, 0o0777)

matched_cutoff = 50
min_length = 50

samfiles = [file for file in os.listdir(samfiles_dir)]
n_conds = len(samfiles)


for n, samfile in enumerate(samfiles):
    condition = samfile.split('_')[0]
    print(f'Processing condition ({n + 1}/{n_conds}): {condition}')

    print(f'\t1.1. Filtering SAM at {matched_cutoff}% percent matched')
    filterSAMbyPercentMatched(f'{samfiles_dir}/{condition}_MD_tags.sam',
                              matched_cutoff=matched_cutoff,
                              output_path=f'{PM_samfiles_dir}/{condition}_PM_{matched_cutoff}.sam')

    print(f'\t1.2. Filtering SAM by read length')
    filterSAMbyReadLength(os.path.join(PM_samfiles_dir, f'{condition}_PM_{matched_cutoff}.sam'),
                          min_length=min_length,
                          output_path=os.path.join(PM_samfiles_dir, f'{condition}_PM_{matched_cutoff}_L_{min_length}.sam'))

# Iterate over identity cutoff
for identity_cutoff in [95, 98, 100]:

    filtered_dir = f'filtered_{identity_cutoff}_sam_files'
    if not os.path.exists(filtered_dir):
            os.makedirs(filtered_dir)
            os.chmod(filtered_dir, 0o0777)

    counts_dir = f'counts_{identity_cutoff}'
    if not os.path.exists(counts_dir):
            os.makedirs(counts_dir)
            os.chmod(counts_dir, 0o0777)

    tpm_dir = f'tpm_{identity_cutoff}'
    if not os.path.exists(tpm_dir):
            os.makedirs(tpm_dir)
            os.chmod(tpm_dir, 0o0777)

    # Iterate over conditions
    failed_conditions = []
    for n, samfile in enumerate(samfiles):

        condition = samfile.split('_')[0]
        print(f'Processing condition ({n + 1}/{n_conds}): {condition}')

        try:

                print(f'\t2.Filtering SAM at {identity_cutoff}% identity')
                filterSAMbyIdentity(os.path.join(PM_samfiles_dir, f'{condition}_PM_{matched_cutoff}_L_{min_length}.sam'),
                                    identity_cutoff=identity_cutoff,
                                    output_path=f'{PM_samfiles_dir}/{condition}_filtered_at_{identity_cutoff}.sam')

                print('\t3.Sorting SAM by name') # required by htseq-count
                sortSAMbyName(f'{PM_samfiles_dir}/{condition}_filtered_at_{identity_cutoff}.sam',
                              output_dir=f'{filtered_dir}/{condition}_filtered_at_{identity_cutoff}_sorted.sam')

                print('\t4.Counting reads')
                htseqCount(f'{filtered_dir}/{condition}_filtered_at_{identity_cutoff}_sorted.sam',
                           gtf_file, feature_type='gene',
                           feature_id='gene_id', output_dir=os.path.join(counts_dir, f'{condition}_counts.tsv'))

                print('\t5.TPM normalizing counts\n')
                tpmNormalizeHtseqOutput(os.path.join(counts_dir, f'{condition}_counts.tsv'), gbk_file,
                                        output_dir=os.path.join(tpm_dir, f'{condition}_tpm.tsv'))

        except Exception:
                failed_conditions.append(condition)
                print(f'Failed condition: {condition} with exception {Exception}')


    aggregateTPMresults(tpm_dir=tpm_dir, output_dir=f'tpm_{identity_cutoff}.tsv')
    aggregateCountsresults(counts_dir=counts_dir, output_dir=f'counts_{identity_cutoff}.tsv')

deleteTemporaryFiles(PM_samfiles_dir)
