#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Semidán Robaina Estévez (srobaina@ull.edu.es)
(day samples, reads >= 50 pb)

conda activate diffexpr
"""

import os
import numpy as np
import pandas as pd
from diffexpr.py_deseq import py_DESeq2
import logging
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
rpy2_logger.setLevel(logging.ERROR)


def deseq2Normalize(counts, coldata, factor):
    """
    Run deseq2
    """
    dds = py_DESeq2(count_matrix=counts,
                    design_matrix=coldata[factor],
                    design_formula=f'~ 1',
                    gene_column='gene_id')
    dds.run_deseq() 
    deseq2_counts = dds.normalized_count()
    return deseq2_counts


numerical_factors = [
    'Latitude', 'Longitude', 'Depth.nominal', 'Temperature', 'Oxygen',
    'ChlorophyllA', 'Carbon.total', 'Salinity',
    'Gradient.Surface.temp(SST)', 'Fluorescence', 'CO3', 'HCO3', 'Density',
    'PO4', 'PAR.PC', 'NO3', 'Si', 'Alkalinity.total', 'Ammonium.5m',
    'Depth.Mixed.Layer', 'Lyapunov', 'NO2', 'Depth.Min.O2', 'NO2NO3',
    'Nitracline', 'Brunt.Väisälä', 'Iron.5m', 'Depth.Max.O2', 'Okubo.Weiss',
    'Residence.time'
]

condition_names = ["95", "98"]
read_length = 50
min_number_of_non_psbA_genes_with_counts = 10

# Initialize directories
for dtype in ['counts', 'deseq2']:
    for condition_name in condition_names:
        dirname = f'deseq2_normalized_counts/sorted_{dtype}_{condition_name}'
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        
if not os.path.exists('deseq2_normalized_counts/sorted_samples'):
    os.makedirs('deseq2_normalized_counts/sorted_samples')
    

def main():

    # Sort count and deseq2 data
    failed_conditions = []
    for condition_name in condition_names:

        day_data = pd.read_csv(os.path.join('raw_data', f'counts_{condition_name}.tsv'), sep='\t').set_index('gene_id')
        day_samples = list(day_data.columns)

        # Read excel
        meta = pd.read_excel('taragenes.xlsx', header=1, engine='openpyxl')[['ENA_ID', 'polar', 'Date/Time'] + numerical_factors]
        meta[meta == 0] = np.nan
        meta = meta.rename(columns={'ENA_ID': 'Sample'})
        meta = meta.set_index('Sample').reindex(index=day_samples)
        
        # Remove non-coding genes (some contain very high count numbers affecting results)
        day_data = day_data.loc[[gene for gene in day_data.index if gene.strip().split('_')[1].isdigit()], :]
        day_data.to_csv(f'deseq2_normalized_counts/raw_data/counts_coding_{condition_name}.tsv', sep='\t')
   
        # Only take samples with 10 or more non-psbA genes with non-zero counts 
        day_data_no_psbA = day_data.loc[[gene for gene in day_data.index if 'P9301_02451' not in gene.strip()], :]
        samples_to_discard = list(day_data.sum()[
           (np.count_nonzero(day_data_no_psbA, axis=0) <= min_number_of_non_psbA_genes_with_counts)
           ].index)

        # Remove 3 samples with low temperatures (< 16ºC) to match experimental set up.
        samples_to_discard.extend(["ERS490187", "ERS494173", "ERS494583"])

        # Drop samples with overall low count number
        meta = meta.drop(samples_to_discard, axis=0)
        meta.to_excel(f'deseq2_normalized_counts/metadata_{condition_name}.xlsx')

        # Get sorted lists of samples by environmental factor gradient
        coldata = {}
        for factor in numerical_factors:
            sorted_df = meta.sort_values(factor)[factor].dropna()
            coldata[factor] = sorted_df.to_frame()
            sorted_df.to_csv(f'deseq2_normalized_counts/sorted_samples/samples_sorted_by_{factor}_{condition_name}.csv')

        # Remove genes with zero counts accross conditions
        counts = pd.read_csv(f'deseq2_normalized_counts/raw_data/counts_coding_{condition_name}.tsv', sep='\t')
        counts = counts[counts.sum(axis=1) > 0]

        # Loop over environmental factors
        for factor in numerical_factors:

            print(f'Doing condition = {condition_name} and factor: {factor}')

            try:
                ordered_samples = coldata[factor].index.tolist()

                sorted_counts = counts.reindex(
                    ['gene_id'] + ordered_samples, axis=1
                )

                deseq2_counts = deseq2Normalize(sorted_counts, coldata, factor)

                sorted_counts.to_excel(f'deseq2_normalized_counts/sorted_counts_{condition_name}/sorted_counts_{condition_name}_{factor}.xlsx')
                deseq2_counts.to_excel(f'deseq2_normalized_counts/sorted_deseq2_{condition_name}/deseq2_norm_{condition_name}_{factor}.xlsx')

            except Exception as e:
                failed_conditions.append(f'condition = {condition_name}, factor: {factor} failed with exception: {e}')
                
    with open('deseq2_normalized_counts/failed_conditions.txt', 'w') as tfile:
        for elem in failed_conditions:
            tfile.write(elem + '\n')
                

if __name__ == '__main__':
    main()