#!/usr/bin/env python3
# coding: utf-8
# conda activate metat2genome

"""
seqkit grep -n -f /usr/storage/semidan/proyecto_laura/98_filtered_fastq_names/ERS488299.txt -o temp_1.fastq.gz /usr/gonzalez/metagenomes/salazar2019/download/pe/fastqfinal/ERS488299_1.fastq.gz
"""


import os
import pandas as pd

from metaT2genome.src.utils import filterFastqByReadNames


names_dir = '/usr/gonzalez/laura/results'
fastq_files_dir = "/usr/gonzalez/metagenomes/salazar2019/download/pe/fastqfinal/"
names_output_dir = "/usr/storage/semidan/proyecto_laura/prochlorococcus/98_filtered_fastq_names/"
fastq_output_dir = "/usr/storage/semidan/proyecto_laura/prochlorococcus/98_filtered_fastq/"


# Preprocess lists of filtered read names to deal with missing pairs
type_1, type_2 = '_1', '_2'
conditions = set([file.split("_")[0] for file in os.listdir(names_dir)])

for condition in conditions:

    names1 = os.path.join(names_dir, f"{condition}{type_1}.txt")
    names2 = os.path.join(names_dir, f"{condition}{type_2}.txt")
    names = []
    if os.path.exists(names1) and os.path.exists(names2):
        df1 = pd.read_csv(names1)
        values1 = df1.values.flatten().tolist()
        df2 = pd.read_csv(names2)
        values2 = df1.values.flatten().tolist()
        names = list(set(values1 + values2))
    elif os.path.exists(names1):
        df = pd.read_csv(names1)
        names = df.values.flatten().tolist()
    elif os.path.exists(names2):
        df = pd.read_csv(names2)
        names = df.values.flatten().tolist()
    
    if names:
        lines = [name + "\n" for name in names]
        with open(os.path.join(names_output_dir, f"{condition}.txt"), "w") as outfile:
            outfile.writelines(lines)

# Filter Fastq files by read names
for n, condition in enumerate(conditions):
    print(f"Running condition: {condition} ({n+1} / {len(conditions)})")
    names = os.path.join(names_output_dir, f"{condition}.txt")
    file1 = os.path.join(fastq_files_dir, f"{condition}{type_1}.fastq.gz")
    file2 = os.path.join(fastq_files_dir, f"{condition}{type_2}.fastq.gz")
    outfile1 = os.path.join(fastq_output_dir, f"{condition}{type_1}.fastq.gz")
    outfile2 = os.path.join(fastq_output_dir, f"{condition}{type_2}.fastq.gz")

    if not os.path.exists(outfile1):
        filterFastqByReadNames(
            file1,
            names,
            output_file=outfile1
        )
    
    if not os.path.exists(outfile2):
        filterFastqByReadNames(
            file2,
            names,
            output_file=outfile2
        )