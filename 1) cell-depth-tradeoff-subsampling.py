import os
import pandas as pd

# for file in ./*/fastqs/full/*S1_L001_R1_001.fastq.gz ; do echo $file; zcat $file | echo $((`wc -l`/4)); 

# to call this snakemake use (replace -j 10 for the number of processes you want)
# snakemake -j 10 -s cell-depth-tradeoff-subsampling.py--keep-going --rerun-incomplete -pn

# path to folder where indices, whitelist and transcript to gene files are located
REF_PATH = '/data/references'
# how to call kallisto
KALLISTO = 'kallisto'
#how to call bustools
BUSTOOLS = 'bustools'

# csv file with metadata information for each dataset
# must have the species information, subsampling levels for reads, total reads,
# technology used and path to one fastq file for R1 and one for R2
# if you have multiple files for each read tehy must be concatenated first
metadatas=pd.read_csv('./cell-depth-tradeoff-metadata.tsv', sep='\t')

def make_t2g_path(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    species = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['species'].values[0]
    if species=='mouse':
        T2G = os.path.join(REF_PATH,'mus_musculus-ensembl-96/transcripts_to_genes.txt')
    if species=='human':
        T2G = os.path.join(REF_PATH,'homo_sapiens-ensembl-96/transcripts_to_genes.txt')
    return T2G

def make_index_path(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    species = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['species'].values[0]
    if species=='mouse':
        INDEX = os.path.join(REF_PATH,'mus_musculus-ensembl-96/transcriptome.idx')
    if species=='human':
        INDEX = os.path.join(REF_PATH,'homo_sapiens-ensembl-96/transcriptome.idx')
    return INDEX

def make_whitelist_path(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    TECHNOLOGY = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['technology'].values[0]
    if TECHNOLOGY=='10xv3':
        WHITELIST = os.path.join(REF_PATH,'10xv3_whitelist.txt')
    if TECHNOLOGY=='10xv2':
        WHITELIST = os.path.join(REF_PATH,'10xv2_whitelist.txt')
    return WHITELIST

def fetch_technology(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    TECHNOLOGY = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['technology'].values[0]
    return TECHNOLOGY

def fetch_read1_filepath(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    DATASET_SAMPLE_PATH = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['dataset_sample_path'].values[0]   
    READ1_FILES = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['concat_read1_file'].values[0].split(',')
    #remove trailing spaces
    READ1_FILES = [ read1_file.strip() for read1_file in READ1_FILES]
    READ1_FILEPATHS = [os.path.join(DATASET_SAMPLE_PATH, read_filename) for read_filename in READ1_FILES]    
    return READ1_FILEPATHS
    
def fetch_read2_filepath(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    DATASET_SAMPLE_PATH = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['dataset_sample_path'].values[0] 
    READ2_FILES = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['concat_read2_file'].values[0].split(',')
    #remove trailing spaces
    READ2_FILES = [ read2_file.strip() for read2_file in READ2_FILES]
    READ2_FILEPATHS = [os.path.join(DATASET_SAMPLE_PATH, read_filename) for read_filename in READ2_FILES]      
    return READ2_FILEPATHS

def fetch_subsampling_depths(wildcards):
    DATASET_SAMPLE_ID = wildcards.dataset_sample_id
    subsampling_depths = metadatas[metadatas['dataset_sample_id']==DATASET_SAMPLE_ID]['subsampling_depths'].values[0]
    subsampling_depths = [int(x) for x in subsampling_depths.split(',')]
    return subsampling_depths

final_subsampled_fastqs = []
final_count_matrices = []
final_genecount_mtx = []

for dataset_sample_id in metadatas[metadatas['process']==1]['dataset_sample_id']:
    for sub_string in metadatas[metadatas['dataset_sample_id']==dataset_sample_id]['subsampling_depths']:
        subsampling_depths = [int(x) for x in sub_string.split(',')]
        for sub in subsampling_depths:
            final_genecount_mtx.append(dataset_sample_id + '/genecounts/genecounts_subsampled_'+str(sub)+'/genecounts.mtx')
print( '==========================================================')
print(final_genecount_mtx)

rule all:
    input:
        final_genecount_mtx,
        
rule subsample:
    input:
        READ1_FILEPATH = fetch_read1_filepath,
        READ2_FILEPATH = fetch_read2_filepath
    params:
        SUBSAMPLING=lambda wildcards: f'{wildcards.subsampling}',

    output:
        R1 = '{dataset_sample_id}/tmp/subsampled_{subsampling}/subsampled_{subsampling}_R1.fq',
        R2 = '{dataset_sample_id}/tmp/subsampled_{subsampling}/subsampled_{subsampling}_R2.fq'
    benchmark:
        "benchmarks/{dataset_sample_id}/{subsampling}/seqtk.txt"
    shell:
        """       
        if (({params.SUBSAMPLING}==0)); then 
        cp {input.READ1_FILEPATH} {output.R1}
        cp {input.READ2_FILEPATH} {output.R2}
        fi
        
        if (({params.SUBSAMPLING}!=0)); then
        seqtk sample -2 -s100 {input.READ1_FILEPATH} {params.SUBSAMPLING} > {output.R1} && \
        seqtk sample -2 -s100 {input.READ2_FILEPATH} {params.SUBSAMPLING} > {output.R2}
        fi
        """  
        
        
rule run_kallisto:
    input:
        R1 = '{dataset_sample_id}/tmp/subsampled_{sub}/subsampled_{sub}_R1.fq',
        R2 = '{dataset_sample_id}/tmp/subsampled_{sub}/subsampled_{sub}_R2.fq',
    params: 
        DATASET_SAMPLE_ID = '{dataset_sample_id}',
        INDEX=make_index_path,
        TECHNOLOGY = fetch_technology
    output:
        KALLISTO_OUT=directory('{dataset_sample_id}/kallisto_out/kallisto_out_subsampled_{sub}')
    benchmark:
        "benchmarks/{dataset_sample_id}/{sub}/kallisto.txt"
    shell:
        """
        {KALLISTO} bus -i {params.INDEX} -x {params.TECHNOLOGY} -o {output} {input.R1} {input.R2} && \
        rm {input.R1} && rm {input.R2}
        """
        
rule run_bustools_correct_sort:
    input: 
        WHITELIST=make_whitelist_path,
        KALLISTO_OUT='{dataset_sample_id}/kallisto_out/kallisto_out_subsampled_{sub}',
    params:
        temp_folder ='{dataset_sample_id}/tmp/',
        genecounts_directory = directory('{dataset_sample_id}/genecounts/genecounts_subsampled_{sub}/')
    output:
        "{dataset_sample_id}/genecounts/genecounts_subsampled_{sub}/output.correct.sort.bus"
    benchmark:
        "benchmarks/{dataset_sample_id}/{sub}/correct.txt"
    shell:
        """
        mkdir -p {params.genecounts_directory} 
        {BUSTOOLS} correct -w {input.WHITELIST} {input.KALLISTO_OUT}/output.bus -p | bustools sort -o {output} - && \
        rm {input.KALLISTO_OUT}/output.bus
        """
        
rule run_bustools_count:
# this rule will run bustools count and then delete all the bus files created 
    input: 
        busfile="{dataset_sample_id}/genecounts/genecounts_subsampled_{sub}/output.correct.sort.bus",
        kallisto_out=lambda wildcards: f'{wildcards.dataset_sample_id}/kallisto_out/kallisto_out_subsampled_{wildcards.sub}',
        T2G=make_t2g_path
    params:
        genecounts_directory = directory('{dataset_sample_id}/genecounts/genecounts_subsampled_{sub}')
    output:
        '{dataset_sample_id}/genecounts/genecounts_subsampled_{sub}/genecounts.barcodes.txt',
        '{dataset_sample_id}/genecounts/genecounts_subsampled_{sub}/genecounts.genes.txt',
        '{dataset_sample_id}/genecounts/genecounts_subsampled_{sub}/genecounts.mtx'
    benchmark:
        "benchmarks/{dataset_sample_id}/{sub}/count.txt"
    shell:
        """
         {BUSTOOLS} count -o {params.genecounts_directory}/genecounts \
         -g {input.T2G} \
         -e {input.kallisto_out}/matrix.ec \
         -t {input.kallisto_out}/transcripts.txt \
         --genecounts {input.busfile}  && \
         rm {input.busfile} 
        """
