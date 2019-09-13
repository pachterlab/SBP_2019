Code for producing the analysis in the "Quantifying the tradeoff between sequencing depth and cell number in single-cell RNA-seq" by Valentine Svensson, Eduardo Beltrame and Lior Pachter



The workflow has 4 steps. The output data after each step can be downloaded from CaltechDATA at https://data.caltech.edu/records/1276

### 1) FASTQ Subsampling and processing with kallisto bus 

Starting from the raw FASTQ files, subsample and process with kallisto bus to generate genecount matrices. This is done with the snakemake file `1) cell-depth-tradeoff-subsampling.py`. It can be called  by doing `snakemake -j 10 -s 1) cell-depth-tradeoff-subsampling.py`. `-j 10` runs 10 parallel jobs. It uses the `.tsv` file `cell-depth-tradeoff-metadata.tsv` to provide information for each run, described below.

To run the snakemake file you need to have `seqkt` v1.3 or greater, `kallisto` v0.46 or greater and `bustools` v0.39.2 or greater. All of them can be installed with conda and the bioconda channel:
```
conda install -c bioconda seqtk
conda install -c bioconda kallisto
conda install -c bioconda bustools
```

The `FASTQ` files for the 3 datasets used in the paper can be downloaded from the 10x website:
```
https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3
https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/heart_10k_v3
https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/neuron_10k_v3
```

Once downloaded, for each dataset the multiple `.fastq.gz` files should be concatenated into a single file for R1 and a single file for R2 for each dataset. For example:
```
cat pbmc_10k_v3_S1_L001_R1_001.fastq.gz pbmc_10k_v3_S1_L002_R1_001.fastq.gz > pbmc_10k_v3_R1_concat_1.fastq.gz
```
It is necessary to have a single file for each read so that seqtk can subsample from it.

The output files produced for this step in the paper in CaltechDATA are in the tar file `snakemake_subsampling_output.tar.gz`

#### Metadata parameters for the subsampling with `cell-depth-tradeoff-metadata.tsv`

The path to each subsampled file should be provided in the `dataset_sample_path` column of the `cell-depth-tradeoff-metadata.tsv`. The name of the R1 file should be in the `concat_read1_file` column, and the same for R2 and the `concat_read2_file` column. 

The column `process` should be 1 for a dataset to be processed by the snakemake
The `species` and `technology` columns should match the lookup for the references in the snakemake file (indices, whitelist and t2g file)

Finally, the depths to be subsampled should be in the `subsampling_depths` column, in a single cell, and separated by commas (this is why it's a tsv file instead of csv, to avoid confusion).


### 2) Create `.H5AD` files with subsampled count matrices

For each dataset, we need to parse all gene count matrices produces (outputted by bustools in .mtx format) into an h5 file using anndata. The short notebook `2) create_H5AD.ipynb` does this and provides a quick check that it worked. 

The output files produced for this step in the paper in CaltechDATA are `pbmc10k_subsamples.h5ad`, `heart10k_subsamples.h5ad` and `neurons10k_subsamples.h5ad`


### 3) Run scVI

The script `3) run_scvi.py` will process the `.h5ad` files defined in the `adata_files` dictionary. Note that even though each sampled point is relatively quick to train in scVI (10-15 minutes), because we sampled each dataset across multiple depths and numbers of cells, we end up with ~200 points to train per dataset. This can take a day or two to finish, since the script does not run parallel jobs by default. 


The default files used to reproduce the workflow are:
```
adata_files= {
'neurons10k' :'./neurons10k_subsamples.h5ad',     
'heart10k' :'./heart10k_subsamples.h5ad',
'pbmc10k' :'./pbmc10k_subsamples.h5ad',    
}
```

We keep 15% of the cells in the validation set and use 85% for training. This can be adjusted by changing the fraction of `total_cells` used for `n_retained_cells`.
```
    total_cells = adata.n_obs
    n_retained_cells = int(0.85*total_cells)
```

The list `cells_sizes` stores the numbers of cells to sample. In our workflow we started at 250 cells, and sampled upwards geometrically, multiplying by a factor of square root of 2 (1.41) until reaching the size of the training set.   

```
    sampling_size = 250
    while sampling_size < n_retained_cells:
        cells_sizes.append(sampling_size)
        sampling_size = int(sampling_size*np.sqrt(2))
```
After scVI training, a t-SNE is computed in the embedding space and the coordinates are also saved. 
The results of each trained point are saved in a csv file named `scvi_output_{ds}<dataset>_c<number of sampled cells>_d<subsampled depth>.csv`, for example `heart10k_c1409_d200000.csv`. This naming convention is used by the final analysis notebook to parse the csv files and perform plotting.

The output files produced for this step in the paper in CaltechDATA are in the tar file `scvi_output_all_datasets.tar.gz`

### 4) Process scVI results and generate figure

The notebook `4) make_all_plots.ipynb` will use the `.csv` and `.h5ad` files to perform analysis and figure creation. Further comments are provided in the notebook itself.
