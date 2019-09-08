Code for producing the analysis in the "Quantifying the tradeoff between sequencing depth and cell number in single-cell RNA-seq" by Valentine Svensson, Eduardo Beltrame and Lior Pachter


The workflow has 4 steps:

### 1) FASTQ Subsampling and processing with kallisto bus 

Starting from the raw FASTQ files, subsample and process with kallisto bus to generate genecount matrices

### 2) Create `.H5AD` files with subsampled count matrices

For each dataset, parse all gene count matrices produces (in  .mtx format) into an h5 file using anndata 

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

### 4) Process scVI results and generate figure
