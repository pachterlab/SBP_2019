Code for producing the analysis in the "Quantifying the tradeoff between sequencing depth and cell number in single-cell RNA-seq" by Valentine Svensson, Eduardo Beltrame and Lior Pachter


The workflow has 4 steps:

### 1) FASTQ Subsampling and processing with kallisto bus 

Starting from the raw FASTQ files, subsample and process with kallisto bus to generate genecount matrices

### 2) Create `.H5AD` files with subsampled count matrices

For each dataset, parse all gene count matrices produces (in  .mtx format) into an h5 file using anndata 

### 3) Run scVI

### 4) Process scVI results and generate figure
