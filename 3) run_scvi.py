import numpy as np
import sklearn
import pandas as pd
import plotnine as p
from glob import glob
from tqdm import tqdm
import anndata
import torch
from scvi.dataset import GeneExpressionDataset
from scvi.models import VAE
from scvi.inference import UnsupervisedTrainer
from scvi.inference.posterior import Posterior
from openTSNE import TSNE
from openTSNE.callbacks import ErrorLogger
import os

adata_files= {
# 'neurons10k' :'./neurons10k_subsamples.h5ad',     
'heart10k' :'./heart10k_subsamples.h5ad',
'pbmc10k' :'./pbmc10k_subsamples.h5ad',    
}

for ds in adata_files:

    print('============== PROCESSING DATASET: ', ds, '===============')

    adata = anndata.read(adata_files[ds])


    total_cells = adata.n_obs
    n_retained_cells = int(0.85*total_cells)
    print('Total cells:', total_cells)
    print('Retained cells:', n_retained_cells)


    cells_sizes = []
    sampling_size = 250
    while sampling_size < n_retained_cells:
        cells_sizes.append(sampling_size)
        sampling_size = int(sampling_size*np.sqrt(2))

    # cells_sizes = np.logspace(np.log2(500), np.log2(n_retained_cells), num=9, base=2).astype(int)
    print('Number of sampled cells for ', ds, cells_sizes)


    cells_dataset = GeneExpressionDataset()
    X_ = adata.layers['0']
    cells_dataset.populate_from_data(X_, gene_names=adata.var.index.values)

    cells_dataset.subsample_genes(1000)
    sel_genes = cells_dataset.gene_names

    n_validation = adata.shape[0] - n_retained_cells
    print(ds, ' n_validation:', n_validation)

    validation_cells = np.random.choice(adata.obs.index, size=n_validation, replace=False)
    learning_cells = adata.obs.index.difference(validation_cells)

    val_adata = adata[validation_cells]
    lea_adata = adata[learning_cells]

    ne_cells = X_.sum(axis=1) > 0
    to_keep = np.where(ne_cells)[0]

    log_counts = np.log(X_[to_keep].sum(axis=1))
    local_means = (np.mean(log_counts) * np.ones((X_.shape[0], 1))).astype(np.float32)
    local_vars = (np.var(log_counts) * np.ones((X_.shape[0], 1))).astype(np.float32)

    batch_indices = np.ones((X_.shape[0], 1))
    labels = np.zeros_like(batch_indices)

    validation_cells_dataset = GeneExpressionDataset()
    X_ = val_adata[:, sel_genes].layers['0']
    validation_cells_dataset.populate_from_data(X_)

    idx = np.random.choice(learning_cells, size=n_retained_cells, replace=False)
    
    # this will store all dataframes we also save as csv
    results_list = []

    for ss_cells in cells_sizes:
        print('                🍞  🥯  🧁  🍰  🥮  🥞  🥧  🎂  🍯  🥪  🥖  🥐  🥟   ')
        print('======== ', ds, ' ========  SUBSAMPLED CELLS:', ss_cells, ' ================== ')
        idx = np.random.choice(learning_cells, size=ss_cells, replace=False)
        tmp_adata = lea_adata[idx]
        n_epochs = int(27 * adata.shape[0] / ss_cells)

        for ss_depth in tmp_adata.layers.keys():
            print(ds, 'now running depth:', ss_depth)


            result_dict = {'ss_depth': ss_depth, 'ss_cells': ss_cells}

            X_ = tmp_adata[:, sel_genes].layers[ss_depth]

            ne_cells = X_.sum(axis=1) > 0
            to_keep = np.where(ne_cells)[0]

            log_counts = np.log(X_[to_keep].sum(axis=1))
            local_means = (np.mean(log_counts) * np.ones((X_.shape[0], 1))).astype(np.float32)
            local_vars = (np.var(log_counts) * np.ones((X_.shape[0], 1))).astype(np.float32)

            batch_indices = np.ones((X_.shape[0], 1))
            labels = np.zeros_like(batch_indices)

            cells_dataset = GeneExpressionDataset()
            cells_dataset.populate_from_data(X_, gene_names=sel_genes)

            vae = VAE(cells_dataset.nb_genes,
                      dispersion='gene',
                      n_layers=2,
                      n_hidden=128,
                      reconstruction_loss='nb')

            trainer = UnsupervisedTrainer(vae, cells_dataset)

            trainer.train(n_epochs=n_epochs)

            validation_posterior = Posterior(vae, validation_cells_dataset, use_cuda=False)
            print(X_.shape)
            result_dict['validation_error'] = validation_posterior.reconstruction_error()

            # Get expression rate parameterization from representation space
            Z_hat = vae.sample_from_posterior_z(torch.from_numpy(cells_dataset.X.toarray()).float())
            Z_hat = np.array(Z_hat.detach()).astype(np.double)

            tsne = TSNE(callbacks=ErrorLogger(),
                        initialization='random',
                        negative_gradient_method='fft',
                        callbacks_every_iters=100,
                        neighbors='approx')

            YY = tsne.fit(Z_hat)

            df = pd.DataFrame(index=tmp_adata.obs.index)
            df['ss_depth'] = result_dict['ss_depth']
            df['ss_cells'] = result_dict['ss_cells']
            df['validation_error'] = result_dict['validation_error']
            df['tsne_0'] = YY[:, 0]
            df['tsne_1'] = YY[:, 1]

            out_file = f'scvi_output_{ds}/{ds}_c{ss_cells}_d{ss_depth}.csv'
            if not os.path.exists(os.path.dirname(out_file)):
                os.makedirs(os.path.dirname(out_file))
                
            df.to_csv(out_file)
            results_list.append(df)
            print(f'Saved: {out_file}')

    # combines all separate depths into a single csv file
    all_results = pd.concat(results_list).reset_index()
    all_results.to_csv(ds + '_all_scvi_results.csv')

