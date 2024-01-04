# RoCell

Notebooks associated with the publication:

## RoCell: An End-to-End Workflow for Single-Cell RNAseq Analysis

Guanya Yang, Jiawei Yan, Qing Zhang

_Abstract:_ 

These notebooks implement the algorithms and demonstrate how to apply to a common dataset.

** Usage

The __RoCell Workflow__ takes an `h5ad` file as the input, and also requires a few hyperparameters which we list below:

* `min_var_pr` is the minimum variance of Pearson residual, a recommended value is 1.3.
* `min_mean_exp` is the minimum value of mean expression (count).
* `batch_label` is the column name for batches, this is dataset-specific.
* `hvg_cutoff` is the minimum value for selecting highly-variable genes. a default value is set to be 2.2.

After installing RoCell, the workflow can be easily performed by a few lines of code:

```python
from RoCell import Workflow
workflow = Workflow(file = 'path_to_data/adata.h5ad')

workflow.compute_Pearson_residual(min_var_pr = 1.3, min_mean_exp = 0.0008)
genetic_variation = workflow.compute_features(batch_label = 'Batch')
adata_hvg = workflow.select_hvg(genetic_variation, hvg_cutoff = 2.2)
adata_latent = workflow.dimension_reduction(adata_hvg)
```

which gives you the embedding of each cell (in the format of `adata`). The users can visualise the results through the standard UMAP algorithm, e.g.,

```python
sc.pl.umap(adata_latent,color=['ae'], legend_loc='on data', frameon=False, add_outline=True,legend_fontsize=8)
sc.pl.umap(adata_latent,color=['Cell_type'], frameon=False, legend_fontoutline=1, add_outline=True, legend_fontsize=6)
```

A runnable example can be found at `demo_RoCell.ipynb`

