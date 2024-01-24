# RoCell

Notebooks associated with the publication:

## RoCell: An End-to-End Workflow for Single-Cell RNAseq Analysis

Guanya Yang, Jiawei Yan, Qing Zhang

_Abstract:_ 

These notebooks implement the algorithms and demonstrate how to apply to a common dataset.

** Usage

The __RoCell Workflow__ takes an `h5ad` file as the input, and also requires a few hyperparameters which we list below:

* `gene_id_col_label` is the name of column containing gene id information in the h5ad file.
* `min_var_pr` is the minimum variance of Pearson residual, a recommended value is 1.3.
* `min_mean_exp` is the minimum value of mean expression (count).
* `batch_label` is the column name for batches, this is dataset-specific.
* `hvg_min` is the minimum value for selecting highly-variable genes. a default value is set to be 2.2.

After installing RoCell, the workflow can be easily performed by a few lines of code:

```python
import RoCell
workflow = RoCell.Workflow(file = 'path_to_data/adata.h5ad')

workflow = RoCell.Workflow(file_name = '/data/autocell/mix_raw_fullinfo.h5ad')
workflow.preprocessing(gene_id_col_label = 'gene_ids')
workflow.compute_Pearson_residual(min_var_pr = 1.3, min_mean_exp = 0.1)
workflow.compute_biological_heterogeneity(batch_label = 'Species Call by SNV')
adata_hvg = workflow.select_hvg(hvg_min = 1.0)
adata_latent = workflow.dimension_reduction(adata_hvg, device='cpu')
```

which gives you the embedding of each cell (in the format of `adata`). The users can visualise the results through the standard UMAP algorithm, e.g.,

```python
sc.pl.umap(adata_latent,color=['ae'], legend_loc='on data', frameon=False, add_outline=True,legend_fontsize=8)
sc.pl.umap(adata_latent,color=['Cell_type'], frameon=False, legend_fontoutline=1, add_outline=True, legend_fontsize=6)
```

A runnable example can be found at `demo.ipynb`

