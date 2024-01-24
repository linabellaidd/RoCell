# RoCell

These notebooks implement the algorithms associated with the publication:

## RoCell: An End-to-End Workflow for Single-Cell RNAseq Analysis

Guanya Yang, Jiawei Yan, Qing Zhang

_Abstract:_ 

Single-cell RNA sequencing plays a crucial role in drug discovery by identifying elusive cell subpopulations, elucidating cell-specific gene expression signatures, and revealing novel cellular phenotypes, thereby enhancing our understanding of pathogenesis of complex diseases and facilitating the development of precision therapies. However, technical noise in scRNA-seq data, stemming from stochastic sampling, amplification bias, dropout events, batch effects, obscure the underlying cellular heterogeneity and introduce biases, all of which then aggravate the reliability and accuracy of downstream analysis. To eliminate the technical noise for better interrogating the intricate biological signals present in genes across distinct cell types, we propose an end-to-end AI solution termed as the Roche single cell RNA-seq platform, a comprehensive scRNA-seq data analysis workflow that seamlessly integrates all steps in scRNA-seq analysis, from quality control, to dimensionality reduction and cell clustering. Particularly, we built an autoencoder-based deep learning model to learn a compressed representation of the data. To reduce the bias from manual feature selection and errors in the batch effect removal step, we proposed to quantify the heterogeneity of gene expressions by a precise negative binomial model-based scoring algorithm, and select genes with higher heterogeneity score as the input to the model. To validate the effectiveness of our method, RoCell is benchmarked on multiple public datasets and shows superior accuracy compared to existing algorithms and remains robust across different types of datasets. Overall, our solution addresses the challenge of technical noise in scRNA-seq data analysis, providing a robust method for uncovering cellular heterogeneity and facilitating downstream analysis.

### Usage:

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

