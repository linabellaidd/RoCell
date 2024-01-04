import os
import pandas as pd
import scanpy as sc
from itertools import combinations
from multiprocessing import Pool

from anndata import read_h5ad
from RoCell import utils

def load_preprocessing(file_name: str, 
                       min_genes: int, 
                       min_counts: int, 
                       min_cells: int):
    
    if file_name.endswith('h5ad') == False:
        raise TypeError('RoCell only supports h5ad file')
    
    adata_raw = read_h5ad(file_name)
    mito_genes = adata_raw.var['gene_ids'].str.startswith('MT-')
    adata_raw.obs['n_genes'] = sum(adata_raw.X > 0, axis=1)
    adata_raw.var['n_cells'] = sum(adata_raw.X > 0, axis=0).T
    if type(adata_raw.X) == np.ndarray:
        adata_raw.obs['percent_mito'] = np.sum(
            adata_raw[:, mito_genes].X, axis=1) / np.sum(adata_raw.X, axis=1)
        adata_raw.obs['n_counts'] = adata_raw.X.sum(axis=1)
    else:
        adata_raw.obs['percent_mito'] = np.sum(adata_raw[:, mito_genes].X, axis=1).A1 / np.sum(adata_raw.X, axis=1).A1
        adata_raw.obs['n_counts'] = adata_raw.X.sum(axis=1).A1

    adata = utils.qc(adata_raw, min_genes, min_counts, min_cells)

    return adata


class Workflow:
    def __init__(self, 
                 file_name: str,
                 min_genes: int = 200, 
                 min_counts: int = 1_000, 
                 min_cells: int = 2):
        self.adata = load_preprocessing(file_name, min_genes, min_counts, min_cells)
        
    def compute_Pearson_residual(self, min_var_pr: float = 1.3, min_mean_exp: float = 0.0008):
        self.adata.var['qc_gene_idx'] = list(range(len(self.adata.var)))
        z = utils.analytic_pearson_residuals(scipy.sparse.csr_matrix(self.adata.X), 100)
        self.adata.var['pr_var'] = np.squeeze(np.asarray(np.var(z, axis = 0))) 
        means_exp = np.mean(adata.X, axis = 0)
        self.adata.var['mean_expression'] = np.squeeze(np.asarray(means_exp))
        self.adata_apr = self.adata.copy()
        self.adata.X = z
        self.aprg_mask = (self.adata.var.pr_var > min_var_pr) & (self.adata.var.mean_expression > min_mean_exp)
        
    
    def _get_batch_expression_descriptors(self, batch):
            adata_X = self.adata.X
            batch_id_a = batch
            if 'qc_cell_idx' not in self.adata.obs.columns:
                print('cell index is not define!')
            else:
                df_a = utils.get_nb_parameter_df(adata_X[self.adata.obs.loc[adata.obs[batch_label] == batch_id_a].qc_cell_idx,:], self.aprg_mask)
                df_a.to_excel(os.path.join(path,'{}.xlsx'.format(batch_id_a)))
                print('excel file of batch {} has saved!'.format(batch_id_a))
                
    def _def get_batch_cos_sim(pair):
        sample_id_a, sample_id_b = pair
        df_a = pd.read_excel(os.path.join(path, sample_id_a + '.xlsx'))
        df_b = pd.read_excel(os.path.join(path, sample_id_b + '.xlsx'))
        df_ab = utils.gene_batch_effect_profiling(df_a, df_b, sample_id_a, sample_id_b)
        df_ab.to_excel(os.path.join(path_,sample_id_a+'_vs_'+ sample_id_b +'.xlsx'))
        print('calculation of cos between {} and {} is done and saved!'.format(sample_id_a, sample_id_b))


                
    def compute_features(self, batch_label: str = 'Batch', num_cores: int = 36):
        self.adata.obs['qc_cell_idx'] = list(range(len(self.adata.obs)))
        
        path = './{}/batch_effect_feature'.format(batch_label)
        if not os.path.exists(path):
            os.makedirs(path)
            
        batch = np.unique(self.adata.obs[batch_label])
        pairs = [[p[0], p[1]] for p in combinations(batch, 2)]

        pool = Pool(num_cores)                        
        pool.map(self.get_batch_expression_descriptors, batch) 
        pool.close()
        pool.join()
        
        
        path_  = './{}/cos_similarity'.format(batch)
        if not os.path.exists(path_):
            os.makedirs(path_)
        
        pool = Pool(num_cores)                      
        pool.map(self._get_batch_cos_sim, pairs) 
        pool.close()
        pool.join()

        cs_dflist = []
        corr_dflist = []
        for i, (sample_a, sample_b) in enumerate(pair):
            filename = sample_a + '_vs_' + sample_b + '.xlsx'
            df_ = pd.read_excel(os.path.join(path_,filename))
            cs_dflist.append(df_.iloc[:,1])
            corr_dflist.append(df_.iloc[:,2])
        df_cs = pd.concat(cs_dflist,axis=1)
        df_cs.index = self.adata[:,self.aprg_mask].var.gene_ids.values
        df_rbg = df_cs.loc[df_cs.index.str.startswith('RPS') | df_cs.index.str.startswith('RPL')]        
        relative_cs = (df_cs  - df_rbg.mean()).fillna(0).abs()
        
        variability = pd.DataFrame(relative_cs.sum(axis=1))
        return variability
        
    def select_hvg(self, variability: pd.DataFrame, hvg_cutoff: float = 2.2):
                
        hvg_df['qc_gene_idx'] = self.adata[:,self.aprg_mask].var.qc_gene_idx.values
        adata_hvg = self.adata_apr[:,hvg_df[hvg_df.bio_intensity > 2.2].qc_gene_idx]
        return adata_hvg
        
    def dimension_reduction(self, adata_hvg):
        adata_hvg = sc.pp.scale(adata_hvg, zero_center=True, max_value=3,copy=True)
        adata_latent = utils.dim_reduction_ae(adata_hvg)
        adata.obs['ae'] = adata_latent.obs.autoencoder.values
        
        return adata_latent
