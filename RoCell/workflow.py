import os, scipy
import pandas as pd
import numpy as np
import scanpy as sc
from itertools import combinations
from multiprocessing import Pool
from tqdm.auto import tqdm
from scipy.stats import kstest, nbinom
from sklearn.metrics.pairwise import cosine_similarity

from anndata import read_h5ad
from RoCell import utils

class Workflow:
    def __init__(self, file_name: str):
        if file_name.endswith('h5ad') == False:
            raise TypeError('RoCell only supports h5ad file')
        else:
            self.adata_raw = read_h5ad(file_name)
        
    def preprocessing(self, gene_id_col_label: str = 'gene_ids', **qc_kwargs):
        self.gene_id_col_label = gene_id_col_label
        
        mito_genes = self.adata_raw.var[gene_id_col_label].str.startswith('MT-')
        self.adata_raw.obs['n_genes'] = (self.adata_raw.X > 0).sum(axis=1)
        self.adata_raw.var['n_cells'] = (self.adata_raw.X > 0).sum(axis=0).T
        if type(self.adata_raw.X) == np.ndarray:
            self.adata_raw.obs['percent_mito'] = np.sum(
                self.adata_raw[:, mito_genes].X, axis=1) / np.sum(self.adata_raw.X, axis=1)
            self.adata_raw.obs['n_counts'] = self.adata_raw.X.sum(axis=1)
        else:
            self.adata_raw.obs['percent_mito'] = np.sum(self.adata_raw[:, mito_genes].X, axis=1).A1 / np.sum(self.adata_raw.X, axis=1).A1
            self.adata_raw.obs['n_counts'] = self.adata_raw.X.sum(axis=1).A1

        self.adata = utils.qc(adata = self.adata_raw, **qc_kwargs)
        
    def compute_Pearson_residual(self, min_var_pr: float = 1.3, min_mean_exp: float = 0.0008):
        self.adata.var['qc_gene_idx'] = list(range(len(self.adata.var)))
        z = utils.analytic_pearson_residuals(scipy.sparse.csr_matrix(self.adata.X), 100)
        self.adata.var['pr_var'] = np.squeeze(np.asarray(np.var(z, axis = 0))) 
        self.adata.var['mean_expression'] = np.squeeze(np.asarray(np.mean(self.adata.X, axis = 0)))
        self.adata_apr = self.adata.copy()
        self.adata_apr.X = np.asarray(z)
        self.aprg_mask = (self.adata.var.pr_var > min_var_pr) & (self.adata.var.mean_expression > min_mean_exp)
    
    def _get_batch_feature_vectors(self, batch_id):
        cell_filter = self.adata.obs[self.adata.obs[self.batch_label] == batch_id].qc_cell_idx
        batch_x = self.adata.X[cell_filter][:,self.aprg_mask].toarray()
        
        var_h = utils.nb_fit(batch_x.mean(0), batch_x.var(0))
        r, p = utils.nb_p(batch_x.mean(0), var_h)

        feature_vectors= pd.DataFrame([kstest(batch_x[:,i], nbinom.cdf, args=(r[i],p[i]), method='exact').statistic for i in range(len(r))],
                                        columns=['loss'])
        feature_vectors['skewness'] = scipy.stats.nbinom.stats(r, p, moments='s')
        feature_vectors['theta'] = r
            
        feature_vectors.to_excel(self.path_batch_effect_feature.format(batch_id), index=False)
        
        return True
                
    def _get_batch_cosine_similarity(self, pair):
        sample_id_a, sample_id_b = pair
        sample_id_a, sample_id_b = sample_id_a.replace('/',':'), sample_id_b.replace('/',':')
        
        file_name = os.path.join(self.path_cos_similarity.format(sample_id_a, sample_id_b))
        if os.path.isfile(file_name) == False:
            df_a = pd.read_excel(self.path_batch_effect_feature.format(sample_id_a))
            df_b = pd.read_excel(self.path_batch_effect_feature.format(sample_id_b))
            
            cs = pd.DataFrame(cosine_similarity(df_a, df_b).diagonal(),
                              columns=['{}_{}_cos_similarity'.format(sample_id_a, sample_id_b)])
            cs.to_excel(file_name, index=False)
        return True
                
    def compute_biological_heterogeneity(self, batch_label: str = 'Batch', num_cores: int = 36):
        self.batch_label = batch_label
        self.path_batch_effect_feature = './{}/batch_effect_feature'.format(batch_label)
        if not os.path.exists(self.path_batch_effect_feature):
            os.makedirs(self.path_batch_effect_feature)
        self.path_batch_effect_feature = self.path_batch_effect_feature + '/{}.xlsx'
        
        self.path_cos_similarity  = './{}/cos_similarity'.format(batch_label)
        if not os.path.exists(self.path_cos_similarity):
            os.makedirs(self.path_cos_similarity)
        self.path_cos_similarity = self.path_cos_similarity + '/{}_vs_{}.xlsx'
        
        self.adata.obs['qc_cell_idx'] = list(range(len(self.adata.obs)))
        
        
            
        batches = sorted(list(self.adata.obs[batch_label].value_counts()[self.adata.obs[batch_label].value_counts() > 0].index))
        pairs = list(combinations(batches, 2))
        
        if 'qc_cell_idx' not in self.adata.obs.columns:
            raise ValueError('cell index is not define!')
        else:
            print('compute batch-specific feature vector')
            with Pool(processes = num_cores) as pool:
                o = list(tqdm(pool.imap_unordered(self._get_batch_feature_vectors, batches), total = len(batches))) 
            
            print('compute cosine similarity')
            with Pool(processes = num_cores) as pool:
                o = list(tqdm(pool.imap_unordered(self._get_batch_cosine_similarity, pairs), total = len(pairs))) 

        cs_dflist = []
        for i, (sample_a, sample_b) in enumerate(pairs):
            cs_dflist.append(pd.read_excel(self.path_cos_similarity.format(sample_a, sample_b)))
        
        df_cs = pd.concat(cs_dflist,axis=1)
        df_cs.index = self.adata[:,self.aprg_mask].var[self.gene_id_col_label].values
        df_rbg = df_cs.loc[df_cs.index.str.startswith('RPS') | df_cs.index.str.startswith('RPL')]        
        relative_cs = (df_cs  - df_rbg.mean()).fillna(0).abs()
        
        self.biological_heterogeneity = pd.DataFrame(relative_cs.sum(axis=1), columns = ['var'])
        self.biological_heterogeneity['qc_gene_idx'] = self.adata[:,self.aprg_mask].var.qc_gene_idx.values
        
    def select_hvg(self, hvg_min: float = 2.2):
        adata_hvg = self.adata_apr[:, self.biological_heterogeneity[self.biological_heterogeneity['var'] > hvg_min].qc_gene_idx]
        return adata_hvg
        
    def dimension_reduction(self, adata_hvg, device:str = 'cpu'):
        adata_hvg = sc.pp.scale(adata_hvg, zero_center=True, max_value=3,copy=True)
        adata_latent = utils.dim_reduction_ae(adata_hvg, device)
        self.adata.obs['ae'] = adata_latent.obs.autoencoder.values
        
        return adata_latent
