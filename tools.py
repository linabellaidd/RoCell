import sys
sys.path.append('/data/singlecell/pipline/')
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
from scipy.optimize import leastsq
from sklearn.metrics import mean_absolute_error 
from scipy.special import comb
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from model import AutoEncoder


def dim_reduction_ae(adata_ae):
    lr = 0.06
    decay_rate = 0.995
    batch_size = 64
    max_iter = 30
    n_neighbors = 20
    resolution = 1.2
    cuda_num = 6
    device = torch.device("cuda: {}".format(cuda_num))
    adata_X = adata_ae.X
    torch.manual_seed(3407)
    AE = AutoEncoder([adata_X.shape[1],64]).to(device)
    train = True
    loss_df = []
    if train:
        for param in AE.parameters():
            if param.dim() > 1:
                nn.init.xavier_uniform_(param)

        optimizer = optim.RMSprop(AE.parameters(), lr=lr)
        lambda_epoch = lambda i_iter: lr * (decay_rate ** i_iter)
        scheduler = optim.lr_scheduler.LambdaLR( optimizer, lr_lambda=lambda_epoch)

        for i_iter in range(max_iter):
            adata_index_set = np.arange(len(adata_X))
            loss_average, k = 0, 0
            while len(adata_index_set):
                batch_index = np.random.choice(adata_index_set, batch_size)
                adata_index_set = list(set(adata_index_set) - set(batch_index))
                input_adata = torch.tensor(adata_X[batch_index], device=device).float()
                recon_adata, latent_space = AE(input_adata)
                loss = F.mse_loss(recon_adata, input_adata)
                k += 1
                loss_average += loss.item()
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
            scheduler.step()
            loss_df.append([i_iter,loss_average/k])
            if i_iter % int(max_iter/10)==0 and i_iter:
                print(i_iter, loss_average/k)
    AE.eval()
    adata_index_set = np.arange(len(adata_X))
    batch_num = int(np.ceil(len(adata_X) / batch_size))
    loss_average, k = 0, 0
    X_latent = []
    for i in range(batch_num):
        s_index = i*batch_size
        e_index = min((i+1)*batch_size, len(adata_X))
        input_adata = torch.tensor(adata_X[s_index: e_index], device=device).float()
        recon_adata, latent_space = AE(input_adata)
        X_latent.extend(latent_space.tolist())
        loss = F.mse_loss(recon_adata, input_adata)
        k += 1
        loss_average += loss.item()
    X_latent = np.array(X_latent)
    adata_latent = sc.AnnData(X_latent)
    adata_latent.obs = adata_ae.obs
    sc.pp.neighbors(adata_latent, n_neighbors=n_neighbors, use_rep="X")
    sc.tl.leiden(adata_latent, resolution=resolution)
    Y_pred_init = adata_latent.obs['leiden']
    df_latent = pd.DataFrame(X_latent, index=np.arange(len(adata_ae)))
    df_latent['cluster'] = Y_pred_init.values
    adata_latent.obsm["X_Embeded_z"] = X_latent
    sc.pp.neighbors(adata_latent, n_neighbors=n_neighbors, use_rep="X_Embeded_z") #ae latend space
    sc.tl.umap(adata_latent,n_components=3)
    adata_latent.obs['autoencoder'] = np.array(Y_pred_init)
    print('------------------------------training is done!----------------------------')
    pd.DataFrame(loss_df).iloc[:,1].plot(figsize=(2,2))
    return adata_latent



# Normalization
def analytic_pearson_residuals(counts, theta):
    if type(counts) == np.ndarray:
        counts = scipy.sparse.csr_matrix(counts)
    counts_sum0 = np.sum(counts, axis=0)
    counts_sum1 = np.sum(counts, axis=1)
    counts_sum  = np.sum(counts) # pg
    mu = counts_sum1 @ counts_sum0 / counts_sum
    z = (counts - mu) / np.sqrt(mu + np.power(mu,2)/theta)
    n = counts.shape[0]
    z[z >  np.sqrt(n)] =  np.sqrt(n)
    z[z < -np.sqrt(n)] = -np.sqrt(n)
    
    return z




def filter (adata,
            max_genes = None,
            min_genes = None,
            max_counts = None,
            min_counts = None,
            min_cells = None,
            max_mito = None,
            annotation_type = None,
            species = 'human'):
    if isinstance(adata, anndata.AnnData):

        ncells = adata.shape[0]
        ngenes = adata.shape[1]
        print('started with ', str(ncells), ' total cells and ', str(ngenes), ' total genes')
        if (max_counts is not None or min_counts is not None):
            if adata.obs.get('n_counts') is None:
                adata.obs['n_counts'] = adata.X.sum(axis=1)        
        if (min_genes is not None or max_genes is not None):
            if adata.obs.get('n_genes') is None:
                adata.obs['n_genes'] = sum(adata.X > 0, axis=1)
        if min_cells is not None:
            if adata.var.get('n_cells') is None:
                adata.var['n_cells'] = sum(adata.X > 0, axis=0).T

        if max_mito is not None:
            if adata.obs.get('percent_mito') is None:
                mito_list = get_mito_genes(species, annotation_type)
                mito_genes = [name for name in adata.var_names if name in mito_list] # ensembl
                # for each cell compute fraction of counts in mito genes vs. all genes
                n_counts = sum(adata.X, axis=1).A1
                n_counts[n_counts == 0] = float('inf')
                adata.obs['percent_mito'] = sum(adata[:, mito_genes].X, axis=1).A1 / n_counts

        #perform actual filtering for all given parameters
        if max_genes is not None:
            curr_cells = adata.shape[0]
            adata = adata[adata.obs.get('n_genes') <= max_genes,:].copy()
            new_cells = adata.shape[0]
            print('removed', str(curr_cells - new_cells), 'cells that expressed more than', str(max_genes), 'genes')        
        if min_genes is not None:
            curr_cells = adata.shape[0]
            adata = adata[adata.obs.get('n_genes') >= min_genes,:].copy()
            new_cells = adata.shape[0]
            print('removed', str(curr_cells - new_cells), 'cells that did not express at least', str(min_genes), ' genes')        
        if max_counts is not None:
            curr_cells = adata.shape[0]
            adata = adata[adata.obs.get('n_counts') <= max_counts,:].copy()
            new_cells = adata.shape[0]
            print('removed', str(curr_cells - new_cells), 'cells that had more than', str(max_counts), ' counts')
        
        if min_counts is not None:
            curr_cells = adata.shape[0]
            adata = adata[adata.obs.get('n_counts') >= min_counts,:].copy()
            new_cells = adata.shape[0]
            print('removed', str(curr_cells - new_cells), 'cells that did not have at least', str(min_counts), 'counts')

        if min_cells is not None:
            curr_genes = adata.shape[1]
            adata = adata[:,adata.var.get('n_cells') >= min_cells].copy()
            new_genes = adata.shape[1]
            print('removed', str(curr_genes - new_genes), 'genes that were not expressed in at least', str(min_cells), 'cells')

        if max_mito is not None:
            curr_cells = adata.shape[0]
            adata = adata[adata.obs.get('percent_mito') < max_mito,:].copy()
            new_cells = adata.shape[0]
            print('removed ', str(curr_cells - new_cells), ' cells that expressed ', str(max_mito * 100), 'percent mitochondrial genes or more')
        
        ncells_final = adata.shape[0]
        ngenes_final = adata.shape[1]
        print('finished with', str(ncells_final), ' total cells and', str(ngenes_final), 'total genes')
        
        return adata


def get_nb_parameter_df(batch_x,rbg):
    
    def fitfunc(p,x):
        P_res = []
        for x_i in x:
            p[1] = int(p[1])
            p1 = comb(x_i + p[1] - 1, p[1] - 1) #p[1] is theta p[0] is p
            p2 = p[0]**p[1]
            p3 = (1-p[0])**x_i
            p4 = p[2]
            P = p1*p2*p3*p4
            P_res.append(P)
        return np.asarray(P_res)


    def errfunc(p,x,y):
        err = y - fitfunc(p,x)
        return err

    def calcualt_mu_lambda(x,freq,theta):
        xdata1, ydata1 = x,freq
        init  = [0.01,theta,1e5]
        out1 = leastsq(errfunc, init, args=(xdata1, ydata1))
        mu1 = (1 - out1[0][0])/out1[0][0] * out1[0][1]
        lambda1 = out1[0][2]
        p1 = out1[0][0]
        theta1 = out1[0][1]
        loss = mean_absolute_error(fitfunc(out1[0], xdata1),ydata1)
        return mu1,lambda1,p1,theta1,loss

    para_lis = []
    for rbg_i in rbg:
        x,freq = np.unique(batch_x[:,rbg_i].reshape(-1,1), return_counts=True)
        if len(x) <= 3:
            best_para = [0,0,0,0,0,rbg_i]
        elif len(x) > 3:
            res_lis = []
            n_theta = 20

            for theta_ in range(1,n_theta):
                mu,lamb,p,theta,loss = calcualt_mu_lambda(x,freq,theta_)
                iter_res = [mu,lamb,p,theta,loss]
                res_lis.append(iter_res)
            df_res = pd.DataFrame(res_lis)
            df_res.columns =["mu","lambda","p","theta","loss"]
            df_res['rbg_id'] = rbg_i
            _mu,_lambda,_p,_theta,_loss, _rbg_id = df_res.iloc[df_res["loss"].sort_values().index[0],:].values   
            best_para = df_res.iloc[df_res["loss"].sort_values().index[0],:].values

        para_lis.append(best_para)
    
    df = pd.DataFrame(para_lis)
    df.columns =["mu","lambda","p","theta","loss",'gene_id']
    df['Kurtosis'] = (6/df.theta) + df.theta*((1-df.p)**2/df.p)
    df['Skewness'] = (1 + df.p) / np.sqrt(df.p * df.theta)
    return df