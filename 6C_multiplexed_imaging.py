import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.spatial import distance

import anndata as ad
import scanpy as sc
import SOAPy_st as sp

# Usage: https://github.com/KaWingLee9/in_house_tools/tree/main/multiplexed_images_pipeline#spatial-interaction-among-three-cell-types
exec(open('https://github.com/KaWingLee9/in_house_tools/blob/main/multiplexed_images_pipeline/multiplexed_image_processing.py').read())

# e.g. for one dataset
cell_info=pd.read_csv(path+'cell_info_ljr.csv',index_col=0)
adata=ad.AnnData(X=csr_matrix(pd.DataFrame([0 for i in range(0,cell_info.shape[0])])),
                 obs=cell_info)
adata.obsm['spatial']=np.array(cell_info[['location_x','location_y']])
x=adata.obs['image_id'].value_counts()

# analysis of multiple samples (parameters designed similar to SOAPy)
adata=relative_distance_analysis(adata,sample_key='image_id',cluster_key='cell_type',CT1='Macrophage',CT2='T cell',CT3='Tumor_cell',
                                 shuffle_type='shuffle_except1',quantatative_method='RN',shuffle_times=100,n_jobs=50)
adata.uns['RelativeDistance']['RN_result']

# subset one sample
adata_x=adata[:,adata.obs['image_id']==1]
adata_x=subtype_by_RelativeDistance(adata_x,CT1='Macrophage',CT2='T cell',CT3='Tumor_cell',cluster_key='cell_type')
adata_x.obs.loc[adata_x.obs['relative_analysis_subtype']=='CT1_CT2','relative_analysis_subtype']='M01'
adata_x.obs.loc[adata_x.obs['relative_analysis_subtype']=='CT1_CT3','relative_analysis_subtype']='M02'

fig1=plot_pixel(adata_x,color_panel,max_quantile=1,show_boundary=False,show_legend=False)

fig2=sc.pl.spatial(adata_x,color='cell_type_other',spot_size=8,
                  palette={'Macrophage':'#FF0000',
                           'T cell':'#0000FF',
                           'Tumor cell':'#00FF00'},
                  na_color='#E8E8E8',
                  show=False)

fig3=sc.pl.spatial(adata_x,color='relative_analysis_subtype',spot_size=8,
              palette={'Macrophage':'#FF000040',
                       'T cell':'#0000FF40',
                       'Tumor cell':'#00FF0040',
                       'M01':'#FFFF00FF',
                       'M02':'#FF00FFFF'},
              na_color='#E8E8E840',show=False)