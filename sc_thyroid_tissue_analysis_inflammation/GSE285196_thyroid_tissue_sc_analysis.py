import scanpy as sc
import squidpy as sq
import numpy as np # 2.0.2
import pandas as pd # 2.3.3
import matplotlib.pyplot as plt  
import os

# Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285196
# Paper: https://pubmed.ncbi.nlm.nih.gov/40710355/
#        Pathogenesis of Graves' Disease Determined Using Single-Cell Sequencing with Thyroid Autoantigen Peptide Stimulation in B Cells


# Base directory
base_path = "/Users/anna_gurina/Desktop/KSE_Bioinformatics/BioLab/datasets/thyroid_datasets_analysis_sandbox/sc_thyroid_tissue_analysis_inflammation/GSE285196_dataset"

# List all sample directories
sample_dirs = [
    "GSM8698011"
    # ,
    # "GSM8698012",
    # "GSM8698013",
    # "GSM8698014",
    # "GSM8698015",
    # "GSM8698016",
    # "GSM8698017",
    # "GSM8698018",
    # "GSM8698019",
    # "GSM8698020"
]

# Create a list to store AnnData objects

adatas = []
adata = []

for sample in sample_dirs:
    sample_path = os.path.join(base_path, sample)

    try:
        # Read the files separately
        matrix_file = f"{sample}_matrix.mtx.gz"
        features_file = f"{sample}_features.tsv.gz"
        barcodes_file = f"{sample}_barcodes.tsv.gz"
        
        # Read matrix
        adata = sc.read_mtx(os.path.join(sample_path, matrix_file)).transpose()  # (cells x genes)
        # print(f"!!!adata.shape: {adata.shape}")
        
        # Read genes (features)
        features = pd.read_csv(os.path.join(sample_path, features_file), 
                             header=None, sep='\t',
                             names=['gene_id', 'gene_symbol', 'feature_type'])
        # print(f"!!!features: {len(features)}")
        
        # Read barcodes and convert to strings
        barcodes = pd.read_csv(os.path.join(sample_path, barcodes_file), 
                              header=None, sep='\t',
                              names=['barcode'])
        barcodes['barcode'] = barcodes['barcode'].astype(str)
        # print(f"!!!barcodes: {len(barcodes)}")


        # Set gene names and barcodes
        adata.var_names = features['gene_symbol']  # Using gene symbols
        # print(f"!!!adata.var_names: {adata.var_names}")
        adata.var['gene_ids'] = features['gene_id'].values
        # print(f"!!!adata.gene_id: { adata.var['gene_ids']}")
        adata.var['feature_type'] = features['feature_type'].values
        # print(f"!!!adata.feature_type: {adata.var['feature_type']}")
        adata.obs_names = barcodes['barcode']
        # print(f"!!!adata.obs_names: {adata.obs_names}")
        
        # Add sample information
        adata.obs['sample'] = sample
        
        print(f"Successfully read {sample}")
        print(f"Shape: {adata.shape}")
        print(f"Number of genes: {len(adata.var_names)}")
        print(f"Number of cells: {len(adata.obs_names)}")
        
        adatas.append(adata)
        
    except Exception as e:
        print(f"Error reading {sample}: {str(e)}")

if adatas:
    # Concatenate all samples
    adata = adatas[0].concatenate(adatas[1:], join='outer')
    print("Combined data shape:", adata.shape)
else:
    print("No data was successfully loaded")



# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-') # mitochondrial genes percentage
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# print("Cell metrics:", adata.obs[['total_counts', 'n_genes_by_counts']])
# print("Gene metrics:", adata.var[['n_cells_by_counts', 'mean_counts']])

print(f"Available columns in adata.var: {adata.var.columns}")
adata.obs[['total_counts_mt','total_counts', 'n_genes_by_counts']].to_csv('figures/qc_metrics_per_cell.csv')
adata.var[['feature_type','mt','n_cells_by_counts', 'mean_counts','pct_dropout_by_counts','total_counts']].to_csv('figures/qc_metrics_per_gene.csv')


# Plot QC metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, show=False)
plt.savefig('figures/qc_violin_plot.png')
# plt.show()

## Filter cells
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# print(f"!!! Filtered adata.shape: {adata.shape}")

# Normalize data
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
print(f"!!! Normalised adata.shape: {adata.shape}")

# Find variable genes
sc.pp.highly_variable_genes(adata)
print(f"Number of highly variable genes: {adata.var['highly_variable'].sum()}")
var_genes_df = adata.var[adata.var['highly_variable']].copy() # Create DataFrame with variable genes information
var_genes_df['gene_name'] = var_genes_df.index # Add gene names as column
var_genes_df = var_genes_df.sort_values('dispersions_norm', ascending=False) # Sort by variance

var_genes_df.to_csv('figures/highly_variable_genes.csv')

print("Top 20 highly variable genes:")
print(var_genes_df[['dispersions_norm', 'means', 'dispersions']].head(20))

# Create visualization
sc.pl.highly_variable_genes(adata, save='.png', show=False)

plt.figure(figsize=(10, 6))
plt.scatter(var_genes_df['means'], var_genes_df['dispersions_norm'], 
           c='gray', alpha=0.5)
plt.scatter(var_genes_df[var_genes_df['highly_variable']]['means'],
           var_genes_df[var_genes_df['highly_variable']]['dispersions_norm'],
           c='red', alpha=0.5)
plt.xlabel('Mean expression')
plt.ylabel('Normalized dispersion')
plt.title('Highly Variable Genes')
plt.savefig('figures/variable_genes_scatter.png')
plt.close()



##################--PCA--##########################################################################################
sc.pp.pca(adata) # specify number of components ???

print("Variance ratio (first 10 PCs):")
print(adata.uns['pca']['variance_ratio'][:10])
print("\nTotal variance explained:", 
      np.sum(adata.uns['pca']['variance_ratio'][:10]))

# PCA coordunates
pca_coords = pd.DataFrame(adata.obsm['X_pca'], 
                         index=adata.obs_names,
                         columns=[f'PC{i+1}' for i in range(adata.obsm['X_pca'].shape[1])])
pca_coords.to_csv('figures/pca_coordinates.csv')

# Plot
sc.pl.pca(adata, 
          color='sample',  # you can color by any column in adata.obs
          components=['1,2', '3,4'],  # show first 4 components
          ncols=2,
          save='_components.png',
          show=False)

# PCA loadings (gene contributions)
loadings = pd.DataFrame(adata.varm['PCs'], 
                       index=adata.var_names,
                       columns=[f'PC{i+1}' for i in range(adata.varm['PCs'].shape[1])])

# Get top contributing genes for each PC
n_top = 10  # number of top genes to show
top_genes = {}
for pc in range(loadings.shape[1]):
    genes = loadings.iloc[:,pc].abs().sort_values(ascending=False)
    top_genes[f'PC{pc+1}'] = genes.head(n_top).index.tolist()

loadings.to_csv('figures/pca_loadings.csv')


##################--UMAP--##########################################################################################

sc.pp.neighbors(adata) # n_neighbors=15, n_pcs=30
sc.tl.umap(adata)

#UMAP coordinates
umap_coords = pd.DataFrame(adata.obsm['X_umap'], 
                          index=adata.obs_names,
                          columns=['UMAP1', 'UMAP2'])
umap_coords.to_csv('figures/umap_coordinates.csv')

# Basic UMAP plot by sample
sc.pl.umap(adata,
           color=['sample'],  # you can add more categories here
           frameon=False,
           save='_by_sample.png',
           show=False)

n_top_genes = 6  # get top variable genes to plot
top_genes = adata.var['dispersions_norm'].sort_values(ascending=False).index[:n_top_genes]

sc.pl.umap(adata,
           color=top_genes,
           ncols=3,
           frameon=False,
           save='_by_genes.png',
           show=False)

# UMAP with QC metrics
sc.pl.umap(adata,
           color=['n_genes_by_counts', 'total_counts'],
           ncols=2,
           frameon=False,
           save='_by_qc.png')

# save neighbor graph info
neighbor_info = pd.DataFrame({
    'n_neighbors': [adata.uns['neighbors']['params']['n_neighbors']],
    'method': [adata.uns['neighbors']['params']['method']],
    'metric': [adata.uns['neighbors']['params']['metric']],
    'random_state': [adata.uns['neighbors']['params']['random_state']]
})
neighbor_info.to_csv('figures/neighbor_graph_params.csv')

# combined plot with multiple features
fig, axes = plt.subplots(2, 2, figsize=(12, 12))
fig.suptitle('UMAP Visualization with Different Features')

sc.pl.umap(adata, color='sample', ax=axes[0,0], show=False, title='By Sample')
sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[0,1], show=False, title='Number of Genes')
sc.pl.umap(adata, color='total_counts', ax=axes[1,0], show=False, title='Total Counts')
sc.pl.umap(adata, color=top_genes[0], ax=axes[1,1], show=False, title=f'Expression of {top_genes[0]}')

plt.tight_layout()
plt.savefig('figures/umap_combined_features.png')
plt.close()

## Print summary statistics
# print("\nUMAP Summary:")
# print(f"Number of cells: {adata.n_obs}")
# print(f"UMAP coordinates shape: {adata.obsm['X_umap'].shape}")
# print("\nNeighbor graph parameters:")
# print(f"n_neighbors: {adata.uns['neighbors']['params']['n_neighbors']}")
# print(f"method: {adata.uns['neighbors']['params']['method']}")
# print(f"metric: {adata.uns['neighbors']['params']['metric']}")
# print(f"random_state: {adata.uns['neighbors']['params']['random_state']}")

##################--Leiden clustering--##########################################################################################
sc.tl.leiden(adata)

# UMAP plot colored by leiden clusters
sc.pl.umap(adata, 
           color='leiden',
           legend_loc='on data',
           frameon=False,
           save='_leiden_clusters.png')

# save cluster info
cluster_info = pd.DataFrame({
    'Cluster': adata.obs['leiden'],
    'Sample': adata.obs['sample']
})
cluster_info.to_csv('figures/leiden_clusters.csv')

# Cluster statistics
cluster_stats = pd.DataFrame({
    'Cluster_size': adata.obs['leiden'].value_counts(),
    'Percentage': (adata.obs['leiden'].value_counts() / len(adata.obs) * 100).round(2)
})
cluster_stats.to_csv('figures/leiden_cluster_stats.csv')

# Summary
print("\nLeiden clustering summary:")
print(f"Number of clusters: {len(adata.obs['leiden'].unique())}")
print("\nCluster sizes:")
print(cluster_stats)

# Cluster sizes
plt.figure(figsize=(10, 5))
cluster_stats['Cluster_size'].plot(kind='bar')
plt.title('Cells per Cluster')
plt.xlabel('Cluster')
plt.ylabel('Number of Cells')
plt.tight_layout()
plt.savefig('figures/leiden_cluster_sizes.png')
plt.close()

# Find marker genes for each cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_markers.png')

# Save marker genes
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
marker_genes = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'scores', 'pvals_adj']})
marker_genes.to_csv('figures/leiden_marker_genes.csv')


