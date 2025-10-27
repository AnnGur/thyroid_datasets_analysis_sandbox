import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests 
from scipy import stats
from sklearn.decomposition import PCA
import scanpy as sc
from typing import Dict

from src.config import AnalysisConfig

class ThyroidAnalysis:
    """Main class for thyroid expression data analysis"""
    
    def __init__(self, config: AnalysisConfig):
        self.config = config
        self._setup_directories()
        self.sample_info = None
        self.expr_matrix = None
        self.adata = None
    
    def _setup_directories(self) -> None:
        """Create necessary directories"""
        os.makedirs(self.config.output_dir, exist_ok=True)
        for subdir in ['qc', 'pca', 'differential_expression']:
            os.makedirs(os.path.join(self.config.output_dir, subdir), exist_ok=True)
    
    def load_data(self, series_matrix_file: str, sample_dir: str) -> None:
        """Load and prepare all necessary data"""
        print(f"Loading series matrix from: {series_matrix_file}")
        if not os.path.exists(series_matrix_file):
            raise FileNotFoundError(f"Series matrix file not found: {series_matrix_file}")
        
        if not os.path.exists(sample_dir):
            raise FileNotFoundError(f"Sample directory not found: {sample_dir}")
        
        try:
            self._parse_series_matrix(series_matrix_file)
            print(f"Parsed series matrix, found {len(self.sample_info)} samples")
            
            self._read_expression_data(sample_dir)
            print("Expression data loaded successfully")
            
            self._create_anndata()
            print("AnnData object created successfully")
        except Exception as e:
            raise Exception(f"Error during data loading: {str(e)}")
    
    def run_analysis(self) -> Dict:
        """Run complete analysis pipeline"""
        results = {
            'qc': self._run_qc(),
            'pca': self._run_pca(),
            'differential_expression': self._run_differential_expression()
        }
        return results
    
    def _parse_series_matrix(self, file_path: str) -> None:
        """Parse series matrix file"""
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        data = {
            'sample_titles': None,
            'series_sample_id': None,
            'braf_status': None,
            'tissue_type': None
        }
        
        for line in lines:
            if line.startswith('!Sample_title'):
                data['sample_titles'] = line.strip().split('\t')[1:]
            elif line.startswith('!Series_sample_id'):
                # Clean and split the Series_sample_id line
                sample_ids_str = line.replace('!Series_sample_id', '').strip().strip('"')
                data['series_sample_id'] = [x.strip() for x in sample_ids_str.split()]
            elif line.startswith('!Sample_characteristics_ch1'):
                if 'braf status:' in line:
                    data['braf_status'] = [x.split(': ')[1].strip('"') 
                                        for x in line.strip().split('\t')[1:]]
                elif 'tissue:' in line:
                    data['tissue_type'] = [x.split(': ')[1].strip('"') 
                                        for x in line.strip().split('\t')[1:]]
        
        # Print debug information
        print(f"Found {len(data['series_sample_id'])} sample IDs")
        print(f"Found {len(data['sample_titles'])} sample titles")
        print(f"Found {len(data['braf_status'])} BRAF status entries")
        print(f"Found {len(data['tissue_type'])} tissue type entries")
        
        lengths = {
            'sample_ids': len(data['series_sample_id']),
            'titles': len(data['sample_titles']),
            'braf_status': len(data['braf_status']),
            'tissue_type': len(data['tissue_type'])
        }
        
        if len(set(lengths.values())) > 1:
            raise ValueError(f"Inconsistent data lengths: {lengths}")
    
        # Create DataFrame using series_sample_id
        self.sample_info = pd.DataFrame({
            'sample_id': data['series_sample_id'],
            'title': data['sample_titles'],
            'braf_status': data['braf_status'],
            'tissue_type': data['tissue_type']
        })
        
        # Verify DataFrame creation
        print("\nSample info shape:", self.sample_info.shape)
        print("First few sample IDs:", self.sample_info['sample_id'].head().tolist())
        
        self._create_analysis_groups()
    
    def _create_analysis_groups(self) -> None:
        """Create analysis groups from tissue types"""
        group_mapping = {
            'PTC with a BRAF mutation': 'PTC_BRAF_mut',
            'PTC without a BRAF mutation': 'PTC_BRAF_wt',
            'normal thyroid tissue': 'Normal',
            'FA': 'FA',
            'FVPTC': 'FVPTC'
        }
        self.sample_info['group'] = self.sample_info['tissue_type'].map(group_mapping)
    
    def _read_expression_data(self, sample_dir: str) -> None:
        """Read expression data from individual files"""
        expression_data = {}
        
        print(f"Reading expression data from: {sample_dir}")
        
        for sample_id in self.sample_info['sample_id']:
            file_path = os.path.join(sample_dir, f"{sample_id}-tbl-1.txt")
            print(f"Processing file: {file_path}")
            
            if os.path.exists(file_path):
                try:
                    # Specify dtypes explicitly
                    data = pd.read_csv(file_path, 
                                    sep='\t',
                                    header=None,
                                    names=['probe_id', sample_id],
                                    dtype={'probe_id': str, sample_id: float})
                    
                    expression_data[sample_id] = data.set_index('probe_id')[sample_id]
                    print(f"Successfully read data for {sample_id}")
                except Exception as e:
                    print(f"Error reading {sample_id}: {str(e)}")
            else:
                print(f"File not found: {file_path}")
        
        if not expression_data:
            raise ValueError("No expression data could be loaded")
        
        self.expr_matrix = pd.DataFrame(expression_data)
        print(f"Final expression matrix shape: {self.expr_matrix.shape}")
    
    def _create_anndata(self) -> None: 
        """Create AnnData object"""
        self.adata = sc.AnnData(X=self.expr_matrix.T)
        self.adata.obs = self.sample_info.set_index('sample_id')
    
    def _run_qc(self) -> Dict:
        """Run quality control analysis"""
        qc_dir = os.path.join(self.config.output_dir, 'qc')
        
        # Calculate QC metrics
        qc_stats = pd.DataFrame({
            'n_genes': (self.expr_matrix > 0).sum(),
            'total_counts': self.expr_matrix.sum(),
            'group': self.sample_info['group']
        })
        
        # Save QC results
        qc_stats.to_csv(f'{qc_dir}/qc_stats.csv')
        self._plot_qc(qc_stats, qc_dir)
        
        return {'stats': qc_stats}
    
    def _plot_qc(self, qc_stats: pd.DataFrame, output_dir: str) -> None:
        """Create QC plots"""
        plt.figure(figsize=(12, 6))
        sns.boxplot(data=self.expr_matrix.melt(), x='variable', y='value')
        plt.xticks(rotation=90)
        plt.title('Expression Distribution by Sample')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/expression_boxplot.pdf')
        plt.close()
    
    def _run_pca(self) -> Dict:
        """Run PCA analysis"""
        pca_dir = os.path.join(self.config.output_dir, 'pca')
        
        # Run PCA
        scaled_data = stats.zscore(self.expr_matrix.T)
        pca = PCA()
        pca_result = pca.fit_transform(scaled_data)
        
        # Create and save PCA plot
        self._plot_pca(pca_result, pca.explained_variance_ratio_, pca_dir)
        
        return {'pca_result': pca_result}
    
    def _plot_pca(self, pca_result: np.ndarray, var_ratio: np.ndarray, 
                  output_dir: str) -> None:
        """Create PCA plot"""
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x=pca_result[:, 0], y=pca_result[:, 1],
                       hue=self.sample_info['group'])
        plt.title('PCA of Expression Data')
        plt.xlabel(f'PC1 ({var_ratio[0]:.2%})')
        plt.ylabel(f'PC2 ({var_ratio[1]:.2%})')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/pca_plot.pdf')
        plt.close()
    
    def _run_differential_expression(self) -> Dict:
        """Run differential expression analysis"""
        de_dir = os.path.join(self.config.output_dir, 'differential_expression')
        results = {}
        
        for group1, group2 in self.config.comparisons:
            result = self._analyze_de(group1, group2)
            self._plot_volcano(result, group1, group2, de_dir)
            results[f'{group1}_vs_{group2}'] = result
        
        return results
    
    def _analyze_de(self, group1: str, group2: str) -> pd.DataFrame:
        """Analyze differential expression between two groups"""
        samples1 = self.sample_info[self.sample_info['group'] == group1]['sample_id']
        samples2 = self.sample_info[self.sample_info['group'] == group2]['sample_id']
        
        results = []
        for gene in self.expr_matrix.index:
            expr1 = self.expr_matrix.loc[gene, samples1]
            expr2 = self.expr_matrix.loc[gene, samples2]
            
            t_stat, p_val = stats.ttest_ind(expr1, expr2)
            fc = np.mean(expr1) - np.mean(expr2)
            
            results.append({
                'gene': gene,
                'log2FC': fc,
                'pvalue': p_val,
                'mean_expr1': np.mean(expr1),
                'mean_expr2': np.mean(expr2)
            })
        
        results_df = pd.DataFrame(results)
        
        # Add error handling for multiple testing correction
        try:
            _, pvals_adj, _, _ = multipletests(results_df['pvalue'], 
                                            method='fdr_bh')
            results_df['padj'] = pvals_adj
        except Exception as e:
            print(f"Warning: Error in multiple testing correction: {str(e)}")
            results_df['padj'] = results_df['pvalue']  # Use uncorrected p-values as fallback
        
        return results_df
    
    def _plot_volcano(self, de_results: pd.DataFrame, group1: str, group2: str, 
                     output_dir: str) -> None:
        """Create volcano plot"""
        plt.figure(figsize=(10, 8))
        plt.scatter(de_results['log2FC'], -np.log10(de_results['padj']),
                   alpha=0.5)
        plt.xlabel('Log2 Fold Change')
        plt.ylabel('-log10(Adjusted P-value)')
        plt.title(f'Volcano Plot: {group1} vs {group2}')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/volcano_{group1}_vs_{group2}.pdf')
        plt.close()