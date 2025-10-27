from pathlib import Path
from src.config import AnalysisConfig
from src.analysis import ThyroidAnalysis
import sys

# Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54958
#         https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse138198
# Paper https://www.nature.com/articles/s41598-025-90499-x.pdf
#       Validation of DNA methylation and transcriptional characteristics in CCL5 and CXCL8 genes in autoimmune thyroiditis with varying iodine levels
#
# base_path = Path("/Users/anna_gurina/Desktop/KSE_Bioinformatics/BioLab/thyroid_datasets_analysis/methilation_data_analysis/methilation_datasets/GSE54958/GSE54958_RAW")

def main():
    # Define base path
    base_path = Path("/Users/anna_gurina/Desktop/KSE_Bioinformatics/BioLab/datasets/thyroid_datasets_analysis_sandbox/methylation_data_analysis/methylation_datasets/GSE54958/GSE54958_RAW")
    
    # Verify paths exist
    print(f"Checking base path: {base_path}")
    if not base_path.exists():
        print(f"Error: Base path does not exist: {base_path}")
        return
    
    # Setup configuration
    config = AnalysisConfig(
        output_dir=str(base_path / 'results'),
        comparisons=[
            ('PTC_BRAF_mut', 'PTC_BRAF_wt'),
            ('PTC_BRAF_mut', 'Normal'),
            ('FA', 'FVPTC')
        ]
    )
    
    # Initialize analysis
    analysis = ThyroidAnalysis(config)
    
    try:
        # Load data
        print("\nLoading data...")
        analysis.load_data(
            series_matrix_file=str(base_path / "GSE54958_series_matrix.txt"),
            sample_dir=str(base_path / "GSE54958_family.xml")
        )
        
        # Run analysis
        print("\nRunning analysis...")
        results = analysis.run_analysis()
        
        # Print summary
        print("\nAnalysis Complete!")
        print(f"Results saved in: {config.output_dir}")
        
    except FileNotFoundError as e:
        print(f"\nError: File not found: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"\nError: Invalid data: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\nError during analysis: {str(e)}")
        print("Full error:", e.__class__.__name__)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()