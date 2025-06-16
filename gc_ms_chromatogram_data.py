import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import signal
from scipy.interpolate import interp1d
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import warnings
warnings.filterwarnings('ignore')

# Step 1: Data Preparation
def prepare_jatropha_data():
    """
    Prepare the Jatropha curcas GC-MS data for chromatogram generation
    """
    
    # Raw data from your PDF
    data = {
        'RT_min': [6.45, 8.89, 11.62, 13.95, 16.73, 19.21],
        'Compound': ['Saponins', 'Tannins', 'Flavonoids', 'Alkaloids', 'Terpenoids', 'Coumarins'],
        'Base_Peak_mz': [457, 333, 289, 300, 411, 179],
        'Key_Fragments': [
            [439, 411, 343],
            [171, 127, 83],
            [155, 137, 109],
            [282, 226, 138],
            [393, 315, 207],
            [151, 123, 95]
        ],
        'Sample_A_Stems': [7.84, 4.56, 3.92, 2.87, 17.45, 1.92],
        'Sample_B_Seeds': [7.69, 5.47, 9.41, 3.11, 12.98, 5.63]
    }
    
    # Create DataFrame
    df = pd.DataFrame(data)
    
    # Add additional calculated columns
    df['Intensity_A_Normalized'] = df['Sample_A_Stems'] / df['Sample_A_Stems'].max() * 100
    df['Intensity_B_Normalized'] = df['Sample_B_Seeds'] / df['Sample_B_Seeds'].max() * 100
    
    # Create molecular weight categories
    df['MW_Category'] = pd.cut(df['Base_Peak_mz'], 
                               bins=[0, 200, 300, 400, 500], 
                               labels=['Low MW', 'Medium MW', 'High MW', 'Very High MW'])
    
    # Add compound class information
    compound_classes = {
        'Saponins': 'Glycosides',
        'Tannins': 'Polyphenols',
        'Flavonoids': 'Polyphenols',
        'Alkaloids': 'Nitrogen Compounds',
        'Terpenoids': 'Terpenes',
        'Coumarins': 'Phenolic Compounds'
    }
    
    df['Compound_Class'] = df['Compound'].map(compound_classes)
    
    print("Data preparation completed!")
    print("\nDataset Overview:")
    print(df.head())
    print(f"\nDataset shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    
    return df

# Step 2: Generate synthetic chromatogram data points
def generate_chromatogram_data(df, sample_type='A', time_points=1000):
    """
    Generate synthetic chromatogram data points for visualization
    """
    
    # Time range from 0 to 25 minutes
    time_range = np.linspace(0, 25, time_points)
    
    # Initialize intensity array
    intensity = np.zeros(time_points)
    
    # Column selection based on sample type
    if sample_type == 'A':
        abundance_col = 'Sample_A_Stems'
        sample_name = 'Stems (Methanol Extract)'
    else:
        abundance_col = 'Sample_B_Seeds'
        sample_name = 'Seeds (Methanol Extract)'
    
    # Generate peaks for each compound
    for idx, row in df.iterrows():
        rt = row['RT_min']
        abundance = row[abundance_col]
        
        # Create Gaussian peak
        peak_width = 0.3  # Standard deviation for peak width
        peak_intensity = abundance * 1000  # Scale up for visualization
        
        # Generate Gaussian peak
        gaussian_peak = peak_intensity * np.exp(-0.5 * ((time_range - rt) / peak_width) ** 2)
        
        # Add to total intensity
        intensity += gaussian_peak
    
    # Add baseline noise
    noise_level = 50
    noise = np.random.normal(0, noise_level, time_points)
    intensity += noise
    
    # Ensure no negative values
    intensity = np.maximum(intensity, 0)
    
    return time_range, intensity, sample_name

# Step 3: Save prepared data
def save_prepared_data(df, filename='jatropha_gcms_data.csv'):
    """
    Save the prepared data to CSV file
    """
    df.to_csv(filename, index=False)
    print(f"\nData saved to {filename}")
    
    # Also save individual sample data
    sample_a_data = df[['RT_min', 'Compound', 'Base_Peak_mz', 'Sample_A_Stems']].copy()
    sample_a_data.columns = ['RT_min', 'Compound', 'Base_Peak_mz', 'Abundance']
    sample_a_data.to_csv('sample_a_stems.csv', index=False)
    
    sample_b_data = df[['RT_min', 'Compound', 'Base_Peak_mz', 'Sample_B_Seeds']].copy()
    sample_b_data.columns = ['RT_min', 'Compound', 'Base_Peak_mz', 'Abundance']
    sample_b_data.to_csv('sample_b_seeds.csv', index=False)
    
    print("Individual sample data saved:")
    print("- sample_a_stems.csv")
    print("- sample_b_seeds.csv")

# Main execution
if __name__ == "__main__":
    print("=== JATROPHA CURCAS GC-MS DATA PREPARATION ===\n")
    
    # Prepare the data
    jatropha_df = prepare_jatropha_data()
    
    # Generate chromatogram data for both samples
    time_a, intensity_a, name_a = generate_chromatogram_data(jatropha_df, 'A')
    time_b, intensity_b, name_b = generate_chromatogram_data(jatropha_df, 'B')
    
    # Save chromatogram data
    chrom_data_a = pd.DataFrame({'Time_min': time_a, 'Intensity': intensity_a})
    chrom_data_b = pd.DataFrame({'Time_min': time_b, 'Intensity': intensity_b})
    
    chrom_data_a.to_csv('chromatogram_stems.csv', index=False)
    chrom_data_b.to_csv('chromatogram_seeds.csv', index=False)
    
    # Save prepared data
    save_prepared_data(jatropha_df)
    
    print("\n=== DATA PREPARATION COMPLETE ===")
    print("\nFiles created:")
    print("1. jatropha_gcms_data.csv - Main dataset")
    print("2. sample_a_stems.csv - Stems sample data")
    print("3. sample_b_seeds.csv - Seeds sample data")
    print("4. chromatogram_stems.csv - Chromatogram data for stems")
    print("5. chromatogram_seeds.csv - Chromatogram data for seeds")
    
    print("\n=== SUMMARY STATISTICS ===")
    print(f"Total compounds identified: {len(jatropha_df)}")
    print(f"Retention time range: {jatropha_df['RT_min'].min():.2f} - {jatropha_df['RT_min'].max():.2f} minutes")
    print(f"Most abundant compound in stems: {jatropha_df.loc[jatropha_df['Sample_A_Stems'].idxmax(), 'Compound']} ({jatropha_df['Sample_A_Stems'].max():.2f}%)")
    print(f"Most abundant compound in seeds: {jatropha_df.loc[jatropha_df['Sample_B_Seeds'].idxmax(), 'Compound']} ({jatropha_df['Sample_B_Seeds'].max():.2f}%)")
