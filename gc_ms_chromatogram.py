import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
import seaborn as sns
from scipy import signal
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.signal import find_peaks, savgol_filter
import warnings
import os
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings('ignore')

class JatrophaChromatogramAnalyzer:
    """
    Ultra-Realistic GC-MS Chromatogram Analyzer for Jatropha curcas 
    with Heavy Noise Patterns Matching Real-World Data
    """
    
    def __init__(self):
        """Initialize with hardcoded Jatropha curcas data"""
        
        print("🔍 Initializing with hardcoded Jatropha curcas data...")
        
        # Hardcoded data - no external file needed
        data = {
            'RT_min': [6.23, 8.71, 11.45, 14.02, 16.88, 19.34,],
            'Compound': ['glycosidetriterpenoid', '2,3,4,5-tetrahydroxybenzoic acid decagalloyl glucose ester', '2-(3,4-dihydroxyphenyl)-3,5,7-trihydroxy-4H-chromen-4-one', 'tropane-lactone', '1-methyl-4-(1-methylethenyl)cyclohexene', '7-hydroxy-6-methoxychromen-2-one',],
            'Base_Peak_mz': [455, 331, 287, 298, 409, 177,],
            'Sample_A_Stems': [12.45, 9.87, 15.32, 5.21, 18.90, 4.33,],
            'Sample_B_Seeds': [8.23, 6.54, 10.76, 3.45, 14.67,2.89, ],
            'Compound_Class': ['Glycosides', 'Polyphenols', 'Polyphenols', 'Nitrogen Compounds', 'Terpenes', 'Phenolic Compounds',]
        }
        
        try:
            self.df = pd.DataFrame(data)
            print(f"✅ Data loaded successfully: {len(self.df)} compounds")
            print(f"Columns available: {list(self.df.columns)}")
            print(f"Sample data preview:")
            print(self.df.head(3))
            
        except Exception as e:
            print(f"❌ ERROR creating DataFrame: {e}")
            self.df = None
    
    def create_output_directories(self):
        """Create output directories for PNG and PDF files"""
        directories = {
            'png': './chromatograms/png',
            'pdf': './chromatograms/pdf',
            'data': './chromatograms/data'
        }
        
        for dir_type, path in directories.items():
            if not os.path.exists(path):
                os.makedirs(path)
                print(f"📁 Created directory: {path}")
        
        return directories
    
    def save_figure(self, fig, filename_base, save_png=True, save_pdf=True, directories=None):
        """Save figure as PNG and/or PDF with high quality"""
        if directories is None:
            directories = self.create_output_directories()
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        saved_files = []
        
        try:
            # Save as PNG
            if save_png:
                png_filename = f"{filename_base}_{timestamp}.png"
                png_path = os.path.join(directories['png'], png_filename)
                fig.savefig(png_path, dpi=300, bbox_inches='tight', 
                           facecolor='white', edgecolor='none', format='png')
                saved_files.append(png_path)
                print(f"💾 PNG saved: {png_path}")
            
            # Save as PDF
            if save_pdf:
                pdf_filename = f"{filename_base}_{timestamp}.pdf"
                pdf_path = os.path.join(directories['pdf'], pdf_filename)
                fig.savefig(pdf_path, bbox_inches='tight', 
                           facecolor='white', edgecolor='none', format='pdf')
                saved_files.append(pdf_path)
                print(f"📄 PDF saved: {pdf_path}")
        
        except Exception as e:
            print(f"❌ Failed to save figure: {e}")
        
        return saved_files
    
    def add_heavy_realistic_noise(self, intensity, time_range, noise_level=0.8):
        """Add HEAVY realistic noise patterns matching real GC-MS data"""
        
        # Start with the original intensity
        noisy_intensity = intensity.copy()
        
        # 1. MASSIVE baseline fluctuations (like real chromatograms)
        baseline_freq1 = 0.05  # Very slow drift
        baseline_freq2 = 0.15  # Medium frequency variations
        baseline_freq3 = 0.3   # Faster variations
        
        baseline_noise = (
            800 * np.sin(2 * np.pi * baseline_freq1 * time_range) + 
            400 * np.sin(2 * np.pi * baseline_freq2 * time_range + 1.5) +
            200 * np.sin(2 * np.pi * baseline_freq3 * time_range + 2.8) +
            300 * np.sin(2 * np.pi * baseline_freq1 * 1.7 * time_range + 0.8) +
            150 * np.cos(2 * np.pi * baseline_freq2 * 2.3 * time_range)
        )
        
        # 2. TONS of random small peaks (chemical noise)
        chemical_noise = np.zeros_like(intensity)
        n_random_peaks = np.random.randint(80, 150)  # MANY more noise peaks
        
        for _ in range(n_random_peaks):
            peak_pos = np.random.uniform(0, 30)
            peak_idx = np.argmin(np.abs(time_range - peak_pos))
            peak_width = np.random.uniform(0.02, 0.15)  # Various widths
            peak_height = np.random.uniform(20, 800)     # Various heights
            
            # Add random peaks with realistic shapes
            width_points = int(peak_width * len(time_range) / 30)
            width_points = max(1, width_points)
            
            for i in range(max(0, peak_idx - width_points*2), min(len(time_range), peak_idx + width_points*3)):
                distance = abs(i - peak_idx)
                if distance <= width_points*2:
                    # Mix of gaussian and exponential decay
                    if i <= peak_idx:
                        gaussian_factor = np.exp(-(distance**2) / (2 * (width_points/2)**2))
                    else:
                        gaussian_factor = np.exp(-distance / width_points)
                    
                    chemical_noise[i] += peak_height * gaussian_factor
        
        # 3. High-frequency electronic noise (much more aggressive)
        electronic_noise = np.random.normal(0, 80, len(intensity))  # Higher amplitude
        
        # Add spiky electronic interference
        spike_positions = np.random.choice(len(intensity), size=int(len(intensity) * 0.1), replace=False)
        for pos in spike_positions:
            electronic_noise[pos] += np.random.normal(0, 200)
        
        # 4. Column bleed and temperature effects (more pronounced)
        column_bleed = np.linspace(0, 500, len(intensity)) * (time_range / 30)**1.5
        
        # Temperature fluctuations causing baseline shifts
        temp_fluctuations = (
            100 * np.sin(2 * np.pi * 0.02 * time_range) +
            80 * np.sin(2 * np.pi * 0.08 * time_range + 1.2) +
            60 * np.sin(2 * np.pi * 0.12 * time_range + 2.1)
        )
        
        # 5. Pump pulsations and flow irregularities
        pump_noise = (
            150 * np.sin(2 * np.pi * 0.7 * time_range) + 
            100 * np.sin(2 * np.pi * 1.2 * time_range + 0.5) +
            80 * np.sin(2 * np.pi * 2.1 * time_range + 1.8) +
            50 * np.sin(2 * np.pi * 3.4 * time_range + 2.3)
        )
        
        # 6. Ion source contamination (random bursts)
        contamination_noise = np.zeros_like(intensity)
        n_contamination_events = np.random.randint(15, 30)
        
        for _ in range(n_contamination_events):
            event_start = np.random.uniform(0, 28)
            event_duration = np.random.uniform(0.1, 1.5)
            event_intensity = np.random.uniform(100, 1000)
            
            start_idx = np.argmin(np.abs(time_range - event_start))
            end_idx = np.argmin(np.abs(time_range - (event_start + event_duration)))
            
            for i in range(start_idx, min(end_idx + 1, len(intensity))):
                decay_factor = np.exp(-(i - start_idx) / (end_idx - start_idx + 1))
                contamination_noise[i] += event_intensity * decay_factor * np.random.uniform(0.3, 1.0)
        
        # 7. Detector saturation and recovery artifacts
        saturation_noise = np.zeros_like(intensity)
        high_peaks = intensity > np.percentile(intensity, 85)
        
        for i in range(len(intensity)):
            if high_peaks[i]:
                # Add tailing and fronting artifacts
                for j in range(max(0, i-10), min(len(intensity), i+20)):
                    distance = abs(j - i)
                    if distance > 0:
                        artifact_intensity = np.random.uniform(50, 300) * np.exp(-distance/5)
                        saturation_noise[j] += artifact_intensity
        
        # 8. Mass spectrometer vacuum fluctuations
        vacuum_noise = 50 * np.random.beta(2, 5, len(intensity)) * np.sin(2 * np.pi * 0.3 * time_range)
        
        # 9. Solvent impurities and background ions
        solvent_background = np.zeros_like(intensity)
        # Common solvent peaks at random positions
        solvent_peaks = [2.1, 3.8, 7.2, 12.5, 18.9, 24.3]
        for peak_rt in solvent_peaks:
            if 0 <= peak_rt <= 30:
                peak_idx = np.argmin(np.abs(time_range - peak_rt))
                peak_height = np.random.uniform(200, 800)
                width = np.random.randint(5, 15)
                
                for i in range(max(0, peak_idx - width), min(len(intensity), peak_idx + width)):
                    distance = abs(i - peak_idx)
                    if distance <= width:
                        solvent_background[i] += peak_height * np.exp(-distance**2 / (2 * (width/3)**2))
        
        # 10. Random walk baseline (very realistic)
        random_walk = np.zeros_like(intensity)
        random_walk[0] = np.random.normal(0, 100)
        for i in range(1, len(intensity)):
            random_walk[i] = random_walk[i-1] + np.random.normal(0, 20)
        
        # Combine ALL noise sources with heavy weighting
        total_noise = (
            baseline_noise * 1.2 +
            chemical_noise * 1.0 +
            electronic_noise * 0.8 +
            column_bleed * 0.6 +
            temp_fluctuations * 0.9 +
            pump_noise * 0.7 +
            contamination_noise * 0.8 +
            saturation_noise * 0.5 +
            vacuum_noise * 1.1 +
            solvent_background * 1.0 +
            random_walk * 0.4
        ) * noise_level
        
        # Apply heavy noise
        noisy_intensity = intensity + total_noise
        
        # Ensure no negative values but allow very low baseline
        noisy_intensity = np.maximum(noisy_intensity, np.min(intensity) * 0.1)
        
        # Add final random spikes throughout
        n_spikes = np.random.randint(20, 50)
        spike_positions = np.random.choice(len(noisy_intensity), size=n_spikes, replace=False)
        for pos in spike_positions:
            spike_height = np.random.uniform(100, 1500)
            noisy_intensity[pos] += spike_height
            
            # Add spike tailing
            for j in range(pos + 1, min(pos + 8, len(noisy_intensity))):
                decay = np.exp(-(j - pos) / 3)
                noisy_intensity[j] += spike_height * decay * 0.3
        
        return noisy_intensity
    
    def generate_realistic_chromatogram(self, sample_type='A', time_points=6000, noise_level=0.8):
        """Generate ultra-realistic chromatogram with heavy noise matching real data"""
        
        if self.df is None:
            print("❌ Cannot generate chromatogram: No data loaded!")
            return None, None, None, None
        
        print(f"🔄 Generating ULTRA-realistic chromatogram for sample type: {sample_type}")
        
        time_range = np.linspace(0, 30, time_points)
        intensity = np.zeros(time_points)
        
        # Sample selection
        if sample_type == 'A':
            abundance_col = 'Sample_A_Stems'
            sample_name = 'Jatropha curcas - Stems (Methanol Extract)'
        else:
            abundance_col = 'Sample_B_Seeds'
            sample_name = 'Jatropha curcas - Seeds (Methanol Extract)'
        
        # Check if required columns exist
        required_cols = ['RT_min', 'Compound', 'Base_Peak_mz', abundance_col]
        missing_cols = [col for col in required_cols if col not in self.df.columns]
        if missing_cols:
            print(f"❌ Missing required columns: {missing_cols}")
            return None, None, None, None
        
        peak_info = []
        
        # Start with noisy baseline (much more variable)
        baseline_level = 5000 + np.random.normal(0, 500, time_points)
        baseline_level += 100 * np.sin(2 * np.pi * 0.1 * time_range) 
        baseline_level += np.linspace(0, 1000, time_points)  # Gradual increase
        intensity = baseline_level.copy()
        
        # Add main compound peaks
        for idx, row in self.df.iterrows():
            rt = row['RT_min']
            abundance = row[abundance_col]
            compound = row['Compound']
            base_peak = row['Base_Peak_mz']
            
            # More realistic peak parameters with more variation
            peak_width = np.random.uniform(0.05, 0.25)  # Wide range of widths
            peak_height = abundance * 50000 + np.random.normal(0, 15000)  # More variation
            peak_height = max(peak_height, 10000)  # Ensure minimum height
            
            # Create realistic peak shapes with heavy tailing
            rt_idx = np.argmin(np.abs(time_range - rt))
            width_points = int(peak_width * time_points / 30)
            
            # Asymmetric peak with heavy tailing (like real data)
            for i in range(max(0, rt_idx - width_points * 3), min(len(time_range), rt_idx + width_points * 8)):
                distance = (i - rt_idx) / width_points
                
                if distance <= 0:  # Leading edge - steeper
                    peak_val = peak_height * np.exp(-0.3 * distance**2)
                else:  # Trailing edge - heavy exponential decay
                    peak_val = peak_height * np.exp(-distance * 0.8)  # Heavy tailing
                
                intensity[i] = max(intensity[i], peak_val + baseline_level[i])
            
            # Add peak shoulders and satellite peaks (very common in real data)
            if abundance > 3:
                # Multiple shoulders
                n_shoulders = np.random.randint(1, 4)
                for s in range(n_shoulders):
                    shoulder_rt = rt + np.random.uniform(-0.4, 0.6)
                    shoulder_idx = np.argmin(np.abs(time_range - shoulder_rt))
                    shoulder_height = peak_height * np.random.uniform(0.05, 0.4)
                    
                    for i in range(max(0, shoulder_idx - 5), min(len(time_range), shoulder_idx + 10)):
                        distance = abs(i - shoulder_idx) / 5
                        if distance <= 2:
                            shoulder_val = shoulder_height * np.exp(-distance**2 / 2)
                            intensity[i] = max(intensity[i], shoulder_val + baseline_level[i])
            
            peak_info.append({
                'RT': rt,
                'Compound': compound,
                'Abundance': abundance,
                'Base_Peak': base_peak,
                'Height': peak_height
            })
        
        # Add MASSIVE noise patterns to match real data
        intensity = self.add_heavy_realistic_noise(intensity, time_range, noise_level)
        
        # Apply minimal smoothing (real data is quite noisy)
        if time_points > 50:
            intensity = savgol_filter(intensity, window_length=min(9, time_points//50), polyorder=1)
        
        return time_range, intensity, sample_name, peak_info
    
    def create_ultra_realistic_chromatogram(self, sample_type='A', figsize=(16, 10), noise_level=0.8, 
                                          save_png=True, save_pdf=True, directories=None):
        """Create ultra-realistic chromatogram matching real-world noisy data"""
        
        if self.df is None:
            print("❌ Cannot create chromatogram: No data loaded!")
            return None, []
        
        print(f"🎨 Creating ULTRA-realistic chromatogram for sample {sample_type}...")
        
        # Generate heavily noisy data
        time_data, intensity_data, sample_name, peaks_data = self.generate_realistic_chromatogram(
            sample_type, time_points=6000, noise_level=noise_level
        )
        
        if any(x is None for x in [time_data, intensity_data, sample_name, peaks_data]):
            print("❌ Failed to generate chromatogram data!")
            return None, []
        
        # Create figure with realistic styling
        fig, ax = plt.subplots(1, 1, figsize=figsize, facecolor='white')
        fig.patch.set_facecolor('white')
        
        # Plot with thin line to show all the noise detail
        ax.plot(time_data, intensity_data, 'k-', linewidth=0.4, alpha=0.9)
        ax.set_facecolor('white')
        
        # Add peak labels for significant peaks only (like real chromatograms)
        labeled_peaks = 0
        for peak in peaks_data:
            if peak['Abundance'] > 4 and labeled_peaks < 8:  # Limit labels to avoid clutter
                peak_idx = np.argmin(np.abs(time_data - peak['RT']))
                
                # Find actual local maximum near the expected retention time
                search_range = int(0.5 * len(time_data) / 30)  # 0.5 min search window
                start_idx = max(0, peak_idx - search_range)
                end_idx = min(len(intensity_data), peak_idx + search_range)
                local_max_idx = start_idx + np.argmax(intensity_data[start_idx:end_idx])
                actual_peak_height = intensity_data[local_max_idx]
                actual_rt = time_data[local_max_idx]
                
                # Add compound label
                ax.text(actual_rt, actual_peak_height + max(intensity_data) * 0.05, 
                       f"{peak['Compound'][:25]}...\n{actual_rt:.2f} min", 
                       ha='center', va='bottom', fontsize=7, fontweight='normal',
                       rotation=90, alpha=0.7)
                labeled_peaks += 1
        
        # Realistic formatting
        sample_name_short = "Stems" if sample_type == 'A' else "Seeds"
        ax.set_title(f'udsm lab 02025-0025: Phytochemicals Screening', 
                    fontsize=12, fontweight='bold', pad=20)
        ax.set_ylabel('Abundance', fontsize=11)
        ax.set_xlabel('Retention Time (min)', fontsize=11)
        ax.set_xlim(0, 30)
        ax.set_ylim(0, max(intensity_data) * 1.2)
        
        # Realistic axis styling
        for spine in ax.spines.values():
            spine.set_color('black')
            spine.set_linewidth(1.0)
        
        ax.grid(True, alpha=0.3, linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=10)
        
        plt.tight_layout()
        
        # Save figure
        saved_files = []
        if save_png or save_pdf:
            sample_name_file = "stems" if sample_type == 'A' else "seeds"
            filename_base = f"jatropha_{sample_name_file}_chromatogram_noise{noise_level}"
            saved_files = self.save_figure(fig, filename_base, save_png, save_pdf, directories)
        
        return fig, saved_files
    
    def create_dual_realistic_chromatogram(self, figsize=(16, 12), noise_level=0.8, 
                                         save_png=True, save_pdf=True, directories=None):
        """Create dual ultra-realistic chromatograms (both samples)"""
        
        if self.df is None:
            print("❌ Cannot create chromatogram: No data loaded!")
            return None, []
        
        print("🎨 Creating dual ULTRA-realistic chromatograms...")
        
        # Generate data for both samples with heavy noise
        time_a, int_a, name_a, peaks_a = self.generate_realistic_chromatogram('A', time_points=6000, noise_level=noise_level)
        time_b, int_b, name_b, peaks_b = self.generate_realistic_chromatogram('B', time_points=6000, noise_level=noise_level)
        
        if any(x is None for x in [time_a, int_a, name_a, peaks_a, time_b, int_b, name_b, peaks_b]):
            print("❌ Failed to generate chromatogram data!")
            return None, []
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, facecolor='white')
        fig.patch.set_facecolor('white')
        
        # Plot Sample A (Stems) - Ultra realistic
        ax1.plot(time_a, int_a, 'k-', linewidth=0.4, alpha=0.9)
        ax1.set_facecolor('white')
        
        # Add labels for Sample A
        labeled_a = 0
        for peak in peaks_a:
            if peak['Abundance'] > 4 and labeled_a < 6:
                peak_idx = np.argmin(np.abs(time_a - peak['RT']))
                search_range = int(0.5 * len(time_a) / 30)
                start_idx = max(0, peak_idx - search_range)
                end_idx = min(len(int_a), peak_idx + search_range)
                local_max_idx = start_idx + np.argmax(int_a[start_idx:end_idx])
                actual_peak_height = int_a[local_max_idx]
                actual_rt = time_a[local_max_idx]
                
                ax1.text(actual_rt, actual_peak_height + max(int_a) * 0.08, 
                        f"{peak['Compound'][:20]}...\n{actual_rt:.2f} min", 
                        ha='center', va='bottom', fontsize=7,
                        rotation=90, alpha=0.7)
                labeled_a += 1
        
        ax1.set_title(f'GC-MS Chromatogram: Jatropha curcas Stems Extract (Noise: {noise_level})', 
                     fontsize=12, fontweight='bold', pad=15)
        ax1.set_ylabel('Abundance', fontsize=11)
        ax1.set_xlim(0, 30)
        ax1.set_ylim(0, max(int_a) * 1.3)
        ax1.grid(True, alpha=0.3, linewidth=0.5)
        
        # Plot Sample B (Seeds) - Ultra realistic
        ax2.plot(time_b, int_b, 'k-', linewidth=0.4, alpha=0.9)
        ax2.set_facecolor('white')
        
        # Add labels for Sample B
        labeled_b = 0
        for peak in peaks_b:
            if peak['Abundance'] > 4 and labeled_b < 6:
                peak_idx = np.argmin(np.abs(time_b - peak['RT']))
                search_range = int(0.5 * len(time_b) / 30)
                start_idx = max(0, peak_idx - search_range)
                end_idx = min(len(int_b), peak_idx + search_range)
                local_max_idx = start_idx + np.argmax(int_b[start_idx:end_idx])
                actual_peak_height = int_b[local_max_idx]
                actual_rt = time_b[local_max_idx]
                
                ax2.text(actual_rt, actual_peak_height + max(int_b) * 0.08, 
                        f"{peak['Compound'][:20]}...\n{actual_rt:.2f} min", 
                        ha='center', va='bottom', fontsize=7,
                        rotation=90, alpha=0.7)
                labeled_b += 1
        
        ax2.set_title(f'GC-MS Chromatogram: Jatropha curcas Seeds Extract (Noise: {noise_level})', 
                     fontsize=12, fontweight='bold', pad=15)
        ax2.set_ylabel('Abundance', fontsize=11)
        ax2.set_xlabel('Retention Time (min)', fontsize=11)
        ax2.set_xlim(0, 30)
        ax2.set_ylim(0, max(int_b) * 1.3)
        ax2.grid(True, alpha=0.3, linewidth=0.5)
        
        # Styling for both axes
        for ax in [ax1, ax2]:
            for spine in ax.spines.values():
                spine.set_color('black')
                spine.set_linewidth(1.0)
            ax.tick_params(axis='both', which='major', labelsize=10)
        
        plt.tight_layout()
        
        # Save figure
        saved_files = []
        if save_png or save_pdf:
            filename_base = f"jatropha_dual_chromatogram_noise{noise_level}"
            saved_files = self.save_figure(fig, filename_base, save_png, save_pdf, directories)
        
        return fig, saved_files
    
    def create_comprehensive_pdf_report(self, noise_level=0.8, directories=None):
        """Create a comprehensive PDF report with all chromatograms and analysis"""
        
        if self.df is None:
            print("❌ Cannot create PDF report: No data loaded!")
            return None
        
        if directories is None:
            directories = self.create_output_directories()
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        pdf_filename = f"jatropha_comprehensive_report_noise{noise_level}_{timestamp}.pdf"
        pdf_path = os.path.join(directories['pdf'], pdf_filename)
        
        print(f"📖 Creating comprehensive PDF report...")
        
        with PdfPages(pdf_path) as pdf:
            # Page 1: Title and Analysis Report
            fig1 = plt.figure(figsize=(8.5, 11))
            fig1.text(0.5, 0.95, 'ULTRA-REALISTIC GC-MS ANALYSIS REPORT', 
                     ha='center', va='top', fontsize=16, fontweight='bold')
            fig1.text(0.5, 0.90, 'Jatropha curcas Methanol Extracts', 
                     ha='center', va='top', fontsize=14)
            fig1.text(0.5, 0.87, f'With Heavy Noise Modeling (Level: {noise_level})', 
                     ha='center', va='top', fontsize=12)
            fig1.text(0.5, 0.83, f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', 
                     ha='center', va='top', fontsize=10)
            
            # Add analysis information
            y_pos = 0.75
            report_text = f"""
INSTRUMENT PARAMETERS:
• GC-MS (HP-5MS column, EI mode)
• Plant Species: Jatropha curcas
• Extraction Method: Methanol extraction
• Samples Analyzed: Stems and Seeds
• Noise Modeling: Ultra-heavy realistic patterns

COMPOUND IDENTIFICATION SUMMARY:
• Total compounds identified: {len(self.df)}
• Retention time range: {self.df['RT_min'].min():.2f} - {self.df['RT_min'].max():.2f} minutes
• Mass range: {self.df['Base_Peak_mz'].min()} - {self.df['Base_Peak_mz'].max()} m/z

SAMPLE A (STEMS) - TOP COMPOUNDS:
"""
            stems_top = self.df.nlargest(3, 'Sample_A_Stems')
            for idx, row in stems_top.iterrows():
                report_text += f"• {row['Compound'][:50]}: {row['Sample_A_Stems']:.2f}%\n  (RT: {row['RT_min']:.2f} min, m/z: {row['Base_Peak_mz']})\n"
            
            report_text += "\nSAMPLE B (SEEDS) - TOP COMPOUNDS:\n"
            seeds_top = self.df.nlargest(3, 'Sample_B_Seeds')
            for idx, row in seeds_top.iterrows():
                report_text += f"• {row['Compound'][:50]}: {row['Sample_B_Seeds']:.2f}%\n  (RT: {row['RT_min']:.2f} min, m/z: {row['Base_Peak_mz']})\n"
            
            report_text += """
REALISTIC NOISE CHARACTERISTICS:
• Heavy baseline fluctuations: ±800 units
• Electronic noise: ±80 units with spikes
• Chemical background: 80-150 random peaks
• Ion source contamination: 15-30 events
• Detector artifacts and saturation effects
• Pump pulsations and flow irregularities
• Solvent impurity peaks
• Random walk baseline drift
"""
            
            fig1.text(0.1, y_pos, report_text, ha='left', va='top', fontsize=9, 
                     fontfamily='monospace', wrap=True)
            pdf.savefig(fig1, bbox_inches='tight')
            plt.close(fig1)
            
            # Page 2: Stems Chromatogram
            fig2, saved_files_a = self.create_ultra_realistic_chromatogram('A', noise_level=noise_level, 
                                                                          save_png=False, save_pdf=False)
            if fig2:
                pdf.savefig(fig2, bbox_inches='tight')
                plt.close(fig2)
            
            # Page 3: Seeds Chromatogram
            fig3, saved_files_b = self.create_ultra_realistic_chromatogram('B', noise_level=noise_level, 
                                                                          save_png=False, save_pdf=False)
            if fig3:
                pdf.savefig(fig3, bbox_inches='tight')
                plt.close(fig3)
            
            # Page 4: Dual Chromatogram
            fig4, saved_files_dual = self.create_dual_realistic_chromatogram(noise_level=noise_level, 
                                                                            save_png=False, save_pdf=False)
            if fig4:
                pdf.savefig(fig4, bbox_inches='tight')
                plt.close(fig4)
            
            # Page 5: Data Table
            fig5 = plt.figure(figsize=(11, 8.5))
            ax = fig5.add_subplot(111)
            ax.axis('tight')
            ax.axis('off')
            
            # Create table data
            table_data = self.df[['Compound', 'RT_min', 'Sample_A_Stems', 'Sample_B_Seeds', 
                                'Base_Peak_mz', 'Compound_Class']].copy()
            table_data['Difference'] = table_data['Sample_B_Seeds'] - table_data['Sample_A_Stems']
            table_data = table_data.round(2)
            
            # Truncate long compound names for table
            table_data['Compound'] = table_data['Compound'].apply(lambda x: x[:30] + '...' if len(x) > 30 else x)
            
            table = ax.table(cellText=table_data.values, colLabels=table_data.columns,
                           cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
            table.auto_set_font_size(False)
            table.set_fontsize(8)
            table.scale(1, 2)
            
            plt.title('Compound Analysis Data Table', fontsize=14, fontweight='bold', pad=20)
            pdf.savefig(fig5, bbox_inches='tight')
            plt.close(fig5)
        
        print(f"📄 Comprehensive PDF report saved: {pdf_path}")
        return pdf_path
    
    def generate_analysis_report(self):
        """Generate comprehensive analysis report"""
        
        if self.df is None:
            print("❌ Cannot generate analysis report: No data loaded!")
            return
        
        print("📊 Generating analysis report...")
        
        print("=" * 80)
        print("         ULTRA-REALISTIC GC-MS ANALYSIS REPORT")
        print("                Jatropha curcas Methanol Extracts")
        print("              (With Heavy Noise Modeling)")
        print("=" * 80)
        
        print(f"\nINSTRUMENT: GC-MS (HP-5MS column, EI mode)")
        print(f"PLANT SPECIES: Jatropha curcas")
        print(f"EXTRACTION METHOD: Methanol extraction")
        print(f"SAMPLES ANALYZED: Stems and Seeds")
        print(f"NOISE MODELING: Ultra-heavy realistic patterns")
        
        print(f"\nCOMPOUND IDENTIFICATION SUMMARY:")
        print(f"Total compounds identified: {len(self.df)}")
        print(f"Retention time range: {self.df['RT_min'].min():.2f} - {self.df['RT_min'].max():.2f} minutes")
        print(f"Mass range: {self.df['Base_Peak_mz'].min()} - {self.df['Base_Peak_mz'].max()} m/z")
        
        print(f"\nSAMPLE A (STEMS) - TOP COMPOUNDS:")
        stems_top = self.df.nlargest(3, 'Sample_A_Stems')
        for idx, row in stems_top.iterrows():
            print(f"  • {row['Compound']}: {row['Sample_A_Stems']:.2f}% (RT: {row['RT_min']:.2f} min, m/z: {row['Base_Peak_mz']})")
        
        print(f"\nSAMPLE B (SEEDS) - TOP COMPOUNDS:")
        seeds_top = self.df.nlargest(3, 'Sample_B_Seeds')
        for idx, row in seeds_top.iterrows():
            print(f"  • {row['Compound']}: {row['Sample_B_Seeds']:.2f}% (RT: {row['RT_min']:.2f} min, m/z: {row['Base_Peak_mz']})")
        
        print(f"\nREALISTIC NOISE CHARACTERISTICS:")
        print(f"  • Heavy baseline fluctuations: ±800 units")
        print(f"  • Electronic noise: ±80 units with spikes")
        print(f"  • Chemical background: 80-150 random peaks")
        print(f"  • Ion source contamination: 15-30 events")
        print(f"  • Detector artifacts and saturation effects")
        print(f"  • Pump pulsations and flow irregularities")
        print(f"  • Solvent impurity peaks")
        print(f"  • Random walk baseline drift")
        
        print("=" * 80)

def main():
    """Main execution function with ultra-realistic noise and file saving"""
    
    print("=" * 60)
    print("    ULTRA-REALISTIC JATROPHA CURCAS GC-MS GENERATOR")
    print("    (Matching Real-World Noisy Data)")
    print("    WITH PNG & PDF EXPORT FUNCTIONALITY")
    print("=" * 60)
    
    # Initialize analyzer
    analyzer = JatrophaChromatogramAnalyzer()
    
    if analyzer.df is None:
        print("\n❌ CRITICAL ERROR: Failed to initialize data!")
        return
    
    # Create output directories
    directories = analyzer.create_output_directories()
    
    # Set heavy noise level for ultra-realistic results
    print("\n🎛️ Ultra-Realistic Noise Levels:")
    print("  1. Heavy noise (0.8) - Matches real-world data")
    print("  2. Extreme noise (1.2) - Very challenging conditions")
    print("  3. Custom level (0.5-2.0)")
    
    try:
        choice = input("Select noise level (1-3) or press Enter for heavy (0.8): ").strip()
        
        if choice == '1' or choice == '':
            noise_level = 0.8
            print("✅ Using heavy realistic noise (0.8)")
        elif choice == '2':
            noise_level = 1.2
            print("✅ Using extreme noise (1.2)")
        elif choice == '3':
            custom = input("Enter custom noise level (0.5-2.0): ").strip()
            try:
                noise_level = float(custom)
                if not (0.5 <= noise_level <= 2.0):
                    noise_level = 0.8
                    print("⚠️ Out of range, using default 0.8")
                else:
                    print(f"✅ Using custom noise level: {noise_level}")
            except:
                noise_level = 0.8
                print("⚠️ Invalid input, using default 0.8")
        else:
            noise_level = 0.8
            print("⚠️ Invalid choice, using default 0.8")
    
    except KeyboardInterrupt:
        print("\n\n⏹️ Process interrupted by user")
        return
    except Exception as e:
        print(f"\n⚠️ Input error: {e}")
        noise_level = 0.8
        print("Using default noise level: 0.8")
    
    print(f"\n🎯 Selected noise level: {noise_level}")
    
    # File export options
    print("\n💾 File Export Options:")
    print("  1. Save PNG + PDF + Show plots")
    print("  2. Save PNG only + Show plots")
    print("  3. Save PDF only + Show plots")
    print("  4. Only show plots (no saving)")
    print("  5. Save comprehensive PDF report")
    
    try:
        export_choice = input("Select export option (1-5) or press Enter for PNG+PDF+Show (1): ").strip()
        
        if export_choice == '2':
            save_png, save_pdf, show_plots = True, False, True
        elif export_choice == '3':
            save_png, save_pdf, show_plots = False, True, True
        elif export_choice == '4':
            save_png, save_pdf, show_plots = False, False, True
        elif export_choice == '5':
            # Create comprehensive PDF report
            pdf_path = analyzer.create_comprehensive_pdf_report(noise_level, directories)
            if pdf_path:
                print(f"✅ Comprehensive PDF report created!")
            return
        else:  # Default or choice 1
            save_png, save_pdf, show_plots = True, True, True
        
        print(f"✅ Export settings: PNG={'ON' if save_png else 'OFF'}, PDF={'ON' if save_pdf else 'OFF'}, Display={'ON' if show_plots else 'OFF'}")
    
    except Exception as e:
        print(f"⚠️ Using defaults: PNG ON, PDF ON, Display ON")
        save_png, save_pdf, show_plots = True, True, True
    
    # Generate analysis report first
    print("\n" + "="*50)
    print("📋 GENERATING ANALYSIS REPORT")
    print("="*50)
    analyzer.generate_analysis_report()
    
    # Visualization options
    print("\n🎨 Chromatogram Generation Options:")
    print("  1. Single chromatogram - Stems (Sample A)")
    print("  2. Single chromatogram - Seeds (Sample B)")
    print("  3. Dual chromatogram - Both samples")
    print("  4. Generate all chromatograms")
    print("  5. Create comprehensive PDF report")
    
    all_saved_files = []
    
    try:
        viz_choice = input("\nSelect visualization (1-5) or press Enter for all (4): ").strip()
        
        if viz_choice == '1':
            print("\n🔬 Generating Stems chromatogram...")
            fig, saved_files = analyzer.create_ultra_realistic_chromatogram('A', noise_level=noise_level, 
                                                                           save_png=save_png, save_pdf=save_pdf, directories=directories)
            if fig and show_plots:
                plt.show()
            if fig:
                all_saved_files.extend(saved_files)
                print("✅ Stems chromatogram completed!")
        
        elif viz_choice == '2':
            print("\n🔬 Generating Seeds chromatogram...")
            fig, saved_files = analyzer.create_ultra_realistic_chromatogram('B', noise_level=noise_level, 
                                                                           save_png=save_png, save_pdf=save_pdf, directories=directories)
            if fig and show_plots:
                plt.show()
            if fig:
                all_saved_files.extend(saved_files)
                print("✅ Seeds chromatogram completed!")
        
        elif viz_choice == '3':
            print("\n🔬 Generating dual chromatogram...")
            fig, saved_files = analyzer.create_dual_realistic_chromatogram(noise_level=noise_level, 
                                                                          save_png=save_png, save_pdf=save_pdf, directories=directories)
            if fig and show_plots:
                plt.show()
            if fig:
                all_saved_files.extend(saved_files)
                print("✅ Dual chromatogram completed!")
        
        elif viz_choice == '5':
            print("\n📖 Creating comprehensive PDF report...")
            pdf_path = analyzer.create_comprehensive_pdf_report(noise_level, directories)
            if pdf_path:
                all_saved_files.append(pdf_path)
                print("✅ Comprehensive PDF report completed!")
        
        else:  # Choice 4 or default
            print("\n🔬 Generating ALL chromatograms...")
            
            # Generate Stems chromatogram
            print("  📊 Creating Stems chromatogram...")
            fig_a, saved_files_a = analyzer.create_ultra_realistic_chromatogram('A', noise_level=noise_level, 
                                                                               save_png=save_png, save_pdf=save_pdf, directories=directories)
            if fig_a and show_plots:
                plt.show()
            if fig_a:
                all_saved_files.extend(saved_files_a)
                print("  ✅ Stems chromatogram completed!")
            
            # Generate Seeds chromatogram
            print("  📊 Creating Seeds chromatogram...")
            fig_b, saved_files_b = analyzer.create_ultra_realistic_chromatogram('B', noise_level=noise_level, 
                                                                               save_png=save_png, save_pdf=save_pdf, directories=directories)
            if fig_b and show_plots:
                plt.show()
            if fig_b:
                all_saved_files.extend(saved_files_b)
                print("  ✅ Seeds chromatogram completed!")
            
            # Generate dual chromatogram
            print("  📊 Creating dual chromatogram...")
            fig_dual, saved_files_dual = analyzer.create_dual_realistic_chromatogram(noise_level=noise_level, 
                                                                                    save_png=save_png, save_pdf=save_pdf, directories=directories)
            if fig_dual and show_plots:
                plt.show()
            if fig_dual:
                all_saved_files.extend(saved_files_dual)
                print("  ✅ Dual chromatogram completed!")
            
            print("\n🎉 All chromatograms generated successfully!")
    
    except KeyboardInterrupt:
        print("\n\n⏹️ Visualization interrupted by user")
        return
    except Exception as e:
        print(f"\n❌ Error during visualization: {e}")
        return
    
    # Additional analysis options
    print("\n" + "="*50)
    print("📈 ADDITIONAL ANALYSIS OPTIONS")
    print("="*50)
    print("  1. Show compound comparison table")
    print("  2. Generate statistical summary")
    print("  3. Export data to CSV")
    print("  4. Create comprehensive PDF report")
    print("  5. Skip additional analysis")
    
    try:
        extra_choice = input("Select option (1-5) or press Enter to skip (5): ").strip()
        
        if extra_choice == '1':
            print("\n📊 COMPOUND COMPARISON TABLE:")
            print("-" * 80)
            comparison_df = analyzer.df[['Compound', 'RT_min', 'Sample_A_Stems', 'Sample_B_Seeds', 'Compound_Class']].copy()
            comparison_df['Difference'] = comparison_df['Sample_B_Seeds'] - comparison_df['Sample_A_Stems']
            comparison_df = comparison_df.round(2)
            print(comparison_df.to_string(index=False))
        
        elif extra_choice == '2':
            print("\n📈 STATISTICAL SUMMARY:")
            print("-" * 50)
            print(f"Stems (Sample A) Statistics:")
            print(f"  Mean abundance: {analyzer.df['Sample_A_Stems'].mean():.2f}%")
            print(f"  Std deviation: {analyzer.df['Sample_A_Stems'].std():.2f}%")
            print(f"  Max compound: {analyzer.df.loc[analyzer.df['Sample_A_Stems'].idxmax(), 'Compound']}")
            print(f"  Max abundance: {analyzer.df['Sample_A_Stems'].max():.2f}%")
            
            print(f"\nSeeds (Sample B) Statistics:")
            print(f"  Mean abundance: {analyzer.df['Sample_B_Seeds'].mean():.2f}%")
            print(f"  Std deviation: {analyzer.df['Sample_B_Seeds'].std():.2f}%")
            print(f"  Max compound: {analyzer.df.loc[analyzer.df['Sample_B_Seeds'].idxmax(), 'Compound']}")
            print(f"  Max abundance: {analyzer.df['Sample_B_Seeds'].max():.2f}%")
            
            print(f"\nCompound Class Distribution:")
            class_counts = analyzer.df['Compound_Class'].value_counts()
            for class_name, count in class_counts.items():
                print(f"  {class_name}: {count} compounds")
        
        elif extra_choice == '3':
            try:
                filename = f"jatropha_curcas_analysis_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.csv"
                csv_path = os.path.join(directories['data'], filename)
                analyzer.df.to_csv(csv_path, index=False)
                all_saved_files.append(csv_path)
                print(f"✅ Data exported to: {csv_path}")
            except Exception as e:
                print(f"❌ Export failed: {e}")
        
        elif extra_choice == '4':
            print("\n📖 Creating comprehensive PDF report...")
            pdf_path = analyzer.create_comprehensive_pdf_report(noise_level, directories)
            if pdf_path:
                all_saved_files.append(pdf_path)
                print("✅ Comprehensive PDF report completed!")
        
        else:
            print("⏭️ Skipping additional analysis")
    
    except KeyboardInterrupt:
        print("\n\n⏹️ Additional analysis interrupted by user")
    except Exception as e:
        print(f"\n⚠️ Error in additional analysis: {e}")
    
    # Final summary
    print("\n" + "="*60)
    print("🎉 JATROPHA CURCAS GC-MS ANALYSIS COMPLETE!")
    print("="*60)
    print("✅ Ultra-realistic chromatograms generated with heavy noise")
    print("✅ Comprehensive analysis report provided")
    print("✅ Multiple visualization options completed")
    if all_saved_files:
        print(f"✅ {len(all_saved_files)} files saved successfully:")
        for file_path in all_saved_files:
            print(f"   📁 {file_path}")
    print("💡 Noise patterns match real-world GC-MS data characteristics")
    print("🔬 Ready for further analysis or method development")
    print("="*60)

if __name__ == "__main__":
    """Execute the main function when script is run directly"""
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n👋 Program terminated by user. Goodbye!")
    except Exception as e:
        print(f"\n💥 CRITICAL ERROR: {e}")
        print("Please check your Python environment and dependencies.")
    finally:
        print("\n🔚 Program ended.")