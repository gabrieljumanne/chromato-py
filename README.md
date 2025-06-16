
# 📊 GC-MS Chromatogram Analyzer

A **Python-based tool** for processing, analyzing, and visualizing gas chromatography–mass spectrometry (GC-MS) chromatogram data. Designed for flexibility and ease of use, this tool empowers researchers, analysts, and students to work efficiently with a wide range of GC-MS datasets. It offers exportable visual reports, customizable data processing, and smooth integration into other scientific workflows.

---

## 🧪 Features

- **Load Chromatogram Data:** Import chromatogram data directly from `.csv` files supporting two-column formats (RetentionTime, Intensity).
- **Data Visualization:** Generate clear, publication-ready chromatogram plots with support for optional noise reduction and smoothing.
- **Export & Reporting:** Save plots and comprehensive reports as PDF documents, timestamped and ready for sharing or archiving.
- **Multi-Sample Comparison:** Compare chromatograms across multiple samples in a single workflow—either by providing several CSV files or by merging data.
- **Customizable Workflow:** Modular code structure allows easy parameter tuning (e.g., smoothing, plot aesthetics) and extension or integration into other pipelines.
- **Lightweight & Dependency-Minimal:** Only relies on a handful of essential Python libraries for efficient performance.

---

## 📁 Project Structure

```
├── chromatograms/                # Output folder for reports/plots
│   └── pdf/                      # Exported chromatogram PDFs
├── gc_ms_chromatogram.py        # Main script for GC-MS analysis
├── gc_ms_chromatogram_data.py   # Helper functions for data handling
├── jatropha_gcms_data.csv       # Example dataset (can be replaced)
├── sample_a_stems.csv           # Sample chromatogram file
├── sample_b_seeds.csv           # Another sample chromatogram file
├── requirements.txt             # Python dependencies
├── .gitignore
└── README.md
```

---

## 🚀 Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/gcms-chromatogram-analyzer.git
cd gcms-chromatogram-analyzer
```

### 2. Set Up Your Python Environment

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

---

## 🔧 Usage

Run the main script to process and visualize chromatogram data:

```bash
python gc_ms_chromatogram.py
```

- **Default Input:** `jatropha_gcms_data.csv`
- **Output:** Timestamped PDF files in `chromatograms/pdf/` (e.g. `chromatogram_report_noise1.4_20250616_022535.pdf`)

### Command-Line Options

- **Custom Input:** Specify a different CSV file or multiple files as input.
- **Noise/Smoothing:** Adjust noise reduction and smoothing parameters via script arguments or by modifying the script.
- **Output Location:** Change or customize the output directory as needed.

---

## 📈 Output Example

After running the script, you will find files such as:

```
chromatograms/pdf/chromatogram_report_noise1.4_20250616_022535.pdf
```

Each report contains:
- The processed chromatogram plot
- Details of the applied smoothing/noise reduction
- Sample metadata (if present)
- Comparison plots (if multiple samples provided)

---

## 📦 Dependencies

All dependencies are managed via `requirements.txt`:

- [`pandas`](https://pandas.pydata.org/): Data manipulation and CSV handling
- [`numpy`](https://numpy.org/): Numerical operations
- [`matplotlib`](https://matplotlib.org/): Plotting and visualization
- [`scipy`](https://scipy.org/): Signal processing (smoothing, noise reduction)

Install them automatically with:

```bash
pip install -r requirements.txt
```

---

## 🧬 Data Format

**Expected CSV Structure:**

```
RetentionTime,Intensity
0.01,12.3
0.02,15.8
...
```

- Each row represents a retention time point and its corresponding intensity.
- Additional columns or sample metadata can be ignored or handled as needed.
- For comparing samples, supply additional CSV files or merge data columns manually.

---

## 🛠️ Customization

The analyzer is designed to be **modular and extensible**. You can easily:

- **Adjust Smoothing/Noise Reduction:** Edit parameters in `gc_ms_chromatogram.py` (e.g., Gaussian kernel width, Savitzky–Golay window).
- **Change Plot Styles:** Modify matplotlib code to customize colors, labels, titles, and export options.
- **Integrate with Pipelines:** Import functions from `gc_ms_chromatogram_data.py` in your own scripts or notebooks.
- **Automate Batch Processing:** Adapt the script for high-throughput analysis or integration into laboratory automation systems.

---

## 💡 Tips & Best Practices

- For best results, pre-process raw GC-MS files to remove obvious artifacts and ensure correct formatting.
- Use virtual environments to isolate dependencies.
- Regularly export and backup your chromatogram reports for reproducibility and compliance.

---

## 📄 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## 👤 Developer

**Gabriel Wambura**  
Python Developer | Data Tools Builder

Feel free to fork, adapt, and contribute!  
Pull requests, feature suggestions, and bug reports are welcome.

---
