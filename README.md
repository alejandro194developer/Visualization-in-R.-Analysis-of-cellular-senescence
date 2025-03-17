# README for Biol5379 Assessment Code

This README provides an overview of the R code provided for the Biol5379 assessment. The code is designed to analyze a bulk RNA-seq dataset comparing proliferating cells (prolif), replicative senescent cells (senes), and senescent cells with depleted mitochondria (senes_MtD). The goal is to explore the transcriptomic changes associated with senescence and mitochondrial depletion in human IMR90 fibroblasts.

## Overview of the Code

The code is divided into two main scripts:
1. **Main Analysis Script (`main.R`)**: This script performs the bulk of the data analysis, including data loading, preprocessing, visualization, and statistical analysis.
2. **Functions Script (`functions.r`)**: This script contains custom functions used in the main analysis script. These functions are sourced at the beginning of the main script.

### Key Steps in the Analysis

1. **Loading and Preprocessing Data**:
   - The script loads several data files, including:
     - `EM.csv`: Expression matrix.
     - `DE_Senes_MtD_vs_Prolif.csv`, `DE_Senes_MtD_vs_Senes.csv`, `DE_Senes_vs_Prolif.csv`: Differential expression tables.
     - `Human_Background_GRCh38.p13.csv`: Background gene information.
   - Missing values are removed using `na.omit()`.

2. **Data Visualization**:
   - **Density Plots**: The script generates density plots to visualize the distribution of gene expression across samples.
   - **MA Plots**: MA plots are created to compare the log2 fold changes versus the mean expression levels for each differential expression comparison.
   - **PCA Plot**: Principal Component Analysis (PCA) is performed to visualize the clustering of samples based on gene expression.
   - **Volcano Plots**: Volcano plots are generated to highlight significantly up- and down-regulated genes in each comparison.
   - **Pathway Analysis**: The script performs pathway enrichment analysis on up- and down-regulated genes using the `clusterProfiler` package.
   - **Heatmaps**: Heatmaps are created to visualize the expression patterns of significant genes across samples.
   - **Box and Jitter Plots**: These plots are used to compare the expression levels of specific gene signatures across different sample groups.

3. **Differential Expression Analysis**:
   - The script identifies significantly up- and down-regulated genes in each comparison using a p-value cutoff of 0.05 and a log2 fold change cutoff of 1.
   - Three gene signatures are defined based on specific patterns of differential expression across the comparisons:
     - **Signature 1**: Genes up-regulated in senescent cells compared to proliferating cells and down-regulated in senescent cells with depleted mitochondria.
     - **Signature 2**: Genes down-regulated in both senescent cells and senescent cells with depleted mitochondria compared to proliferating cells.
     - **Signature 3**: Genes down-regulated in senescent cells compared to proliferating cells but up-regulated in senescent cells with depleted mitochondria.

4. **Output**:
   - The script generates a series of plots that are saved as PDF files. These plots are intended to be included in the final report to illustrate the findings of the analysis.

### Key Functions Used

- **`facet_density_plot()`**: Generates density plots for gene expression data.
- **`create_master_table()`**: Combines expression data with background gene information.
- **`add_significant_flag_to_de_table()`**: Adds a significance flag to differential expression tables based on p-value and log2 fold change cutoffs.
- **`format_de_tables()`**: Formats differential expression tables for further analysis.
- **`plot_ma()`**: Generates MA plots.
- **`plot_pca()`**: Performs PCA and generates PCA plots.
- **`get_up_regulated_genes_table_from_de()`**: Extracts up-regulated genes from differential expression tables.
- **`get_down_regulated_genes_table_from_de()`**: Extracts down-regulated genes from differential expression tables.
- **`plot_volcano()`**: Generates volcano plots.
- **`pathway_plot()`**: Performs pathway enrichment analysis and generates pathway plots.
- **`plot_heatmap()`**: Generates heatmaps for gene expression data.
- **`plot_box_gitter_plot()`**: Generates box and jitter plots for gene expression data.

### Requirements

- **R Packages**: The script requires several R packages, including `ggplot2`, `reshape2`, `dplyr`, `tibble`, `ggrepel`, `amap`, `clusterProfiler`, `org.Hs.eg.db`, `tidyr`, and `stringr`.
- **Data Files**: The script expects the following input files:
  - `EM.csv`
  - `DE_Senes_MtD_vs_Prolif.csv`
  - `DE_Senes_MtD_vs_Senes.csv`
  - `DE_Senes_vs_Prolif.csv`
  - `Human_Background_GRCh38.p13.csv`

### How to Run the Code

1. Ensure all required R packages are installed.
2. Place the data files in a `data` directory within the working directory.
3. Source the `functions.r` script at the beginning of the main analysis script.
4. Run the main analysis script (`main.R`).

### Output

The script generates several plots, which should be saved as PDF files. These plots are intended to be included in the final report to illustrate the findings of the analysis.

### Notes

- The code is heavily commented to explain each step of the analysis.
- The script uses a p-value cutoff of 0.05 and a log2 fold change cutoff of 1 for identifying significant genes. These cutoffs can be adjusted if needed.
- The pathway analysis uses the `clusterProfiler` package and the `org.Hs.eg.db` database for human gene annotations.

### Conclusion

This script provides a comprehensive analysis of the RNA-seq dataset, focusing on the transcriptomic changes associated with senescence and mitochondrial depletion. The generated plots and analysis results should help in understanding the biological implications of these changes.

---

This README provides a high-level overview of the code. For more detailed information, please refer to the comments within the R scripts.
