## HiCHub2: Hi-C Interaction Hub Analysis Toolkit
HiCHub2 is a Python-based CLI tool designed to identify, compare, and cluster high-frequency interaction hubs from Hi-C data. It utilizes community detection algorithms (Leiden) to locate dense interaction regions and K-means clustering to reveal structural dynamics across multiple samples.

## Core Features
singlehub: Identify static interaction hubs in a single sample.

pairhub: Perform differential analysis between two samples to find significantly changing hubs.

varhub: Identify hubs with high variance (dynamic fluctuations) across multiple samples.

cluster_hubs: Perform standard Z-score normalization and K-means clustering on existing hub results.

cluster_pairhubs: Full-auto pipeline. Merges multiple pairwise BED files, extracts intensities from .hic files for all samples, and performs final clustering.

## Installation
Ensure you have the bedtools suite installed on your system, then install the following Python dependencies:

```
pip install hicstraw igraph leidenalg pybedtools qnorm scikit-learn seaborn pandas numpy scipy
```

## Command Details

1. singlehub (Single Sample Identification)Used for initial exploration of interaction hubs in a single cell line or tissue.
```
python hichub2.py singlehub \
    --path ./data/ \
    --name CellA \
    --resolution 10000 \
    --output CellA_hubs.bed
```
2. pairhub (Two-Sample Differential)Uses Wilcoxon rank-sum tests to find differential interaction regions between two samples.
```
python hichub2.py pairhub \
    --path ./data/ \
    --name CellA CellB \
    --qval 0.001 \
    --output AB_diff.bed
```
3. varhub (Multi-Sample Variance)Identifies hubs that show high variance across all input samples.
```
python hichub2.py varhub \
    --path ./data/ \
    --name CellA CellB CellC \
    --pval 0.05 \
    --output multi_var.bed
```
4. cluster_hubs (Simple Clustering)Clusters an existing output file and generates a heatmap.
```
python hichub2.py cluster_hubs \
    --input hubs.bed \
    --name CellA CellB CellC \
    --n_clusters 8 \
    --plot_out heatmap.png
```
5. cluster_pairhubs (Integrated Merge & Cluster)A powerful sub-command that automatically handles input BED files located in the --path.
```
python hichub2.py cluster_pairhubs \
    --path /home/manager/hic_folder/ \
    --name GM12878 K562 HepG2 \
    --input_beds GM_K562_diff.bed K562_HepG2_diff.bed \
    --n_clusters 6 \
    --threads 23 \
    --output final_combined_results.bed \
    --plot_out cluster_heatmap.png
```
## Parameter List

Argument	Description	Default
--path	Base path for .hic files and input .bed files	Required
--name	List of sample names (must match .hic filenames)	Required
--resolution	Hi-C resolution in base pairs (BP)	10000
--norm	Hi-C normalization method (NONE, VC, KR, SCALE, etc.)	NONE
--data_type	Data extraction type (oe: observed/expected, observed)	oe
--distance	Max interaction span in bins (e.g., 200 = 2Mb at 10kb)	200
--score_thres	Threshold for interaction score filtering	2.0
--res_community	Resolution for Leiden algorithm (higher = more fragmented hubs)	0.6
--gap_size	Max gap (in bins) allowed when stitching adjacent regions	2
--threads	Number of CPU cores for parallel processing	23
--n_clusters	Number of clusters for K-means	6


## Output Description
BED File:

chr1, start1, end1, chr2, start2, end2: Genomic coordinates of the Hub.

SampleName: The average OE interaction intensity for that sample within the Hub.

cluster: The assigned K-means cluster ID.

Heatmap (PNG):

Visualizes the clusters. The Y-axis represents Hubs, and the X-axis represents Samples. Colors indicate relative intensity changes after Z-score normalization.

## Important Notes
Memory Management: When running cluster_pairhubs, reading multiple large .hic files simultaneously consumes significant RAM. Adjust --threads based on your system capacity.

Pathing: cluster_pairhubs automatically looks for the files listed in --input_beds within the directory specified by --path. Ensure your BED files are placed there before running.

