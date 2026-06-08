# Reproducing SVGbench Analyses

This guide maps the major analyses in the SVGbench paper to repository scripts.
The full benchmark requires the downloaded datasets, intermediate result files
and method-specific software environments.

## Data Preparation

| Analysis | Scripts |
| --- | --- |
| General benchmark dataset preprocessing | `01_preprocessing/preprocess_benchmarking_datasets.R` |
| Slide-seqV1 preprocessing | `01_preprocessing/preprocess_slideseqV1.R` |
| Slide-seqV2 preprocessing | `01_preprocessing/preprocess_slideseqV2.R` |
| DBiT-seq preprocessing | `01_preprocessing/preprocess_DBiT-seq.R` |
| SM-Omics preprocessing | `01_preprocessing/preprocess_SM-Omics.R` |

## Running SVG Detection Methods

| Method | Benchmark datasets | Simulations | Gene/spot subsampling | MOSTA |
| --- | --- | --- | --- | --- |
| Giotto | `02_run methods/run_giotto/callSVG_giotto.R` | `02_run methods/run_giotto/simulations_giotto.R` | `02_run methods/run_giotto/subsample_genes_giotto.R`, `02_run methods/run_giotto/subsample_spots_giotto.R` | `02_run methods/run_giotto/mosta_giotto.R` |
| MERINGUE | `02_run methods/run_meringue/callSVG_meringue.R` | `02_run methods/run_meringue/simulations_meringue.R` | `02_run methods/run_meringue/subsample_genes_meringue.R`, `02_run methods/run_meringue/subsample_spots_meringue.R` | `02_run methods/run_meringue/mosta_meringue.R` |
| Moran's I | `02_run methods/run_moransi/callSVG_moransi.R` | `02_run methods/run_moransi/simulations_moransi.R` | `02_run methods/run_moransi/subsample_genes_moransi.R`, `02_run methods/run_moransi/subsample_spots_moransi.R` | `02_run methods/run_moransi/mosta_moransi.R` |
| nnSVG | `02_run methods/run_nnSVG/callSVG_nnSVG.R` | `02_run methods/run_nnSVG/simulations_nnSVG.R` | `02_run methods/run_nnSVG/subsample_genes_nnSVG.R`, `02_run methods/run_nnSVG/subsample_spots_nnSVG.R` | `02_run methods/run_nnSVG/mosta_nnSVG.R` |
| SOMDE | `02_run methods/run_somde/callSVG_somde.py` | `02_run methods/run_somde/simulations_somde.py` | `02_run methods/run_somde/subsample_genes_somde.py`, `02_run methods/run_somde/subsample_spots_somde.py` | `02_run methods/run_somde/mosta_somde.py` |
| SpatialDE | `02_run methods/run_spatialde/callSVG_spatialDE.py` | `02_run methods/run_spatialde/simulations_spatialDE.py` | `02_run methods/run_spatialde/subsample_genes_spatialde.py`, `02_run methods/run_spatialde/subsample_spots_spatialde.py` | `02_run methods/run_spatialde/mosta_spatialDE.py` |
| SPARK-X | `02_run methods/run_sparkx/callSVG_sparkx.R` | `02_run methods/run_sparkx/simulations_sparkx.R` | `02_run methods/run_sparkx/subsample_genes_sparkx.R`, `02_run methods/run_sparkx/subsample_spots_sparkx.R` | `02_run methods/run_sparkx/mosta_sparkx.R` |

## Evaluation And Plotting

| Paper analysis | Scripts |
| --- | --- |
| Pairwise agreement between SVG rankings | `03_SVG detection/plot_pairwise_spearmans.R`, `09_additional figures/plot_pairwise_spearmans_heatmap.R` |
| P-value distribution and significant SVG counts | `03_SVG detection/plot_p_value_proportions.R`, `09_additional figures/plot_number_sigSVGs.R` |
| Overlap between detected SVGs | `03_SVG detection/plot_prop_overlapping_sigSVGs.R`, `09_additional figures/plot_average_jaccard_heatmaps.R` |
| Gene statistic interdependency | `04_gene interdependency/plot_gene_stats_expr_correlation.R`, `04_gene interdependency/plot_gene_sampling_spearmans.R` |
| Robustness to sparsity | `05_robustness to sparsity/plot_prop_of_overlap_and_unique_SVGs.R`, `05_robustness to sparsity/plot_spot_sampling_spearmans.R` |
| Downstream clustering | `06_downstream clustering/perform_kmeans_hc_clustering.R`, `06_downstream clustering/perform_bayesspace_clustering.R`, `06_downstream clustering/perform_spaGCN_clustering.py`, `06_downstream clustering/perform_SINFONIA_clustering.py`, `06_downstream clustering/plot_clustering_concordance.R` |
| Simulation accuracy | `07_simulations/simulate_masks.R`, `07_simulations/calculate_accuracy.R`, `07_simulations/plot_prop_groundtruth_sigSVGs.R`, `09_additional figures/plot_simulations_auroc_curve.R` |
| Computational resources | `08_computational resource/plot_time.R`, `08_computational resource/plot_memory.R` |
| Dataset QC and supplementary panels | scripts in `09_additional figures/` |

## Practical Notes

- Start from a small subset of datasets before running the full benchmark.
- Create output directories before running method scripts, or adapt scripts to
  create them automatically.
- Keep R and Python method environments separate when package dependencies
  conflict.
- Use `config.yml.example` as the single source of truth for local path layout.
