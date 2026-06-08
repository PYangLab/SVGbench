[![DOI](https://zenodo.org/badge/672695737.svg)](https://zenodo.org/doi/10.5281/zenodo.10295502)

# SVGbench

SVGbench contains the source code used to benchmark spatially variable gene
(SVG) detection methods for spatial transcriptomics data.

The repository accompanies:

> Chen, C., Kim, H. J. & Yang, P. Evaluating spatially variable gene detection
> methods for spatial transcriptomics data. *Genome Biology* 25, 18 (2024).
> <https://doi.org/10.1186/s13059-023-03145-y>

![Fig1_evaluation workflow_v4](https://github.com/PYangLab/SVGbench/assets/10148940/9e39d918-5a7c-4b3f-8068-24cee617c533)

## What Is Included

This is a manuscript reproducibility repository. It provides scripts for:

- preprocessing spatial transcriptomics benchmark datasets;
- running SVG detection methods including Giotto, MERINGUE, Moran's I, nnSVG,
  SOMDE, SpatialDE and SPARK-X;
- comparing SVG calls across methods and datasets;
- testing robustness to gene and spot subsampling;
- simulating spatial expression patterns and evaluating accuracy;
- measuring computational time and memory;
- generating paper and supplementary figures.

## Data

Simulation data and raw gene expression matrices with spatial coordinates are
available from the project Google Drive folder:

<https://drive.google.com/drive/folders/1ECInc9LNjTjK-8Oa8CR8eyvnc4pDepXl?usp=sharing>

Most scripts expect each dataset to contain:

- `Spatial_count.csv`: spatial locations or spots by genes before transposition
  in most R scripts;
- `Locations.csv`: spatial coordinates and dataset-specific metadata.

The original scripts use relative paths such as
`../../../SpatialData/SpatialBenchmark/DataUploadSubset/`. For a new machine,
copy `config.yml.example` to `config.yml` and use it as a guide for arranging
or symlinking the downloaded data.

## Repository Layout

| Folder | Purpose |
| --- | --- |
| `01_preprocessing/` | Prepare benchmark datasets from raw count and coordinate files. |
| `02_run methods/` | Run each SVG detection method on benchmark, simulation, subsampling and MOSTA datasets. |
| `03_SVG detection/` | Compare method outputs, p-value behaviour and overlap between significant SVGs. |
| `04_gene interdependency/` | Assess dependence between gene-level statistics and method rankings. |
| `05_robustness to sparsity/` | Evaluate robustness to gene and spot subsampling. |
| `06_downstream clustering/` | Test SVG sets in downstream spatial-domain clustering workflows. |
| `07_simulations/` | Generate simulations and calculate accuracy against ground truth. |
| `08_computational resource/` | Plot computational time and memory usage. |
| `09_additional figures/` | Generate supporting and supplementary figures. |

## Setup

Install R dependencies with:

```r
source("install_R_dependencies.R")
```

Create the Python environment with:

```bash
conda env create -f environment.yml
conda activate svgbench
```

Some benchmarked methods have complex or version-sensitive installation
requirements. The setup files install the packages used most commonly across the
scripts, but method-specific packages such as Giotto, SOMDE, SpatialDE and
SpaGCN may require following the upstream installation instructions for the
platform being used.

## Running The Benchmark

The full benchmark is computationally intensive. A typical workflow is:

1. Download the data and arrange it according to `config.yml.example`.
2. Run preprocessing scripts in `01_preprocessing/`.
3. Run method-specific scripts in `02_run methods/`.
4. Use scripts in `03_SVG detection/` to compare SVG calls across methods.
5. Run robustness, simulation, clustering and resource analyses as needed.
6. Use `docs/reproduce_figures.md` to map manuscript figures to scripts.

Example for one R-based method:

```bash
Rscript "02_run methods/run_nnSVG/callSVG_nnSVG.R"
```

Example for one Python-based method:

```bash
python "02_run methods/run_spatialde/callSVG_spatialDE.py"
```

## Citation

If you use SVGbench, please cite:

```bibtex
@article{Chen2024SVGbench,
  title = {Evaluating spatially variable gene detection methods for spatial transcriptomics data},
  author = {Chen, Carissa and Kim, Hani Jieun and Yang, Pengyi},
  journal = {Genome Biology},
  year = {2024},
  volume = {25},
  number = {18},
  doi = {10.1186/s13059-023-03145-y},
  url = {https://doi.org/10.1186/s13059-023-03145-y}
}
```

## License

This repository is distributed under the CC-BY-4.0 license.
