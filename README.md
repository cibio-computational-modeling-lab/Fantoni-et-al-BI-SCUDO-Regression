# BI-SCUDO Regression - A novel algorithm to infer clinical variables from omics datasets
**Authors:** Giacomo Fantoni, Mario Lauria, Carlo Alberto Rossi, Roberto Bizzotto and Luca Marchetti

Code for the paper _TODO: add title and DOI once available._

This repository contains the complete Python implementation of **BI-SCUDO regression**, an extension of the SCUDO (Signature-based Clustering for Diagnostic Purposes) framework from classification to the prediction of **continuous variables**. Starting from an expression matrix (e.g. transcriptomics or proteomics) and a continuous metadata variable, the pipeline:

1. Filters features by their association with the extremes of the variable of interest (Wilcoxon rank-sum, t-test or data-driven penalties).
2. Builds a per-sample **rank-based signature** made of the top `n1` up-regulated and bottom `n2` down-regulated features.
3. Computes a **weighted enrichment-score distance** between every pair of samples.
4. Predicts the variable of interest of each sample as the distance-weighted mean of its `k` nearest neighbours.
5. Optionally tunes all the hyperparameters (filtering threshold, `n1`, `n2`, percentile split, `k`) with a **genetic algorithm** that minimizes cross-validated MSE while favouring compact biomarkers.

The resulting set of features retained in the signatures is returned as the candidate **biomarker**.

---

## ⚙️ INSTALLATION

### With pip

Pre-compiled wheels are available for Linux (x86_64, aarch64), macOS (Intel, Apple Silicon) and Windows (x86_64), for Python 3.10–3.14. Download the one matching your system and Python version (e.g. `cp312` for Python 3.12) and install it with:

```bash
git clone https://github.com/cibio-computational-modeling-lab/Fantoni-et-al-BI-SCUDO-Regression
cd Fantoni-et-al-BI-SCUDO-Regression
pip install <platform>/scudo_regression-<version>-<python>-<platform>.whl
```


## 📥 INPUT FORMAT

- **Expression data** (`data_file`): CSV with features (probes/genes/proteins) as **rows** and samples as **columns**. The first column contains feature identifiers.
- **Metadata** (`metadata_file`): CSV with one row per sample. The first column contains sample identifiers matching the expression data columns; the column named as `variable_of_interest` holds the continuous target. Samples with a missing target are dropped.

Only the samples present in both files are used.

---

## 🚀 HOW TO RUN AN ANALYSIS

```bash
scudo-regression <output_dir> <data_file> <metadata_file> <variable_of_interest> [options]
```

or equivalently `python -m scudo_regression ...`, or with Docker:

```bash
docker run --rm -v "$PWD:/data" scudo-regression \
    /data/results /data/expression.csv /data/metadata.csv age
```

### Example: reproducing `test_data/results`

`test_data/` contains a synthetic training cohort (`training_data.csv`, `training_metadata.csv`) and an independent validation cohort (`validation_data.csv`, `validation_metadata.csv`), each with 1000 features and 150 samples. The continuous target is the `to_be_predicted` column. The output of a complete run is in `test_data/results/`.

#### Without re-running the optimization

Pass the optimized parameters from `test_data/results/parameters.csv` as fixed values. This reproduces the predictions, signatures, biomarker and metrics in `test_data/results/`:

```bash
scudo-regression results/ \
    test_data/training_data.csv test_data/training_metadata.csv to_be_predicted \
    --percentile 0.1 \
    --filtering-threshold 0.04118778943423801 --n1 23 --n2 25 --n-neighbors 3 \
    --validation-data test_data/validation_data.csv \
    --validation-metadata test_data/validation_metadata.csv \
    --permutations 1000
```

#### Re-running the optimization

The results were obtained with the genetic algorithm searching `n1`/`n2` in `[0, 150]` (`min_n=0`, `max_n=150`). These bounds are not exposed by the CLI, so set them through `ScudoConfig`:

```python
from scudo_regression.data import ScudoConfig
from scudo_regression.io import build_expression_data, save_scudo_results
from scudo_regression.pipeline import ScudoPipeline
from scudo_regression.validation.independent import validate_independent
from scudo_regression.validation.permutation import permutation_test

# The guard is required: the genetic algorithm runs in parallel with multiprocessing.
if __name__ == "__main__":
    config = ScudoConfig(variable_of_interest="to_be_predicted", min_n=0, max_n=150)

    training = build_expression_data(
        "test_data/training_data.csv", "test_data/training_metadata.csv", "to_be_predicted"
    )
    results = ScudoPipeline(config).run(training)

    validation_data = build_expression_data(
        "test_data/validation_data.csv", "test_data/validation_metadata.csv", "to_be_predicted"
    )
    validation = validate_independent(
        results.reduced_data,
        results.starting_data.expression_data.metadata,
        results.starting_data.expression_data.subjects,
        results.optimization_results,
        validation_data,
        results.optimization_results.n_neighbors,
    )

    p_value = permutation_test(results, n_permutations=1000, n_jobs=config.n_jobs)

    save_scudo_results("results/", results, validation, p_value=p_value)
```

The genetic algorithm is not seeded, so a new optimization can end at slightly different parameters and metrics than those in `test_data/results/`.

The same optimization with the CLI's default bounds (`n1`/`n2` in `[10, 250]`) is simply:

```bash
scudo-regression results/ \
    test_data/training_data.csv test_data/training_metadata.csv to_be_predicted \
    --validation-data test_data/validation_data.csv \
    --validation-metadata test_data/validation_metadata.csv \
    --permutations 1000
```

### Options

| Option | Default | Description |
|---|---|---|
| `--logarithmic {true,false}` | `false` | Whether data values are already log-transformed |
| `--percentile` | `0.1` | Percentile used to define the `High`/`Low` classes |
| `--cv` | `10` | Number of cross-validation folds (`0` to skip) |
| `--filtering-test` | `Wilcoxon` | One of `Wilcoxon`, `t-test`, `DataDriven`, `none` |
| `--filtering-threshold` | _optimized_ | Fixed p-value threshold for feature filtering |
| `--n1` / `--n2` | _optimized_ | Fixed sizes of the up-/down-regulated signature |
| `--n-neighbors` | `5` | Number of neighbours used for prediction |
| `--validation-data` / `--validation-metadata` | – | Independent validation cohort |
| `--permutations` | `0` | Number of permutations for significance testing |
| `--n-jobs` | `-1` | Parallel workers (`-1` = all cores) |

If any of `--filtering-threshold`, `--n1` or `--n2` is omitted, the genetic algorithm optimizes the filtering threshold, `n1`, `n2`, the percentile (among `0.1, 0.15, 0.2, 0.25, 0.3`) and the number of neighbours (`3`–`15`). Genetic algorithm settings (population size, generations, elitism, …) can be changed through `ScudoConfig` in `data.py`.

### Using it as a library

```python
from scudo_regression.data import ScudoConfig
from scudo_regression.io import build_expression_data, save_scudo_results
from scudo_regression.pipeline import ScudoPipeline

if __name__ == "__main__":
    expr = build_expression_data("expression.csv", "metadata.csv", "age")
    results = ScudoPipeline(ScudoConfig(variable_of_interest="age")).run(expr)
    save_scudo_results("results/", results)
```

---

## 📤 OUTPUT

All results are written to `output_dir`:

- `parameters.csv`: Final (fixed or optimized) parameters.
- `biomarker.csv`: Features included in the signatures, with their p-value and weight.
- `signatures.csv`: Per-sample ranked signatures.
- `distance_matrix.csv`: Sample × sample distance matrix.
- `predictions.csv`: Observed vs. predicted values on the training set.
- `metrics.csv`: Cross-validated R², correlation, MSE, RMSE, RPD, elapsed time and (if requested) permutation p-value.
- `cv_fold_<k>/`: Predictions and metrics for each cross-validation fold.
- `validation/`: Predictions, metrics and signatures on the independent validation cohort (if provided).

---


For any questions or issues, feel free to open an [Issue](https://github.com/cibio-computational-modeling-lab/Fantoni-et-al-BI-SCUDO-Regression/issues) on this repository! 🎯
