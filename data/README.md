# data/

Shipped example datasets, ported from `../../ValidMLInference-main/src/ValidMLInference/data/`. All data is public.

## Contents

| File | Source | Used by |
|---|---|---|
| `remote_work_data.dta` | upstream `remote_work_data.csv` | `examples/remote_work.do`, fixture `remote_work_ols` |
| `topic_model/estimation.dta` | NPZ keys `estimation_data_ly` + `covars` | `examples/topic_model_example.do` |
| `topic_model/theta_full.csv` | NPZ `theta_est_full` (916 × 2) | document-topic share matrix W (full sample) |
| `topic_model/theta_samp.csv` | NPZ `theta_est_samp` (916 × 2) | document-topic share matrix W (10% subsample) |
| `topic_model/beta_full.csv` | NPZ `beta_est_full` (2 × V) | topic-word distribution matrix B (full sample) |
| `topic_model/beta_samp.csv` | NPZ `beta_est_samp` (2 × V) | topic-word distribution matrix B (10% subsample) |
| `topic_model/lda_data.csv` | NPZ `lda_data` | feature-count vector, used to build κ |
| `topic_model/gamma_draws.csv` | NPZ `gamma_draws` | MCMC draws for the joint-estimation comparison row |
| `topic_model/load.do` | (helper) | reads the four matrix CSVs into Stata matrices `W_full`/`W_samp`/`B_full`/`B_samp` + selection matrix `S` and the LDA-count matrix `LDA` |

The Fed-sentiment data shipped with the Python package (`fed_sentiment_meetings.csv`, `fed_sentiment_B.csv`, `fed_sentiment_kappa.csv`) is not unpacked here because none of the three replication notebooks touch it.

## Regenerating from source

The conversion script is `build_data.py`. It uses the upstream Python interfaces, so the `vmli` conda env (Python ≥ 3.10 with `ValidMLInference 1.4.0` installed) is required:

```
C:\Users\loren\.conda\envs\vmli\python.exe ValidMLInference-stata\data\build_data.py
```

The script is a one-shot — outputs are committed to the repo so end users do not need Python to run the do-files.
