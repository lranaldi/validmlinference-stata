"""build_data.py — one-shot script that converts the upstream Python-package
datasets into Stata-friendly files for the ValidMLInference Stata port.

Run from the repository root:

    C:\\Users\\loren\\.conda\\envs\\vmli\\python.exe ValidMLInference-stata\\data\\build_data.py

Outputs (relative to this script's directory):
    remote_work_data.dta              full remote-work job-postings sample
    topic_model/estimation.dta         916 obs of (ly, q1..q11)
    topic_model/theta_full.csv         916 x 2 document-topic shares (full)
    topic_model/theta_samp.csv         916 x 2 document-topic shares (10% subsample)
    topic_model/beta_full.csv          2 x V topic-word distributions (full)
    topic_model/beta_samp.csv          2 x V topic-word distributions (10% subsample)
    topic_model/lda_data.csv           n_lda x 2 (C_i for full and subsample)
    topic_model/gamma_draws.csv        MCMC draws for the "Joint" comparison row

All upstream data is public.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ValidMLInference import remote_work_data, topic_model_data


def _ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _write_remote_work(out_dir: Path) -> None:
    df = remote_work_data()
    df = df.drop(columns=[c for c in df.columns if c.startswith("Unnamed")])
    # pandas writes Stata strL for object dtypes by default; the explicit
    # convert_dates={} suppresses datetime guesses. version=118 -> Stata 14+
    df.to_stata(
        out_dir / "remote_work_data.dta",
        write_index=False,
        version=118,
        variable_labels={
            "city_name": "City name",
            "naics_2022_2": "NAICS-2022 industry code (2-digit)",
            "salary": "Annual salary (USD)",
            "wfh_rwo": "AI/ML-generated WFH indicator",
            "soc_2021_2": "SOC-2021 occupation code (2-digit)",
            "employment_type_name": "Employment type (full/part-time)",
        },
    )


def _write_topic_model(out_dir: Path) -> None:
    sub = _ensure_dir(out_dir / "topic_model")
    data = topic_model_data()

    Y = np.asarray(data["estimation_data"]["ly"]).reshape(-1)
    Z = np.asarray(data["covars"])
    if Z.shape[0] != Y.shape[0]:
        raise RuntimeError(f"Y has {Y.shape[0]} rows but Z has {Z.shape[0]}")

    cols = {"ly": Y}
    for j in range(Z.shape[1]):
        cols[f"q{j + 1}"] = Z[:, j]
    df = pd.DataFrame(cols)
    df.to_stata(sub / "estimation.dta", write_index=False, version=118)

    np.savetxt(sub / "theta_full.csv", np.asarray(data["theta_est_full"]), delimiter=",")
    np.savetxt(sub / "theta_samp.csv", np.asarray(data["theta_est_samp"]), delimiter=",")
    np.savetxt(sub / "beta_full.csv", np.asarray(data["beta_est_full"]), delimiter=",")
    np.savetxt(sub / "beta_samp.csv", np.asarray(data["beta_est_samp"]), delimiter=",")
    np.savetxt(sub / "lda_data.csv", np.asarray(data["lda_data"]), delimiter=",")
    np.savetxt(sub / "gamma_draws.csv", np.asarray(data["gamma_draws"]), delimiter=",")


def main() -> None:
    here = Path(__file__).resolve().parent
    _ensure_dir(here)
    print(f"[build_data] writing remote-work data to {here}")
    _write_remote_work(here)
    print(f"[build_data] writing topic-model data to {here / 'topic_model'}")
    _write_topic_model(here)
    print("[build_data] done.")


if __name__ == "__main__":
    main()
