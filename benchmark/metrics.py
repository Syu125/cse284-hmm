from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd


@dataclass
class ComparisonSummary:
    n_sites: int
    concordance: float
    cohen_kappa: float
    switches_per_mb_model: float
    switches_per_mb_ref: float
    median_tract_bp_model: float
    median_tract_bp_ref: float
    yri_pct_model: float
    yri_pct_ref: float


def normalize_labels(series: pd.Series, label_map: Dict[str, str] | None = None) -> pd.Series:
    if label_map is None:
        label_map = {
            "0": "YRI",
            "1": "CEU",
            "AFR": "YRI",
            "EUR": "CEU",
            "Y": "YRI",
            "C": "CEU",
        }

    normalized = series.astype(str).str.upper().map(lambda value: label_map.get(value, value))
    return normalized


def cohen_kappa(y_true: Sequence[str], y_pred: Sequence[str]) -> float:
    if len(y_true) == 0:
        return float("nan")

    labels = sorted(set(y_true) | set(y_pred))
    if len(labels) == 1:
        return 1.0

    label_to_idx = {label: i for i, label in enumerate(labels)}
    matrix = np.zeros((len(labels), len(labels)), dtype=float)

    for truth, pred in zip(y_true, y_pred):
        matrix[label_to_idx[truth], label_to_idx[pred]] += 1

    total = matrix.sum()
    p0 = np.trace(matrix) / total
    row_marginal = matrix.sum(axis=1) / total
    col_marginal = matrix.sum(axis=0) / total
    pe = float(np.sum(row_marginal * col_marginal))

    if np.isclose(1.0 - pe, 0.0):
        return 0.0

    return (p0 - pe) / (1.0 - pe)


def switches_per_mb(positions: Sequence[int], labels: Sequence[str]) -> float:
    if len(positions) < 2:
        return float("nan")

    switches = 0
    for i in range(1, len(labels)):
        if labels[i] != labels[i - 1]:
            switches += 1

    span_bp = int(positions[-1]) - int(positions[0])
    if span_bp <= 0:
        return float("nan")

    span_mb = span_bp / 1_000_000
    return switches / span_mb


def tract_lengths_bp(positions: Sequence[int], labels: Sequence[str]) -> List[int]:
    if len(positions) == 0:
        return []

    tracts: List[int] = []
    start_idx = 0

    for i in range(1, len(labels)):
        if labels[i] != labels[i - 1]:
            tracts.append(int(positions[i - 1]) - int(positions[start_idx]))
            start_idx = i

    tracts.append(int(positions[-1]) - int(positions[start_idx]))
    return tracts


def global_percentage(labels: Sequence[str], target: str = "YRI") -> float:
    if len(labels) == 0:
        return float("nan")
    return 100.0 * sum(label == target for label in labels) / len(labels)


def compare_sample(aligned_df: pd.DataFrame, model_col: str, ref_col: str) -> ComparisonSummary:
    positions = aligned_df["position"].to_numpy()
    model = aligned_df[model_col].to_numpy()
    ref = aligned_df[ref_col].to_numpy()

    concordance = float(np.mean(model == ref)) if len(aligned_df) else float("nan")
    kappa = cohen_kappa(ref.tolist(), model.tolist())

    model_tracts = tract_lengths_bp(positions, model)
    ref_tracts = tract_lengths_bp(positions, ref)

    return ComparisonSummary(
        n_sites=len(aligned_df),
        concordance=concordance,
        cohen_kappa=kappa,
        switches_per_mb_model=switches_per_mb(positions, model),
        switches_per_mb_ref=switches_per_mb(positions, ref),
        median_tract_bp_model=float(np.median(model_tracts)) if model_tracts else float("nan"),
        median_tract_bp_ref=float(np.median(ref_tracts)) if ref_tracts else float("nan"),
        yri_pct_model=global_percentage(model, "YRI"),
        yri_pct_ref=global_percentage(ref, "YRI"),
    )


def align_predictions(
    model_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    sample_col: str,
    pos_col: str,
    model_label_col: str,
    ref_label_col: str,
    valid_labels: Tuple[str, ...] = ("YRI", "CEU", "HET"),
) -> pd.DataFrame:
    model = model_df[[sample_col, pos_col, model_label_col]].rename(
        columns={sample_col: "sample_id", pos_col: "position", model_label_col: "model_label"}
    )
    ref = ref_df[[sample_col, pos_col, ref_label_col]].rename(
        columns={sample_col: "sample_id", pos_col: "position", ref_label_col: "ref_label"}
    )

    model["model_label"] = normalize_labels(model["model_label"])
    ref["ref_label"] = normalize_labels(ref["ref_label"])

    aligned = model.merge(ref, on=["sample_id", "position"], how="inner")
    if valid_labels:
        valid = set(valid_labels)
        aligned = aligned[
            aligned["model_label"].isin(valid) & aligned["ref_label"].isin(valid)
        ]
    return aligned.sort_values(["sample_id", "position"]).reset_index(drop=True)
