"""Shared statistical utilities for the manuscript figure scripts.

Centralises the bootstrap helpers so every figure uses the same stratified-
bootstrap implementation, the same paired ΔAUC bootstrap, and the same
AUC<0.5 sanity warning:
  - stratified bootstrap (resample positives + negatives separately)
  - paired bootstrap for ΔAUC inference
  - AUC<0.5 sanity warning
"""

from __future__ import annotations
import warnings
import numpy as np
from sklearn.metrics import roc_auc_score


def stratified_bootstrap_indices(y, rng):
    """Resample positives and negatives separately, preserving class
    proportions in every iteration. Returns the concatenated index vector."""
    y = np.asarray(y).astype(int)
    pos = np.flatnonzero(y == 1)
    neg = np.flatnonzero(y == 0)
    return np.concatenate([
        rng.choice(pos, size=len(pos), replace=True),
        rng.choice(neg, size=len(neg), replace=True),
    ])


def auc_with_ci(z, y, n_boot=2000, seed=20260427, label=""):
    """AUC + stratified-bootstrap 95% CI.

    Warn if AUC < 0.5 (signals predictor direction is flipped for this
    endpoint and the score may be inversely associated with it).
    """
    z = np.asarray(z, dtype=float)
    y = np.asarray(y).astype(int)
    rng = np.random.default_rng(seed)
    auc = roc_auc_score(y, z)
    if auc < 0.5:
        warnings.warn(
            f"[auc_with_ci{(' ' + label) if label else ''}] "
            f"AUC = {auc:.3f} < 0.5 — predictor may be inversely associated "
            f"with the endpoint. Inspect direction.")
    aucs = np.full(n_boot, np.nan)
    for b in range(n_boot):
        idx = stratified_bootstrap_indices(y, rng)
        if len(np.unique(y[idx])) < 2:
            continue
        try:
            aucs[b] = roc_auc_score(y[idx], z[idx])
        except Exception:
            pass
    aucs = aucs[~np.isnan(aucs)]
    if len(aucs) < 30:
        return auc, np.nan, np.nan
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return auc, lo, hi


def paired_delta_auc_ci(z1, z2, y, n_boot=2000, seed=20260427, label=""):
    """Paired stratified bootstrap of ΔAUC = AUC(z1) - AUC(z2) using the
    same resampled indices for both predictors.

    Returns (delta_obs, lo, hi, p_two_sided, n_valid).
    """
    z1 = np.asarray(z1, dtype=float)
    z2 = np.asarray(z2, dtype=float)
    y = np.asarray(y).astype(int)
    rng = np.random.default_rng(seed)
    delta_obs = roc_auc_score(y, z1) - roc_auc_score(y, z2)
    deltas = np.full(n_boot, np.nan)
    for b in range(n_boot):
        idx = stratified_bootstrap_indices(y, rng)
        if len(np.unique(y[idx])) < 2:
            continue
        try:
            deltas[b] = (roc_auc_score(y[idx], z1[idx])
                          - roc_auc_score(y[idx], z2[idx]))
        except Exception:
            pass
    deltas = deltas[~np.isnan(deltas)]
    if len(deltas) < 30:
        return delta_obs, np.nan, np.nan, np.nan, int(len(deltas))
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    if delta_obs >= 0:
        p_two = 2 * np.mean(deltas <= 0)
    else:
        p_two = 2 * np.mean(deltas >= 0)
    p_two = min(max(p_two, 1.0 / max(len(deltas), 1)), 1.0)
    return delta_obs, lo, hi, p_two, int(len(deltas))
