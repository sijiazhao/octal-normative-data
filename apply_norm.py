"""Drop-in helper for applying the OCTAL age-stratified norms.

Recommended use: `age_percentile(metric, value, age)`. Returns the participant's
percentile rank among age-matched healthy controls, with the sign flipped for
lower-is-better metrics so a *low* percentile always means impairment regardless
of which metric you're scoring. This single rule works for every OCTAL metric
in the norm file, including the skewed ones (most of them).

Example
-------
    from apply_norm import age_percentile
    pct = age_percentile("DSST_RT", value=4.2, age=65)    # → ~20 (slower than 80% of HC)
    pct = age_percentile("OIS_LTM_Acc", value=55, age=65) # → ~37 (lower than ~63% of HC)

Vectorised over a DataFrame:
    import pandas as pd
    df["pct_DSST_RT"] = age_percentile("DSST_RT", df["DSST_RT"], df["age"])

Advanced (z-score)
------------------
For the small subset of metrics whose age-detrended distribution is approximately
normal (see `Octal_Normality_Summary.csv`), `age_z_score(metric, value, age)` is
also available. It uses the raw or log scale automatically per metric, and
sign-flips so higher z = better. For most metrics, the percentile approach is
the safer recommendation — z-score on a skewed metric undercovers the tails.

Reads `Octal_MetricsByAge_Raw.csv` from the same folder by default.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd

ArrayLike = Union[float, int, np.ndarray, pd.Series]

NORMS_PATH = Path(__file__).parent / "Octal_MetricsByAge_Raw.csv"

# Recommended approach per metric: ("raw", "log", or "percentile")
RECOMMENDED: dict[str, str] = {
    # Accuracy / score metrics — distribution-free is safest (ceiling-skewed)
    "DSST_pHit":            "percentile",
    "OMT_Acc":              "percentile",
    "OMT_TargetDetection":  "percentile",
    "OMT_Guessing":         "percentile",
    "ROCF_CopyScore":       "percentile",
    "ROCF_RecallScore":     "percentile",
    "ROCF_Remember":        "percentile",
    "OIS_STM_Acc":          "percentile",
    "OIS_STM_SemanticAcc":  "percentile",
    "OIS_LTM_SemanticAcc":  "percentile",
    # Approximately normal on raw scale — raw z is fine
    "OIS_LTM_Acc":          "raw",
    "OMT_Misbinding":       "raw",
    "OIS_STM_LocErr":       "raw",
    "OIS_LTM_LocErr":       "raw",
    # RT / TMT — log-transformed z (approximately normal after log)
    "DSST_RT":              "log",
    "OMT_IdeRT":            "log",
    "OMT_LocRT":            "log",
    "TMT_A":                "log",
    "TMT_B":                "log",
    "TMT_BdA":              "log",
    "TMT_average":          "log",
    # Spatial errors — log-normal (raw is right-skewed)
    "OMT_LocErr":           "log",
    "OMT_Imprecision":      "log",
    # Mild-skew count — raw z is acceptable
    "DSST_nCorrectResponse": "raw",
}

# Direction: lower-is-better metrics (sign-flip z so higher z = better)
LOWER_IS_BETTER = {
    "DSST_RT", "OMT_IdeRT", "OMT_LocRT", "TMT_A", "TMT_B", "TMT_BdA",
    "TMT_average", "OIS_STM_LocErr", "OIS_LTM_LocErr", "OMT_LocErr",
    "OMT_Imprecision", "OMT_Misbinding", "OMT_Guessing",
}

# Percentile column points stored in the file
PERCENTILES = [5, 10, 25, 50, 75, 90, 95]


def _load_norms(path: Path = NORMS_PATH) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.set_index("age")
    return df


_NORMS = None


def _norms() -> pd.DataFrame:
    global _NORMS
    if _NORMS is None:
        _NORMS = _load_norms()
    return _NORMS


def _row_for_age(age: ArrayLike) -> pd.DataFrame:
    """Return a DataFrame of norm rows for the given age(s), nearest-neighbour."""
    norms = _norms()
    ages = np.atleast_1d(np.asarray(age, dtype=float))
    ages_clipped = np.clip(np.round(ages), norms.index.min(), norms.index.max())
    return norms.loc[ages_clipped]


def _ensure_metric(metric: str) -> None:
    if metric not in RECOMMENDED:
        raise KeyError(
            f"Unknown metric '{metric}'. Known: {sorted(RECOMMENDED)}"
        )


def age_z_score(metric: str, value: ArrayLike, age: ArrayLike) -> ArrayLike:
    """Return the z-score using the recommended scale (raw or log).

    Sign convention: higher z = better. For lower-is-better metrics (RTs,
    errors), the z is sign-flipped automatically.

    For metrics where the recommended approach is `percentile`, this falls
    back to a raw z-score with a printed warning; prefer `age_percentile`
    for those.
    """
    _ensure_metric(metric)
    approach = RECOMMENDED[metric]
    rows = _row_for_age(age)

    if approach == "log":
        mean = rows[f"{metric}_log_mean"].to_numpy()
        sd = rows[f"{metric}_log_sd"].to_numpy()
        v = np.log(np.asarray(value, dtype=float))
    else:
        if approach == "percentile":
            import warnings
            warnings.warn(
                f"{metric}: distribution is not normal; raw z is approximate."
                " Use age_percentile() for a robust standardised score."
            )
        mean = rows[f"{metric}_mean"].to_numpy()
        sd = rows[f"{metric}_sd"].to_numpy()
        v = np.asarray(value, dtype=float)

    with np.errstate(divide="ignore", invalid="ignore"):
        z = (v - mean) / sd
    if metric in LOWER_IS_BETTER:
        z = -z
    if isinstance(value, pd.Series):
        return pd.Series(np.asarray(z), index=value.index)
    if np.isscalar(value):
        return float(np.asarray(z).item())
    return np.asarray(z)


def age_percentile(metric: str, value: ArrayLike, age: ArrayLike) -> ArrayLike:
    """Return the percentile rank within the age-matched HC distribution.

    Higher percentile = better performance (for both higher-is-better and
    lower-is-better metrics — for the latter, a low raw value (e.g. fast RT)
    corresponds to a high percentile).
    """
    _ensure_metric(metric)
    rows = _row_for_age(age)
    cols = [f"{metric}_p{p}" for p in PERCENTILES]
    grid_p = np.array(PERCENTILES, dtype=float)
    grid_vals = rows[cols].to_numpy()

    v = np.atleast_1d(np.asarray(value, dtype=float))
    out = np.empty_like(v, dtype=float)
    for i, x in enumerate(v):
        row_vals = grid_vals[i] if grid_vals.ndim == 2 else grid_vals
        finite = ~np.isnan(row_vals)
        if not finite.any() or np.isnan(x):
            out[i] = np.nan
            continue
        rv = row_vals[finite]
        gp = grid_p[finite]
        # Interpolate in the metric→percentile mapping.
        # For lower-is-better metrics (RT, errors), invert so the percentile
        # always means "fraction of HC scoring worse than the patient".
        order = np.argsort(rv)
        rv_sorted = rv[order]
        gp_sorted = gp[order]
        if x <= rv_sorted[0]:
            pct = gp_sorted[0] / 2  # below smallest tabulated bracket
        elif x >= rv_sorted[-1]:
            pct = (gp_sorted[-1] + 100) / 2  # above largest
        else:
            pct = float(np.interp(x, rv_sorted, gp_sorted))
        if metric in LOWER_IS_BETTER:
            pct = 100 - pct
        out[i] = pct

    if isinstance(value, pd.Series):
        return pd.Series(out, index=value.index)
    if np.isscalar(value):
        return float(out[0])
    return out


# Convenience: vectorised wrappers for a DataFrame
def add_z_scores(df: pd.DataFrame, age_col: str = "age",
                 metrics: list[str] | None = None) -> pd.DataFrame:
    """Add `z_<metric>` columns to a DataFrame for every recognised metric
    present in the columns. Returns a new DataFrame (does not modify in place).
    """
    out = df.copy()
    targets = metrics or [m for m in RECOMMENDED if m in df.columns]
    for m in targets:
        if m not in df.columns:
            continue
        out[f"z_{m}"] = age_z_score(m, df[m], df[age_col])
    return out


def add_percentiles(df: pd.DataFrame, age_col: str = "age",
                    metrics: list[str] | None = None) -> pd.DataFrame:
    out = df.copy()
    targets = metrics or [m for m in RECOMMENDED if m in df.columns]
    for m in targets:
        if m not in df.columns:
            continue
        out[f"pct_{m}"] = age_percentile(m, df[m], df[age_col])
    return out


if __name__ == "__main__":
    # Small self-test
    print("Example: patient age 65")
    print("  DSST_RT = 4.2s:", age_z_score("DSST_RT", 4.2, 65),
          "(log z; higher z = better)")
    print("  DSST_RT pct:   ", age_percentile("DSST_RT", 4.2, 65), "(higher = better)")
    print("  OIS_LTM_Acc=55:", age_z_score("OIS_LTM_Acc", 55, 65),
          "(raw z; expect deprecation warning)")
    print("  OIS_LTM_Acc pct:", age_percentile("OIS_LTM_Acc", 55, 65))
