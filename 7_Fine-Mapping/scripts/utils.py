"""Shared utility functions for the meta-locus fine-mapping workflow."""

from __future__ import annotations

import pandas as pd
import numpy as np


def clean_chr_value(x):
    """Return chromosome labels as strings without chr prefix or .0 suffix."""
    if pd.isna(x):
        return np.nan
    x = str(x).strip().replace("chr", "").replace("CHR", "")
    try:
        f = float(x)
        if f.is_integer():
            return str(int(f))
    except Exception:
        pass
    return x


def clean_chr_series(s: pd.Series) -> pd.Series:
    return s.apply(clean_chr_value)


def to_bool(s: pd.Series) -> pd.Series:
    """Convert common string/number truthy values to boolean."""
    if s.dtype == bool:
        return s
    return s.astype(str).str.lower().isin(["true", "1", "t", "yes"])


def collapse_unique(x) -> str:
    vals = sorted(set(x.dropna().astype(str)))
    vals = [v for v in vals if v not in ["", ".", "NA", "N/A", "nan", "None"]]
    return ";".join(vals)


def ensure_cols(df: pd.DataFrame, cols, label="dataframe"):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def read_table(path: str, sep: str | None = None) -> pd.DataFrame:
    """Read comma/tab/whitespace table. If sep is None, pandas infers."""
    if sep is None:
        return pd.read_csv(path, sep=None, engine="python")
    return pd.read_csv(path, sep=sep)


def priority_sex(cat: str) -> str:
    if cat in ["sex_moderated only", "sex_main + sex_moderated"]:
        return "Sex moderated"
    if cat == "sex_main only":
        return "Sex main"
    return "No sex locus"


def priority_generation(cat: str) -> str:
    if cat in ["generation_moderated only", "generation_main + generation_moderated"]:
        return "Generation moderated"
    if cat == "generation_main only":
        return "Generation main"
    return "No generation locus"
